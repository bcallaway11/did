# edid-adaptive.R
# Adaptive event-study estimation under uncertain parallel trends
# (Section 5.2, Proposition 5.1 of Chen, Sant'Anna & Xie 2025), building on the
# minimax shrinkage function delta*(.; rho^2) of Armstrong, Kline & Sun
# (Econometrica 2025) as tabulated in their MissAdapt replication package.

# Load the vendored MissAdapt lookup tables. `aks_lookup.rds` is a faithful
# conversion (via R.matlab::readMat) of policy.mat / thresholds.mat /
# emse_corr.mat from the MissAdapt replication package (Armstrong, Kline & Sun,
# Econometrica 2025; github.com/lsun20/MissAdapt commit 98d823a, also archived
# as Zenodo record 16890198; MIT license, see inst/COPYRIGHTS); the original
# .mat files are shipped alongside it in inst/extdata/aks_lookup for
# provenance, and the conversion is built/checked by data-raw/aks_lookup.R
# (which embeds repo/commit/license/grid-convention attributes in the .rds).
# Contents: y_grid (481-point t_O grid on [-12, 12]), psi_mat
# (481 x 60 delta* policy values; ROWS = y-grid points, COLS = corr-grid
# points), st / ht (soft/hard-threshold lookups), mse_lambda (ERM lambda
# lookup), corr_grid (the 60-point |corr| grid abs(tanh(seq(-3, -0.05, 0.05))),
# decreasing, indexing the columns).
.edid_aks_lookup <- function() {
  path <- system.file("extdata", "aks_lookup", "aks_lookup.rds", package = "did")
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("The MissAdapt lookup tables (inst/extdata/aks_lookup/aks_lookup.rds) are not ",
         "available in this installation of the did package.", call. = FALSE)
  }
  tab <- readRDS(path)
  needed <- c("y_grid", "psi_mat", "st", "ht", "mse_lambda", "corr_grid")
  if (!all(needed %in% names(tab))) {
    stop("The MissAdapt lookup table file is malformed (missing components: ",
         paste(setdiff(needed, names(tab)), collapse = ", "), ").", call. = FALSE)
  }
  tab
}

# Core adaptive computation of Proposition 5.1 / eqn (5.4) on the bivariate
# scalar pair: YR / VR restricted (efficient) estimate and variance, YU / VU
# unrestricted (conservative), VUR their covariance. Ports the reference
# implementation used in the paper's empirical application exactly (same
# over-identification algebra, lookup keys, spline methods, and clamped
# off-grid interpolation). Returns the adaptive estimate and all components.
#
# `assume_efficient = TRUE` imposes the Hausman covariance identity VUR = VR
# (the Proposition 5.1 efficient-restricted special case, the same convention
# the MissAdapt README describes: "If one assumes that, in the absence of
# bias, the restricted estimator YR is efficient, then VUR can be set equal
# to VR"). Under the identity sigma_UO = -sigma_O^2, so corr =
# -sqrt(1 - VR/VU) and the GMM combination collapses to YR; both identities
# are imposed exactly (not via the generic floating-point algebra).
.edid_aks_core <- function(YR, VR, YU, VU, VUR, tables = .edid_aks_lookup(),
                           assume_efficient = FALSE) {
  stopifnot(is.numeric(YR), is.numeric(VR), is.numeric(YU), is.numeric(VU), is.numeric(VUR))
  if (isTRUE(assume_efficient)) VUR <- VR

  # ---- Over-identification direction (Proposition 5.1) ----
  YO  <- YR - YU                       # Y_O-hat = ES_avg-hat - ES_avg-check
  VO  <- VR - 2 * VUR + VU             # sigma_O^2 = sigma_R^2 - 2 sigma_UR + sigma_U^2 > 0
  VUO <- VUR - VU                      # sigma_UO = sigma_UR - sigma_U^2
  if (!is.finite(VO) || VO <= 0) {
    stop(sprintf(paste0("edid_adaptive: VO = VR - 2*VUR + VU = %.6g is not positive; the ",
                        "restricted/unrestricted variances admit no valid over-identification ",
                        "direction (Proposition 5.1 requires sigma_O^2 > 0)."), VO), call. = FALSE)
  }
  tO  <- YO / sqrt(VO)                 # over-identification statistic t_O

  # Efficient GMM (CUE) combination and its variance. Under assume_efficient
  # VUO = -VO exactly, so GMM = YU + YO = YR and V_GMM = VU - VO = VR; impose
  # the identities exactly rather than leaving them to FP rounding.
  GMM   <- YU - VUO / VO * YO
  V_GMM <- VU - VUO^2 / VO
  if (isTRUE(assume_efficient)) {
    GMM   <- YR
    V_GMM <- VR
  }

  # Correlation between the unrestricted estimator and the over-id direction;
  # delta*'s second argument is rho_AKS^2 = corr^2. The tables are indexed by
  # |corr| = |rho_AKS|, which is equivalent.
  corr       <- VUO / sqrt(VO) / sqrt(VU)
  rho_aks_sq <- corr^2

  corr_grid <- tables$corr_grid
  y_grid    <- tables$y_grid
  psi_mat   <- tables$psi_mat
  Ky        <- length(y_grid)

  # Off-grid guards: the lookup tables are tabulated only on these ranges.
  # Clamp to the grid (with a warning) so out-of-range inputs are loud rather
  # than a silent spline extrapolation; interior inputs are exact no-ops.
  corr_lo <- min(corr_grid); corr_hi <- max(corr_grid)
  acorr <- abs(corr)
  if (acorr < corr_lo || acorr > corr_hi) {
    warning(sprintf("edid_adaptive: |corr| = %.4f outside the tabulated grid [%.4f, %.4f]; clamping (extrapolation suppressed).",
                    acorr, corr_lo, corr_hi), call. = FALSE)
  }
  acorr <- min(max(acorr, corr_lo), corr_hi)
  tO_lo <- min(y_grid); tO_hi <- max(y_grid)
  if (tO < tO_lo || tO > tO_hi) {
    warning(sprintf("edid_adaptive: tO = %.4f outside the tabulated y-grid [%.4f, %.4f]; clamping for the nonlinear estimate.",
                    tO, tO_lo, tO_hi), call. = FALSE)
  }
  tO_eval <- min(max(tO, tO_lo), tO_hi)

  # ---- Nonlinear adaptive estimate (delta* of AKS Theorem 1(ii)) ----
  # For each tO grid point, interpolate psi across correlation values, then
  # across tO (the MissAdapt two-stage spline interpolation).
  psi_grid <- numeric(Ky)
  for (i in seq_len(Ky)) {
    psi_fun <- stats::splinefun(corr_grid, psi_mat[i, ], method = "fmm", ties = mean)
    psi_grid[i] <- psi_fun(acorr)
  }
  psi_extrap <- stats::splinefun(y_grid, psi_grid, method = "natural")
  t_tilde  <- psi_extrap(tO_eval)                       # delta*(t_O; rho_AKS^2)
  adaptive <- VUO / sqrt(VO) * t_tilde + GMM            # eqn (5.4)

  # ---- Soft-threshold variant ----
  st_fun <- stats::splinefun(corr_grid, tables$st, method = "fmm", ties = mean)
  st <- st_fun(acorr)
  adaptive_st <- VUO / sqrt(VO) * ((tO > st) * (tO - st) + (tO < -st) * (tO + st)) + GMM

  # ---- Hard-threshold variant ----
  ht_fun <- stats::splinefun(corr_grid, tables$ht, method = "fmm", ties = mean)
  ht <- ht_fun(acorr)
  adaptive_ht <- VUO / sqrt(VO) * ((tO > ht) * tO + (tO < -ht) * tO) + GMM

  # ---- Pre-test estimate (hard threshold at 1.96) ----
  pretest <- VUO / sqrt(VO) * ((tO > 1.96) * tO + (tO < -1.96) * tO) + GMM

  # ---- ERM estimates ----
  erm <- VUO / sqrt(VO) * tO * (tO^2 / (tO^2 + 1)) + GMM
  lambda_fun <- stats::splinefun(corr_grid, tables$mse_lambda, method = "fmm", ties = mean)
  erm_lambda <- lambda_fun(acorr)
  adaptive_erm <- VUO / sqrt(VO) * tO * (tO^2 / (tO^2 + erm_lambda)) + GMM

  list(
    # Inputs (VUR is the value actually used: VR when assume_efficient)
    YR = YR, YU = YU, VR = VR, VU = VU, VUR = VUR,
    assume_efficient = isTRUE(assume_efficient),
    # Over-identification components (Proposition 5.1)
    YO = YO, VO = VO, VUO = VUO, tO = tO, corr = corr, rho_aks_sq = rho_aks_sq,
    # Efficient GMM combination
    GMM = GMM, V_GMM = V_GMM, se_GMM = sqrt(V_GMM),
    # Adaptive estimates
    adaptive = adaptive, adaptive_nonlinear = adaptive,
    adaptive_st = adaptive_st, adaptive_ht = adaptive_ht,
    pretest = pretest, erm = erm, adaptive_erm = adaptive_erm,
    # Thresholds
    soft_threshold = st, hard_threshold = ht, erm_lambda = erm_lambda
  )
}

#' Adaptive event-study estimator under uncertain parallel trends
#'
#' Implements the adaptive event-study estimator of Proposition 5.1 in Chen,
#' Sant'Anna & Xie (2025), eqn (5.4), which adapts the construction of
#' Armstrong, Kline & Sun (Econometrica 2025) to the bivariate pair
#' \eqn{(\widecheck{ES}_{avg}, \widehat{ES}_{avg})}: the efficient PT-All
#' estimator plays the role of the \emph{restricted} estimator (efficient under
#' the additional restrictions; biased if they fail) and the conservative
#' PT-Post estimator the role of the \emph{unrestricted} estimator. With
#' \eqn{\widehat\sigma_O^2 = \widehat\sigma_R^2 - 2\widehat\sigma_{UR} +
#' \widehat\sigma_U^2}, \eqn{\widehat\sigma_{UO} = \widehat\sigma_{UR} -
#' \widehat\sigma_U^2}, \eqn{\widehat{t}_O = (\widehat{ES}_{avg} -
#' \widecheck{ES}_{avg})/\widehat\sigma_O}, and \eqn{\widehat\rho^2_{AKS} =
#' \widehat\sigma_{UO}^2 / (\widehat\sigma_U^2 \widehat\sigma_O^2)}, the
#' adaptive estimator is
#' \deqn{\widehat{ES}_{avg}^{AKS} = \widecheck{ES}_{avg} +
#'   \frac{\widehat\sigma_{UO}}{\widehat\sigma_O}\left[\delta^*(\widehat{t}_O;
#'   \widehat\rho^2_{AKS}) - \widehat{t}_O\right],}
#' where \eqn{\delta^*} is the smooth minimax shrinkage function of Armstrong,
#' Kline & Sun, interpolated from the lookup tables of their MissAdapt
#' replication package. It minimizes the worst-case ratio of actual to oracle
#' mean-squared error over all bias bounds simultaneously, avoiding the
#' variance discontinuity that hard pre-testing pays.
#'
#' @inheritParams edid_hausman
#' @param parameter \code{"overall"} (default) applies the construction to the
#'   paper's headline scalar pair \eqn{(\widecheck{ES}_{avg},
#'   \widehat{ES}_{avg})}; \code{"event_study"} applies the same bivariate
#'   construction to each post-treatment \eqn{ES(e)} separately (the per-\eqn{e}
#'   analogue mentioned in the paper, with the usual multiple-comparison
#'   caveats).
#' @param assume_efficient logical or \code{NULL} (default \code{NULL} = AUTO).
#'   Selects how the covariance \eqn{\widehat\sigma_{UR}} between the two
#'   estimators is obtained. \code{FALSE}: estimate it empirically as the
#'   (cluster-robust) sample covariance of the two aggregations' per-unit
#'   influence functions, alongside the two variances. \code{TRUE}: impose the
#'   Hausman covariance identity \eqn{\widehat\sigma_{UR} =
#'   \widehat\sigma^2_R} -- the Proposition 5.1 reduction, justified exactly
#'   when the restricted estimator attains the efficiency bound (equivalently
#'   the MissAdapt README's "if one assumes that, in the absence of bias, the
#'   restricted estimator is efficient, then \code{VUR} can be set equal
#'   to \code{VR}"). Under the identity \eqn{\widehat\sigma^2_O =
#'   \widehat\sigma^2_U - \widehat\sigma^2_R}, \eqn{\widehat\rho_{AKS} =
#'   -\sqrt{1 - \widehat\sigma^2_R/\widehat\sigma^2_U}}, and the GMM
#'   combination collapses to the restricted estimate exactly
#'   (\code{GMM == YR}, \code{V_GMM == VR}).
#'
#'   \code{NULL} (AUTO, the default) resolves scheme-aware from the restricted
#'   fit: \code{TRUE} when \code{fit_restricted} is bound-attaining --
#'   \code{weight_scheme = "efficient"}, OR a no-covariate fit with any
#'   \code{weight_scheme} other than \code{"uniform"} (with no covariates the
#'   efficient, averaged, and gmm weights coincide and attain the bound) --
#'   and \code{FALSE} otherwise. The rationale: the empirical
#'   influence-function recombination is the unconditional-GMM-style
#'   construction, which coincides with the efficient-influence-function
#'   covariance only without covariates; when the restricted estimator
#'   attains the bound, the Proposition 5.1 identity is the asymptotically
#'   exact convention and is imposed, while a non-bound-attaining restricted
#'   fit (uniform weights; or constant-weight schemes with covariates) has no
#'   such identity and keeps the internally consistent empirical covariance.
#'   An explicit \code{TRUE}/\code{FALSE} always overrides the AUTO rule. See
#'   Details.
#'
#' @details
#' All variances and covariances are computed from the per-unit influence
#' functions of the two aggregations (cluster-robust when the fits carry
#' cluster assignments). The construction requires
#' \eqn{\widehat\sigma_O^2 > 0}; when the two estimators coincide (no
#' over-identification direction, e.g. a just-identified design) the function
#' stops with an informative error. Inputs whose \eqn{|corr|} or
#' \eqn{\widehat{t}_O} fall outside the tabulated grids are clamped to the grid
#' boundary with a warning (no silent spline extrapolation).
#'
#' \strong{The two covariance conventions.} When the restricted fit attains
#' the efficiency bound under its maintained restrictions, the Hausman
#' covariance identity \eqn{\widehat\sigma_{UR} = \widehat\sigma^2_R} holds
#' asymptotically (the Proposition 5.1 reduction), so the empirical-covariance
#' convention (\code{assume_efficient = FALSE}) and the imposed-identity
#' convention (\code{assume_efficient = TRUE}) coincide in the limit and the
#' two adaptive estimates converge to the same value. In finite samples the
#' empirical influence-function covariance deviates from
#' \eqn{\widehat\sigma^2_R} by sampling noise, so the two conventions give
#' (slightly) different \eqn{\widehat{t}_O}, \eqn{\widehat\rho^2_{AKS}}, and
#' adaptive estimates. \code{assume_efficient = TRUE} reproduces the
#' convention used in the MissAdapt \code{example.R} (\code{VUR <- VR}) and
#' requires \eqn{\widehat\sigma^2_R < \widehat\sigma^2_U} (otherwise
#' \eqn{\widehat\sigma_O^2 \le 0} and the function stops). Note that
#' \code{assume_efficient = TRUE} is an \emph{assumption}, not an estimate: it
#' is justified exactly when the restricted estimator attains the efficiency
#' bound, and the AUTO default (\code{NULL}) imposes it only then -- for a
#' bound-attaining restricted fit (\code{weight_scheme = "efficient"}, or any
#' non-uniform scheme without covariates, where the empirical
#' unconditional-GMM-style recombination coincides with the
#' efficient-influence-function construction anyway). If the restricted fit
#' is not bound-attaining (e.g. comparing two conservative fits, or
#' constant-weight schemes with covariates), the identity is wrong and AUTO
#' keeps the empirical covariance.
#'
#' \strong{No inference accompanies the adaptive estimate.} Proposition 5.1 is
#' a point-estimation (risk) result; the paper makes no coverage claims for
#' \eqn{\widehat{ES}_{avg}^{AKS}}, and no standard error is reported for it.
#' Use the efficient or conservative fits (or the robustness frontier,
#' \code{\link{edid_frontier}}) for inference.
#'
#' \strong{Lookup-table provenance.} The shipped tables
#' (\code{inst/extdata/aks_lookup/}) are the \code{policy.mat},
#' \code{thresholds.mat}, and \code{emse_corr.mat} lookup tables of the
#' MissAdapt replication package of Armstrong, Kline & Sun (Econometrica 2025),
#' vendored byte-identically from
#' \url{https://github.com/lsun20/MissAdapt} (commit \code{98d823a}; also
#' archived as Zenodo record 16890198) and distributed under the MIT license
#' (Copyright (c) 2023 Sophie Sun; see \code{inst/COPYRIGHTS}). The
#' \code{aks_lookup.rds} conversion the function reads (so no MATLAB-file
#' reader is required at runtime) is built by \code{data-raw/aks_lookup.R},
#' which documents the grid conventions and runs orientation/monotonicity/
#' symmetry sanity checks, including a regression against the published
#' MissAdapt vignette example; the provenance, commit, license, and grid
#' conventions are embedded as attributes of the \code{.rds}.
#'
#' @return An object of class \code{edid_adaptive}. For
#'   \code{parameter = "overall"}: a list with the adaptive estimate
#'   (\code{adaptive}, the eqn (5.4) nonlinear estimator), the components
#'   (\code{YU}, \code{YR}, \code{VU}, \code{VR}, \code{VUR}, \code{YO},
#'   \code{VO}, \code{VUO}, \code{tO}, \code{corr}, \code{rho_aks_sq}), the
#'   efficient GMM combination (\code{GMM}, \code{V_GMM}, \code{se_GMM}), the
#'   soft-threshold / hard-threshold / pre-test / ERM variants, and the
#'   interpolated thresholds. For \code{parameter = "event_study"}: the same
#'   quantities as a per-\eqn{e} data.frame in \code{$table}. \code{VUR} is
#'   the covariance actually used (equal to \code{VR} under
#'   \code{assume_efficient = TRUE}, in which case \code{GMM} equals \code{YR}
#'   exactly); \code{$assume_efficient} records the RESOLVED convention and
#'   \code{$assume_efficient_auto} whether it came from the AUTO rule.
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025). Efficient
#'   Difference-in-Differences and Event Study Estimators. Section 5.2,
#'   Proposition 5.1. \cr
#'   Armstrong, T. B., Kline, P., & Sun, L. (2025). Adapting to
#'   Misspecification. \emph{Econometrica}, 93(6), 1981-2005. Replication
#'   package: MissAdapt, Zenodo 16890198.
#'
#' @seealso \code{\link{edid}}, \code{\link{edid_hausman}},
#'   \code{\link{edid_frontier}}
#'
#' @examples
#' \donttest{
#' df <- data.frame(
#'   id   = rep(1:120, each = 6),
#'   time = rep(1:6, 120),
#'   g    = rep(sample(c(3, 5, Inf), 120, replace = TRUE), each = 6)
#' )
#' df$y <- rnorm(120)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
#'   rnorm(nrow(df), 0, 0.5)
#' fit_R <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
#'               aggregate = "event_study", cband = FALSE)
#' fit_U <- edid(df, "y", "id", "time", "g", pt_assumption = "post",
#'               aggregate = "event_study", cband = FALSE)
#' edid_adaptive(fit_U, fit_R)
#' }
#'
#' @export
edid_adaptive <- function(fit_unrestricted, fit_restricted,
                          parameter = c("overall", "event_study"),
                          e_set = NULL, assume_efficient = NULL) {
  parameter <- match.arg(parameter)
  assume_efficient_auto <- is.null(assume_efficient)
  if (!assume_efficient_auto) {
    stopifnot(is.logical(assume_efficient), length(assume_efficient) == 1L,
              !is.na(assume_efficient))
  }
  .edid_toolkit_check_fits(fit_unrestricted, fit_restricted)
  if (assume_efficient_auto) {
    # AUTO (scheme-aware): impose the Proposition 5.1 identity VUR = VR exactly when the restricted fit is
    # bound-attaining -- weight_scheme = "efficient", OR a no-covariate fit with any non-uniform scheme
    # (without covariates the efficient/averaged/gmm weights coincide and attain the bound). Otherwise keep
    # the empirical influence-function covariance. Explicit TRUE/FALSE always bypasses this rule.
    ws      <- fit_restricted$weight_scheme
    has_cov <- !is.null(fit_restricted$xformla) && inherits(fit_restricted$xformla, "formula") &&
               length(all.vars(fit_restricted$xformla)) > 0L
    assume_efficient <- identical(ws, "efficient") ||
      (!has_cov && !is.null(ws) && !identical(ws, "uniform"))
  }

  n      <- fit_restricted$n
  ci     <- fit_restricted$cluster_indices
  tables <- .edid_aks_lookup()

  # Bivariate variance of the (unrestricted, restricted) ESTIMATES:
  # cluster-robust sandwich of the two aggregation IFs, order 1/n (the same
  # scale as the reference implementation's mean(if^2)/n).
  .biv <- function(psi_U, psi_R) {
    S <- cluster_cov_edid(cbind(psi_U, psi_R), ci, n)
    list(VU = S[1L, 1L], VR = S[2L, 2L], VUR = S[1L, 2L])
  }

  if (parameter == "overall") {
    oU <- .edid_param_ifs(fit_unrestricted, "overall")
    oR <- .edid_param_ifs(fit_restricted,  "overall")
    v  <- .biv(oU$IF[, 1L], oR$IF[, 1L])
    out <- .edid_aks_core(YR = oR$est, VR = v$VR, YU = oU$est, VU = v$VU, VUR = v$VUR,
                          tables = tables, assume_efficient = assume_efficient)
    out$parameter <- "overall"
  } else {
    e_set <- .edid_shared_e_set(fit_unrestricted, fit_restricted, e_set)
    pU <- .edid_param_ifs(fit_unrestricted, "event_study", e_set)
    pR <- .edid_param_ifs(fit_restricted,  "event_study", e_set)
    cols <- c("YU", "YR", "VU", "VR", "VUR", "YO", "VO", "VUO", "tO", "corr", "rho_aks_sq",
              "GMM", "se_GMM", "adaptive", "adaptive_st", "adaptive_ht", "pretest",
              "adaptive_erm", "soft_threshold", "hard_threshold")
    rows <- vector("list", length(e_set))
    for (j in seq_along(e_set)) {
      v <- .biv(pU$IF[, j], pR$IF[, j])
      cj <- .edid_aks_core(YR = pR$est[j], VR = v$VR, YU = pU$est[j], VU = v$VU, VUR = v$VUR,
                           tables = tables, assume_efficient = assume_efficient)
      rows[[j]] <- cbind(data.frame(e = e_set[j]), as.data.frame(cj[cols]))
    }
    out <- list(parameter = "event_study", e_set = e_set,
                table = do.call(rbind, rows))
    rownames(out$table) <- NULL
  }

  out$assume_efficient <- assume_efficient
  out$assume_efficient_auto <- assume_efficient_auto
  out$n <- n
  out$clustered <- !is.null(ci)
  class(out) <- c("edid_adaptive", "list")
  out
}

#' @describeIn edid_adaptive Print method.
#' @param x an \code{edid_adaptive} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_adaptive <- function(x, digits = 4, ...) {
  cat("\nAdaptive event-study estimator (Chen, Sant'Anna & Xie 2025, Proposition 5.1;\n")
  cat("Armstrong, Kline & Sun 2025)\n")
  if (isTRUE(x$assume_efficient)) {
    cat("  Covariance convention: imposed Hausman identity sigma_UR = sigma_R^2\n")
    cat(sprintf("  (assume_efficient = TRUE%s; GMM = YR by construction)\n",
                if (isTRUE(x$assume_efficient_auto)) " [AUTO: restricted fit is bound-attaining]" else ""))
  } else if (isTRUE(x$assume_efficient_auto)) {
    cat("  Covariance convention: empirical influence-function covariance\n")
    cat("  (assume_efficient = FALSE [AUTO: restricted fit is not bound-attaining])\n")
  }
  if (identical(x$parameter, "overall")) {
    cat(sprintf("  Parameter: ES_avg%s\n", if (isTRUE(x$clustered)) " (cluster-robust)" else ""))
    cat(sprintf("  Conservative (unrestricted): %s   Efficient (restricted): %s\n",
                format(x$YU, digits = digits), format(x$YR, digits = digits)))
    cat(sprintf("  t_O = %s, rho_AKS^2 = %s\n",
                format(x$tO, digits = digits), format(x$rho_aks_sq, digits = digits)))
    cat(sprintf("  Adaptive estimate (eqn 5.4): %s\n", format(x$adaptive, digits = digits)))
    cat(sprintf("  [GMM: %s; soft-threshold: %s; hard-threshold: %s; pre-test: %s]\n",
                format(x$GMM, digits = digits), format(x$adaptive_st, digits = digits),
                format(x$adaptive_ht, digits = digits), format(x$pretest, digits = digits)))
  } else {
    cat(sprintf("  Parameter: ES(e) per event time%s\n",
                if (isTRUE(x$clustered)) " (cluster-robust)" else ""))
    tab <- x$table[, c("e", "YU", "YR", "tO", "rho_aks_sq", "adaptive", "GMM")]
    num <- vapply(tab, is.numeric, logical(1L))
    tab[num] <- lapply(tab[num], function(z) signif(z, digits))
    print(tab, row.names = FALSE)
  }
  cat("\nNote: the adaptive estimate is a point-estimation (risk) construction; no\n")
  cat("standard error or coverage claim accompanies it (see ?edid_adaptive).\n")
  invisible(x)
}
