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
# decreasing, indexing the columns), and the B-FLCI critical-value tables of
# AKS Section 4.2.2: flci_B_grid (the 91-point B-tilde = B/sigma_O grid
# c(0.01, seq(0.1, 9, 0.1)) indexing the rows), flci_cv_adaptive /
# flci_cv_st (91 x 60 c_.05 tables for the adaptive and soft-threshold
# estimators; columns = corr_grid points), flci_minimax_B_grid /
# flci_cv_minimax (90-row analogue for the B-minimax estimator -- vendored
# for completeness but consumed by no exported interface, because
# edid_adaptive() computes no B-minimax point estimate; internal hook only).
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
    stop(sprintf(paste0("edid_adaptive: VO = VR - 2*VUR + VU = %.6g is not positive, so there is ",
                        "no over-identification direction to adapt over (Proposition 5.1 requires ",
                        "sigma_O^2 > 0). This typically means the efficient (restricted) and ",
                        "conservative (unrestricted) fits coincide -- common when the thin-cohort ",
                        "guard has pinned every cell to its just-identified moment (check the fits' ",
                        "$thin_cohorts and warnings, or edid_weights()) -- or, under ",
                        "assume_efficient = TRUE, that the restricted fit is not empirically more ",
                        "precise than the unrestricted one (VR >= VU)."), VO), call. = FALSE)
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
    soft_threshold = st, hard_threshold = ht, erm_lambda = erm_lambda,
    # Quantities the B-FLCI layer (.edid_aks_ci) reuses so that the CIs
    # inherit exactly the clamping conventions and interpolants of the
    # estimates they are centered at: the clamped |corr| actually used for
    # every corr-grid lookup, and the two-stage psi interpolant
    # delta*(.; rho_AKS^2) on the y grid.
    acorr_eval = acorr, psi_fun = psi_extrap
  )
}

# ============================================================================
# B-FLCI inference layer: the fixed-length confidence intervals of Armstrong,
# Kline & Sun (Econometrica 2025), Section 4.2.2 (arXiv v6 numbering; their
# eqs. (7)-(8)), centered at the adaptive (and soft-threshold) estimates.
# Every interval is {center +- c_.05(B/sigma_O; rho, delta) * sigma_U} with
# sigma_U = sqrt(VU), the UNRESTRICTED estimator's standard error -- never
# se_GMM, never sqrt(VR): because the local bias b cannot be consistently
# estimated, the adaptive estimator has no consistently estimable asymptotic
# distribution and hence no own standard error (AKS Section 4.2.1).
# ============================================================================

# The b-tilde (= b/sigma_O) grid on which coverage is evaluated and on which
# the inner sup of AKS eq. (8) is taken: seq(-9, 9, 0.025), exactly the
# `b.grid` of MissAdapt's risk.mat (721 points; generated inline rather than
# vendored).
.edid_aks_b_grid <- function() seq(-9, 9, by = 0.025)

# Deterministic Gaussian-quadrature coverage of {estimate +- cval * sigma_U}
# at scaled bias b-tilde, from the AKS distributional representation (their
# eq. (7)):
#   (theta_hat - theta)/sigma_U = rho [delta(Z1 + b) - b] + sqrt(1 - rho^2) Z2,
# so conditional on Z1 = z the pivot is normal with mean
# m(z) = rho (delta(z + b) - b) and sd s = sqrt(1 - rho^2), giving
#   coverage(b) = E_Z1[ Phi((cval - m)/s) - Phi((-cval - m)/s) ].
# The Z1 expectation uses the fixed grid z = seq(-8.5, 8.5, 0.01) with
# dnorm(z) * 0.01 weights. This deterministic quadrature replaces MissAdapt's
# set.seed(1), 100000-draw Monte Carlo (agreement within ~4e-3, the MC noise
# level) and is the package's single coverage convention. delta_fun is the
# shrinkage policy: the psi interpolant for the adaptive estimator, the
# soft-threshold map for the soft-threshold estimator.
.edid_aks_flci_coverage <- function(delta_fun, corr, cval, b_grid) {
  z  <- seq(-8.5, 8.5, by = 0.01)
  wz <- stats::dnorm(z) * 0.01
  s  <- sqrt(1 - corr^2)
  vapply(b_grid, function(b) {
    m <- corr * (delta_fun(z + b) - b)
    sum((stats::pnorm((cval - m) / s) - stats::pnorm((-cval - m) / s)) * wz)
  }, numeric(1L))
}

# Solve AKS eq. (8) for the soft-threshold estimator at the CORRECT adaptive
# soft threshold lambda = lambda*(rho) -- the thresholds.mat value that
# defines the soft-threshold ESTIMATE the interval is centered at: the
# smallest critical value c such that
#   min over |b-tilde| <= B_tilde (risk.mat b grid) of coverage(b-tilde; c)
# is at least `level`. Coverage is nondecreasing in c for every b-tilde, so
# bisection applies; the returned upper bracket guarantees min coverage >=
# level up to the bisection tolerance. The b-tilde sup uses the nonnegative
# half of the grid only: coverage(b) = coverage(-b) exactly, because delta_S
# is odd and the z grid is symmetric.
#
# This runtime solve exists because the SHIPPED flci_adaptive_st_cv.mat is
# calibrated to a different (smaller, off-grid-extrapolated) soft threshold
# than the one the soft-threshold estimate uses -- a signed-vs-absolute grid
# mispairing in MissAdapt's calculate_B_FLCI.R (its line 36 reassigns the
# interpolation grid to the signed tanh grid while evaluating at abs(corr)).
# At the correct lambda*(rho), the shipped critical values can undercover
# within |b| <= B (e.g. min coverage 0.74 at rho = -0.995, B-tilde = 9); the
# st_cv = "missadapt" option keeps the shipped table for exact replication.
.edid_aks_st_cv_exact <- function(corr, lambda, B_tilde, level = 0.95) {
  z  <- seq(-8.5, 8.5, by = 0.01)
  wz <- stats::dnorm(z) * 0.01
  s  <- sqrt(1 - corr^2)
  bg <- .edid_aks_b_grid()
  b  <- bg[bg >= 0 & bg <= B_tilde + 1e-12]
  # delta_S(z + b) - b is independent of c: precompute the conditional means
  TB <- outer(z, b, "+")
  D  <- (TB > lambda) * (TB - lambda) + (TB < -lambda) * (TB + lambda)
  M  <- corr * (D - matrix(b, nrow = length(z), ncol = length(b), byrow = TRUE))
  cov_min <- function(cc) {
    min(colSums((stats::pnorm((cc - M) / s) - stats::pnorm((-cc - M) / s)) * wz))
  }
  lo <- 0; hi <- 6
  if (cov_min(hi) < level) {
    stop("edid_adaptive: the soft-threshold eq.-(8) critical value exceeds 6; ",
         "this cannot occur on the tabulated (corr, B) grids.", call. = FALSE)
  }
  while (hi - lo > 1e-7) {
    mid <- (lo + hi) / 2
    if (cov_min(mid) >= level) hi <- mid else lo <- mid
  }
  hi
}

# Critical-value lookup c_.05(B_tilde; |corr|) from a vendored 91 x 60 cv
# table: fmm spline of row `row` against corr_grid (the decreasing |corr|
# grid), evaluated at the clamped |corr|. This is numerically IDENTICAL
# (exact fmm mirror symmetry, asserted in data-raw/aks_lookup.R) to
# MissAdapt's spline of the same row against the signed increasing grid
# tanh(seq(-3, -0.05, 0.05)) evaluated at the signed (negative) correlation,
# while remaining on-grid for corr > 0 inputs (the signed-grid convention
# would silently extrapolate far off-grid there).
.edid_aks_flci_cv <- function(cv_mat, row, corr_grid, acorr) {
  stats::splinefun(corr_grid, cv_mat[row, ], method = "fmm", ties = mean)(acorr)
}

# Resolve requested bias bounds B (in sigma_O units, i.e. B-tilde) to rows of
# the FLCI tables. Conventions (documented in ?edid_adaptive):
#   B = 0   -> row 1 (B-tilde = 0.01), the authors' own `B == 0` mapping;
#   B = Inf -> the B-tilde = 9 row, AKS's finite approximation to the
#              infinity-FLCI ("we approximate an infinity-FLCI by setting
#              B = 9 sigma_O");
#   otherwise snap to the grid within 1e-8. MissAdapt's own lookup uses an
#   exact floating-point match, which crashes for requests like B = 0.3 or
#   0.7 (the Matlab-written grid doubles differ from R literals in the last
#   bit); the tolerance match accepts those without changing any on-grid row.
.edid_aks_flci_rows <- function(B, B_grid) {
  if (!is.numeric(B) || length(B) == 0L || anyNA(B)) {
    stop("edid_adaptive: B must be a numeric vector of bias bounds (sigma_O units) ",
         "with no missing values.", call. = FALSE)
  }
  if (any(B < 0)) {
    stop("edid_adaptive: bias bounds B must be nonnegative (B is the bound on ",
         "|b|/sigma_O).", call. = FALSE)
  }
  idx <- vapply(B, function(b) {
    if (b == 0) return(1L)
    if (is.infinite(b)) return(length(B_grid))
    j <- which.min(abs(B_grid - b))
    if (abs(B_grid[j] - b) > 1e-8) {
      stop(sprintf(paste0("edid_adaptive: B = %g is not on the tabulated B-tilde grid; ",
                          "valid values are 0, Inf, 0.01, and 0.1, 0.2, ..., 9 ",
                          "(within 1e-8)."), b), call. = FALSE)
    }
    j
  }, integer(1L))
  data.frame(B = B, row = idx, B_tilde = B_grid[idx])
}

# Assemble the B-FLCI table for one .edid_aks_core() result: one row per
# (variant, B), every interval centered at that variant's point estimate and
# scaled by sigma_U = sqrt(VU).
.edid_aks_ci <- function(core, tables, B, st_cv, level = 0.95) {
  rows    <- .edid_aks_flci_rows(B, tables$flci_B_grid)
  sigma_U <- sqrt(core$VU)
  acorr   <- core$acorr_eval
  corr_eval <- if (core$corr < 0) -acorr else acorr  # clamped, signed

  cv_ad <- vapply(rows$row, function(r) {
    .edid_aks_flci_cv(tables$flci_cv_adaptive, r, tables$corr_grid, acorr)
  }, numeric(1L))
  cv_st <- if (identical(st_cv, "exact")) {
    vapply(rows$B_tilde, function(bt) {
      .edid_aks_st_cv_exact(corr_eval, core$soft_threshold, bt, level)
    }, numeric(1L))
  } else {
    vapply(rows$row, function(r) {
      .edid_aks_flci_cv(tables$flci_cv_st, r, tables$corr_grid, acorr)
    }, numeric(1L))
  }

  out <- rbind(
    data.frame(variant = "adaptive", B = rows$B, B_tilde = rows$B_tilde,
               center = core$adaptive, sigma_U = sigma_U, cv = cv_ad,
               lower = core$adaptive - cv_ad * sigma_U,
               upper = core$adaptive + cv_ad * sigma_U,
               cv_source = "table"),
    data.frame(variant = "soft_threshold", B = rows$B, B_tilde = rows$B_tilde,
               center = core$adaptive_st, sigma_U = sigma_U, cv = cv_st,
               lower = core$adaptive_st - cv_st * sigma_U,
               upper = core$adaptive_st + cv_st * sigma_U,
               cv_source = if (identical(st_cv, "exact")) "exact" else "missadapt_table")
  )
  rownames(out) <- NULL
  out
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
#' @param ci logical (default \code{TRUE}): attach the fixed-length confidence
#'   intervals (B-FLCIs) of Armstrong, Kline & Sun (2025, Section 4.2 of the
#'   published version; Section 4.2.2 and eqs. (7)-(8) in the arXiv-v6
#'   numbering used below) to the adaptive and soft-threshold estimates. See
#'   Details ("Inference: the AKS B-FLCIs").
#' @param B numeric vector or \code{NULL}: additional bias bounds, in
#'   \eqn{\widehat\sigma_O} units (\eqn{\widetilde{B} = B/\sigma_O}, the units
#'   of the AKS tables), at which to report B-FLCIs. The 0-FLCI
#'   (\code{B = 0}) and the \eqn{\infty}-FLCI (\code{B = Inf}) are always
#'   reported, per the AKS recommendation to "report alongside an adaptive
#'   estimate the critical values for a 0-FLCI and \eqn{\infty}-FLCI, thereby
#'   summarizing the range of critical values needed to guarantee coverage
#'   under different assumptions". Values must lie on the tabulated grid
#'   \{0.01, 0.1, 0.2, ..., 9\} (matched within 1e-8, so \code{B = 0.3}
#'   works; MissAdapt's own exact floating-point match crashes there);
#'   \code{B = 0} maps to the \eqn{\widetilde{B} = 0.01} row (the authors'
#'   convention) and \code{B = Inf} to the \eqn{\widetilde{B} = 9} row (AKS:
#'   "we approximate an \eqn{\infty}-FLCI by setting \eqn{B = 9\sigma_O}").
#' @param level confidence level; must be \code{0.95}. The MissAdapt
#'   critical-value tables are tabulated for the 95\% level only (no other
#'   level is tabulated anywhere in their package), so any other value is an
#'   error.
#' @param st_cv \code{"exact"} (default) or \code{"missadapt"}: source of the
#'   critical value for the \emph{soft-threshold} B-FLCI. \code{"exact"}
#'   solves AKS eq. (8) at runtime (deterministic quadrature + bisection) at
#'   the same soft threshold \eqn{\lambda^*(\rho)} that defines the
#'   soft-threshold estimate the interval is centered at. \code{"missadapt"}
#'   reproduces the shipped \code{flci_adaptive_st_cv.mat} critical values
#'   exactly, for comparability with MissAdapt output. The two differ because
#'   the shipped table is calibrated to a different threshold than the
#'   estimate: \code{calculate_B_FLCI.R} in MissAdapt interpolates the soft
#'   threshold for its coverage simulation against the \emph{signed}
#'   correlation grid while evaluating at \code{abs(corr)}, extrapolating
#'   beyond the grid and yielding a threshold of about 0.45-0.54 for every
#'   \eqn{\rho} instead of \eqn{\lambda^*(\rho)} from \code{thresholds.mat}
#'   (which reaches 1.12 at \eqn{|\rho| = 0.97}). At the correct
#'   \eqn{\lambda^*(\rho)}, the shipped critical values can undercover within
#'   \eqn{|b| \le B} (quadrature minimum coverage 0.74 at \eqn{\rho = -0.995},
#'   \eqn{\widetilde{B} = 9}, versus the nominal 0.95). The \emph{adaptive}
#'   (nonlinear) cv table has no such issue and is always used as shipped.
#'
#' @details
#' All variances and covariances are computed from the per-unit influence
#' functions of the two aggregations (cluster-robust when the fits carry
#' cluster assignments). The construction requires
#' \eqn{\widehat\sigma_O^2 > 0}; when the two estimators coincide (no
#' over-identification direction --- e.g. a just-identified design, or every
#' cell pinned to its just-identified moment by the thin-cohort guard, see
#' \code{min_pair_units} in \code{\link{edid}}) the function stops with an
#' informative error. Inputs whose \eqn{|corr|} or
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
#' \strong{Inference: the AKS B-FLCIs.} Proposition 5.1 itself is a
#' point-estimation (risk) result, and no conventional standard error attaches
#' to the adaptive estimate: in the AKS normal limit experiment the local bias
#' \eqn{b} of the restricted estimator cannot be consistently estimated, so
#' neither can the asymptotic distribution of the adaptive estimator
#' (Armstrong, Kline & Sun 2025, Section 4.2). Armstrong, Kline & Sun instead
#' construct \emph{fixed-length confidence intervals} (B-FLCIs)
#' \deqn{\{\hat\theta \pm c_{\alpha}(B/\sigma_O;\, \rho, \delta)\,\sigma_U\},}
#' centered at the adaptive (or soft-threshold) estimate and scaled by
#' \eqn{\sigma_U = \sqrt{VU}}, the \emph{unrestricted} (conservative)
#' estimator's standard error. The critical value solves their eq. (8)
#' (arXiv-v6 numbering): the smallest \eqn{\chi} such that
#' \eqn{\sup_{|\tilde b| \le \tilde B} P(|\rho[\delta(Z_1+\tilde b)-\tilde b]
#' + \sqrt{1-\rho^2} Z_2| > \chi) \le \alpha}, using the distributional
#' representation of their eq. (7). The exact guarantee is: \emph{in the
#' normal limit experiment with known covariance matrix, the B-FLCI covers the
#' target with probability at least \eqn{1-\alpha} for every \eqn{(\theta, b)}
#' with \eqn{|b| \le B}} -- i.e., uniformly over the bias \eqn{b} in
#' \eqn{[-B, B]} (with \eqn{B} in absolute units; \eqn{\widetilde{B} =
#' B/\sigma_O} in the tables' units), for all \eqn{\theta}, at the plugged-in
#' correlation \eqn{\rho}. It is not conditional coverage and not uniform over
#' \eqn{B}; the feasible version plugs in consistent estimates of
#' \eqn{(\rho, \sigma_U, \sigma_O)}, justified by the local-asymptotic
#' framework in which those are consistently estimable while \eqn{b} is not.
#' Coverage degrades smoothly for \eqn{|b| > B}; setting \eqn{B = \infty}
#' recovers the usual interval centered at the unrestricted estimator, and no
#' interval centered at the adaptive estimate can be both short and uniformly
#' valid over all biases (Armstrong & Kolesar 2021, Section 4). Reporting the
#' \eqn{B = 0} and \eqn{B = \infty} intervals together -- the default here --
#' brackets the critical values needed under any bias bound, which is the AKS
#' recommendation. Coverage diagnostics in this implementation (and the
#' eq.-(8) solve under \code{st_cv = "exact"}) use deterministic Gaussian
#' quadrature over the representation (7) rather than MissAdapt's seeded
#' Monte Carlo; the two agree to the MC noise level (~4e-3).
#'
#' \strong{Lookup-table provenance.} The shipped tables
#' (\code{inst/extdata/aks_lookup/}) are the \code{policy.mat},
#' \code{thresholds.mat}, \code{emse_corr.mat}, \code{flci_adaptive_cv.mat},
#' \code{flci_adaptive_st_cv.mat}, and \code{flci_minimax_cv.mat} lookup
#' tables of the MissAdapt replication package of Armstrong, Kline & Sun
#' (Econometrica 2025), vendored byte-identically from
#' \url{https://github.com/lsun20/MissAdapt} (commit \code{98d823a}; also
#' archived as Zenodo record 16890198) and distributed under the MIT license
#' (Copyright (c) 2023 Sophie Sun; see \code{inst/COPYRIGHTS}). The
#' \code{aks_lookup.rds} conversion the function reads (so no MATLAB-file
#' reader is required at runtime) is built by \code{data-raw/aks_lookup.R},
#' which documents the grid conventions and runs orientation/monotonicity/
#' symmetry sanity checks, including a regression against the published
#' MissAdapt vignette example; the provenance, commit, license, and grid
#' conventions are embedded as attributes of the \code{.rds}. The FLCI
#' critical-value tables are 95\%-only and tabulated to two decimals on the
#' \eqn{\widetilde{B}} grid \{0.01, 0.1, ..., 9\} by the signed correlation
#' grid \code{tanh(seq(-3, -0.05, 0.05))}; the lookup splines each
#' \eqn{\widetilde{B}} row across the \eqn{|\rho|} grid and evaluates at the
#' clamped \eqn{|\widehat\rho|} (exactly equivalent, by spline mirror
#' symmetry, to the authors' signed-grid lookup for \eqn{\widehat\rho < 0},
#' and well-defined -- not an off-grid extrapolation -- for
#' \eqn{\widehat\rho > 0}, where the critical value is symmetric in
#' \eqn{\rho}). The \code{flci_minimax_cv.mat} table (critical values for the
#' B-minimax estimator) is vendored for completeness but not exposed:
#' \code{edid_adaptive} computes no B-minimax point estimate, so there is no
#' estimate for that interval to be centered at; the table is available
#' internally as \code{.edid_aks_lookup()$flci_cv_minimax} for a future
#' \code{B}-minimax estimator.
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
#'   \code{$assume_efficient_fallback} is \code{TRUE} when the AUTO rule
#'   selected \code{assume_efficient = TRUE} from the \code{"efficient"} weight
#'   scheme label but the restricted fit was \emph{not} empirically tighter
#'   than the unrestricted one (\eqn{\widehat\sigma_R^2 \ge
#'   \widehat\sigma_U^2}, so the imposed identity would give a non-positive
#'   over-identification variance), in which case it fell back to the empirical
#'   covariance (\code{assume_efficient = FALSE}) with a message rather than
#'   erroring -- the internally-consistent choice when the efficient leg is not
#'   tighter (an \emph{explicit} \code{assume_efficient = TRUE} still errors).
#'
#'   When \code{ci = TRUE} (the default), the object additionally carries:
#'   \code{$ci}, a data.frame with one row per (variant, B) -- for
#'   \code{parameter = "event_study"} also per \code{e} -- with columns
#'   \code{variant} (\code{"adaptive"}, the headline interval, or
#'   \code{"soft_threshold"}), \code{B} (the requested bound, \code{0} /
#'   \code{Inf} / user-supplied), \code{B_tilde} (the table row actually used:
#'   0.01 for \code{B = 0}, 9 for \code{B = Inf}), \code{center} (that
#'   variant's point estimate), \code{sigma_U} (\eqn{= \sqrt{VU}}, the scale
#'   of every interval), \code{cv} (the 95\% critical value
#'   \eqn{c_{.05}(\widetilde{B}; \hat\rho)}), \code{lower}/\code{upper}
#'   (\code{center} \eqn{\pm} \code{cv * sigma_U}), and \code{cv_source}
#'   (\code{"table"} for the shipped MissAdapt tables, \code{"exact"} for the
#'   runtime eq.-(8) solve -- the corrected-versus-shipped flag for the
#'   soft-threshold rows); plus \code{$sigma_U} (overall parameter only),
#'   \code{$ci_level} (\code{0.95}), \code{$st_cv} (the resolved soft-threshold
#'   cv source), and \code{$ci_note} (the one-line validity statement).
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025). Efficient
#'   Difference-in-Differences and Event Study Estimators. Section 5.2,
#'   Proposition 5.1. \cr
#'   Armstrong, T. B., Kline, P., & Sun, L. (2025). Adapting to
#'   Misspecification. \emph{Econometrica}, 93(6), 1981-2005. Replication
#'   package: MissAdapt, Zenodo 16890198. (B-FLCIs: Section 4.2; equation
#'   numbers (7)-(8) cited here follow the arXiv v6 manuscript.) \cr
#'   Armstrong, T. B., & Kolesar, M. (2021). Sensitivity analysis using
#'   approximate moment condition models. \emph{Quantitative Economics},
#'   12(1), 77-108. (Impossibility of short CIs valid over all biases.)
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
                          e_set = NULL, assume_efficient = NULL,
                          ci = TRUE, B = NULL, level = 0.95,
                          st_cv = c("exact", "missadapt")) {
  parameter <- match.arg(parameter)
  st_cv <- match.arg(st_cv)
  stopifnot(is.logical(ci), length(ci) == 1L, !is.na(ci))
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      abs(level - 0.95) > 1e-12) {
    stop("edid_adaptive: level must be 0.95. The MissAdapt B-FLCI critical-value ",
         "tables are tabulated for the 95% level only (no other level is ",
         "tabulated anywhere in their package).", call. = FALSE)
  }
  # 0- and Inf-FLCIs always reported (AKS recommendation); user B values added
  # in between (sort() places Inf last)
  B_report <- sort(unique(c(0, Inf, B)))
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

  n        <- fit_restricted$n
  clus_idx <- fit_restricted$cluster_indices
  tables   <- .edid_aks_lookup()
  if (ci) {
    flci_needed <- c("flci_B_grid", "flci_cv_adaptive", "flci_cv_st")
    if (!all(flci_needed %in% names(tables))) {
      stop("edid_adaptive: the installed aks_lookup.rds predates the B-FLCI layer ",
           "(missing components: ",
           paste(setdiff(flci_needed, names(tables)), collapse = ", "),
           "); rebuild it with data-raw/aks_lookup.R or reinstall the package.",
           call. = FALSE)
    }
  }

  # Bivariate variance of the (unrestricted, restricted) ESTIMATES:
  # cluster-robust sandwich of the two aggregation IFs, order 1/n (the same
  # scale as the reference implementation's mean(if^2)/n).
  .biv <- function(psi_U, psi_R) {
    S <- cluster_cov_edid(cbind(psi_U, psi_R), clus_idx, n)
    list(VU = S[1L, 1L], VR = S[2L, 2L], VUR = S[1L, 2L])
  }

  # AUTO convention fallback (Brazil-gate footgun): when AUTO selected
  # assume_efficient = TRUE from the scheme LABEL (weight_scheme = "efficient") but the
  # efficient (restricted) leg is NOT empirically tighter than the conservative one
  # (VR >= VU), the imposed identity VUR = VR forces sigma_O^2 = VU - VR <= 0 and the
  # core stops with "VO ... is not positive". On big-N covariate fits the efficient
  # weights can be empirically NOISIER than the conservative ones (the kernel collapses /
  # the over-identified moments add variance), so the LABEL is wrong about bound-attainment
  # there. Rather than erroring out, fall back to the empirical-covariance convention
  # (assume_efficient = FALSE) with a one-time message -- the internally-consistent choice
  # when the identity does not hold. An EXPLICIT assume_efficient = TRUE still errors (the
  # user asserted the identity; honor the documented contract). `$assume_efficient_fallback`
  # records that the fallback fired. Returns the core list (always assume_efficient = FALSE
  # on the fallback path).
  auto_fallback_fired <- FALSE
  .core_auto <- function(YR, VR, YU, VU, VUR) {
    ae <- assume_efficient
    if (assume_efficient_auto && isTRUE(assume_efficient)) {
      VO_try <- VR - 2 * VR + VU            # VO under the imposed identity VUR = VR (= VU - VR)
      if (!is.finite(VO_try) || VO_try <= 0) {
        if (!auto_fallback_fired) {
          message("edid_adaptive [AUTO]: weight_scheme = \"efficient\" labels the restricted fit as ",
                  "bound-attaining, but it is NOT empirically more precise than the conservative fit ",
                  "(VR >= VU), so the Hausman identity VUR = VR would give a non-positive ",
                  "over-identification variance. Falling back to assume_efficient = FALSE (empirical ",
                  "influence-function covariance) -- the internally-consistent convention when the ",
                  "efficient leg is not tighter. Pass assume_efficient = TRUE explicitly to override ",
                  "and error instead.")
          auto_fallback_fired <<- TRUE
        }
        ae <- FALSE
      }
    }
    .edid_aks_core(YR = YR, VR = VR, YU = YU, VU = VU, VUR = VUR,
                   tables = tables, assume_efficient = ae)
  }

  if (parameter == "overall") {
    oU <- .edid_param_ifs(fit_unrestricted, "overall")
    oR <- .edid_param_ifs(fit_restricted,  "overall")
    v  <- .biv(oU$IF[, 1L], oR$IF[, 1L])
    core <- .core_auto(YR = oR$est, VR = v$VR, YU = oU$est, VU = v$VU, VUR = v$VUR)
    out <- core
    out$psi_fun <- NULL   # internal interpolant; not part of the user-facing object
    out$acorr_eval <- NULL
    out$parameter <- "overall"
    if (ci) {
      out$ci <- .edid_aks_ci(core, tables, B_report, st_cv, level)
      out$sigma_U <- sqrt(core$VU)
    }
  } else {
    e_set <- .edid_shared_e_set(fit_unrestricted, fit_restricted, e_set)
    pU <- .edid_param_ifs(fit_unrestricted, "event_study", e_set)
    pR <- .edid_param_ifs(fit_restricted,  "event_study", e_set)
    cols <- c("YU", "YR", "VU", "VR", "VUR", "YO", "VO", "VUO", "tO", "corr", "rho_aks_sq",
              "GMM", "se_GMM", "adaptive", "adaptive_st", "adaptive_ht", "pretest",
              "adaptive_erm", "soft_threshold", "hard_threshold")
    rows <- vector("list", length(e_set))
    ci_rows <- if (ci) vector("list", length(e_set)) else NULL
    for (j in seq_along(e_set)) {
      v <- .biv(pU$IF[, j], pR$IF[, j])
      cj <- .core_auto(YR = pR$est[j], VR = v$VR, YU = pU$est[j], VU = v$VU, VUR = v$VUR)
      rows[[j]] <- cbind(data.frame(e = e_set[j]), as.data.frame(cj[cols]))
      if (ci) {
        ci_rows[[j]] <- cbind(data.frame(e = e_set[j]),
                              .edid_aks_ci(cj, tables, B_report, st_cv, level))
      }
    }
    out <- list(parameter = "event_study", e_set = e_set,
                table = do.call(rbind, rows))
    rownames(out$table) <- NULL
    if (ci) {
      out$ci <- do.call(rbind, ci_rows)
      rownames(out$ci) <- NULL
    }
  }

  if (ci) {
    out$ci_level <- level
    out$st_cv <- st_cv
    out$ci_note <- paste(
      "Each B-FLCI covers with probability >= 0.95 uniformly over biases",
      "|b| <= B*sigma_O at the plug-in rho, in the AKS normal limit experiment",
      "(Armstrong, Kline & Sun 2025, Sec. 4.2); B = Inf uses their B = 9*sigma_O",
      "approximation. All intervals are center +- cv*sigma_U with sigma_U = sqrt(VU).")
  }
  # Report the EFFECTIVE convention: the AUTO fallback (above) demotes a label-driven
  # assume_efficient = TRUE to FALSE when the efficient leg is not empirically tighter.
  out$assume_efficient <- if (auto_fallback_fired) FALSE else assume_efficient
  out$assume_efficient_auto <- assume_efficient_auto
  out$assume_efficient_fallback <- isTRUE(auto_fallback_fired)
  out$n <- n
  out$clustered <- !is.null(clus_idx)
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
    if (isTRUE(x$assume_efficient_fallback)) {
      cat("  (assume_efficient = FALSE [AUTO fallback: the \"efficient\"-labelled restricted\n")
      cat("   fit is not empirically tighter (VR >= VU), so the Hausman identity was dropped])\n")
    } else {
      cat("  (assume_efficient = FALSE [AUTO: restricted fit is not bound-attaining])\n")
    }
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
    if (!is.null(x$ci)) {
      ad <- x$ci[x$ci$variant == "adaptive", , drop = FALSE]
      cat(sprintf("\n  95%% adaptive FLCIs (estimate +- cv * sigma_U, sigma_U = sqrt(VU) = %s):\n",
                  format(x$sigma_U, digits = digits)))
      for (k in seq_len(nrow(ad))) {
        cat(sprintf("    B = %-4s (B~ = %-4s): cv = %s, CI = [%s, %s]\n",
                    format(ad$B[k]), format(ad$B_tilde[k]),
                    format(ad$cv[k], digits = digits),
                    format(ad$lower[k], digits = digits),
                    format(ad$upper[k], digits = digits)))
      }
    }
  } else {
    cat(sprintf("  Parameter: ES(e) per event time%s\n",
                if (isTRUE(x$clustered)) " (cluster-robust)" else ""))
    tab <- x$table[, c("e", "YU", "YR", "tO", "rho_aks_sq", "adaptive", "GMM")]
    num <- vapply(tab, is.numeric, logical(1L))
    tab[num] <- lapply(tab[num], function(z) signif(z, digits))
    print(tab, row.names = FALSE)
    if (!is.null(x$ci)) {
      cat("\n  95% adaptive FLCIs (estimate +- cv * sigma_U, sigma_U = sqrt(VU)):\n")
      ad <- x$ci[x$ci$variant == "adaptive",
                 c("e", "B", "B_tilde", "center", "cv", "lower", "upper"), drop = FALSE]
      num <- vapply(ad, is.numeric, logical(1L))
      ad[num] <- lapply(ad[num], function(z) signif(z, digits))
      print(ad, row.names = FALSE)
    }
  }
  if (!is.null(x$ci)) {
    cat("\nEach FLCI covers with prob >= 0.95 uniformly over PT violations |b| <= B*sigma_O\n")
    cat("(plug-in rho; AKS 2025, Sec. 4.2; B = Inf via their B = 9*sigma_O approximation).\n")
    cat("The soft-threshold FLCI (cv source: ", x$st_cv,
        ") is in $ci; see ?edid_adaptive.\n", sep = "")
  } else {
    cat("\nNote: no conventional standard error attaches to the adaptive estimate (its\n")
    cat("asymptotic distribution is not consistently estimable); call with ci = TRUE\n")
    cat("for the AKS fixed-length confidence intervals (see ?edid_adaptive).\n")
  }
  invisible(x)
}
