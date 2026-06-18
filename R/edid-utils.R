# edid-utils.R
# Internal constants and small shared helpers for the EDiD estimator.
# Do NOT modify this file to change the estimator logic -- see the relevant
# edid-*.R file for each component.

#' @keywords internal
EDID_COND_THRESH <- 1e12   # condition number above which pseudoinverse is used

#' @keywords internal
EDID_DENOM_EPS <- 1e-12    # denominator threshold below which uniform weights are used

#' @keywords internal
EDID_CLIP_LO <- 1 / 20     # ratio clipping lower bound (deferred: covariate path)

#' @keywords internal
EDID_CLIP_HI <- 20         # ratio clipping upper bound (deferred: covariate path)

#' @keywords internal
EDID_SE_EPS <- sqrt(.Machine$double.eps) * 10  # SE below which is treated as zero/NA

#' @keywords internal
# Comfort threshold for the thin-cohort radar note (informational only; NO behavior
# change). Finite treated cohorts at or above min_pair_units but below this size are
# flagged as small enough that the COVARIATE-path over-identified efficient weights may
# be unreliable (the audited Nguyen failure: 14- and 33-unit cohorts at 0.6-1.4% share,
# fatal with d=4 X, while the no-covariate path is fine). 36 sits just above the
# Nguyen 33-unit failing cohort and below the >= 75-unit cohorts the Dias-Fontes gate
# ran cleanly. The no-covariate path never trips on it (radar is PT-All only and the
# message names the covariate path explicitly).
EDID_THIN_COHORT_COMFORT <- 36L

#' @keywords internal
# Minimum number of clusters below which the Section-5 toolkit's cluster-robust statistics
# are flagged as unreliable (message + field; the statistic is still returned). With G < 5
# the cluster-level moment covariance is too noisy / degenerate for the chi-square (or
# G/(G-1)) reference to hold -- the ACA gate's 2-3-state cohorts produce Sargan "rejections"
# (H ~ 350, p ~ 1e-73) that are few-cluster artifacts, not parallel-trends evidence.
EDID_FEWCLUSTER_MIN <- 5L

#' @keywords internal
# Net cross-cohort hedge-mass red-flag threshold for the fit diagnostics (informational
# only). Calibrated from the gate evidence: healthy efficient fits carry net cross-cohort
# hedge mass ~0.01-0.43 (the "hedges" hedge -- gross negative mass offsets gross
# positive), while the audited broken with-X fits show net ~= gross >= 0.6 with zero
# negative mass (the cross-cohort control variates stop hedging and CARRY the estimand).
# 0.55 sits above the healthy band and below every broken sighting (0.628, 0.681, 0.857,
# 0.878).
EDID_NET_HEDGE_FLAG <- 0.55

#' @keywords internal
# Estimability auto-guard (OPT-IN; edid_auto_excise_unstable_pairs) post-trim thresholds.
# A cross-cohort comparison cohort is excised when, on the units surviving overlap trimming,
# its fitted propensity ratio still exceeds EDID_RATIO_EXCISE_THRESH (the |r| > 100 scale the
# extreme-ratio diagnostic uses -- a ratio this large AFTER trimming is a poisoned,
# unestimable cross moment), or when fewer than EDID_RATIO_EXCISE_MINKEEP units survive for
# that comparison (the trim removed essentially all its mass).
EDID_RATIO_EXCISE_THRESH <- 100
#' @keywords internal
EDID_RATIO_EXCISE_MINKEEP <- 5L

#' Is fork-based parallelism unsafe on this platform's BLAS?
#'
#' macOS Apple Accelerate (vecLib) BLAS is not safe to call from a process forked by
#' \code{parallel::mclapply}: a forked worker that enters Accelerate (the covariate-path
#' cell loop's \code{crossprod} / kernel solves) can crash, which \code{mclapply}
#' surfaces only as a missing/NULL result -- corrupting or aborting the fit with no R
#' error. Returns \code{TRUE} on Darwin when the linked BLAS reports as an
#' Accelerate/vecLib library, \code{FALSE} otherwise (Linux/Windows, or macOS linked
#' against a fork-safe BLAS such as OpenBLAS). Used by \code{\link{edid}} to default
#' \code{cores > 1} back to serial on the unsafe configuration (override:
#' \code{options(edid_allow_fork_blas = TRUE)}). Cheap and dependency-free
#' (\code{extSoftVersion()} string match); not exported.
#' @keywords internal
.edid_fork_blas_unsafe <- function() {
  if (!identical(Sys.info()[["sysname"]], "Darwin")) return(FALSE)
  blas <- tryCatch(tolower(extSoftVersion()[["BLAS"]] %||% ""),
                   error = function(e) "")
  # Accelerate.framework / vecLib.framework / libBLAS.dylib are the Apple BLAS markers;
  # a fork-safe replacement (openblas, libRblas, mkl, ...) does not match.
  grepl("accelerate|veclib", blas, fixed = FALSE) ||
    (grepl("libblas\\.dylib$", blas) && grepl("/system/library/frameworks/", blas))
}

#' @keywords internal
# Variance-inflation ceiling for the misspec_robust weight-estimation channel: the largest factor by which
# folding psi_Omega may inflate a cell's EIF variance (=> SE inflation <= sqrt of this). A genuine first-order
# correction inflates the cell SE by at most ~2-3x; this 100 (SE <= 10x) is far above that yet far below the
# catastrophic blowups (SE ~1e14) the sieve psi can produce in poor-overlap / placebo cells, where huge
# inverse-propensity prefactors meet a near-singular series basis Gram. Beyond it the channel is not a credible
# influence function (it is no longer mean-zero), so that cell falls back to the plug-in efficient SE.
EDID_PSI_VAR_RATIO <- 100

#' Is the weight-estimation channel a credible influence function for this cell?
#'
#' A valid \eqn{\psi_\Omega} is finite and mean-zero, and -- being a first-order, root-n-vanishing correction --
#' inflates the cell EIF variance only modestly. In poor-overlap / placebo cells the sieve OLS-projection IF can
#' instead explode (the eigen-floor bounds the coupling but not its product with large inverse-propensity
#' prefactors and a near-singular basis Gram), giving a non-mean-zero \eqn{\psi} and an absurd SE. This gate
#' rejects such \eqn{\psi} so the caller can fall back to the (finite) plug-in efficient SE for that cell --
#' the same per-cell skip convention used elsewhere on this path. Tested against \code{eif_base}, the
#' ACH-corrected, mean-zero baseline EIF the channel would be folded into.
#'
#' The variance ratio is computed on the SAME quantity the reported SE uses: when \code{cluster_indices} is
#' supplied (the cluster-robust SE sums the EIF within clusters first), the inflation is measured on the
#' cluster-summed EIF, so the SE-inflation ceiling (the square root of \code{EDID_PSI_VAR_RATIO}) holds for the
#' clustered SE too (a \eqn{\psi} that is modest per unit but strongly within-cluster correlated -- which
#' inflates the clustered SE far more than the i.i.d. SE -- is then judged on the metric that governs the
#' reported number). Without clustering it is the plain i.i.d. sum of squares.
#'
#' @param psi numeric length-n weight-estimation influence function for the cell
#' @param eif_base numeric length-n baseline EIF (mean-zero) that \code{psi} would be added to
#' @param cluster_indices optional length-n cluster id vector; when non-NULL the ratio is computed on the
#'   cluster-summed EIF, matching the cluster-robust SE. \code{NULL} (default) => i.i.d. sum of squares.
#' @return \code{TRUE} if \code{psi} is finite and its (clustered, if applicable) variance inflation is within
#'   \code{EDID_PSI_VAR_RATIO}
#' @keywords internal
psi_channel_credible_edid <- function(psi, eif_base, cluster_indices = NULL) {
  if (is.null(psi) || !all(is.finite(psi))) return(FALSE)
  if (is.null(cluster_indices)) {                                  # i.i.d.: per-unit sum of squares
    base <- eif_base; fold <- eif_base + psi
  } else {                                                         # cluster-robust: sum within clusters first
    base <- rowsum(eif_base, cluster_indices); fold <- rowsum(eif_base + psi, cluster_indices)
  }
  v0 <- sum(base^2); v1 <- sum(fold^2)
  is.finite(v1) && v1 <= EDID_PSI_VAR_RATIO * max(v0, .Machine$double.eps)
}

# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------

#' Biased sample covariance (divide by n, not n-1)
#'
#' @param x numeric vector
#' @param y numeric vector, same length as x
#' @return scalar
#' @keywords internal
cov_nn_edid <- function(x, y) {
  mean((x - mean(x)) * (y - mean(y)))
}

# ---------------------------------------------------------------------------
# Weighted (Hajek) versions of the group-mean / group-covariance primitives.
# These power the observation-weights (`weightsname`) path. When `w` is NULL --
# the unweighted default -- each dispatches to the EXACT unweighted expression
# above, so the no-weights path is BYTE-IDENTICAL (same floating-point ops).
# Weights need NOT sum to anything in particular here: the Hajek mean and the
# group-share-normalized covariance are scale-invariant in `w`.
# ---------------------------------------------------------------------------

#' Weighted (Hajek) mean
#'
#' \code{wmean_edid(x, NULL)} is \code{mean(x)} bit-for-bit; otherwise
#' \eqn{\sum_i w_i x_i / \sum_i w_i}.
#' @param x numeric vector
#' @param w numeric vector of nonnegative weights (same length), or NULL
#' @return scalar
#' @keywords internal
wmean_edid <- function(x, w = NULL) {
  if (is.null(w)) return(mean(x))
  sw <- sum(w)
  if (sw == 0) return(NA_real_)
  sum(w * x) / sw
}

#' Weighted biased covariance, group-share normalized
#'
#' \code{wcov_nn_edid(x, y, NULL)} is \code{cov_nn_edid(x, y)} bit-for-bit.
#' With weights it returns the Hajek-weighted second moment of the demeaned
#' vectors, \eqn{\sum_i w_i (x_i-\bar x_w)(y_i-\bar y_w)/\sum_i w_i}. This is
#' the weighted analog of the biased (divide-by-n) covariance: it is the
#' covariance under the reweighted empirical measure with masses
#' \eqn{w_i/\sum w}.
#' @param x,y numeric vectors of equal length
#' @param w numeric vector of nonnegative weights (same length), or NULL
#' @return scalar
#' @keywords internal
wcov_nn_edid <- function(x, y, w = NULL) {
  if (is.null(w)) return(mean((x - mean(x)) * (y - mean(y))))
  sw <- sum(w)
  if (sw == 0) return(NA_real_)
  mx <- sum(w * x) / sw
  my <- sum(w * y) / sw
  sum(w * (x - mx) * (y - my)) / sw
}

#' Sampling-variance term of a group mean (\eqn{\mathrm{Cov}(\bar x_g, \bar y_g)})
#'
#' Returns the (cross-)sampling-variance contribution of two group means built
#' on the SAME group of units, which the no-covariate \eqn{\Omega^*} builder adds
#' as \code{cov_nn_edid(x, y) / n_g} in the unweighted case. The weighted (Hajek)
#' generalization is the design-based variance of the ratio (Hajek) mean,
#' \deqn{\sum_{i\in g} w_i^2 (x_i - \bar x_w)(y_i - \bar y_w) / W_g^2,\quad
#'   W_g = \sum_{i\in g} w_i,}
#' which the EIF identity \eqn{\Omega^* = \Psi'\Psi/n^2} requires (the per-unit
#' moment influence carries an explicit \eqn{w_i} factor and a \eqn{1/\pi_g}, with
#' \eqn{\pi_g = W_g/n}). With \code{w = NULL} (or all-equal weights after the
#' mean-1 normalization) it is \code{cov_nn_edid(x, y) / n_g} bit-for-bit:
#' \eqn{n_g \,\mathrm{cov}_{nn}/n_g^2 = \mathrm{cov}_{nn}/n_g}.
#'
#' @param x,y numeric vectors of equal length (the group's difference vectors)
#' @param w numeric vector of the group's nonnegative weights, or NULL
#' @return scalar
#' @keywords internal
wvar_term_edid <- function(x, y, w = NULL) {
  if (is.null(w)) return(mean((x - mean(x)) * (y - mean(y))) / length(x))
  sw <- sum(w)
  if (sw == 0) return(NA_real_)
  mx <- sum(w * x) / sw
  my <- sum(w * y) / sw
  sum(w * w * (x - mx) * (y - my)) / (sw * sw)
}

#' Effective sample size (Kish ESS) for the ridge / Ledoit-Wolf intensity
#'
#' The vanishing ridge / Ledoit-Wolf intensities that stabilize the WEIGHTS in
#' \code{\link{edid}} scale as the reciprocal of the number of independent
#' contributions backing the estimated moment covariance \eqn{\widehat\Omega^*}.
#' Unweighted that count is the active-unit count; under dispersed observation
#' weights (\code{weightsname}) the heavily-weighted units dominate
#' \eqn{\widehat\Omega^*}, so the right count is the Kish effective sample size
#' \deqn{n_{\mathrm{eff}} = \frac{(\sum_i w_i)^2}{\sum_i w_i^2} \le m,}
#' the design-based ESS of the units active in that cell's weighted
#' \eqn{\widehat\Omega^*}. Using the raw count \eqn{n} instead under-regularizes
#' the scale-invariant weighted \eqn{\widehat\Omega^*} by the factor
#' \eqn{n / n_{\mathrm{eff}} \ge 1}.
#'
#' \strong{Byte-identity (non-negotiable).} With \code{w = NULL} (the unweighted
#' default) this returns \code{n_full} (\code{panel_obj$n}) UNCHANGED -- exactly
#' the count every legacy ridge/LW intensity used -- so all unweighted intensities
#' are bit-for-bit identical regardless of which units are active in the cell.
#' Only the WEIGHTED branch uses the active-unit Kish ESS, matching the mandate
#' \code{n_eff = if (no weights) panel_obj$n else (sum w)^2/sum(w^2)} over the
#' cell's active units. After the mean-1 normalization in
#' \code{prepare_edid_panel()} a CONSTANT weight column is all-ones on its active
#' units, so \eqn{n_{\mathrm{eff}} = (\sum 1)^2/\sum 1 = n_{\mathrm{act}}}; the
#' headline weighted-application designs activate all units
#' (\eqn{n_{\mathrm{act}} = n}), so a constant column is a no-op there too.
#'
#' \strong{Asymptotics.} For a fixed weight distribution \eqn{n_{\mathrm{eff}}}
#' grows proportionally to the active count, so any intensity of the form
#' \eqn{c / n_{\mathrm{eff}}} still vanishes as the sample grows -- the
#' semiparametric-efficiency limit is preserved.
#'
#' @param w numeric vector of nonnegative unit weights (\code{panel_obj$unit_weights}),
#'   or \code{NULL} on the unweighted path
#' @param active_mask logical vector over ALL units (length \code{panel_obj$n})
#'   selecting the units that enter this cell's weighted \eqn{\widehat\Omega^*}
#'   (the nonzero rows of \eqn{\Psi}). Used only on the weighted branch.
#' @param n_full scalar full-sample size (\code{panel_obj$n}); the value returned
#'   verbatim on the unweighted path (\code{w = NULL}). Defaults to
#'   \code{length(active_mask)}.
#' @return scalar effective sample size (\code{>= 1}); \code{n_full} exactly when
#'   \code{w} is \code{NULL}
#' @keywords internal
n_eff_edid <- function(w, active_mask, n_full = length(active_mask)) {
  if (is.null(w)) return(n_full)                            # unweighted: legacy full n, byte-identical
  ww <- w[active_mask]
  sw <- sum(ww)
  if (!is.finite(sw) || sw <= 0) return(n_full)            # degenerate: fall back to full count
  ne <- (sw * sw) / sum(ww * ww)
  if (!is.finite(ne) || ne < 1) 1 else ne
}

#' Active-unit mask for a no-covariate (g, t) cell's weighted Omega*
#'
#' The nonzero rows of \eqn{\Psi} (\code{compute_psi_moments_nocov_edid()}) -- the
#' units that actually enter the cell's \eqn{\widehat\Omega^*} -- are the treated
#' cohort \eqn{g}, the never-treated group, and every comparison cohort
#' \eqn{g'_j} appearing in \code{pairs$gp}. Returns their union as a logical mask
#' over all units, for \code{\link{n_eff_edid}}.
#'
#' @param target_g scalar cohort value
#' @param pairs data.frame with column \code{gp} (the comparison cohorts), H rows
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @return logical vector length \code{panel_obj$n}
#' @keywords internal
active_mask_nocov_edid <- function(target_g, pairs, panel_obj) {
  m <- panel_obj$cohort_masks[[as.character(target_g)]] | panel_obj$never_treated_mask
  for (gp in unique(pairs$gp)) {
    cm <- panel_obj$cohort_masks[[as.character(gp)]]
    if (!is.null(cm)) m <- m | cm
  }
  m
}

#' Safe mean: returns NA on empty vector instead of NaN
#'
#' @param x numeric vector
#' @return scalar
#' @keywords internal
safe_mean_edid <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

# ---------------------------------------------------------------------------
# Package imports (merged from edid-imports.R). stats::pnorm/qnorm/quantile/sd/
# setNames are declared in imports.R; only the additional symbols are added here.
# ---------------------------------------------------------------------------
#' @importFrom stats sd quantile as.formula
NULL

# ---------------------------------------------------------------------------
# Linear algebra helpers (merged from edid-linalg.R). Base-R svd(), no MASS dep.
# ---------------------------------------------------------------------------

#' SVD-based Moore-Penrose pseudoinverse
#'
#' @param mat numeric matrix
#' @param tol tolerance for zero singular values; defaults to
#'   \code{max(dim(mat)) * max(svd$d) * .Machine$double.eps}
#' @return matrix of same dimensions as \code{t(mat)}
#' @keywords internal
compute_pseudoinverse_edid <- function(mat, tol = NULL) {
  s <- svd(mat)
  d <- s$d
  if (is.null(tol)) {
    tol <- max(dim(mat)) * max(c(d, 0)) * .Machine$double.eps
  }
  # Zero out singular values below tolerance
  d_inv <- ifelse(d > tol, 1 / d, 0)
  # Pseudoinverse: V diag(d_inv) U'
  s$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(s$u)
}

#' Condition number of a matrix via SVD
#'
#' Singular values at or below \code{tol * max(d)} (\code{tol = 100 * .Machine$double.eps}) are treated as
#' structural zeros: an exactly (or numerically) singular matrix returns \code{Inf}, not the large-but-finite
#' ratio \code{max(d) / min(d[d > 0])} of its FP-noise smallest singular value -- which would let a caller
#' compare a rank-deficient matrix against a finite condition threshold and wrongly take the \code{solve()} path.
#'
#' @param mat numeric matrix
#' @return scalar: max singular value / min singular value above the relative tolerance.
#'   Returns \code{Inf} if any singular value is a structural zero (or the matrix is all zero).
#' @keywords internal
check_condition_edid <- function(mat) {
  d <- svd(mat, nu = 0L, nv = 0L)$d
  if (length(d) == 0L || max(d) == 0) return(Inf)
  tol <- 100 * .Machine$double.eps * max(d)
  if (any(d <= tol)) return(Inf)
  max(d) / min(d)
}

#' Weighted OLS helper
#'
#' Computes \eqn{\hat\beta = (X'WX)^{-1} X'Wy} using \code{.lm.fit()}.
#' Falls back to SVD-based pseudoinverse if the normal equations are
#' numerically singular.
#'
#' @param X numeric matrix (n x p)
#' @param y numeric vector length n
#' @param weights numeric vector length n (NULL = uniform)
#' @return named list with elements \code{coef}, \code{fitted}, \code{residuals}
#' @keywords internal
solve_ols_edid <- function(X, y, weights = NULL) {
  n <- nrow(X)
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  W <- sqrt(weights)
  Xw <- X * W
  yw <- y * W
  fit <- tryCatch(
    stats::.lm.fit(Xw, yw),
    error = function(e) NULL
  )
  if (!is.null(fit) && all(is.finite(fit$coefficients))) {
    beta  <- fit$coefficients
    yhat  <- drop(X %*% beta)
    resid <- y - yhat
    return(list(coef = beta, fitted = yhat, residuals = resid))
  }
  # Fallback: pseudoinverse
  XtWX <- t(Xw) %*% Xw
  XtWy <- t(Xw) %*% yw
  beta  <- drop(compute_pseudoinverse_edid(XtWX) %*% XtWy)
  yhat  <- drop(X %*% beta)
  resid <- y - yhat
  list(coef = beta, fitted = yhat, residuals = resid)
}

# Evaluated edid() arguments for an internal refit of a fitted model: the
# argument list (everything except `data`) with which `fit` was estimated, for
# the refit tools (edid_sargan()'s moment-set refits, edid_refit_bootstrap() /
# edid_perturbation_bootstrap()). Fits store this snapshot at fit time
# ($args), so refits reuse the materials captured when the fit was made. For
# fits created before $args existed (e.g. loaded from disk), falls back to the
# legacy idiom: normalize the stored call to named arguments and re-evaluate
# them in `envir` -- which breaks if a variable referenced in the call has
# changed, is no longer reachable, or is a `..N` promise from a
# programmatically built call.
.edid_refit_args <- function(fit, envir = parent.frame()) {
  if (!is.null(fit$args)) return(fit$args)
  mc <- match.call(definition = edid, call = fit$call)
  args <- as.list(mc)[-1L]
  args$data <- NULL
  lapply(args, function(a) eval(a, envir = envir))
}

# Recover the estimation data for a refit: the supplied `data`, else the data expression stored in the
# fit's call, re-evaluated in `envir` (the update() idiom; same recovery edid_sargan() uses).
.edid_recover_data <- function(fit, data = NULL, envir = parent.frame()) {
  if (!is.null(data)) return(as.data.frame(data))
  d <- tryCatch(as.data.frame(eval(fit$call$data, envir = envir)), error = function(e) NULL)
  if (is.null(d))
    stop("Could not recover the estimation data from the fit's call; pass `data` explicitly.", call. = FALSE)
  d
}

# Memoization cache for plug-in refits (.edid_plugin_refit). The over-identification toolkit calls the
# plug-in refit on the SAME legs many times in one operation (edid_hausman event_study + overall,
# edid_sargan, and once per window-grow step), and the plug-in refit is a PURE FUNCTION of (fit, data) --
# corrections are forced off, options come from fit$args, data is fixed -- so caching it is bit-identical to
# recomputing. Keyed by a stable fingerprint of the fit (att-vector moments) + the refit args + nrow(data).
# Bounded (clear-all on overflow) to cap memory. Disable via options(edid_plugin_cache = FALSE).
.edid_plugin_cache <- new.env(parent = emptyenv())

.edid_plugin_key <- function(fit, args) {
  a  <- fit$att_gt$att
  af <- if (is.null(a) || !length(a)) "na" else
        sprintf("%d:%.12g:%.12g:%.12g", length(a),
                sum(a, na.rm = TRUE), sum(a^2, na.rm = TRUE), max(abs(a), na.rm = TRUE))
  xf <- if (is.null(args$xformla)) "" else paste(deparse(args$xformla), collapse = "")
  nd <- if (is.data.frame(args$data)) nrow(args$data) else 0L
  paste(af, args$weight_scheme, args$ratio_method, args$pt_assumption,
        as.character(if (is.null(args$omega_cov_shrink)) "" else args$omega_cov_shrink), xf,
        paste(if (is.null(args$weightsname)) "" else args$weightsname, collapse = ","),
        paste(if (is.null(args$clustervars)) "" else args$clustervars, collapse = ","),
        args$idname, args$tname, args$gname, args$yname,
        if (is.null(args$min_pair_units)) "" else args$min_pair_units,
        if (is.null(args$anticipation)) "" else args$anticipation, nd, sep = "|")
}

#' Clear the plug-in-refit memoization cache
#'
#' The over-identification toolkit memoizes its internal plug-in refits within a session (see
#' \code{options(edid_plugin_cache)}). This clears that cache; rarely needed (entries are keyed by a fit
#' fingerprint, so they never collide), but available for long-running sessions or benchmarking.
#' @return Invisibly \code{NULL}.
#' @keywords internal
#' @export
edid_clear_plugin_cache <- function() {
  rm(list = ls(.edid_plugin_cache, all.names = TRUE), envir = .edid_plugin_cache)
  invisible(NULL)
}

# Refit a fit in the BARE PLUG-IN configuration -- all three estimation-effect channels off
# (misspec_robust / estimation_effect / higher_order) -- so its influence functions are the EFFICIENT
# plug-in EIF. This is the influence function the over-identification toolkit (edid_hausman / edid_sargan /
# edid_frontier) must use: the over-identification (J / Hausman) object lives on the efficient
# inverse-variance variance, NOT the misspecification-robust variance (Andrews, Chen & Tecchio 2025,
# Sec 5 / Prop 5.2; the misspecification-robust SE is for INFERENCE ON THE ESTIMAND, their Sec 4 -- a
# separate object). Point estimates are unchanged (the channels are variance-only), so only the IFs revert
# to the efficient ones, and the toolkit is invariant to how the leg was originally fit. Estimation options
# (pt_assumption, moment_set, weight_scheme, smoother, shrinkage, clustering, weights, ...) come from the
# fit's stored snapshot ($args); `data` is recovered from the call when NULL. A plug-in refit is also
# cheaper than the original misspec_robust fit (it skips the expensive psi_Omega weight-estimation channel).
.edid_plugin_refit <- function(fit, data = NULL, envir = parent.frame()) {
  args <- .edid_refit_args(fit, envir)
  args$data              <- .edid_recover_data(fit, data, envir)
  args$misspec_robust    <- FALSE
  args$estimation_effect <- FALSE
  args$higher_order      <- FALSE
  args[["cband_method"]] <- NULL     # analytic default applies (no bootstrap below)
  args$cband             <- FALSE
  args$bstrap            <- FALSE
  args$aggregate         <- "none"   # the toolkit aggregates on demand via .edid_param_ifs
  ## MEMOIZE: the plug-in refit is a pure function of (fit, data); the toolkit calls this on the SAME legs
  ## repeatedly. Return the cached refit on a hit (bit-identical -- the do.call AND the att-reproduction
  ## guard already ran on the miss). options(edid_plugin_cache = FALSE) forces the always-refit path.
  use_cache <- isTRUE(getOption("edid_plugin_cache", TRUE))
  key <- if (use_cache) tryCatch(.edid_plugin_key(fit, args), error = function(e) NULL) else NULL
  if (!is.null(key) && exists(key, envir = .edid_plugin_cache, inherits = FALSE))
    return(get(key, envir = .edid_plugin_cache, inherits = FALSE))
  refit <- do.call(edid, args)
  # Safety net against a data mismatch. The dropped channels are variance-only, so the plug-in refit MUST
  # reproduce the fit's point estimates; if it does not, the supplied/recovered `data` is not the data this
  # fit was made from (e.g. a data-symbol collision when the fit was built inside a function, where the
  # call's data expression resolves to a different object in the caller's frame). Error LOUDLY rather than
  # return a silently-wrong over-identification statistic; passing `data` explicitly resolves it.
  a0 <- if (!is.null(fit$att_gt)) fit$att_gt$att else NULL
  a1 <- if (!is.null(refit$att_gt)) refit$att_gt$att else NULL
  bad <- is.null(a0) || is.null(a1) || length(a0) != length(a1) ||
    { ok <- is.finite(a0) & is.finite(a1)
      sum(ok) == 0L || max(abs(a1[ok] - a0[ok])) > 1e-6 * (1 + max(abs(a0[ok]))) }
  if (isTRUE(bad))
    stop("edid over-identification toolkit: the supplied/recovered `data` does not reproduce this fit ",
         "(the plug-in refit's point estimates differ). The data could not be recovered unambiguously ",
         "from the fit's call -- e.g. the fit was built inside a function or another object shares the ",
         "data variable's name. Pass `data` explicitly.", call. = FALSE)
  if (!is.null(key)) {
    if (length(ls(.edid_plugin_cache, all.names = TRUE)) >= 16L)
      rm(list = ls(.edid_plugin_cache, all.names = TRUE), envir = .edid_plugin_cache)   # bound memory
    assign(key, refit, envir = .edid_plugin_cache)
  }
  refit
}
