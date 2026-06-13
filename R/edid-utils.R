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

# Refit helpers need to preserve no-covariate shrinkage targets that are
# configured through global options rather than formal edid() arguments. Legacy
# fits without these fields fall through to the caller's current options, as
# before.
.edid_with_nocov_shrink_options <- function(fit, expr) {
  target <- fit$nocov_shrink_target
  if (is.null(target) || length(target) != 1L || is.na(target)) {
    return(eval.parent(substitute(expr)))
  }
  target <- as.character(target)
  rho <- fit$nocov_shrink_rho
  if (is.null(rho) || length(rho) != 1L || !is.finite(rho)) rho <- 0
  old <- options(edid_nocov_shrink_target = target,
                 edid_nocov_ar1_rho = as.numeric(rho))
  on.exit(options(old), add = TRUE)
  eval.parent(substitute(expr))
}
