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
#' @param mat numeric matrix
#' @return scalar: max singular value / min positive singular value.
#'   Returns \code{Inf} if min singular value is 0.
#' @keywords internal
check_condition_edid <- function(mat) {
  d <- svd(mat, nu = 0L, nv = 0L)$d
  if (length(d) == 0L || max(d) == 0) return(Inf)
  min_pos <- min(d[d > 0])
  if (length(min_pos) == 0L) return(Inf)
  max(d) / min_pos
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
