# edid-linalg.R
# Linear algebra helpers for the EDiD estimator.
# Uses base R svd() for pseudoinverse -- no MASS dependency.

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
    .lm.fit(Xw, yw),
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
