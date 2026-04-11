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
