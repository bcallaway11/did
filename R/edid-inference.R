# edid-inference.R
# Analytical standard error and inference helpers for the EDiD estimator.

#' Safely compute SE, CI, and p-value from an EIF vector
#'
#' Dispatches to \code{compute_eif_se_edid()} with optional cluster aggregation.
#' If the resulting SE is not valid (zero, NA, or non-finite), all inference
#' results are set to \code{NA}.
#'
#' @param eif numeric vector length n (or NULL, for NA cells)
#' @param cluster_indices integer vector length n (1..G) or NULL
#' @param alpha significance level in (0, 1)
#' @param att scalar ATT estimate (used for t-stat; may be NA for inference check)
#'
#' @return named list:
#'   \code{se}, \code{ci_lower}, \code{ci_upper}, \code{t_stat},
#'   \code{p_value}, \code{inference_valid}
#' @keywords internal
safe_inference_edid <- function(eif, cluster_indices = NULL, alpha = 0.05,
                                att = NA_real_) {
  na_result <- list(
    se             = NA_real_,
    ci_lower       = NA_real_,
    ci_upper       = NA_real_,
    t_stat         = NA_real_,
    p_value        = NA_real_,
    inference_valid = FALSE
  )

  if (is.null(eif) || !is.numeric(eif) || length(eif) == 0L) {
    return(na_result)
  }

  if (is.null(cluster_indices)) {
    n  <- length(eif)
    se <- compute_eif_se_edid(eif, n)
  } else {
    G            <- length(unique(cluster_indices))
    if (G <= 1L) return(na_result)                       # cluster-robust SE undefined with < 2 clusters
    cluster_sums <- drop(rowsum(eif, cluster_indices))
    se           <- sqrt((G / (G - 1)) * sum(cluster_sums^2) / (length(eif)^2))
  }

  valid_se <- is.finite(se) && se > EDID_SE_EPS

  if (!valid_se) return(na_result)

  # att must be finite for CIs and p-value to be meaningful.
  # If att is NA/non-finite, return the SE but mark inference as invalid
  # (inference_valid = FALSE signals that CIs and p-value are not available).
  if (!is.finite(att)) {
    return(list(
      se              = se,
      ci_lower        = NA_real_,
      ci_upper        = NA_real_,
      t_stat          = NA_real_,
      p_value         = NA_real_,
      inference_valid = FALSE
    ))
  }

  z_crit  <- qnorm(1 - alpha / 2)
  t_stat  <- att / se
  p_value <- 2 * pnorm(-abs(t_stat))

  list(
    se              = se,
    ci_lower        = att - z_crit * se,
    ci_upper        = att + z_crit * se,
    t_stat          = t_stat,
    p_value         = p_value,
    inference_valid = TRUE
  )
}

#' Compute SE from EIF vector
#'
#' \deqn{SE = \sqrt{\sum_i \text{eif}_i^2 / n^2}}
#'
#' @param eif_vec numeric vector (may be cluster-aggregated sums)
#' @param n integer denominator (number of units or clusters)
#'
#' @return scalar SE
#' @keywords internal
compute_eif_se_edid <- function(eif_vec, n) {
  sqrt(sum(eif_vec^2) / (n^2))
}

#' Aggregate EIF to cluster level (centered)
#'
#' Returns the vector of cluster sums of \code{eif}, mean-subtracted.
#' The small-sample correction \eqn{G/(G-1)} is applied in the SE formula
#' (in \code{safe_inference_edid}), not here.
#'
#' @param eif numeric vector length n
#' @param cluster_indices integer vector length n (values 1..G)
#'
#' @return numeric vector length G (cluster sums, centered)
#' @keywords internal
cluster_aggregate_edid <- function(eif, cluster_indices) {
  cluster_sums <- drop(rowsum(eif, cluster_indices))  # length G
  cluster_sums - mean(cluster_sums)
}
