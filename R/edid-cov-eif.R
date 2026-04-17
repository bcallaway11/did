# edid-cov-eif.R
# Generated outcome and EIF computation for the EDiD covariate path.
# Implements Chen, Sant'Anna & Xie (2025) Eq. (3.7)-(3.12).

# ---------------------------------------------------------------------------
# Generated outcomes (doubly-robust, n x H matrix)
# ---------------------------------------------------------------------------

#' Compute doubly-robust generated outcomes for a (g, t) cell
#'
#' Returns the n x H matrix \eqn{\Phi} where column j corresponds to pair
#' \eqn{j = (g'_j, t_{pre,j})} and row i to unit i.  The formula implements
#' Eq. (3.9) of Chen, Sant'Anna & Xie (2025).
#'
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param g scalar: treatment cohort
#' @param t scalar: target time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param prop_ratios named list of n-vectors keyed by \code{as.character(gp)}:
#'   the cross-fitted \eqn{\hat r_{g,g'}(X_i)}
#' @param cond_means named list of n-vectors keyed by
#'   \code{paste0(gp, "_", period)}: the cross-fitted conditional means
#' @param pt_assumption \code{"all"} or \code{"post"}
#'
#' @return numeric matrix n x H; entries may be NA if nuisances are NA
#' @keywords internal
compute_generated_outcomes_cov_edid <- function(
  panel_obj,
  g,
  t,
  pairs,
  prop_ratios,
  cond_means,
  pt_assumption
) {
  H    <- nrow(pairs)
  n    <- panel_obj$n
  ow   <- panel_obj$outcome_wide

  mask_g  <- panel_obj$cohort_masks[[as.character(g)]]
  pi_g    <- panel_obj$cohort_fractions[[as.character(g)]]
  Ig      <- as.numeric(mask_g)

  col_t   <- panel_obj$period_to_col[[as.character(t)]]
  col_1   <- panel_obj$period_to_col[[as.character(panel_obj$period_1)]]

  gen_out_mat <- matrix(NA_real_, nrow = n, ncol = H)

  for (j in seq_len(H)) {
    gp_j    <- pairs$gp[j]
    tpre_j  <- pairs$tpre[j]

    Igp_j <- if (is.infinite(gp_j)) {
      as.numeric(panel_obj$never_treated_mask)
    } else {
      mask_gp <- panel_obj$cohort_masks[[as.character(gp_j)]]
      if (is.null(mask_gp)) {
        warning(sprintf("compute_generated_outcomes_cov_edid: no mask for gp=%g", gp_j))
        next
      }
      as.numeric(mask_gp)
    }

    r_j  <- prop_ratios[[as.character(gp_j)]]
    m_t  <- cond_means[[paste0(gp_j, "_", t)]]
    m_tp <- cond_means[[paste0(gp_j, "_", tpre_j)]]

    if (is.null(r_j) || is.null(m_t) || is.null(m_tp)) {
      warning(sprintf(
        "compute_generated_outcomes_cov_edid: missing nuisance for pair (gp=%g, tpre=%g).",
        gp_j, tpre_j
      ))
      next
    }

    col_tp <- panel_obj$period_to_col[[as.character(tpre_j)]]

    if (pt_assumption == "post") {
      # PT-Post: single pair (Inf, t_base); tpre is the base period
      # phi_j = (Ig/pi_g - Igp*r_j) * [(Y_t - Y_tpre) - (m_t - m_tp)]
      #       + (m_t - m_tp)
      m_diff <- m_t - m_tp
      y_diff <- (ow[, col_t] - ow[, col_tp])
      phi_j  <- (Ig / pi_g - Igp_j * r_j) * (y_diff - m_diff) + m_diff

    } else {
      # PT-All: Eq. (3.9), using period_1 as baseline
      # phi_j =
      #   (Ig/pi_g) * [(Y_t - Y_1) - (Y_tpre - Y_1)]
      #   - (Igp*r_j) * [(Y_t - Y_1) - (Y_tpre - Y_1)]
      #   + (m_t - m_tp)
      #   - (Ig/pi_g) * (m_t - m_tp)
      #   + (Igp*r_j) * (m_t - m_tp)
      # which simplifies to:
      #   (Ig/pi_g - Igp*r_j) * [(Y_t - Y_tpre) - (m_t - m_tp)] + (m_t - m_tp)
      # (same form as PT-Post; the difference is which tpre/col_tp is used)
      m_diff <- m_t - m_tp
      y_diff <- (ow[, col_t] - ow[, col_tp])
      phi_j  <- (Ig / pi_g - Igp_j * r_j) * (y_diff - m_diff) + m_diff
    }

    gen_out_mat[, j] <- phi_j
  }

  gen_out_mat
}

# ---------------------------------------------------------------------------
# Conditional Omega* (H x H) via Nadaraya-Watson
# ---------------------------------------------------------------------------

#' Compute the averaged conditional covariance matrix Omega*(X)
#'
#' Estimates \eqn{\Omega^* = n^{-1} \sum_i \hat\Omega^*(X_i)} where each
#' \eqn{\hat\Omega^*(X_i)} is a Nadaraya-Watson kernel smoother of the outer
#' product of EIF residuals.
#'
#' Residuals \eqn{e_{i,j}} are obtained by projecting each column of
#' \code{gen_out_mat} onto the B-spline basis and taking the residual.
#'
#' The kernel is Gaussian with per-dimension bandwidths from
#' \code{stats::bw.nrd0()}.
#'
#' \strong{Computational complexity}: O(n^2 * H^2). Emits a warning for
#' n > 1000.
#'
#' @param panel_obj panel object (needs \code{covariate_matrix})
#' @param gen_out_mat numeric matrix n x H
#' @param bw numeric vector length d or NULL (auto from \code{bw.nrd0})
#'
#' @return numeric matrix H x H (positive semi-definite)
#' @keywords internal
compute_omega_star_cov_edid <- function(panel_obj, gen_out_mat, bw = NULL) {
  X_mat <- panel_obj$covariate_matrix
  n     <- nrow(X_mat)
  d     <- ncol(X_mat)
  H     <- ncol(gen_out_mat)

  if (n > 1000L) {
    warning(sprintf(
      "compute_omega_star_cov_edid: n=%d > 1000; O(n^2) kernel loop may be slow.", n
    ))
  }

  # Step 1: compute residual matrix (n x H)
  # Project each column of gen_out_mat onto the full basis; take residuals.
  B_all_obj <- build_basis_matrix_edid(X_mat, bs_df = 4L)
  B_all     <- unclass(B_all_obj)
  attr(B_all, "bs_objects") <- NULL

  e_mat <- matrix(NA_real_, nrow = n, ncol = H)
  for (j in seq_len(H)) {
    phi_j <- gen_out_mat[, j]
    if (anyNA(phi_j)) {
      e_mat[, j] <- NA_real_
      next
    }
    fit_j       <- solve_ols_edid(B_all, phi_j)
    e_mat[, j]  <- fit_j$residuals
  }

  # Step 2: bandwidths
  if (is.null(bw)) {
    bw <- numeric(d)
    for (k in seq_len(d)) {
      h_k <- tryCatch(stats::bw.nrd0(X_mat[, k]), error = function(e) 0)
      if (!is.finite(h_k) || h_k < .Machine$double.eps) {
        warning(sprintf(
          "compute_omega_star_cov_edid: bandwidth for covariate %d is 0 or NA; using h=1.", k
        ))
        h_k <- 1
      }
      bw[k] <- h_k
    }
  }

  # Step 3: kernel-weighted Omega* at each unit i, then average
  # Use vectorised outer computation for efficiency
  # K_mat[i, ell] = product_k dnorm((X[ell,k] - X[i,k]) / h_k) / h_k
  # We compute this column-by-column over covariates.

  # Initialise: K_mat is n x n  (columns = ell, rows = i)
  K_mat <- matrix(1, nrow = n, ncol = n)
  for (k in seq_len(d)) {
    diff_k  <- outer(X_mat[, k], X_mat[, k], "-")  # [i, ell] = X[i,k] - X[ell,k]
    K_mat   <- K_mat * stats::dnorm(diff_k / bw[k]) / bw[k]
  }

  # Accumulate Omega*
  Omega_hat <- matrix(0, nrow = H, ncol = H)
  n_valid   <- 0L

  for (i in seq_len(n)) {
    K_i   <- K_mat[i, ]
    K_sum <- sum(K_i)
    if (!is.finite(K_sum) || K_sum < 1e-15) next

    # e_i: H-vector of residuals at unit i
    e_i <- e_mat[i, ]
    if (anyNA(e_i)) next

    # Omega_i[j,k] = sum_ell K_i[ell] * e_mat[ell, j] * e_mat[ell, k] / K_sum
    # = (e_mat' * diag(K_i) * e_mat)[j,k] / K_sum
    # but e_mat rows may have NAs; handle via masking
    valid_ell <- which(!apply(e_mat, 1L, anyNA))
    K_v       <- K_i[valid_ell]
    e_v       <- e_mat[valid_ell, , drop = FALSE]
    K_sum_v   <- sum(K_v)

    if (!is.finite(K_sum_v) || K_sum_v < 1e-15) next

    # Weighted outer sum: t(e_v) %*% diag(K_v) %*% e_v
    Omega_i <- t(e_v * K_v) %*% e_v / K_sum_v

    Omega_hat <- Omega_hat + Omega_i
    n_valid   <- n_valid + 1L
  }

  if (n_valid == 0L) {
    warning("compute_omega_star_cov_edid: no valid kernel evaluations; returning zero matrix.")
    return(matrix(0, nrow = H, ncol = H))
  }

  Omega_hat / n_valid
}

# ---------------------------------------------------------------------------
# EIF with covariate adjustment
# ---------------------------------------------------------------------------

#' Compute the efficient influence function for a cell with covariates
#'
#' Given the n x H matrix of doubly-robust generated outcomes and the
#' H-vector of efficient weights, returns the zero-mean EIF for this cell.
#'
#' @param panel_obj panel object (used for n)
#' @param gen_out_mat numeric matrix n x H
#' @param weights numeric vector length H summing to 1
#' @param att_gt scalar point estimate
#'
#' @return numeric vector length n, zero-mean
#' @keywords internal
compute_eif_cov_edid <- function(panel_obj, gen_out_mat, weights, att_gt) {
  eif <- drop(gen_out_mat %*% weights) - att_gt
  # Numerical centering for exact zero mean
  eif <- eif - mean(eif)
  eif
}
