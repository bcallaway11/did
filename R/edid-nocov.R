# edid-nocov.R
# No-covariate path for the EDiD estimator:
#   compute_omega_star_nocov_edid()
#   compute_efficient_weights_edid()
#   compute_generated_outcomes_nocov_edid()
#   compute_eif_nocov_edid()

# ---------------------------------------------------------------------------
# Helper: get column index from panel_obj
# ---------------------------------------------------------------------------
.col <- function(panel_obj, period_val) {
  panel_obj$period_to_col[[as.character(period_val)]]
}

# ---------------------------------------------------------------------------
# Omega* covariance matrix (H x H)
# ---------------------------------------------------------------------------

#' Compute the Omega* covariance matrix for the no-covariate EDiD path
#'
#' Builds the \eqn{H \times H} sample covariance matrix of the identifying
#' moments for cell \code{(target_g, target_t)}.
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param pt_assumption \code{"all"} or \code{"post"}
#'
#' @return numeric matrix H x H
#' @keywords internal
compute_omega_star_nocov_edid <- function(
  target_g, target_t, pairs, panel_obj, pt_assumption
) {
  H   <- nrow(pairs)
  n   <- panel_obj$n
  ow  <- panel_obj$outcome_wide

  mask_g   <- panel_obj$cohort_masks[[as.character(target_g)]]
  mask_inf <- panel_obj$never_treated_mask
  n_g      <- sum(mask_g)
  n_inf    <- sum(mask_inf)

  col_t  <- .col(panel_obj, target_t)
  col_1  <- .col(panel_obj, panel_obj$period_1)

  if (pt_assumption == "post") {
    # PT-Post: 1x1 matrix = var of standard DiD moment
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)

    delta_g   <- ow[mask_g,   col_t] - ow[mask_g,   col_base]
    delta_inf <- ow[mask_inf, col_t] - ow[mask_inf, col_base]

    omega <- matrix(
      cov_nn_edid(delta_g, delta_g) / n_g +
        cov_nn_edid(delta_inf, delta_inf) / n_inf,
      nrow = 1L, ncol = 1L
    )
    return(omega)
  }

  # ---------------------------------------------------------------------------
  # PT-All: H x H matrix, entry-by-entry
  # ---------------------------------------------------------------------------
  # Pre-compute treated-group change (same for all j, k)
  delta_g_t_1 <- ow[mask_g, col_t] - ow[mask_g, col_1]

  # Pre-compute never-treated changes for each unique tpre
  unique_tpre <- unique(pairs$tpre)
  delta_inf_cache <- vector("list", length(unique_tpre))
  names(delta_inf_cache) <- as.character(unique_tpre)
  for (tp in unique_tpre) {
    col_pre <- .col(panel_obj, tp)
    delta_inf_cache[[as.character(tp)]] <-
      ow[mask_inf, col_t] - ow[mask_inf, col_pre]
  }

  # Pre-compute comparison-cohort changes for each unique (gp, tpre)
  unique_gp_tpre <- unique(pairs[, c("gp", "tpre")])
  delta_gp_cache <- list()
  for (rr in seq_len(nrow(unique_gp_tpre))) {
    gp_val  <- unique_gp_tpre$gp[rr]
    tp_val  <- unique_gp_tpre$tpre[rr]
    key     <- paste0(gp_val, "_", tp_val)
    mask_gp <- panel_obj$cohort_masks[[as.character(gp_val)]]
    col_pre <- .col(panel_obj, tp_val)
    delta_gp_cache[[key]] <- ow[mask_gp, col_pre] - ow[mask_gp, col_1]
  }

  omega <- matrix(0, nrow = H, ncol = H)

  for (j in seq_len(H)) {
    gp_j   <- pairs$gp[j]
    tpre_j <- pairs$tpre[j]
    key_j  <- paste0(gp_j, "_", tpre_j)
    n_gp_j <- sum(panel_obj$cohort_masks[[as.character(gp_j)]])
    delta_inf_j <- delta_inf_cache[[as.character(tpre_j)]]
    delta_gp_j  <- delta_gp_cache[[key_j]]

    for (k in seq_len(H)) {
      if (k < j) {
        omega[j, k] <- omega[k, j]  # symmetric
        next
      }
      gp_k   <- pairs$gp[k]
      tpre_k <- pairs$tpre[k]
      key_k  <- paste0(gp_k, "_", tpre_k)
      n_gp_k <- sum(panel_obj$cohort_masks[[as.character(gp_k)]])
      delta_inf_k <- delta_inf_cache[[as.character(tpre_k)]]
      delta_gp_k  <- delta_gp_cache[[key_k]]

      # Term A: treated group variance (always present; same for all j, k)
      term_a <- cov_nn_edid(delta_g_t_1, delta_g_t_1) / n_g

      # Term B: never-treated cross-covariance
      term_b <- cov_nn_edid(delta_inf_j, delta_inf_k) / n_inf

      # Term C_j: non-zero only if gp_j == target_g
      term_cj <- 0
      if (is.finite(gp_j) && gp_j == target_g) {
        term_cj <- cov_nn_edid(delta_g_t_1, delta_gp_j) / n_g
      }

      # Term C_k: non-zero only if gp_k == target_g
      term_ck <- 0
      if (is.finite(gp_k) && gp_k == target_g) {
        term_ck <- cov_nn_edid(delta_g_t_1, delta_gp_k) / n_g
      }

      # Term D: non-zero only if gp_j == gp_k
      term_d <- 0
      if (gp_j == gp_k) {  # works for both finite and Inf
        term_d <- cov_nn_edid(delta_gp_j, delta_gp_k) / n_gp_j
      }

      omega[j, k] <- term_a + term_b - term_cj - term_ck + term_d
    }
  }

  omega
}

# ---------------------------------------------------------------------------
# Efficient weights
# ---------------------------------------------------------------------------

#' Compute efficient inverse-covariance weights
#'
#' Implements \eqn{w = (\Omega^{*-1} \mathbf{1}) / (\mathbf{1}' \Omega^{*-1} \mathbf{1})}
#' with fallback to uniform weights when the matrix is degenerate.
#'
#' @param omega_star numeric matrix H x H
#'
#' @return numeric vector length H, summing to 1
#' @keywords internal
compute_efficient_weights_edid <- function(omega_star) {
  H       <- nrow(omega_star)
  ones_H  <- rep(1, H)
  unif    <- ones_H / H

  # Degenerate: all zeros -> uniform
  if (all(omega_star == 0)) return(unif)

  # Check condition number
  kappa <- check_condition_edid(omega_star)
  if (!is.finite(kappa) || kappa > EDID_COND_THRESH) {
    inv_omega <- compute_pseudoinverse_edid(omega_star)
  } else {
    inv_omega <- tryCatch(
      solve(omega_star),
      error = function(e) compute_pseudoinverse_edid(omega_star)
    )
  }

  num   <- drop(inv_omega %*% ones_H)
  denom <- sum(num)
  if (!is.finite(denom) || abs(denom) < EDID_DENOM_EPS) return(unif)
  num / denom
}

# ---------------------------------------------------------------------------
# Generated outcomes (scalar moments)
# ---------------------------------------------------------------------------

#' Compute generated-outcome scalars for each valid pair
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param pt_assumption \code{"all"} or \code{"post"}
#'
#' @return numeric vector length H
#' @keywords internal
compute_generated_outcomes_nocov_edid <- function(
  target_g, target_t, pairs, panel_obj, pt_assumption
) {
  H   <- nrow(pairs)
  ow  <- panel_obj$outcome_wide

  mask_g   <- panel_obj$cohort_masks[[as.character(target_g)]]
  mask_inf <- panel_obj$never_treated_mask

  col_t  <- .col(panel_obj, target_t)

  if (pt_assumption == "post") {
    # PT-Post: one pair (Inf, base)
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)
    y_hat <- mean(ow[mask_g, col_t] - ow[mask_g, col_base]) -
             mean(ow[mask_inf, col_t] - ow[mask_inf, col_base])
    return(y_hat)  # scalar; will be treated as length-1 vector
  }

  # PT-All
  col_1 <- .col(panel_obj, panel_obj$period_1)
  term_g <- mean(ow[mask_g, col_t] - ow[mask_g, col_1])

  y_hat <- numeric(H)
  for (j in seq_len(H)) {
    gp_j    <- pairs$gp[j]
    tpre_j  <- pairs$tpre[j]
    col_pre <- .col(panel_obj, tpre_j)

    term_inf <- mean(ow[mask_inf, col_t] - ow[mask_inf, col_pre])

    mask_gp  <- panel_obj$cohort_masks[[as.character(gp_j)]]
    term_gp  <- mean(ow[mask_gp, col_pre] - ow[mask_gp, col_1])

    y_hat[j] <- term_g - term_inf - term_gp
  }
  y_hat
}

# ---------------------------------------------------------------------------
# Efficient Influence Function (per unit, length n)
# ---------------------------------------------------------------------------

#' Compute the no-covariate efficient influence function for a (g, t) cell
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param weights numeric vector length H (efficient weights)
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param att_gt scalar ATT estimate for this cell
#' @param pt_assumption \code{"all"} or \code{"post"}
#'
#' @return numeric vector length n (zero-mean by construction)
#' @keywords internal
compute_eif_nocov_edid <- function(
  target_g, target_t, pairs, weights, panel_obj, att_gt, pt_assumption
) {
  n    <- panel_obj$n
  H    <- nrow(pairs)
  ow   <- panel_obj$outcome_wide

  mask_g   <- panel_obj$cohort_masks[[as.character(target_g)]]
  mask_inf <- panel_obj$never_treated_mask
  pi_g     <- panel_obj$cohort_fractions[[as.character(target_g)]]
  pi_inf   <- sum(mask_inf) / n

  col_t <- .col(panel_obj, target_t)

  eif <- numeric(n)

  if (pt_assumption == "post") {
    # PT-Post: one pair; base = g - 1 - anticipation
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)
    w_1      <- weights[1L]

    delta_g   <- ow[mask_g,   col_t] - ow[mask_g,   col_base]
    mean_g    <- mean(delta_g)
    eif[mask_g] <- eif[mask_g] + w_1 * (delta_g - mean_g) / pi_g

    delta_inf  <- ow[mask_inf, col_t] - ow[mask_inf, col_base]
    mean_inf   <- mean(delta_inf)
    eif[mask_inf] <- eif[mask_inf] - w_1 * (delta_inf - mean_inf) / pi_inf

    # gp_j == Inf: no comparison cohort term
  } else {
    # PT-All
    col_1    <- .col(panel_obj, panel_obj$period_1)
    col_base <- col_1   # base for treated group

    delta_g_t_base <- ow[mask_g, col_t] - ow[mask_g, col_base]
    mean_g_t_base  <- mean(delta_g_t_base)

    for (j in seq_len(H)) {
      w_j     <- weights[j]
      gp_j    <- pairs$gp[j]
      tpre_j  <- pairs$tpre[j]
      col_pre <- .col(panel_obj, tpre_j)

      # Treated group contribution (always present)
      # Y_hat_j enters via centering: phi_{ij} = I(G=g)/pi_g * (delta - Y_hat_j)
      # But the EIF is: score - att_gt, and score = sum w_j phi_{ij}
      # phi_{ij} for treated group = I(G=g)/pi_g * (delta_{g,t,base} - Y_hat_j)
      # We compute: sum_j w_j * I(G=g)/pi_g * (delta - Y_hat_j)
      # = I(G=g)/pi_g * [ (delta - mean_g) + (mean_g - sum_j w_j Y_hat_j) ]
      # But it's simpler to accumulate directly:
      eif[mask_g] <- eif[mask_g] +
        w_j * (delta_g_t_base - mean_g_t_base) / pi_g

      # Never-treated contribution (subtract)
      delta_inf_t_pre <- ow[mask_inf, col_t] - ow[mask_inf, col_pre]
      mean_inf_t_pre  <- mean(delta_inf_t_pre)
      eif[mask_inf] <- eif[mask_inf] -
        w_j * (delta_inf_t_pre - mean_inf_t_pre) / pi_inf

      # Comparison cohort contribution (subtract)
      mask_gp          <- panel_obj$cohort_masks[[as.character(gp_j)]]
      pi_gp            <- panel_obj$cohort_fractions[[as.character(gp_j)]]
      delta_gp_pre_base <- ow[mask_gp, col_pre] - ow[mask_gp, col_base]
      mean_gp_pre_base  <- mean(delta_gp_pre_base)
      eif[mask_gp] <- eif[mask_gp] -
        w_j * (delta_gp_pre_base - mean_gp_pre_base) / pi_gp
    }
  }

  # The score accumulated above already has mean 0 by construction:
  # each group's contribution is demeaned (delta_i - mean(delta)).
  # The EIF is the score itself; do NOT subtract att_gt again.
  eif
}
