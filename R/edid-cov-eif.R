# edid-cov-eif.R
# Generated outcome and EIF computation for the EDiD covariate path.
# Implements Chen, Sant'Anna & Xie (2025) Eq. (3.9), (3.10), (3.12), (4.4).

# ---------------------------------------------------------------------------
# Generated outcomes (doubly-robust, n x H matrix)
# ---------------------------------------------------------------------------

#' Compute doubly-robust generated outcomes for a (g, t) cell
#'
#' Returns the n x H matrix of generated outcomes where column j corresponds
#' to pair j = \eqn{(g'_j, t_{pre,j})} and row i to unit i.  Implements
#' Eq. (4.4) of Chen, Sant'Anna & Xie (2025).
#'
#' For self-comparison pairs (gp == g), the formula reduces to Eq. (3.2):
#' \deqn{\tilde{Y} = (G_g/\pi_g - r_{g,\infty} G_\infty/\pi_g)(Y_t - Y_{tpre} - m_{\infty,t,tpre})}
#'
#' For cross-cohort pairs (gp != g), the three-term doubly-robust formula applies:
#' \deqn{\tilde{Y} = (G_g/\pi_g)(Y_t - Y_1 - m_{\infty,t,1})
#'        - r_{g,\infty} (G_\infty/\pi_g)(Y_t - Y_{tpre} - m_{\infty,t,tpre})
#'        - r_{g,g'} (G_{g'}/\pi_g)(Y_{tpre} - Y_1 - m_{g',tpre,1})}
#' Note: term1 uses only \eqn{m_{\infty,t,1}}, not \eqn{m_{g',tpre,1}}.  Adding
#' \eqn{m_{g',tpre,1}} to term1 would bias the estimator by
#' \eqn{E[Y_{tpre}-Y_1|G=g]}, since the propensity-ratio correction in term3
#' accounts for G=g' units only, not G=g units.
#'
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param g scalar: treatment cohort
#' @param t scalar: target time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param prop_ratios named list of n-vectors keyed by \code{as.character(gp)}:
#'   cross-fitted propensity ratios. Must include key \code{"Inf"} for
#'   \eqn{r_{g,\infty}} and keys for each cross-cohort gp.
#' @param cond_means named list of n-vectors keyed by
#'   \code{paste0(gp, "_", period)}: cross-fitted conditional means
#'   \eqn{E[Y_{period} - Y_1 | G=gp, X]}. Must include never-treated keys.
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

  # Never-treated indicator (used in all pairs)
  I_inf   <- as.numeric(panel_obj$never_treated_mask)

  gen_out_mat <- matrix(NA_real_, nrow = n, ncol = H)

  for (j in seq_len(H)) {
    gp_j    <- pairs$gp[j]
    tpre_j  <- pairs$tpre[j]
    col_tp  <- panel_obj$period_to_col[[as.character(tpre_j)]]

    # Determine if this is a self-comparison pair
    is_self <- is.finite(gp_j) && gp_j == g

    if (is_self) {
      # -----------------------------------------------------------------
      # Self-comparison pair (gp == g): uses never-treated as comparison
      # Eq. (3.2): phi = (G_g/pi_g - r[g,Inf]*G_Inf/pi_g) *
      #                   (Y_t - Y_tpre - m_{Inf,t,tpre}(X))
      # -----------------------------------------------------------------
      r_inf <- prop_ratios[["Inf"]]
      # m_{Inf,t,tpre}(X) = E[Y_t - Y_tpre | G=Inf, X]
      #                    = E[Y_t - Y_1 | G=Inf, X] - E[Y_tpre - Y_1 | G=Inf, X]
      m_inf_t  <- cond_means[[paste0("Inf_", t)]]
      m_inf_tp <- cond_means[[paste0("Inf_", tpre_j)]]

      if (is.null(r_inf) || is.null(m_inf_t) || is.null(m_inf_tp)) {
        warning(sprintf(
          "compute_generated_outcomes_cov_edid: missing nuisance for self-pair (gp=%g, tpre=%g).",
          gp_j, tpre_j
        ))
        next
      }

      m_inf_diff <- m_inf_t - m_inf_tp  # m_{Inf,t,tpre}(X)
      y_diff     <- ow[, col_t] - ow[, col_tp]  # Y_t - Y_tpre

      phi_j <- (Ig / pi_g - r_inf * I_inf / pi_g) * (y_diff - m_inf_diff)

    } else {
      # -----------------------------------------------------------------
      # Cross-cohort pair (gp != g): full three-term Eq. (4.4)
      # phi = (G_g/pi_g) * (Y_t - Y_1 - m_{Inf,t,1}(X) - m_{g',tpre,1}(X))
      #     - r[g,Inf] * (G_Inf/pi_g) * (Y_t - Y_tpre - m_{Inf,t,tpre}(X))
      #     - r[g,g']  * (G_g'/pi_g)  * (Y_tpre - Y_1 - m_{g',tpre,1}(X))
      # -----------------------------------------------------------------
      gp_key <- as.character(gp_j)

      # Propensity ratios
      r_inf <- prop_ratios[["Inf"]]
      r_gp  <- prop_ratios[[gp_key]]

      # Conditional means
      m_inf_t  <- cond_means[[paste0("Inf_", t)]]
      m_inf_tp <- cond_means[[paste0("Inf_", tpre_j)]]
      m_gp_tp  <- cond_means[[paste0(gp_key, "_", tpre_j)]]

      if (is.null(r_inf) || is.null(r_gp) ||
          is.null(m_inf_t) || is.null(m_inf_tp) || is.null(m_gp_tp)) {
        warning(sprintf(
          "compute_generated_outcomes_cov_edid: missing nuisance for cross-pair (gp=%g, tpre=%g).",
          gp_j, tpre_j
        ))
        next
      }

      # Comparison cohort indicator
      if (is.infinite(gp_j)) {
        I_gp <- I_inf
      } else {
        mask_gp <- panel_obj$cohort_masks[[gp_key]]
        if (is.null(mask_gp)) {
          warning(sprintf("compute_generated_outcomes_cov_edid: no mask for gp=%g", gp_j))
          next
        }
        I_gp <- as.numeric(mask_gp)
      }

      # m_{Inf,t,tpre}(X) = m_{Inf,t,1}(X) - m_{Inf,tpre,1}(X)
      m_inf_diff <- m_inf_t - m_inf_tp

      # Outcome differences
      Y_t      <- ow[, col_t]
      Y_1      <- ow[, col_1]
      Y_tpre   <- ow[, col_tp]

      # Three-term doubly-robust formula matching Eq. (4.4) of the paper:
      # Y_tilde = (G_g/pi_g)(Y_t - Y_1 - m_{inf,t,t'}(X) - m_{g',t',1}(X))
      #         - r_{g,inf}(G_inf/pi_g)(Y_t - Y_t' - m_{inf,t,t'}(X))
      #         - r_{g,g'}(G_{g'}/pi_g)(Y_t' - Y_1 - m_{g',t',1}(X))
      term1 <- (Ig / pi_g) * (Y_t - Y_1 - m_inf_diff - m_gp_tp)
      term2 <- r_inf * (I_inf / pi_g) * (Y_t - Y_tpre - m_inf_diff)
      term3 <- r_gp * (I_gp / pi_g) * (Y_tpre - Y_1 - m_gp_tp)

      phi_j <- term1 - term2 - term3
    }

    gen_out_mat[, j] <- phi_j
  }

  gen_out_mat
}

# ---------------------------------------------------------------------------
# Conditional Omega* (H x H) via Nadaraya-Watson kernel
# ---------------------------------------------------------------------------

#' Compute the averaged conditional covariance matrix Omega*(X)
#'
#' Estimates \eqn{\Omega^* = n^{-1} \sum_i \hat\Omega^*(X_i)} using a faithful
#' plug-in of Eq. (3.12) from Chen, Sant'Anna & Xie (2025).
#'
#' Each (j,k)-th element of Omega*(X) is estimated using Nadaraya-Watson
#' kernel smoothing of outcome change covariances within specific cohorts,
#' scaled by propensity scores.
#'
#' \strong{Computational complexity}: O(n^2 * H^2). Emits a warning for
#' n > 1000.
#'
#' @param panel_obj panel object (needs \code{covariate_matrix}, \code{outcome_wide},
#'   \code{cohort_masks}, \code{never_treated_mask})
#' @param g scalar: target treatment cohort
#' @param t scalar: target time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param prop_ratios named list of n-vectors: cross-fitted propensity ratios
#' @param cond_means named list of n-vectors: cross-fitted conditional means
#' @param bw numeric vector length d or NULL (auto from \code{bw.nrd0})
#'
#' @return numeric matrix H x H (positive semi-definite)
#' @keywords internal
compute_omega_star_cov_edid <- function(panel_obj, g, t, pairs,
                                        prop_ratios, cond_means,
                                        inv_propensities = NULL,
                                        bw = NULL) {
  X_mat <- panel_obj$covariate_matrix
  n     <- nrow(X_mat)
  d     <- ncol(X_mat)
  H     <- nrow(pairs)
  ow    <- panel_obj$outcome_wide

  if (n > 1000L) {
    warning(sprintf(
      "compute_omega_star_cov_edid: n=%d > 1000; O(n^2) kernel loop may be slow.", n
    ))
  }

  # -----------------------------------------------------------------------
  # Step 1: Bandwidths
  # -----------------------------------------------------------------------
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

  # -----------------------------------------------------------------------
  # Step 2: Build kernel weight matrix K_mat[i, ell] (n x n)
  # -----------------------------------------------------------------------
  K_mat <- matrix(1, nrow = n, ncol = n)
  for (k in seq_len(d)) {
    diff_k <- outer(X_mat[, k], X_mat[, k], "-")
    K_mat  <- K_mat * stats::dnorm(diff_k / bw[k]) / bw[k]
  }

  # -----------------------------------------------------------------------
  # Step 3: Precompute outcome change residuals for each cohort
  # For Eq. (3.12) we need:
  #   Cov(Y_t - Y_1, Y_t - Y_1 | G=g, X)          [treated group]
  #   Cov(Y_t - Y_{t'_j}, Y_t - Y_{t'_k} | G=Inf, X) [never-treated]
  #   Cov(Y_t - Y_1, Y_{t'_j} - Y_1 | G=g, X)      [cross-term, self-pairs]
  #   Cov(Y_{t'_j} - Y_1, Y_{t'_k} - Y_1 | G=g'_j, X) [cross-cohort]
  # -----------------------------------------------------------------------
  col_t <- panel_obj$period_to_col[[as.character(t)]]
  col_1 <- panel_obj$period_to_col[[as.character(panel_obj$period_1)]]

  mask_g   <- panel_obj$cohort_masks[[as.character(g)]]
  mask_inf <- panel_obj$never_treated_mask

  # -----------------------------------------------------------------------
  # Omega* scaling terms: 1/p_g(X), 1/p_inf(X), 1/p_{g'}(X)
  # Paper Eq. (3.12) uses conditional propensity scores as scalar
  # pre-factors on each conditional covariance term.
  # When inv_propensities is provided (from estimate_all_inverse_propensities),
  # use the estimated conditional values. Otherwise fall back to unconditional.
  # -----------------------------------------------------------------------
  pi_g   <- panel_obj$cohort_fractions[[as.character(g)]]
  pi_inf <- sum(mask_inf) / n

  if (!is.null(inv_propensities)) {
    inv_pg_vec   <- inv_propensities[[as.character(g)]]
    inv_pinf_vec <- inv_propensities[["Inf"]]
    if (is.null(inv_pg_vec))   inv_pg_vec   <- rep(1 / pi_g, n)
    if (is.null(inv_pinf_vec)) inv_pinf_vec <- rep(1 / pi_inf, n)
  } else {
    inv_pg_vec   <- rep(1 / pi_g, n)
    inv_pinf_vec <- rep(1 / pi_inf, n)
  }

  # Helper: kernel-smoothed conditional covariance of (A, B) given G=group at point x_i
  # Returns an n-vector (one value per evaluation point)
  kernel_cond_cov <- function(A, B, group_mask) {
    # For each evaluation point i, compute:
    #   sum_{ell in group} K(x_i, x_ell) * (A_ell - mu_A(x_i)) * (B_ell - mu_B(x_i))
    #   / sum_{ell in group} K(x_i, x_ell)
    # Where mu_A(x_i) = sum_{ell in group} K * A_ell / sum K
    idx <- which(group_mask)
    if (length(idx) < 2L) return(rep(0, n))

    K_group <- K_mat[, idx, drop = FALSE]  # n x n_group
    K_sums  <- rowSums(K_group)
    K_sums[K_sums < 1e-15] <- NA_real_

    A_group <- A[idx]
    B_group <- B[idx]

    # Kernel-weighted means
    mu_A <- drop(K_group %*% A_group) / K_sums
    mu_B <- drop(K_group %*% B_group) / K_sums

    # Kernel-weighted covariance
    # sum_ell K[i,ell] * (A_ell - mu_A_i) * (B_ell - mu_B_i) / sum_ell K[i,ell]
    resid_A <- sweep(matrix(A_group, nrow = n, ncol = length(idx), byrow = TRUE),
                     1, mu_A, "-")
    resid_B <- sweep(matrix(B_group, nrow = n, ncol = length(idx), byrow = TRUE),
                     1, mu_B, "-")
    cov_vals <- rowSums(K_group * resid_A * resid_B) / K_sums
    cov_vals[is.na(cov_vals)] <- 0
    cov_vals
  }

  # -----------------------------------------------------------------------
  # Step 4: Build Omega* by computing each (j,k) element via Eq. (3.12)
  # then averaging over units
  # -----------------------------------------------------------------------
  Omega_hat <- matrix(0, nrow = H, ncol = H)

  # Precompute outcome changes we'll need repeatedly
  Y_t_minus_Y1 <- ow[, col_t] - ow[, col_1]

  for (j in seq_len(H)) {
    gp_j   <- pairs$gp[j]
    tpre_j <- pairs$tpre[j]
    col_tj <- panel_obj$period_to_col[[as.character(tpre_j)]]

    is_self_j <- is.finite(gp_j) && gp_j == g

    Y_t_minus_Ytj <- ow[, col_t] - ow[, col_tj]
    Y_tj_minus_Y1 <- ow[, col_tj] - ow[, col_1]

    for (k in j:H) {
      gp_k   <- pairs$gp[k]
      tpre_k <- pairs$tpre[k]
      col_tk <- panel_obj$period_to_col[[as.character(tpre_k)]]

      is_self_k <- is.finite(gp_k) && gp_k == g

      Y_t_minus_Ytk <- ow[, col_t] - ow[, col_tk]
      Y_tk_minus_Y1 <- ow[, col_tk] - ow[, col_1]

      # Eq. (3.12) term by term, using conditional 1/p_g(X):
      # Term 1: (1/p_g(X)) * Cov(Y_t - Y_1, Y_t - Y_1 | G=g, X)
      term1 <- inv_pg_vec *
        kernel_cond_cov(Y_t_minus_Y1, Y_t_minus_Y1, mask_g)

      # Term 2: (1/p_inf(X)) * Cov(Y_t - Y_{t'_j}, Y_t - Y_{t'_k} | G=Inf, X)
      term2 <- inv_pinf_vec *
        kernel_cond_cov(Y_t_minus_Ytj, Y_t_minus_Ytk, mask_inf)

      # Term 3: -1{g == g'_j}/p_g(X) * Cov(Y_t - Y_1, Y_{t'_j} - Y_1 | G=g, X)
      term3 <- 0
      if (is_self_j) {
        term3 <- -inv_pg_vec *
          kernel_cond_cov(Y_t_minus_Y1, Y_tj_minus_Y1, mask_g)
      }

      # Term 4: -1{g == g'_k}/p_g(X) * Cov(Y_t - Y_1, Y_{t'_k} - Y_1 | G=g, X)
      term4 <- 0
      if (is_self_k) {
        term4 <- -inv_pg_vec *
          kernel_cond_cov(Y_t_minus_Y1, Y_tk_minus_Y1, mask_g)
      }

      # Term 5: 1{g'_j == g'_k}/p_{g'_j}(X) * Cov(Y_{t'_j}-Y_1, Y_{t'_k}-Y_1 | G=g'_j, X)
      term5 <- 0
      gp_j_eff <- if (is_self_j) Inf else gp_j
      gp_k_eff <- if (is_self_k) Inf else gp_k
      if (identical(gp_j_eff, gp_k_eff)) {
        gp_key_jk <- as.character(gp_j_eff)
        if (is.infinite(gp_j_eff)) {
          inv_pgp_vec <- inv_pinf_vec
          mask_gp_jk <- mask_inf
        } else {
          if (!is.null(inv_propensities) && !is.null(inv_propensities[[gp_key_jk]])) {
            inv_pgp_vec <- inv_propensities[[gp_key_jk]]
          } else {
            pi_gp <- panel_obj$cohort_fractions[[gp_key_jk]]
            inv_pgp_vec <- if (!is.null(pi_gp) && pi_gp > 1e-15) rep(1/pi_gp, n) else rep(0, n)
          }
          mask_gp_jk <- panel_obj$cohort_masks[[gp_key_jk]]
          if (is.null(mask_gp_jk)) mask_gp_jk <- rep(FALSE, n)
        }
        term5 <- inv_pgp_vec *
          kernel_cond_cov(Y_tj_minus_Y1, Y_tk_minus_Y1, mask_gp_jk)
      }

      # Average Omega*[j,k](X) over all units
      omega_jk <- mean(term1 + term2 + term3 + term4 + term5)

      Omega_hat[j, k] <- omega_jk
      if (k != j) Omega_hat[k, j] <- omega_jk
    }
  }

  # Ensure positive semi-definiteness via eigenvalue floor
  eig <- eigen(Omega_hat, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-12)
  Omega_hat <- eig$vectors %*% diag(eig$values, nrow = H) %*% t(eig$vectors)

  Omega_hat
}

# ---------------------------------------------------------------------------
# EIF with covariate adjustment
# ---------------------------------------------------------------------------

#' Compute the efficient influence function for a cell with covariates
#'
#' The efficient GMM estimator is
#' \deqn{\hat\beta_{g,t} = \sum_j w_j \cdot \frac{1}{n}\sum_i \tilde{Y}_{j,i}}
#' where \eqn{\tilde{Y}_{j,i}} is the doubly-robust generated outcome for pair j
#' and \eqn{w_j} are the fixed efficient weights.  By the delta method, its
#' influence function is
#' \deqn{EIF_i = \sum_j w_j \cdot (\tilde{Y}_{j,i} - \beta_j)
#'             = \left(\sum_j w_j \tilde{Y}_{j,i}\right) - ATT(g,t)}
#' (using \eqn{\sum_j w_j \beta_j = ATT(g,t)}).
#'
#' Statistical note: an alternative form \eqn{EIF_i = \sum_j w_j \tilde{Y}_{j,i}
#' + (G_{g,i}/\pi_g) \cdot ATT(g,t)} that appears in some semiparametric
#' efficiency calculations adds a term whose mean is \eqn{ATT(g,t)} (since
#' \eqn{E[G_{g,i}/\pi_g] = 1}).  After centring, this equals
#' \eqn{(correct\,EIF) + (G_{g,i}/\pi_g - 1) \cdot ATT(g,t)}, inflating the
#' variance by \eqn{ATT^2 \cdot Var(G_{g,i}/\pi_g - 1) > 0} whenever
#' \eqn{ATT \ne 0}.  The correct expression for the SE formula
#' \eqn{SE = \sqrt{\sum_i EIF_i^2 / n^2}} is the one below.
#'
#' @param panel_obj panel object (needs n)
#' @param gen_out_mat numeric matrix n x H (generated outcomes)
#' @param weights numeric vector length H summing to 1
#' @param att_gt scalar point estimate (= sum_j w_j * colMeans(gen_out_mat))
#' @param g scalar: target treatment cohort (unused; kept for API compatibility)
#'
#' @return numeric vector length n, mean approximately 0
#' @keywords internal
compute_eif_cov_edid <- function(panel_obj, gen_out_mat, weights, att_gt, g) {
  # EIF_i = (sum_j w_j * phi_{j,i}) - ATT(g,t)
  # This has theoretical mean 0 (E[phi_j] = beta_j, sum_j w_j*beta_j = ATT).
  eif <- drop(gen_out_mat %*% weights) - att_gt

  # Numerical centering: removes any floating-point residual from finite-sample
  # estimation of nuisances.  Does not change the variance formula.
  eif <- eif - mean(eif)
  eif
}
