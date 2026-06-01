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
#' For cross-cohort pairs (gp != g), the doubly-robust generated outcome of
#' Eq. (4.4) applies:
#' \deqn{\tilde{Y} = (G_g/\pi_g)\,(Y_t - Y_1 - m_{\infty,t,tpre} - m_{g',tpre,1})
#'        - \frac{p_g}{p_\infty}\,\frac{G_\infty}{\pi_g}\,(Y_t - Y_{tpre} - m_{\infty,t,tpre})
#'        - \frac{p_g}{p_{g'}}\,\frac{G_{g'}}{\pi_g}\,(Y_{tpre} - Y_1 - m_{g',tpre,1}),}
#' where \eqn{m_{\infty,t,tpre} = m_{\infty,t,1} - m_{\infty,tpre,1}}. The
#' treated-cohort (G=g) term subtracts \emph{both} conditional-mean adjustments,
#' \eqn{m_{\infty,t,tpre}} and \eqn{m_{g',tpre,1}}: this is what makes the moment
#' doubly robust, identifying \eqn{ATT(g,t)} when either the outcome models or
#' the propensity ratios (but not necessarily both) are correctly specified.
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

    # Determine if this is a self-comparison / two-period DiD pair.
    # Under PT-Post the single moment pair is (gp = Inf, tpre = g-1-anticipation): it is the
    # two-period DiD (treated vs never-treated) on Y_t - Y_tpre, i.e. the SAME structure as a
    # self-comparison pair (Eq 3.2), NOT the three-term cross-cohort formula. Route it to the
    # self/two-period branch so the base period is tpre (= g-1), not period_1. (Without this the
    # cross branch algebraically collapses to the PT-All moment (Y_t - Y_1 - m_{Inf,t,1}) for later
    # cohorts g >= 3, biasing PT-Post; g = 2 coincides since g-1 = period_1. H = 1 here, so the
    # pointwise weight is trivially 1 and Omega needs no PT-Post special case.)
    is_self <- (is.finite(gp_j) && gp_j == g) || identical(pt_assumption, "post")

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
#' \strong{Computational complexity}: O(n^2 * H^2). For large n (> 5000) an
#' informational message is emitted in interactive sessions; silence it with
#' \code{options(edid.quiet = TRUE)}.
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
                                        bw = NULL,
                                        return_pointwise = FALSE) {
  X_mat <- panel_obj$covariate_matrix
  n     <- nrow(X_mat)
  d     <- ncol(X_mat)
  H     <- nrow(pairs)
  ow    <- panel_obj$outcome_wide

  # Performance note (not a correctness condition): the kernel loop is
  # O(n^2 * H^2), which is only a concern for large n. Surface it as an
  # informational message in interactive sessions, silenceable via
  # options(edid.quiet = TRUE); never as a warning (n in the thousands is
  # ordinary for DiD and does not indicate anything wrong with the results).
  if (n > 5000L && interactive() && !isTRUE(getOption("edid.quiet"))) {
    message(sprintf(
      "compute_omega_star_cov_edid: n=%d; the O(n^2) kernel loop may be slow.", n
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
    # Kernel-weighted conditional covariance over the group, for every evaluation point i:
    #   Cov_K(A,B | X_i) = E_K[AB | X_i] - E_K[A | X_i] E_K[B | X_i],   E_K[. | X_i] = (K_i .)/(K_i 1).
    # This weighted-sums form uses only matrix-VECTOR products, so it is algebraically identical to the
    # two-pass residual form  sum_ell K (A_ell-mu_A_i)(B_ell-mu_B_i)/sum_ell K  but never materializes the
    # n x n_group residual matrices (the memory/time bottleneck at large n). A and B are centered by their
    # group means first -- covariance is shift-invariant, so this is exact, and keeping the products small
    # avoids the catastrophic cancellation that the naive E_K[AB]-E_K[A]E_K[B] would suffer near zero cov.
    idx <- which(group_mask)
    if (length(idx) < 2L) return(rep(0, n))

    K_group <- K_mat[, idx, drop = FALSE]  # n x n_group
    K_sums  <- rowSums(K_group)
    K_sums[K_sums < 1e-15] <- NA_real_

    A_c <- A[idx] - mean(A[idx])           # center (shift-invariant) for numerical stability
    B_c <- B[idx] - mean(B[idx])

    mu_A  <- drop(K_group %*% A_c) / K_sums          # E_K[A - Abar | X_i]
    mu_B  <- drop(K_group %*% B_c) / K_sums          # E_K[B - Bbar | X_i]
    mu_AB <- drop(K_group %*% (A_c * B_c)) / K_sums  # E_K[(A-Abar)(B-Bbar) | X_i]
    cov_vals <- mu_AB - mu_A * mu_B                  # = Cov_K(A,B | X_i)
    cov_vals[is.na(cov_vals)] <- 0
    cov_vals
  }

  # -----------------------------------------------------------------------
  # Step 4: Build Omega* by computing each (j,k) element via Eq. (3.12)
  # then averaging over units
  # -----------------------------------------------------------------------
  Omega_hat <- matrix(0, nrow = H, ncol = H)
  # Per-unit Omega*(X_i) array (n x H x H), built only when requested (paper's
  # pointwise efficient weights use Omega*(X_i)^{-1} per observation).
  Omega_array <- if (return_pointwise) array(0, dim = c(n, H, H)) else NULL

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
      # g'_j is the TRUE comparison-cohort label (= g for self-pairs). Conditioning must be on
      # G=g'_j with prefactor 1/p_{g'_j}, exactly as printed in Eq (3.12) and as the no-covariate
      # path does (edid-nocov.R term_d at gp_j==gp_k). Do NOT remap self-pairs to G=Inf: that
      # conditions term5 on the never-treated pre-period covariance instead of the treated
      # cohort's own, corrupting Omega* (and the efficient weights, e.g. negative weights) whenever
      # the cohorts have different pre-period covariance.
      term5 <- 0
      gp_j_eff <- gp_j
      gp_k_eff <- gp_k
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

      # Per-unit Omega*[j,k](X_i), then its average over units.
      omega_jk_i <- term1 + term2 + term3 + term4 + term5
      omega_jk   <- mean(omega_jk_i)

      Omega_hat[j, k] <- omega_jk
      if (k != j) Omega_hat[k, j] <- omega_jk
      if (return_pointwise) {
        Omega_array[, j, k] <- omega_jk_i
        if (k != j) Omega_array[, k, j] <- omega_jk_i
      }
    }
  }

  # Per-unit array path: shrink each pointwise Omega*(X_i) toward the pooled Omega-bar
  # (= Omega_hat) before returning; stabilization/inversion is done downstream by
  # compute_pointwise_weights_edid().
  #
  # Why shrink. The pointwise estimator must estimate an H x H conditional covariance LOCALLY
  # (kernel), which is far noisier than the single pooled Omega-bar used by the constant-weight
  # ("averaged") scheme. When Omega*(X) varies little in X (the common case under good overlap),
  # that extra noise inflates the variance of the efficient estimator BELOW the efficiency the
  # bound promises -- it can do worse than "averaged" in finite samples. Shrinking toward Omega-bar
  # with a data-driven intensity lambda removes that noise: lambda -> 1 (revert to the stable pooled
  # weight) when the across-unit spread of Omega*(X_i) is mostly sampling noise, and lambda -> 0
  # (keep the pointwise weights) when the spread reflects genuine shape variation. Because the
  # kernel estimate sharpens as n grows, lambda -> 0 asymptotically and the estimator coincides with
  # the paper's pointwise-efficient estimator in the limit (the shrinkage is an asymptotically
  # negligible finite-sample regularization, like the eigenvalue floor). Ledoit-Wolf-style rule:
  # lambda = (within-unit sampling variance) / (across-unit variance of Omega*(X_i)), capped to [0,1].
  if (return_pointwise) {
    Hh      <- dim(Omega_array)[2]
    lam_opt <- suppressWarnings(as.numeric(getOption("edid_shrink_lambda", NA_real_)))  # NA = data-driven; 0 disables
    if (length(lam_opt) == 1L && is.finite(lam_opt)) {
      lam <- min(1, max(0, lam_opt))
    } else {
      ksum      <- rowSums(K_mat); ksq <- rowSums(K_mat^2)        # Kish effective local sample size
      m_eff     <- stats::median(ksum^2 / pmax(ksq, .Machine$double.eps))
      shape_var <- mean(apply(Omega_array, c(2, 3), stats::var))  # across-unit spread (signal + noise)
      dg        <- diag(Omega_hat)
      samp_var  <- mean(outer(dg, dg) + Omega_hat^2) / max(m_eff, 1)  # within-unit kernel sampling noise
      lam       <- min(1, max(0, samp_var / max(shape_var, .Machine$double.eps)))
    }
    if (isTRUE(getOption("edid_diag_lambda")))
      options(edid_lambda_acc = c(getOption("edid_lambda_acc", numeric(0)), lam))  # diagnostic accumulator
    if (lam > 0)
      for (jj in seq_len(Hh)) for (kk in seq_len(Hh))
        Omega_array[, jj, kk] <- (1 - lam) * Omega_array[, jj, kk] + lam * Omega_hat[jj, kk]
    attr(Omega_array, "shrink_lambda") <- lam
    return(Omega_array)
  }

  # Ensure positive semi-definiteness via eigenvalue floor (averaged matrix)
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
#' The estimator is the ratio \eqn{\widehat{ATT}_{g,t} = \mathbb{E}_n[w' \tilde{Y}] /
#' \mathbb{E}_n[G_g]} (the \eqn{G_g/\pi_g} factors inside \eqn{\tilde{Y}} make it a
#' ratio in \eqn{\widehat\pi_g}). Its first-order influence function is
#' \deqn{EIF_i = w(X_i)' \tilde{Y}_i - \frac{G_{g,i}}{\pi_g} ATT(g,t),}
#' i.e. the centering is \eqn{-(G_{g,i}/\pi_g)\,ATT}, NOT the constant \eqn{-ATT}.
#' The constant centering omits the first-order contribution of the estimated
#' treated-cohort share \eqn{\widehat\pi_g = \mathbb{E}_n[G_g]} and inflates the
#' variance by \eqn{ATT^2 (1/\pi_g - 1)} with no asymptotic shrinkage. The
#' standard error is \eqn{\widehat{SE} = \sqrt{\sum_i EIF_i^2}/n}.
#'
#' @param panel_obj panel object (needs cohort_masks, cohort_fractions)
#' @param gen_out_mat numeric matrix n x H (generated outcomes)
#' @param weights either a length-H vector (constant weights) or an n x H matrix
#'   of per-observation pointwise weights \eqn{w(X_i)}
#' @param att_gt scalar point estimate (= sum_j w_j * colMeans(gen_out_mat))
#' @param g scalar: target treatment cohort (unused; kept for API compatibility)
#'
#' @return numeric vector length n, mean approximately 0
#' @keywords internal
compute_eif_cov_edid <- function(panel_obj, gen_out_mat, weights, att_gt, g) {
  # Correct first-order influence function for the ratio estimator
  #   ATT_hat = E_n[w' Ytilde] / E_n[G_g]:
  #     EIF_i = w' Ytilde_i - (G_{g,i} / pi_g) * ATT.
  # The constant centering (w' Ytilde_i - ATT) omits the first-order influence
  # of the estimated treated-cohort share pi_hat_g = E_n[G_g]; it inflates the
  # variance by ATT^2 (1/pi_g - 1) with no asymptotic shrinkage. The mean of
  # EIF below is 0 by construction (E_n[G_g] = pi_g), so no de-meaning is used.
  Gg   <- as.numeric(panel_obj$cohort_masks[[as.character(g)]])
  pi_g <- panel_obj$cohort_fractions[[as.character(g)]]
  # w' Ytilde_i : constant weights (length-H vector) or pointwise weights (n x H matrix)
  wY   <- if (is.matrix(weights)) rowSums(gen_out_mat * weights) else drop(gen_out_mat %*% weights)
  wY - (Gg / pi_g) * att_gt
}

#' ACH (Ackerberg, Chen & Hahn 2012) first-step nuisance-estimation correction
#'
#' Returns the length-n vector to SUBTRACT from the plug-in EIF so the influence function
#' accounts for estimation of the first-step sieve nuisances entering the generated outcomes
#' --- the conditional means \eqn{m} and propensity ratios \eqn{r}. The corrected EIF is
#' \eqn{\psi_i - \sum_k [\,\text{score}_k \, H_k^{-1} \Gamma_k\,]_i}, where
#' \eqn{\Gamma_k = \partial E_n[w'\tilde Y]/\partial\theta_k} is the pathwise derivative of the
#' UNCENTERED weighted moment (the centered \eqn{\psi} is mean-zero, so its derivative is the
#' wrong, ~0 object). \eqn{\Gamma_k} is computed numerically by perturbing the fitted prediction
#' along each basis direction and recomputing \eqn{\tilde Y}, with the WEIGHTS HELD FIXED so the
#' \eqn{\Omega}/weight-estimation channel is not re-introduced or double-counted (production keeps
#' \eqn{\Omega} fixed). This is a practical (numerical) form of the ACH two-step variance estimator;
#' \eqn{\tilde Y} is linear in each prediction, so the finite difference is exact up to roundoff.
#' Valid for the plug-in (K = 1, train = test = full) regime; \code{fit_edid_cells} enforces this.
#'
#' @param panel_obj,g,t,pairs,pt_assumption as in \code{compute_generated_outcomes_cov_edid}
#' @param prop_ratios,cond_means named lists of fitted nuisance prediction vectors
#' @param weights frozen weights: length-H vector or n x H matrix (NOT recomputed here)
#' @param m_aux,r_aux named lists (keyed as \code{cond_means}/\code{prop_ratios}) of per-nuisance
#'   pieces \code{list(B_test, score_mat, H_inv, is_fallback)} from the \code{return_aux} path
#' @param eps_rel relative finite-difference step
#' @return numeric vector length n (the term to subtract from the plug-in EIF)
#' @keywords internal
compute_ach_correction_cov_edid <- function(panel_obj, g, t, pairs, prop_ratios,
                                            cond_means, weights, m_aux, r_aux,
                                            pt_assumption = "all", eps_rel = 1e-6) {
  wmoment <- function(pr, cm) {
    go <- compute_generated_outcomes_cov_edid(panel_obj, g, t, pairs, pr, cm, pt_assumption)
    if (is.matrix(weights)) rowSums(go * weights) else drop(go %*% weights)
  }
  m0         <- mean(wmoment(prop_ratios, cond_means))   # uncentered moment at theta_hat
  n          <- panel_obj$n
  correction <- numeric(n)

  # one nuisance key: numerical Gamma along each basis column, then score %*% (H_inv %*% Gamma)
  add_term <- function(correction, a, base, recompute) {
    if (is.null(a) || isTRUE(a$is_fallback) || is.null(a$B_test)) return(correction)
    B   <- a$B_test; p <- ncol(B)
    eps <- eps_rel * (1 + max(abs(base)))
    Gamma <- vapply(seq_len(p), function(j) (mean(recompute(base + eps * B[, j])) - m0) / eps, numeric(1))
    correction + as.vector(a$score_mat %*% drop(a$H_inv %*% Gamma))
  }

  for (key in names(r_aux)) {                            # propensity ratios r_{g,gp}
    correction <- add_term(correction, r_aux[[key]], prop_ratios[[key]],
                           function(newp) { pr <- prop_ratios; pr[[key]] <- newp; wmoment(pr, cond_means) })
  }
  for (key in names(m_aux)) {                            # conditional means m_{gp,period,1}
    correction <- add_term(correction, m_aux[[key]], cond_means[[key]],
                           function(newp) { cm <- cond_means; cm[[key]] <- newp; wmoment(prop_ratios, cm) })
  }
  correction
}

#' Pointwise efficient weights w(X_i) = Omega*(X_i)^(-1) 1 / (1' Omega*(X_i)^(-1) 1)
#'
#' Per-observation semiparametric-efficient weights from the conditional-covariance array. Each
#' Omega*(X_i) is regularized by a DIMENSION-AWARE relative eigenvalue floor.
#' The kernel Omega*(X) is estimated at the uniform Nadaraya-Watson rate rho_n,
#' whose variance exponent is (5-d)/10 (product Gaussian kernel, per-covariate
#' bw.nrd0 ~ n^(-1/5); d = number of covariates). Asymptotic negligibility
#' requires the floor TOL = f(n) to vanish but DOMINATE rho_n, i.e. TOL = n^(-a)
#' with 0 < a < (5-d)/10. We take a = c*(5-min(d,4))/10 with c = 0.7 (strictly
#' interior to the admissible band for d <= 4; for d >= 5 the band (0,(5-d)/10) is
#' EMPTY -- d is clamped to 4 giving fallback a = 0.07, where efficiency is no
#' longer claimed, see the d >= 5 warning in fit_edid_cells): the floor is
#' asymptotically negligible
#' (estimator stays pointwise-efficient and the plug-in SE is consistent in the
#' limit) yet stays above the NW eigenvalue noise for finite-sample stability.
#' Condition number is capped at ~n^(a). Pure per-unit inversion (no floor) is
#' unstable: a few near-singular Omega*(X_i) produce enormous weights. Degenerate
#' units fall back to uniform (1/H).
#'
#' @param omega_array numeric array n x H x H of per-unit Omega*(X_i), from
#'   \code{compute_omega_star_cov_edid(..., return_pointwise = TRUE)}
#' @param d integer, number of covariates entering the kernel (sets the floor rate)
#' @return numeric matrix n x H, each row summing to 1
#' @keywords internal
compute_pointwise_weights_edid <- function(omega_array, d = 1L) {
  n   <- dim(omega_array)[1]
  H   <- dim(omega_array)[2]
  one <- rep(1, H)
  W   <- matrix(NA_real_, n, H)
  # Dimension-aware floor exponent a = c*(5-d)/10, c = 0.7. clamp d to <=4 so a>0
  # (the band (0,(5-d)/10) is empty for d>=5, where the NW conditional covariance
  # is not uniformly consistent; a=0.07 is a conservative fallback there).
  a_floor <- 0.7 * (5 - min(as.integer(d), 4L)) / 10
  tol     <- n^(-a_floor)
  # Diagnostic override (default behavior unchanged): getOption("edid_eig_tol") sets the
  # relative eigenvalue floor directly (condition-number cap = 1/tol). Used to study how
  # regularization strength affects pointwise efficiency; NA/unset keeps the rate above.
  tol_ov <- suppressWarnings(as.numeric(getOption("edid_eig_tol", NA_real_)))
  if (length(tol_ov) == 1L && is.finite(tol_ov) && tol_ov > 0) tol <- tol_ov
  for (i in seq_len(n)) {
    Mi <- omega_array[i, , ]
    Mi <- 0.5 * (Mi + t(Mi))
    e  <- eigen(Mi, symmetric = TRUE)
    mx <- max(e$values)
    if (!is.finite(mx) || mx <= 0) { W[i, ] <- one / H; next }
    ev_floored <- pmax(e$values, mx * tol)
    Minv <- e$vectors %*% diag(1 / ev_floored, H) %*% t(e$vectors)
    v    <- drop(Minv %*% one)
    den  <- sum(v)
    W[i, ] <- if (is.finite(den) && abs(den) > 1e-12) v / den else one / H
  }
  W
}
