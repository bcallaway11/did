# edid-fit.R
# Outer (g, t) cell loop for the EDiD estimator.

#' Fit all (g, t) cells for the EDiD estimator
#'
#' Iterates over all treatment cohorts and all time periods (excluding
#' \code{period_1}), computes point estimates, EIFs, and analytical SEs for
#' each cell.
#'
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param pt_assumption character: \code{"all"} or \code{"post"}
#' @param alpha significance level in (0, 1)
#' @param store_eif logical: if TRUE, include EIF vectors in returned cells
#' @param xformla one-sided formula or NULL: covariate formula (routed to
#'   covariate path when non-trivial and \code{panel_obj$covariate_matrix}
#'   is non-NULL)
#' @param need_eif logical: if TRUE, always store EIF regardless of store_eif
#'   (used internally when \code{n_bootstrap > 0})
#'
#' @return list with elements:
#'   \describe{
#'     \item{\code{cells}}{list of \code{edid_cell_result} objects}
#'     \item{\code{eif_matrix}}{n x n_valid_cells numeric matrix, or NULL}
#'     \item{\code{cell_index}}{data.frame: group, time, cell_id, is_pre}
#'   }
#' @keywords internal
fit_edid_cells <- function(
  panel_obj, pt_assumption, alpha, store_eif, xformla = NULL, seed = NULL,
  need_eif = FALSE, weight_method = c("efficient", "averaged", "gmm", "uniform"),
  estimation_effect = FALSE
) {
  weight_method <- match.arg(weight_method)
  # Determine if covariate path is active
  is_trivial_xformla <- is.null(xformla) ||
    identical(deparse(xformla, width.cutoff = 500L), "~1")
  use_cov_path <- !is_trivial_xformla && !is.null(panel_obj$covariate_matrix)

  # Curse of dimensionality guard for the pointwise efficient weights. The kernel
  # Omega*(X) is consistent only while its local effective sample size n * prod(h_k) ~
  # n^{(5-d)/5} grows, i.e. d < 5. At d >= 5 the admissible floor band (0, (5-d)/10) is
  # empty, Omega_hat*(X) is not consistent, and the regularized weights collapse toward
  # uniform (the estimator stays CONSISTENT via the valid DR moments, but is no longer
  # semiparametrically efficient). Recommend the constant-weight schemes, whose pooled
  # Omega-bar is sqrt(n)-consistent with no curse of dimensionality.
  if (weight_method == "efficient" && use_cov_path) {
    d_cov <- ncol(panel_obj$covariate_matrix)
    if (d_cov >= 5L) {
      warning(sprintf(paste0(
        "weights='efficient' with %d continuous covariates: the pointwise kernel ",
        "Omega*(X) is not consistently estimable (curse of dimensionality; admissible ",
        "floor band is empty for d>=5), so the efficient weights collapse toward uniform ",
        "and the efficiency guarantee is lost (ATT stays consistent). Use ",
        "weights='averaged' (sqrt(n)-consistent, no curse) or condition on a low-",
        "dimensional index (e.g. the propensity score)."), d_cov), call. = FALSE)
    }
  }

  # The "gmm" scheme inverts the unconditional covariance of the generated outcomes, which is
  # estimated from the same data as the moments. This induces a finite-sample two-step bias of
  # order O(1/n) that can be sizeable with strong treatment effects, many moments, or small
  # cohorts. Asymptotically it coincides with "averaged" and is never more efficient, so the
  # constant-weight default "averaged" (or the asymptotically efficient "efficient") is preferred.
  if (weight_method == "gmm") {
    warning(paste0(
      "weights = 'gmm' inverts the unconditional moment covariance and carries a finite-sample ",
      "O(1/n) bias that can be large under strong effects, many periods, or small cohorts; it is ",
      "asymptotically equal to 'averaged' but never more efficient. Prefer 'averaged' or 'efficient'."),
      call. = FALSE)
  }

  # Nuisance functions are estimated by plug-in (train = test = full sample), as in the paper's
  # remark that the efficient-influence-function moments are Neyman orthogonal, so the first-step
  # nuisance estimates do not affect the first-order asymptotic variance. Cross-fitting (K > 1) is
  # not used: the estimated weights and nuisances are already first-order negligible, while
  # sample-splitting can inflate the finite-sample variance of the just-identified DR estimator.
  K_use <- 1L
  if (use_cov_path) {
    fold_id <- rep(1L, panel_obj$n)
  } else {
    fold_id <- NULL
  }

  # The ACH (Ackerberg, Chen & Hahn 2012) first-step correction requires the plug-in
  # (train = test = full sample) M-estimator pieces; it is derived for K = 1. Guard
  # defensively in case cross-fitting is ever enabled upstream.
  if (isTRUE(estimation_effect)) {
    if (!use_cov_path) {
      warning("estimation_effect has no effect without covariates (no first-step nuisances).", call. = FALSE)
      estimation_effect <- FALSE
    } else if (K_use > 1L) {
      stop("estimation_effect = TRUE is only supported with plug-in nuisances (K = 1).", call. = FALSE)
    }
  }

  tgroups   <- panel_obj$treatment_groups
  tperiods  <- panel_obj$time_periods
  period_1  <- panel_obj$period_1
  n         <- panel_obj$n

  # Periods to iterate over: all except period_1
  iter_periods <- tperiods[tperiods != period_1]

  # Pre-allocate cell list
  n_cells <- length(tgroups) * length(iter_periods)
  cells   <- vector("list", n_cells)

  # cell_index tracking
  ci_group  <- numeric(n_cells)
  ci_time   <- numeric(n_cells)
  ci_cell_id <- integer(n_cells)
  ci_is_pre  <- logical(n_cells)

  keep_eif <- store_eif || need_eif
  eif_list <- if (keep_eif) vector("list", n_cells) else NULL

  cell_id <- 0L
  n_extreme_ratio_instances <- 0L  # accumulate extreme-ratio warnings; emit once at end

  for (g in tgroups) {
    for (t in iter_periods) {
      cell_id <- cell_id + 1L
      is_pre  <- (t < g)

      # Step 1: enumerate valid pairs
      pairs <- enumerate_valid_pairs_edid(
        target_g          = g,
        treatment_groups  = tgroups,
        time_periods      = tperiods,
        period_1          = period_1,
        pt_assumption     = pt_assumption,
        anticipation      = panel_obj$anticipation
      )

      # NA cell if no valid pairs
      if (nrow(pairs) == 0L) {
        cells[[cell_id]] <- list(
          group           = g,
          time            = t,
          att             = NA_real_,
          se              = NA_real_,
          ci_lower        = NA_real_,
          ci_upper        = NA_real_,
          t_stat          = NA_real_,
          p_value         = NA_real_,
          n_pairs         = 0L,
          weights         = NULL,
          condition_num   = NA_real_,
          is_pre          = is_pre,
          inference_valid = FALSE,
          eif             = NULL
        )
        ci_group[cell_id]   <- g
        ci_time[cell_id]    <- t
        ci_cell_id[cell_id] <- cell_id
        ci_is_pre[cell_id]  <- is_pre
        if (keep_eif) eif_list[[cell_id]] <- rep(NA_real_, n)
        next
      }

      # Covariate nuisance estimates (per cell)
      prop_ratios  <- NULL
      cond_means   <- NULL
      inv_propensities <- NULL

      if (use_cov_path) {
        # Build nuisance estimation pairs.
        # For Eq. (4.4), we need nuisances for:
        #   - r[g, Inf] always (never-treated propensity ratio)
        #   - r[g, gp] for each cross-cohort gp
        #   - m_{Inf, t} and m_{Inf, tpre} for never-treated conditional means
        #   - m_{gp, tpre} for each cross-cohort gp's pretrend
        # Build pairs_for_nuisance that includes Inf + all cross-cohort gps
        pairs_for_nuisance <- pairs
        # Self-comparison pairs use Inf as comparison
        self_cmp <- is.finite(pairs_for_nuisance$gp) & (pairs_for_nuisance$gp == g)
        if (any(self_cmp)) {
          pairs_for_nuisance$gp[self_cmp] <- Inf
        }
        # Ensure Inf is always present for never-treated nuisances
        # (needed even for cross-cohort pairs)
        cross_pairs <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
        if (nrow(cross_pairs) > 0L) {
          # Add Inf-based pairs for the same tpre values (for m_{Inf,tpre})
          inf_pairs <- data.frame(gp = Inf, tpre = unique(cross_pairs$tpre))
          pairs_for_nuisance <- rbind(pairs_for_nuisance, inf_pairs)
          pairs_for_nuisance <- unique(pairs_for_nuisance)
        }

        prop_ratios <- withCallingHandlers(
          estimate_all_propensity_ratios(
            panel_obj = panel_obj,
            g         = g,
            pairs     = pairs_for_nuisance,
            bs_df    = 4L,
            K_folds  = K_use,
            fold_id  = fold_id,
            return_aux = estimation_effect
          ),
          warning = function(w) {
            if (grepl("Extreme propensity ratios", conditionMessage(w), fixed = TRUE)) {
              n_extreme_ratio_instances <<- n_extreme_ratio_instances + 1L
              invokeRestart("muffleWarning")
            }
          }
        )
        cond_means <- estimate_all_conditional_means(
          panel_obj = panel_obj,
          pairs     = pairs_for_nuisance,
          t_val    = t,
          bs_df    = 4L,
          K_folds  = K_use,
          fold_id  = fold_id,
          return_aux = estimation_effect
        )
        # Split predictions from the ACH aux pieces (the *_aux return path wraps both).
        r_aux <- NULL; m_aux <- NULL
        if (isTRUE(estimation_effect)) {
          r_aux <- prop_ratios$aux; prop_ratios <- prop_ratios$predictions
          m_aux <- cond_means$aux;  cond_means  <- cond_means$predictions
        }
        inv_propensities <- estimate_all_inverse_propensities(
          panel_obj = panel_obj,
          g         = g,
          pairs     = pairs,
          bs_df     = 4L,
          K_folds   = K_use,
          fold_id   = fold_id
        )
      }

      # Steps 2-6: dispatch on covariate vs. no-covariate path
      if (use_cov_path) {
        # --- Covariate path ---
        gen_out_mat <- compute_generated_outcomes_cov_edid(
          panel_obj     = panel_obj,
          g             = g,
          t             = t,
          pairs         = pairs,
          prop_ratios   = prop_ratios,
          cond_means    = cond_means,
          pt_assumption = pt_assumption
        )
        if (isTRUE(getOption("edid_store_genout"))) {            # diagnostic: expose per-cell generated outcomes
          acc <- getOption("edid_genout_acc", list()); acc[[paste0(g, "_", t)]] <- gen_out_mat
          options(edid_genout_acc = acc)
        }
        if (weight_method == "efficient") {
          # Paper's pointwise efficient weights w(X_i)=Omega*(X_i)^{-1}1/(1'Omega*(X_i)^{-1}1),
          # with a dimension-aware eigenvalue-floor regularization (a=0.7*(5-d)/10) that is
          # asymptotically negligible yet dominates the NW estimation noise for stability.
          omega_arr <- compute_omega_star_cov_edid(panel_obj, g, t, pairs,
                                                   prop_ratios, cond_means,
                                                   inv_propensities, return_pointwise = TRUE)
          W_pw     <- compute_pointwise_weights_edid(omega_arr,
                                                     d = ncol(panel_obj$covariate_matrix))  # n x H
          cond_num <- tryCatch(check_condition_edid(apply(omega_arr, c(2, 3), mean)),
                               error = function(e) NA_real_)
          wY_i     <- rowSums(gen_out_mat * W_pw)
          if (anyNA(wY_i))                                                    # NA would split the support:
            stop(sprintf(paste0("edid cell (%s,%s): NA in weighted generated outcomes (typically a ",
              "non-invertible local covariance or a missing covariate prediction). The point estimate ",
              "would use complete cases while the EIF spans all units, breaking the empirical mean-zero ",
              "identity and invalidating the SE. edid requires complete cases: drop incomplete units ",
              "before calling, or inspect cell (%s,%s)."), g, t, g, t), call. = FALSE)
          att_gt   <- mean(wY_i, na.rm = TRUE)                               # E_n[w(X_i)' Ytilde_i]
          eif_gt   <- compute_eif_cov_edid(panel_obj, gen_out_mat, W_pw, att_gt, g)
          if (isTRUE(estimation_effect))                                    # ACH first-step correction (W frozen)
            eif_gt <- eif_gt - compute_ach_correction_cov_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, W_pw, m_aux, r_aux, pt_assumption)
          weights  <- colMeans(W_pw, na.rm = TRUE)                            # store mean weight per pair
          # Diagnostic only (default off): per-component cross-unit SD of the pointwise weights
          # quantifies how much Omega*(X) shape-varies; ~0 => weights ~constant => efficient ~ averaged.
          if (isTRUE(getOption("edid_diag_wpw")))
            message(sprintf("WPWDIAG %d_%d wsd=[%s] meanw=[%s] maxw=%.3f",
              g, t, paste(round(apply(W_pw, 2, stats::sd, na.rm = TRUE), 4), collapse = ","),
              paste(round(colMeans(W_pw, na.rm = TRUE), 3), collapse = ","), max(abs(W_pw), na.rm = TRUE)))
        } else {
          # Constant-weight schemes (valid but not pointwise-efficient):
          #   averaged = invert kernel Omega-bar; gmm = invert unconditional S_hat; uniform = 1/H.
          omega    <- compute_omega_star_cov_edid(panel_obj, g, t, pairs,
                                                  prop_ratios, cond_means,
                                                  inv_propensities)
          cond_num <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)
          H_local  <- nrow(omega)
          weights  <- switch(weight_method,
            uniform = rep(1 / H_local, H_local),
            # pairwise.complete.obs so a single NA moment does not NA-poison the whole
            # covariance (compute_efficient_weights_edid then guards any residual non-finite).
            gmm     = compute_efficient_weights_edid(stats::cov(gen_out_mat, use = "pairwise.complete.obs")),
            compute_efficient_weights_edid(omega))   # "averaged"
          att_gt   <- sum(weights * colMeans(gen_out_mat, na.rm = TRUE))
          eif_gt   <- compute_eif_cov_edid(panel_obj, gen_out_mat, weights, att_gt, g)
          if (isTRUE(estimation_effect))                                    # ACH first-step correction (w frozen)
            eif_gt <- eif_gt - compute_ach_correction_cov_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, weights, m_aux, r_aux, pt_assumption)
        }
      } else {
        # --- No-covariate path ---
        y_hat    <- compute_generated_outcomes_nocov_edid(g, t, pairs, panel_obj, pt_assumption)
        omega    <- compute_omega_star_nocov_edid(g, t, pairs, panel_obj, pt_assumption)
        cond_num <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)
        # Honor weight_method: with no covariates there is no X-variation, so efficient,
        # averaged, and gmm all invert the same unconditional Omega* and coincide; only
        # uniform differs. (Previously this path silently ignored weight_method.)
        H_local  <- length(y_hat)
        weights  <- if (weight_method == "uniform") rep(1 / H_local, H_local) else compute_efficient_weights_edid(omega)
        att_gt   <- sum(weights * y_hat)
        eif_gt   <- compute_eif_nocov_edid(g, t, pairs, weights, panel_obj, att_gt, pt_assumption)
      }

      # Step 7: SE and inference
      inf_res <- safe_inference_edid(eif_gt, panel_obj$cluster_indices, alpha, att_gt)

      # Step 8: store
      eif_store <- if (keep_eif) eif_gt else NULL
      cells[[cell_id]] <- list(
        group           = g,
        time            = t,
        att             = att_gt,
        se              = inf_res$se,
        ci_lower        = inf_res$ci_lower,
        ci_upper        = inf_res$ci_upper,
        t_stat          = inf_res$t_stat,
        p_value         = inf_res$p_value,
        n_pairs         = nrow(pairs),
        weights         = weights,
        condition_num   = cond_num,
        is_pre          = is_pre,
        inference_valid = inf_res$inference_valid,
        eif             = eif_store
      )

      ci_group[cell_id]   <- g
      ci_time[cell_id]    <- t
      ci_cell_id[cell_id] <- cell_id
      ci_is_pre[cell_id]  <- is_pre

      if (keep_eif) eif_list[[cell_id]] <- eif_gt
    }
  }

  if (n_extreme_ratio_instances > 0L) {
    warning(sprintf(
      "Extreme propensity ratios detected (max > 100) in %d estimation step(s). Results may be unstable.",
      n_extreme_ratio_instances
    ))
  }

  # Build EIF matrix if needed
  eif_matrix <- NULL
  if (keep_eif) {
    # Stack EIF vectors as columns: n x n_cells matrix
    eif_matrix <- do.call(cbind, eif_list)
    if (is.null(dim(eif_matrix))) {
      eif_matrix <- matrix(eif_matrix, nrow = n)
    }
  }

  cell_index <- data.frame(
    group   = ci_group,
    time    = ci_time,
    cell_id = ci_cell_id,
    is_pre  = ci_is_pre,
    stringsAsFactors = FALSE
  )

  # Diagnostic: if every post-treatment cell is NA (no admissible pairs), warn instead
  # of silently returning NA. A common cause under pt_assumption='post' is a
  # non-consecutive/irregular time grid where the literal pre-period (g-1-anticipation)
  # is unobserved, so every (g,t) cell yields zero pairs.
  post_idx <- which(!ci_is_pre)
  if (length(post_idx)) {
    post_atts <- vapply(post_idx, function(k) {
      a <- cells[[ci_cell_id[k]]]$att; if (is.null(a)) NA_real_ else a
    }, numeric(1L))
    if (all(!is.finite(post_atts))) {
      warning(sprintf(paste0(
        "All post-treatment ATT(g,t) cells are NA (no admissible pairs). Under ",
        "pt_assumption='%s' the comparison base period is the period immediately before ",
        "treatment; on a non-consecutive/irregular time grid that period may be ",
        "unobserved. Check time spacing and anticipation."), pt_assumption), call. = FALSE)
    }
  }

  list(
    cells      = cells,
    eif_matrix = eif_matrix,
    cell_index = cell_index
  )
}
