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
  need_eif = FALSE
) {
  # Determine if covariate path is active
  is_trivial_xformla <- is.null(xformla) ||
    identical(deparse(xformla, width.cutoff = 500L), "~1")
  use_cov_path <- !is_trivial_xformla && !is.null(panel_obj$covariate_matrix)

  # Nuisance estimation uses the plug-in (no sample-splitting) approach
  # described in the paper's main text (Sec. 5.2, footnote 601), with
  # train = test = full sample. Cross-fitting (K>1) is intentionally
  # disabled in this release: benchmark experiments showed that K=5
  # cross-fitting substantially inflated the point estimator's variance
  # at moderate sample sizes (~35% larger mc_sd vs plug-in at n=500),
  # erasing the efficiency gain over the just-identified DR estimator.
  # The cross-fitting machinery in estimate_all_* (the K_folds>1 branch)
  # is preserved for future re-enablement once the calibration of
  # cross-fitted standard errors is better understood. See
  # benchmark/edid_cov_fix_verification.R for the supporting evidence.
  K_use <- 1L
  if (use_cov_path) {
    fold_id <- rep(1L, panel_obj$n)
  } else {
    fold_id <- NULL
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
            fold_id  = fold_id
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
          fold_id  = fold_id
        )
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
        omega   <- compute_omega_star_cov_edid(panel_obj, g, t, pairs,
                                                prop_ratios, cond_means,
                                                inv_propensities)
        cond_num <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)
        weights <- compute_efficient_weights_edid(omega)
        # ATT: weighted mean of column means of gen_out_mat
        col_means <- colMeans(gen_out_mat, na.rm = TRUE)
        att_gt    <- sum(weights * col_means)
        eif_gt    <- compute_eif_cov_edid(panel_obj, gen_out_mat, weights, att_gt, g)
      } else {
        # --- No-covariate path ---
        y_hat    <- compute_generated_outcomes_nocov_edid(g, t, pairs, panel_obj, pt_assumption)
        omega    <- compute_omega_star_nocov_edid(g, t, pairs, panel_obj, pt_assumption)
        cond_num <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)
        weights  <- compute_efficient_weights_edid(omega)
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

  list(
    cells      = cells,
    eif_matrix = eif_matrix,
    cell_index = cell_index
  )
}
