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
#' @param covariates NULL (covariate path is a stub)
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
  panel_obj, pt_assumption, alpha, store_eif, covariates,
  need_eif = FALSE
) {
  if (!is.null(covariates)) {
    stop("covariate path not yet implemented in edid(). ",
         "Pass covariates = NULL or omit the covariates argument.")
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

      # Step 2: generated outcomes
      y_hat <- compute_generated_outcomes_nocov_edid(g, t, pairs, panel_obj, pt_assumption)

      # Step 3: Omega*
      omega <- compute_omega_star_nocov_edid(g, t, pairs, panel_obj, pt_assumption)

      # Condition number diagnostics
      cond_num <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)

      # Step 4: efficient weights
      weights <- compute_efficient_weights_edid(omega)

      # Step 5: point estimate
      att_gt <- sum(weights * y_hat)

      # Step 6: EIF
      eif_gt <- compute_eif_nocov_edid(g, t, pairs, weights, panel_obj, att_gt, pt_assumption)

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
