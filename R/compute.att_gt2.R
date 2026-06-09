#' @title Get the cohort for the current g,t values
#' @description A utility function to get vector with index of units that will be included in the DiD estimation for the current g,t pair.
#'
#' @param group Group.
#' @param time Time period.
#' @param tfac Time factor, which is 1 for varying base period and 0 for universal base period.
#' @param pret Pre treatment period.
#' @param dp2 DiD parameters v2.0.
#'
#' @return A vector of 1s, 0s and NAs indicating the cohort for the current g,t pair.
#' @noRd
get_did_cohort_index <- function(group, time, tfac, pret, dp2){
  # return a vector of dimension id_size with 1, 0 or NA values
  time_periods <- dp2$time_periods

  if(dp2$panel){
    # --- Balanced panel: use pre-computed cumulative sizes ---
    cohort_vec <- dp2$cohort_counts$cohort
    cum_sizes <- dp2$.cohort_cum_sizes  # pre-computed in compute.att_gt2

    # find the first cohort after the cutoff for not-yet-treated
    if (dp2$control_group == "notyettreated") {
      cutoff <- time_periods[max(time, pret) + tfac] + dp2$anticipation
      min_idx <- match(TRUE, cohort_vec > cutoff)
      min_control_group <- if (is.na(min_idx)) Inf else cohort_vec[min_idx]
    } else {
      min_control_group <- Inf
    }

    # control group boundaries
    ctrl_start_idx <- match(TRUE, cohort_vec >= min_control_group)
    if (is.na(ctrl_start_idx)) {
      start_control <- dp2$id_count + 1L  # no control group found
    } else {
      start_control <- if (ctrl_start_idx == 1L) 1L else cum_sizes[ctrl_start_idx - 1L] + 1L
    }
    end_control <- cum_sizes[length(cum_sizes)]

    # treated group boundaries
    treat_idx <- match(dp2$treated_groups[group], cohort_vec)
    if (is.na(treat_idx)) {
      stop("Internal error: treated group ", dp2$treated_groups[group],
           " not found in cohort_vec. Please report this issue.")
    }
    start_treat <- if (treat_idx == 1L) 1L else cum_sizes[treat_idx - 1L] + 1L
    end_treat <- cum_sizes[treat_idx]

    # fill the cohort index vector
    did_cohort_index <- rep(NA_integer_, dp2$id_count)
    if (start_control <= end_control) did_cohort_index[start_control:end_control] <- 0L
    did_cohort_index[start_treat:end_treat] <- 1L

  } else {
    # getting the index to get units who will participate in the estimation for the (g,t) cell.
    # Note: This works because the data is already ordered in a specific way. Changing that order will break this.

    dat <- dp2$time_invariant_data          # shortcut

    # flags for treated and control units
    Gflag <- dat[[dp2$gname]] == dp2$treated_groups[group]

    if (dp2$control_group == "nevertreated") {
      Cflag <- dat[[dp2$gname]] == Inf
    } else {  # not-yet-treated
      Cflag <- (dat[[dp2$gname]] == Inf) |
        (dat[[dp2$gname]] > dp2$time_periods[max(time, pret) + tfac] + dp2$anticipation &
           dat[[dp2$gname]] != dp2$treated_groups[group])
    }

    # keep only rows observed in pret or t
    keep <- dat[[dp2$tname]] %in% c(dp2$time_periods[time + tfac], dp2$time_periods[pret])

    did_cohort_index <- rep(NA_integer_, nrow(dat))
    did_cohort_index[keep & Cflag] <- 0L
    did_cohort_index[keep & Gflag] <- 1L

  }

  return(did_cohort_index)
}

#' @title Wrapper to run DRDID package
#' @description A utility function to run the DRDID package for the current g,t pair.
#'
#' @param cohort_data Data table with the cohort data for the current g,t pair
#' @param covariates Matrix of covariates to be used in the estimation
#' @param dp2 DiD parameters v2.0.
#'
#' @return A list containing the estimated ATT and the influence function vector.
#' @noRd
run_DRDID <- function(cohort_data, covariates, dp2, g_val = NULL, t_val = NULL, force_rc = FALSE){

  extra_args <- if (is.null(dp2$extra_args)) list() else dp2$extra_args
  gt_label <- if (!is.null(g_val) && !is.null(t_val)) paste0(" for group ", g_val, " in time period ", t_val) else ""
  # whether to compute influence functions (default TRUE; FALSE = point estimates only)
  do_inf <- is.null(dp2$compute_inffunc) || isTRUE(dp2$compute_inffunc)

  if(dp2$panel && !force_rc){
    # --------------------------------------
    # Panel Data
    # --------------------------------------
    # cohort_data is a plain named list of vectors (D, y1, y0, i.weights), built in
    # run_att_gt_estimation(). Working on the vectors directly avoids per-cell
    # data.table construction/extraction; the estimator inputs are identical.

    # still total number of units (not just included in G or C)
    D_all <- cohort_data$D
    n <- length(D_all)

    # pick up the indices for units that will be used to compute ATT(g,t)
    valid_obs <- which(!is.na(D_all))
    D    <- D_all[valid_obs]
    y1   <- cohort_data$y1[valid_obs]
    y0   <- cohort_data$y0[valid_obs]
    i.weights <- cohort_data$i.weights[valid_obs]

    if(dp2$xformla != ~1){
      covariates <- covariates[valid_obs,]
    } else {
      covariates <- rep(1, length(valid_obs))
    }
    covariates <- as.matrix(covariates)

    # num obs. for computing ATT(g,t)
    n1 <- length(D)

    #-----------------------------------------------------------------------------
    # check for overlap and regression problems
    #-----------------------------------------------------------------------------
    custom_est_method <- inherits(dp2$est_method, "function")

    if (!custom_est_method) {
      D_vec <- D

      # checks for pscore based methods
      if (dp2$est_method %in% c("dr", "ipw")) {
        preliminary_logit <- overlap_logit_fit(covariates, D_vec)
        preliminary_pscores <- preliminary_logit$fitted.values
        if (max(preliminary_pscores) >= 0.999) {
          warning(paste0("overlap condition violated", gt_label))
          return(list(att = NA, inf_func = if (do_inf) rep(NA_real_, n) else NULL))
        }
      }

      # check if can run regression using control units
      if (dp2$est_method %in% c("dr", "reg")) {
        control_covs <- covariates[D_vec == 0, , drop = FALSE]
        if (rcond(t(control_covs) %*% control_covs) < .Machine$double.eps) {
          warning(paste0("Not enough control units", gt_label, " to run specified regression"))
          return(list(att = NA, inf_func = if (do_inf) rep(NA_real_, n) else NULL))
        }
      }
    }

    #-----------------------------------------------------------------------------
    # code for actually computing ATT(g,t)
    #-----------------------------------------------------------------------------


    if (inherits(dp2$est_method, "function")) {
      # user-specified function
      attgt <- do.call(dp2$est_method, c(list(
                          y1=y1,
                          y0=y0,
                          D=D,
                          covariates=covariates,
                          i.weights=i.weights,
                          inffunc=do_inf), extra_args))
    } else if (dp2$est_method == "ipw") {
      # inverse-probability weights
      attgt <- std_ipw_did_panel(y1=y1,
                                        y0=y0,
                                        D=D,
                                        covariates=covariates,
                                        i.weights=i.weights,
                                        boot=FALSE, inffunc=do_inf)
    } else if (dp2$est_method == "reg") {
      # regression
      attgt <- reg_did_panel(y1=y1,
                                    y0=y0,
                                    D=D,
                                    covariates=covariates,
                                    i.weights=i.weights,
                                    boot=FALSE, inffunc=do_inf)
    } else {
      # doubly robust, this is default
      attgt <- drdid_panel(y1=y1,
                                  y0=y0,
                                  D=D,
                                  covariates=covariates,
                                  i.weights=i.weights,
                                  boot=FALSE, inffunc=do_inf)
    }

    # adjust influence function to account for only using
    # subgroup to estimate att(g,t)
    if (do_inf) {
      inf_func_vector <- rep(0, n)
      inf_func_not_na <- (n/n1)*attgt$att.inf.func
      inf_func_vector[valid_obs] <- inf_func_not_na
    } else {
      inf_func_vector <- NULL   # point estimates only
    }

  } else {
    # --------------------------------------
    # Repeated Cross-Section
    # --------------------------------------
    # cohort_data is a plain named list of vectors (D, y, post, i.weights, .rowid),
    # built in run_att_gt_estimation(). Working on the vectors directly avoids per-cell
    # data.table construction/extraction; the estimator inputs are identical.
    D_all <- cohort_data$D
    rowid_all <- cohort_data$.rowid   # all rows (needed for the unbalanced aggregation)

    # still total number of obs. (not just included in G or C)
    n <- length(D_all)

    # pick up the indices for units that will be used to compute ATT(g,t)
    valid_obs <- which(!is.na(D_all))
    D    <- D_all[valid_obs]
    y    <- cohort_data$y[valid_obs]
    post <- cohort_data$post[valid_obs]
    i.weights <- cohort_data$i.weights[valid_obs]

    # num obs. for computing ATT(g,t)
    n1 <- length(D)

    if(dp2$xformla != ~1){
      covariates <- covariates[valid_obs,]
    } else {
      covariates <- covariates[valid_obs]
    }

    covariates <- as.matrix(covariates)

    #-----------------------------------------------------------------------------
    # check for overlap and regression problems
    #-----------------------------------------------------------------------------
    custom_est_method <- inherits(dp2$est_method, "function")

    if (!custom_est_method) {
      D_vec <- D

      # determine correct inf_func length for early returns
      n_inf <- if (dp2$allow_unbalanced_panel) dp2$id_count else n

      # checks for pscore based methods
      if (dp2$est_method %in% c("dr", "ipw")) {
        preliminary_logit <- overlap_logit_fit(covariates, D_vec)
        preliminary_pscores <- preliminary_logit$fitted.values
        if (max(preliminary_pscores) >= 0.999) {
          warning(paste0("overlap condition violated", gt_label))
          return(list(att = NA, inf_func = if (do_inf) rep(NA_real_, n_inf) else NULL))
        }
      }

      # check if can run regression using control units
      if (dp2$est_method %in% c("dr", "reg")) {
        control_covs <- covariates[D_vec == 0, , drop = FALSE]
        if (rcond(t(control_covs) %*% control_covs) < .Machine$double.eps) {
          warning(paste0("Not enough control units", gt_label, " to run specified regression"))
          return(list(att = NA, inf_func = if (do_inf) rep(NA_real_, n_inf) else NULL))
        }
      }
    }

    #-----------------------------------------------------------------------------
    # code for actually computing ATT(g,t)
    #-----------------------------------------------------------------------------

    if (inherits(dp2$est_method, "function")) {
      # user-specified function
      attgt <- do.call(dp2$est_method, c(list(
                          y=y,
                          post=post,
                          D=D,
                          covariates=covariates,
                          i.weights=i.weights,
                          inffunc=do_inf), extra_args))
    } else if (dp2$est_method == "ipw") {
      # inverse-probability weights
      attgt <- std_ipw_did_rc(y=y,
                             post=post,
                             D=D,
                             covariates=covariates,
                             i.weights=i.weights,
                             boot=FALSE, inffunc=do_inf)
    } else if (dp2$est_method == "reg") {
      # regression
      attgt <- reg_did_rc(y=y,
                         post=post,
                         D=D,
                         covariates=covariates,
                         i.weights=i.weights,
                         boot=FALSE, inffunc=do_inf)
    } else {
      # doubly robust, this is default
      attgt <- drdid_rc(y=y,
                       post=post,
                       D=D,
                       covariates=covariates,
                       i.weights=i.weights,
                       boot=FALSE, inffunc=do_inf)
    }

    # n/n1 adjusts for estimating the
    # att_gt only using observations from groups
    # G and C
    # adjust influence function to account for only using
    # subgroup to estimate att(g,t)
    if (!do_inf) {
      inf_func_vector <- NULL   # point estimates only
    } else if(dp2$allow_unbalanced_panel){
      # since this is technically a panel data but ran as RCS, we need to adjust the
      # influence function by aggregating influence value by .rowid (several obs of one
      # unit can be used to estimate ATT in each 2x2). rowsum(..., reorder = FALSE)
      # reproduces the data.table `by = .rowid` first-appearance group order exactly.
      inf_func_long <- numeric(n)
      inf_func_long[valid_obs] <- (dp2$id_count / n1) * attgt$att.inf.func
      inf_func_vector <- rowsum(inf_func_long, rowid_all, reorder = FALSE)[, 1L]

    } else {
      inf_func_vector <- rep(0, n)
      inf_func_not_na <- (n/n1)*attgt$att.inf.func
      inf_func_vector[valid_obs] <- inf_func_not_na
    }

  }

  return(list(att = attgt$ATT, inf_func = inf_func_vector))

}


#' @title Run ATT estimation for a given group-time pair
#'
#' @description `run_att_gt_estimation` does the main work for computing
#'  multiperiod group-time average treatment effects
#' @param g group of interest (treated group at time t)
#' @param t time period
#' @param dp2 A DIDparams object v2.0
#'
#' @return a list with the gt cell and the results after performing estimation
#'
#' @keywords internal
run_att_gt_estimation <- function(g, t, dp2){

  if(dp2$print_details){cat("\n", paste0("Evaluating (g,t) = (",dp2$treated_groups[g],",",dp2$time_periods[t],")"))}
  tfac <- if (dp2$base_period != "universal") 1L else 0L
  # set pret
  # varying base period
  pret <- t

  # Use pre-computed pret for universal base period and post-treatment periods
  pret_g <- dp2$.pret_by_group[g]  # pre-computed in compute.att_gt2

  # universal base period
  if (dp2$base_period == "universal") {
    pret <- pret_g
  }

  # check if in post-treatment period
  if ((dp2$treated_groups[g] <= dp2$time_periods[(t+tfac)])) {

    # update pre-period if in post-treatment period to
    # be  period (g-delta-1)
    pret <- pret_g

    # print a warning message if there are no pre-treatment periods
    if (is.na(pret)) {
      warning(paste0("There are no pre-treatment periods for the group first treated at ", dp2$treated_groups[g], "\nUnits from this group are dropped"))

      # if there are no pre-treatment periods, code will
      # jump out of this loop
      return(NULL)
    }
  }

  #-----------------------------------------------------------------------------
  # if we are in period (g-1), normalize results to be equal to 0
  # and break without computing anything
  if (dp2$base_period == "universal") {
    if (dp2$time_periods[pret] == dp2$time_periods[(t+tfac)]) {
      return(list(base_period_norm = TRUE))
    }
  }

  # get units in treatment and control group
  did_cohort_index <- get_did_cohort_index(group = g, time = t, tfac = tfac, pret = pret, dp2)
  # In case of no treatment or control group in the cohort, return NULL
  valid_did_cohort <- any(did_cohort_index == 1) & any(did_cohort_index == 0)
  if(!isTRUE(valid_did_cohort)){
    if(dp2$print_details){cat("\n Skipping (g,t) as no treatment group or control group found")}
    return(NULL)
  }



  if(dp2$panel){
    # Determine which weight period to use based on fix_weights
    use_rc_for_weights <- (!is.null(dp2$fix_weights) && dp2$fix_weights == "varying")

    if (use_rc_for_weights) {
      # fix_weights = "varying": stack into RC format with per-period weights
      n_units <- length(did_cohort_index)
      cohort_data <- data.table(
        D = rep(did_cohort_index, 2),
        y = c(dp2$outcomes_tensor[[pret]], dp2$outcomes_tensor[[t+tfac]]),
        post = rep(c(0L, 1L), each = n_units),
        i.weights = c(dp2$weights_tensor[[pret]], dp2$weights_tensor[[t+tfac]])
      )
      # Use earlier-period covariates for both halves — fix_weights only
      # changes weights, not the covariate conditioning set.
      # Use min(pret, t) to match the panel estimator's convention:
      # with base_period="universal", pret can be later than t for placebo cells.
      cov_early <- dp2$covariates_tensor[[base::min(pret, t)]]
      if (is.matrix(cov_early)) {
        covariates <- rbind(cov_early, cov_early)
      } else {
        covariates <- c(cov_early, cov_early)
      }
    } else {
      # Default or fixed weight options: use panel estimator with single weight vector
      if (is.null(dp2$fix_weights)) {
        # Default: weight from earlier of the two periods
        w_idx <- base::min(pret, t)
      } else if (dp2$fix_weights == "base_period") {
        w_idx <- dp2$.pret_by_group[g]
      } else if (dp2$fix_weights == "first_period") {
        w_idx <- 1L
      }
      # Plain named list of the cohort vectors (the tensors already hold vectors);
      # this avoids constructing/extracting a data.table per (g,t) cell, which the
      # profiler showed is a large share of the run time. run_DRDID()'s panel branch
      # consumes these as vectors, so the result is identical.
      cohort_data <- list(D = did_cohort_index,
                          y1 = dp2$outcomes_tensor[[t + tfac]],
                          y0 = dp2$outcomes_tensor[[pret]],
                          i.weights = dp2$weights_tensor[[w_idx]])
      covariates <- dp2$covariates_tensor[[base::min(pret, t)]]
    }
  } else {

    log_vec <- dp2$time_invariant_data[[ dp2$tname ]] == dp2$time_periods[t+tfac]
    # convert TRUE/FALSE to 1/0 in place (fastest)
    set(dp2$time_invariant_data, j = "post", value = as.integer(log_vec))

    # Handle fix_weights for RC/unbalanced panel
    if (!is.null(dp2$fix_weights) && dp2$fix_weights %in% c("base_period", "first_period")) {
      if (dp2$fix_weights == "base_period") {
        target_period <- dp2$time_periods[dp2$.pret_by_group[g]]
      } else {
        target_period <- dp2$time_periods[1]
      }
      # Build weight lookup from target period
      tid <- dp2$time_invariant_data
      target_mask <- tid[[dp2$tname]] == target_period
      target_ids <- tid[[dp2$idname]][target_mask]
      target_ws <- tid[["weights"]][target_mask]
      # Map each observation's weight from the target period by integer id matching.
      # Equivalent to the previous named-character lookup (first match wins on any
      # duplicate id; unmatched ids -> NA) but avoids coercing every id to character.
      fixed_w <- target_ws[match(tid[[dp2$idname]], target_ids)]
      # Exclude units not observed in target period by setting D to NA
      # (run_DRDID filters on !is.na(D))
      na_w <- is.na(fixed_w)
      if (any(na_w & !is.na(did_cohort_index))) {
        did_cohort_index[na_w] <- NA_integer_
        warning(paste0("Some units not observed in ", dp2$fix_weights,
                       " (period ", target_period, ") for group ",
                       dp2$treated_groups[g], " in time period ",
                       dp2$time_periods[t+tfac], ". These units are excluded."))
      }
      # Plain named list of the cohort vectors (avoids per-cell data.table
      # construction/extraction; run_DRDID()'s RC branch consumes these as vectors).
      cohort_data <- list(D = did_cohort_index, y = tid[[dp2$yname]], post = tid$post,
                          i.weights = fixed_w, .rowid = tid$.rowid)
    } else {
      cohort_data <- list(D = did_cohort_index, y = dp2$time_invariant_data[[dp2$yname]],
                          post = dp2$time_invariant_data$post,
                          i.weights = dp2$time_invariant_data$weights,
                          .rowid = dp2$time_invariant_data$.rowid)
    }
    covariates <- dp2$covariates_matrix
  }

  # run estimation
  force_rc <- !is.null(dp2$fix_weights) && dp2$fix_weights == "varying" && dp2$panel
  did_result <- tryCatch(run_DRDID(cohort_data, covariates, dp2, g_val = dp2$treated_groups[g], t_val = dp2$time_periods[t+tfac], force_rc = force_rc),
                         error = function(e) {
                           warning("Error computing internal 2x2 DiD for (g, t) = (", dp2$treated_groups[g], ", ", dp2$time_periods[t+tfac], "): ", e$message, ". The ATT for this cell will be set to NA.")
                           return(NULL)
                         })

  # When force_rc on balanced panel, the influence function has 2*n_units rows.
  # Half-split is safe here: cohort_data is explicitly stacked as
  # [all pre, all post] via rep(c(0L, 1L), each = n_units) in construction above.
  if (force_rc && !is.null(did_result) && dp2$panel && !is.null(did_result$inf_func)) {
    inf <- did_result$inf_func
    n_half <- length(inf) %/% 2L
    did_result$inf_func <- inf[1:n_half] + inf[(n_half + 1):(2L * n_half)]
  }

  return(did_result)

}

#' @title Compute Group-Time Average Treatment Effects
#'
#' @description `compute.att_gt2` does the (g,t) cell computation and sends it to estimation,
#' then does all the post-processing after estimation
#'
#' @param dp2 A DIDparams object v2.0
#'
#' @return a list with length equal to the number of groups times the
#'  number of time periods; each element of the list contains an
#'  object that contains group-time average treatment effect as well
#'  as which group it is for and which time period it is for. It also exports
#'  the influence function which is used externally to compute
#'  standard errors.
#'
#' @keywords internal
#' @export
compute.att_gt2 <- function(dp2) {

  n <- dp2$id_count  # Total number of units
  # whether to compute influence functions (default TRUE; FALSE = point estimates only)
  do_inf <- is.null(dp2$compute_inffunc) || isTRUE(dp2$compute_inffunc)
  time_periods <- dp2$time_periods # tlist

  tlist.length <- if (dp2$base_period != "universal") length(time_periods) - 1L else length(time_periods)
  tfac <- if (dp2$base_period != "universal") 1L else 0L

  # Pre-compute cumulative cohort sizes for fast indexing in get_did_cohort_index
  if (dp2$panel) {
    dp2$.cohort_cum_sizes <- cumsum(dp2$cohort_counts$cohort_size)
  }


  # Pre-compute pret for each group (used for post-treatment periods and universal base)
  # This avoids calling tail(which(...)) repeatedly inside run_att_gt_estimation
  pret_by_group <- vapply(seq_len(dp2$treated_groups_count), function(g) {
    idx <- which((time_periods + dp2$anticipation) < dp2$treated_groups[g])
    if (length(idx) == 0L) NA_integer_ else idx[length(idx)]
  }, integer(1))
  dp2$.pret_by_group <- pret_by_group

  # in terms of indexes, not calendar times.
  # Build the (g,t) index pairs directly in g-major / t-minor order. This is
  # identical to expand.grid(g, t) followed by order(g, t) but avoids both the
  # data.frame construction and the per-iteration single-row subset below.
  g_vec <- rep.int(seq_len(dp2$treated_groups_count), rep.int(tlist.length, dp2$treated_groups_count))
  t_vec <- rep.int(seq_len(tlist.length), dp2$treated_groups_count)
  total_gt_iterations <- length(g_vec)


  # Running estimation using run_att_gt_estimation() function
  # Helper function for processing each (g,t) pair
  process_gt <- function(g, t) {

    # Run estimation
    gt_result <- run_att_gt_estimation(g, t, dp2)

      # Compute post-treatment indicator
    post.treat <- as.integer(dp2$treated_groups[g] <= dp2$time_periods[t+tfac])

    # Check for NULL first (estimation failed or was skipped). When do_inf = FALSE
    # (point estimates only) no influence-function vector is stored (saving memory).
    if (is.null(gt_result)) {
      inffunc_updates <- if (do_inf) rep(NA_real_, n) else NULL
      gt_result <- list(att = NA, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat, inffunc_updates = inffunc_updates)
      return(gt_result)
    }

    # Base period normalization: ATT is 0 by construction
    if (!is.null(gt_result$base_period_norm)) {
      inffunc_updates <- if (do_inf) rep(0, n) else NULL
      gt_result <- list(att = 0, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat, inffunc_updates = inffunc_updates)
      return(gt_result)
    }

    if (is.null(gt_result$att)) {
      # Estimation returned a result but without an ATT
      inffunc_updates <- if (do_inf) rep(NA_real_, n) else NULL
      gt_result <- list(att = NA, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat, inffunc_updates = inffunc_updates)
      return(gt_result)
    } else {
      att <- gt_result$att
      inf_func <- gt_result$inf_func

      # Handle NaN ATT: treat as estimation failure
      if (is.nan(att)) {
        att <- NA
        if (do_inf) inf_func <- rep(NA_real_, n)
      }

      # Save ATT and influence function (inf_func is NULL when do_inf = FALSE)
      inffunc_updates <- inf_func
      gt_result <- list(att = att, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat, inffunc_updates = inffunc_updates)
      return(gt_result)
    }
  }

  # run the estimation for each (g,t) pair with process_gt
  gt_results <- lapply(seq_len(total_gt_iterations), function(idx) {
    process_gt(g_vec[idx], t_vec[idx])
  })

  # Filter out NULL results
  gt_results <- Filter(Negate(is.null), gt_results)

  if (length(gt_results) == 0L) {
    stop("No valid (g, t) cells found for att_gt() estimation. ",
         "Check treatment timing, control group definition, and anticipation settings.")
  }

  # Post processing: Build sparse influence function matrix directly from triplets.
  # Skipped entirely for point estimates only (do_inf = FALSE) -- no n x k matrix is
  # ever formed, which is the memory saving of compute_inffunc = FALSE.
  if (!do_inf) {
    inffunc <- NULL
  } else {
    n_rows <- length(gt_results[[1]]$inffunc_updates)
    n_cols <- length(gt_results)

    # Collect non-zero (and NA) entries as sparse triplets. Gather per-column
    # results into lists and concatenate once, rather than growing the triplet
    # vectors with c() on every iteration (which reallocates O(n_cols^2) in total).
    # The triplet order is unchanged (column-major, original within-column order),
    # so the resulting sparse matrix is identical.
    nz_list <- vector("list", n_cols)
    val_list <- vector("list", n_cols)
    for (j in seq_len(n_cols)) {
      vec <- gt_results[[j]]$inffunc_updates
      nz <- which(is.na(vec) | vec != 0)
      nz_list[[j]] <- nz
      val_list[[j]] <- vec[nz]
    }
    trip_i <- unlist(nz_list, use.names = FALSE)
    trip_x <- unlist(val_list, use.names = FALSE)
    trip_j <- rep.int(seq_len(n_cols), lengths(nz_list))
    inffunc <- Matrix::sparseMatrix(i = trip_i, j = trip_j, x = trip_x,
                                    dims = c(n_rows, n_cols))
  }

  return(list(attgt.list=gt_results, inffunc=inffunc))
}
