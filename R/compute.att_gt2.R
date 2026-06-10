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

    # flags for treated and control units. The per-group flag, the
    # never-treated flag and the per-period masks are loop-invariant across
    # (g,t) cells, so they are pre-computed once in compute.att_gt2() instead
    # of rescanning every row per cell (complete.cases() ran upstream in
    # pre_process_did2, so == has no NA hazard relative to %in%).
    Gflag <- dp2$.gflag_by_group[[group]]

    if (dp2$control_group == "nevertreated") {
      Cflag <- dp2$.never_treated
    } else {  # not-yet-treated: the cutoff varies with t, so this comparison stays per cell
      Cflag <- dp2$.never_treated |
        (dat[[dp2$gname]] > dp2$time_periods[max(time, pret) + tfac] + dp2$anticipation &
           !Gflag)
    }

    # keep only rows observed in pret or t
    keep <- dp2$.period_masks[[time + tfac]] | dp2$.period_masks[[pret]]

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
#' @return A list containing the estimated ATT and the influence function as
#'  sparse triplets (`if_i` row indices, `if_x` values), except on the force_rc
#'  path which returns the dense stacked vector (`inf_func`) for the caller's
#'  half-fold. `if_i`/`if_x` (or `inf_func`) are NULL for point estimates only.
#' @noRd
run_DRDID <- function(cohort_data, covariates, dp2, g_val = NULL, t_val = NULL,
                      pret_val = NULL, force_rc = FALSE, check_cache_key = NULL){

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
      intercept_only <- FALSE
    } else {
      covariates <- rep(1, length(valid_obs))
      intercept_only <- TRUE
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

      # Guard booleans: reuse the cached result when run_att_gt_estimation()
      # established that this cell's (covariates, D) inputs are bit-identical to
      # an earlier cell's (panel + nevertreated + non-varying weights; see
      # compute.att_gt2). Warnings are still emitted per cell below, so cached
      # failures keep the per-cell warning count and text unchanged.
      guard <- if (!is.null(check_cache_key)) dp2$.check_cache[[check_cache_key]] else NULL
      if (is.null(guard)) {
        # overlap check for pscore based methods; regression-feasibility check
        # for outcome-regression methods, run only if overlap passed (matching
        # the early return on overlap failure below)
        overlap_fail <- (dp2$est_method %in% c("dr", "ipw")) &&
          overlap_check_fail(covariates, D_vec, intercept_only)
        rcond_fail <- !overlap_fail && (dp2$est_method %in% c("dr", "reg")) &&
          rcond_check_fail(covariates, D_vec, intercept_only)
        guard <- list(overlap_fail = overlap_fail, rcond_fail = rcond_fail)
        if (!is.null(check_cache_key)) dp2$.check_cache[[check_cache_key]] <- guard
      }

      if (guard$overlap_fail) {
        warning(paste0("overlap condition violated", gt_label))
        return(list(att = NA, if_i = if (do_inf) seq_len(n) else NULL,
                    if_x = if (do_inf) rep(NA_real_, n) else NULL))
      }
      if (guard$rcond_fail) {
        warning(paste0("Covariate matrix for control units is singular or numerically ill-conditioned", gt_label, "; consider centering/rescaling covariates or removing collinear terms"))
        return(list(att = NA, if_i = if (do_inf) seq_len(n) else NULL,
                    if_x = if (do_inf) rep(NA_real_, n) else NULL))
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
    # subgroup to estimate att(g,t). Stored as sparse triplets (if_i row
    # indices / if_x values) rather than a dense length-n vector: the assembly
    # in compute.att_gt2() concatenates the triplets directly, mirroring the
    # slow path, instead of rescanning dense columns with which().
    if (do_inf) {
      if_i <- valid_obs
      if_x <- (n/n1)*attgt$att.inf.func
    } else {
      if_i <- NULL   # point estimates only
      if_x <- NULL
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
      intercept_only <- FALSE
    } else {
      covariates <- covariates[valid_obs]
      intercept_only <- TRUE
    }

    covariates <- as.matrix(covariates)

    # number of rows of the final influence-function column for early returns:
    # the unbalanced-panel branch aggregates observations to units; force_rc
    # failure triplets are already at unit level (the caller's half-fold only
    # applies to the dense success vector); plain RC has one row per observation
    n_inf <- if (dp2$allow_unbalanced_panel || force_rc) dp2$id_count else n

    skip_this_att_gt <- FALSE
    if (sum(D * post) == 0) {
      warning(paste0("No units in group ", g_val, " in time period ", t_val))
      skip_this_att_gt <- TRUE
    }
    if (sum(D * (1 - post)) == 0) {
      warning(paste0("No units in group ", g_val, " in time period ", pret_val))
      skip_this_att_gt <- TRUE
    }
    if (sum((1 - D) * post) == 0) {
      warning(paste0("No available control units for group ", g_val,
                     " in time period ", t_val))
      skip_this_att_gt <- TRUE
    }
    if (sum((1 - D) * (1 - post)) == 0) {
      warning(paste0("No available control units for group ", g_val,
                     " in time period ", pret_val))
      skip_this_att_gt <- TRUE
    }
    if (skip_this_att_gt) {
      return(list(att = NA, if_i = if (do_inf) seq_len(n_inf) else NULL,
                  if_x = if (do_inf) rep(NA_real_, n_inf) else NULL))
    }

    #-----------------------------------------------------------------------------
    # check for overlap and regression problems
    #-----------------------------------------------------------------------------
    custom_est_method <- inherits(dp2$est_method, "function")

    if (!custom_est_method) {
      D_vec <- D

      # checks for pscore based methods (no caching on the RC branch: the
      # cohort index varies cell by cell here)
      if (dp2$est_method %in% c("dr", "ipw")) {
        if (overlap_check_fail(covariates, D_vec, intercept_only)) {
          warning(paste0("overlap condition violated", gt_label))
          return(list(att = NA, if_i = if (do_inf) seq_len(n_inf) else NULL,
                      if_x = if (do_inf) rep(NA_real_, n_inf) else NULL))
        }
      }

      # check if can run regression using control units
      if (dp2$est_method %in% c("dr", "reg")) {
        if (rcond_check_fail(covariates, D_vec, intercept_only)) {
          warning(paste0("Covariate matrix for control units is singular or numerically ill-conditioned", gt_label, "; consider centering/rescaling covariates or removing collinear terms"))
          return(list(att = NA, if_i = if (do_inf) seq_len(n_inf) else NULL,
                      if_x = if (do_inf) rep(NA_real_, n_inf) else NULL))
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
    # subgroup to estimate att(g,t). Stored as sparse triplets (if_i/if_x),
    # mirroring the slow path; the dense vector is kept only where it is still
    # needed (unbalanced-panel rowsum aggregation, force_rc half-fold).
    if (!do_inf) {
      if_i <- NULL   # point estimates only
      if_x <- NULL
    } else if(dp2$allow_unbalanced_panel){
      # since this is technically a panel data but ran as RCS, we need to adjust the
      # influence function by aggregating influence value by .rowid (several obs of one
      # unit can be used to estimate ATT in each 2x2). rowsum(..., reorder = FALSE)
      # reproduces the data.table `by = .rowid` first-appearance group order exactly.
      # The aggregated dense vector is then converted to triplets with the same
      # is.na | != 0 scan the assembly previously applied to dense columns.
      inf_func_long <- numeric(n)
      inf_func_long[valid_obs] <- (dp2$id_count / n1) * attgt$att.inf.func
      inf_func_vector <- rowsum(inf_func_long, rowid_all, reorder = FALSE)[, 1L]
      if_i <- which(is.na(inf_func_vector) | inf_func_vector != 0)
      if_x <- inf_func_vector[if_i]

    } else if (force_rc) {
      # balanced panel run as RC (fix_weights = "varying"): return the dense
      # stacked [all pre, all post] vector; run_att_gt_estimation() half-folds
      # it to unit level and converts it to triplets there.
      inf_func_vector <- rep(0, n)
      inf_func_vector[valid_obs] <- (n/n1)*attgt$att.inf.func
      return(list(att = attgt$ATT, inf_func = inf_func_vector))
    } else {
      if_i <- valid_obs
      if_x <- (n/n1)*attgt$att.inf.func
    }

  }

  return(list(att = attgt$ATT, if_i = if_i, if_x = if_x))

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



  # key into the overlap/rcond guard cache (set on the panel non-varying branch
  # below when the cache is active; NULL = no caching for this cell)
  check_cache_key <- NULL

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
      # Key for the overlap/rcond guard cache (see compute.att_gt2): under
      # nevertreated the cohort index depends only on g, so the guards'
      # (covariates, D) inputs are bit-identical across all cells of a group
      # sharing a covariate (and weight) period.
      if (!is.null(dp2$.check_cache)) {
        check_cache_key <- paste(g, base::min(pret, t), w_idx, sep = "_")
      }
    }
  } else {

    # post indicator reuses the pre-computed period mask (see compute.att_gt2)
    log_vec <- dp2$.period_masks[[t + tfac]]
    # convert TRUE/FALSE to 1/0 in place (fastest)
    set(dp2$time_invariant_data, j = "post", value = as.integer(log_vec))

    # Handle fix_weights for RC/unbalanced panel
    if (!is.null(dp2$fix_weights) && dp2$fix_weights %in% c("base_period", "first_period")) {
      if (dp2$fix_weights == "base_period") {
        target_period <- dp2$time_periods[dp2$.pret_by_group[g]]
      } else {
        target_period <- dp2$time_periods[1]
      }
      tid <- dp2$time_invariant_data
      # Memoize the weight lookup per target period (env set in compute.att_gt2):
      # target_period depends only on g for "base_period" and is constant for
      # "first_period", and tid's id/time/.w columns are never mutated in the
      # (g,t) loop, so fixed_w/na_w are identical across all cells sharing a
      # target period. Keying by as.character(target_period) memoizes the NA
      # case (group with no pre-period, varying base, pre-treatment cell)
      # under "NA" with unchanged semantics. The warning and the
      # did_cohort_index mutation below stay per cell.
      memo <- dp2$.fixed_w_memo
      memo_key <- as.character(target_period)
      cached <- if (!is.null(memo)) memo[[memo_key]] else NULL
      if (is.null(cached)) {
        # Build weight lookup from target period
        target_mask <- tid[[dp2$tname]] == target_period
        target_ids <- tid[[dp2$idname]][target_mask]
        target_ws <- tid[[".w"]][target_mask]
        # Map each observation's weight from the target period by integer id matching.
        # Equivalent to the previous named-character lookup (first match wins on any
        # duplicate id; unmatched ids -> NA) but avoids coercing every id to character.
        fixed_w <- target_ws[match(tid[[dp2$idname]], target_ids)]
        cached <- list(fixed_w = fixed_w, na_w = is.na(fixed_w))
        if (!is.null(memo)) memo[[memo_key]] <- cached
      }
      fixed_w <- cached$fixed_w
      # Exclude units not observed in target period by setting D to NA
      # (run_DRDID filters on !is.na(D))
      na_w <- cached$na_w
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
                          i.weights = dp2$time_invariant_data[[".w"]],
                          .rowid = dp2$time_invariant_data$.rowid)
    }
    covariates <- dp2$covariates_matrix
  }

  # run estimation
  force_rc <- !is.null(dp2$fix_weights) && dp2$fix_weights == "varying" && dp2$panel
  did_result <- tryCatch(run_DRDID(cohort_data, covariates, dp2, g_val = dp2$treated_groups[g],
                                   t_val = dp2$time_periods[t + tfac],
                                   pret_val = dp2$time_periods[pret],
                                   force_rc = force_rc,
                                   check_cache_key = check_cache_key),
                         error = function(e) {
                           warning("Error computing internal 2x2 DiD for (g, t) = (", dp2$treated_groups[g], ", ", dp2$time_periods[t+tfac], "): ", e$message, ". The ATT for this cell will be set to NA.")
                           return(NULL)
                         })

  # When force_rc on balanced panel, the influence function has 2*n_units rows.
  # Half-split is safe here: cohort_data is explicitly stacked as
  # [all pre, all post] via rep(c(0L, 1L), each = n_units) in construction above.
  # The folded dense vector is converted to sparse triplets (if_i/if_x) with the
  # same is.na | != 0 scan the assembly previously applied to dense columns.
  if (force_rc && !is.null(did_result) && dp2$panel && !is.null(did_result$inf_func)) {
    inf <- did_result$inf_func
    n_half <- length(inf) %/% 2L
    inf_folded <- inf[1:n_half] + inf[(n_half + 1):(2L * n_half)]
    did_result$if_i <- which(is.na(inf_folded) | inf_folded != 0)
    did_result$if_x <- inf_folded[did_result$if_i]
    did_result$inf_func <- NULL
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
  } else {
    # Pre-compute the loop-invariant row masks used by get_did_cohort_index()
    # and the post indicator in run_att_gt_estimation(): per-period membership,
    # the never-treated flag, and per-group treated flags. Each was previously
    # recomputed with full O(rows) scans for every (g,t) cell on the RC /
    # unbalanced-panel path. Memory: (T + n_groups + 1) logical vectors of
    # length nrow(time_invariant_data).
    dat <- dp2$time_invariant_data
    dp2$.period_masks <- lapply(time_periods, function(tp) dat[[dp2$tname]] == tp)
    dp2$.never_treated <- dat[[dp2$gname]] == Inf
    dp2$.gflag_by_group <- lapply(dp2$treated_groups, function(gval) dat[[dp2$gname]] == gval)
    # Memo for the fix_weights = "base_period"/"first_period" weight lookup in
    # run_att_gt_estimation(): fixed_w/na_w depend only on the target period
    # (one per group for "base_period", a single constant for "first_period"),
    # so the full-table mask + match() is computed once per distinct period
    # instead of once per (g,t) cell.
    if (!is.null(dp2$fix_weights) && dp2$fix_weights %in% c("base_period", "first_period")) {
      dp2$.fixed_w_memo <- new.env(parent = emptyenv())
    }
  }

  # Cache for the per-cell overlap/rcond guard booleans. For panel data with
  # control_group = "nevertreated" and non-varying weights, the guards'
  # (covariates, D) inputs are bit-identical across all cells of a group that
  # share a covariate (and weight) period, so each boolean is computed once per
  # (g, covariate-period, weight-period) key (see run_DRDID). Not used for
  # notyettreated (the control cohort varies with t) or fix_weights = "varying"
  # (RC estimator path). Only the booleans are stored; failure warnings are
  # still emitted per cell. options(did.disable_check_cache = TRUE) forces the
  # per-cell checks (escape hatch / used to verify bit-identical equivalence).
  if (dp2$panel && dp2$control_group[1] == "nevertreated" &&
      (is.null(dp2$fix_weights) || dp2$fix_weights != "varying") &&
      !isTRUE(getOption("did.disable_check_cache"))) {
    dp2$.check_cache <- new.env(parent = emptyenv())
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
    # (point estimates only) no influence-function triplets are stored (saving memory).
    # Failed cells propagate NA to every row of their column: if_i = seq_len(n),
    # if_x = NA, exactly the triplets the dense-column scan used to produce.
    if (is.null(gt_result)) {
      gt_result <- list(att = NA, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat,
                        if_i = if (do_inf) seq_len(n) else NULL,
                        if_x = if (do_inf) rep(NA_real_, n) else NULL)
      return(gt_result)
    }

    # Base period normalization: ATT is 0 by construction (empty triplet column)
    if (!is.null(gt_result$base_period_norm)) {
      gt_result <- list(att = 0, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat,
                        if_i = if (do_inf) integer(0) else NULL,
                        if_x = if (do_inf) numeric(0) else NULL)
      return(gt_result)
    }

    if (is.null(gt_result$att)) {
      # Estimation returned a result but without an ATT
      gt_result <- list(att = NA, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat,
                        if_i = if (do_inf) seq_len(n) else NULL,
                        if_x = if (do_inf) rep(NA_real_, n) else NULL)
      return(gt_result)
    } else {
      att <- gt_result$att
      if_i <- gt_result$if_i
      if_x <- gt_result$if_x

      # Handle NaN ATT: treat as estimation failure
      if (is.nan(att)) {
        att <- NA
        if (do_inf) {
          if_i <- seq_len(n)
          if_x <- rep(NA_real_, n)
        }
      }

      # Save ATT and influence triplets (if_i/if_x are NULL when do_inf = FALSE)
      gt_result <- list(att = att, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat, if_i = if_i, if_x = if_x)
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
    n_cols <- length(gt_results)

    # Each cell already carries its influence column as sparse triplets
    # (if_i row indices / if_x values), mirroring the slow path; concatenate
    # them once instead of rescanning dense per-cell vectors with which().
    # The triplet order is unchanged (column-major, original within-column
    # order), so the resulting sparse matrix is identical. Every code path
    # produces columns with dp2$id_count rows (panel units, RC observations,
    # unbalanced-panel units after the rowsum aggregation, force_rc units
    # after the half-fold).
    nz_list <- lapply(gt_results, `[[`, "if_i")
    val_list <- lapply(gt_results, `[[`, "if_x")
    trip_i <- unlist(nz_list, use.names = FALSE)
    trip_x <- unlist(val_list, use.names = FALSE)
    trip_j <- rep.int(seq_len(n_cols), lengths(nz_list))
    inffunc <- Matrix::sparseMatrix(i = trip_i, j = trip_j, x = trip_x,
                                    dims = c(n, n_cols))
  }

  return(list(attgt.list=gt_results, inffunc=inffunc))
}
