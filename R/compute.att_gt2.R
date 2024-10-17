#' @title Get pre treatment period
#' @description A utility function to get the pre treatment period for a given g,t pair.
#'
#' @param group Group.
#' @param time Time period.
#' @param base_period Whether to use a "varying" base period or a
#'  "universal" base period.
#' @param anticipation  The number of time periods before participating
#'  in the treatment where units can anticipate participating in the
#'  treatment and therefore it can affect their untreated potential outcomes
#'
#' @return Time period indicating the pre treatment period.
#' @noRd
get_pret <- function(group, time, base_period, anticipation){
  if(base_period == "universal"){
    pret <- group - 1 - anticipation
  } else {
    pret <- ifelse(time >= group, group - 1 - anticipation, time - 1)
  }

  return(pret)
}

#' @title Get the cohort for the current g,t values
#' @description A utility function to get vector with index of units that will be included in the DiD estimation for the current g,t pair.
#'
#' @param group Group.
#' @param time Time period.
#' @param pret Pre treatment period.
#' @param dp2 DiD parameters v2.0.
#'
#' @return Time period indicating the pre treatment period.
#' @noRd
get_did_cohort_index <- function(group, time, pret, dp2){
  # return a vector of dimension id_size with 1, 0 or NA values
  treated_groups <- dp2$treated_groups
  time_periods <- dp2$time_periods
  # based on control_group option
  min_control_group_index <- ifelse(dp2$control_group == "notyettreated", max(time, pret) + dp2$anticipation + 1, Inf)
  min_control_group <- ifelse(min_control_group_index == Inf | min_control_group_index > dp2$time_periods_count,
                              Inf,
                              time_periods[min_control_group_index])
  max_control_group <- Inf # always include the never treated units as the maximum.


  # select the DiD cohort
  did_cohort_index <- rep(NA, dp2$id_count)

  if(dp2$panel){

    # getting the index to get units who will participate in the estimation for the (g,t) cell.
    start_control <- dp2$cohort_counts[cohort < min_control_group, sum(cohort_size)]+1
    end_control <- dp2$cohort_counts[cohort <= max_control_group, sum(cohort_size)]
    index <- which(dp2$cohort_counts[, cohort] == time_periods[group])
    start_treat <- ifelse(index == 1, 1, dp2$cohort_counts[1:(index-1), sum(cohort_size)]+1)
    end_treat <- dp2$cohort_counts[1:index, sum(cohort_size)]
    # set the cohort index; .C = 0 and .G = 1
    did_cohort_index[start_control:end_control] <- 0
    did_cohort_index[start_treat:end_treat] <- 1
  } else {
    # getting the index to get units who will participate in the estimation for the (g,t) cell.
    # Note: This works because the data is already ordered in a specific way. Changing that order will break this.

    # Pre-extract column names for use within the loop
    col_names <- names(dp2$crosstable_counts)[-1]  # Exclude 'T' column
    min_idx <- match(as.character(min_control_group), col_names)
    max_idx <- match(as.character(max_control_group), col_names)
    # correction in case there is not never-treated units! Pick the very last time period...
    ifelse(is.na(min_idx), min_idx <- match(tail(col_names,1), col_names), min_idx <- min_idx)
    ifelse(is.na(max_idx), max_idx <- match(tail(col_names,1), col_names), max_idx <- max_idx)

    # Start the loop
    for (i in c(pret, time)) {

      # ------------------------------------------------
      # SET UP: Identify index_time, start_time, and precompute relevant counts
      # ------------------------------------------------

      # Identify the index_time for 'i' based on period_counts and compute the starting point
      index_time <- which(dp2$period_counts[, period] == dp2$reverse_mapping[i])
      start_time <- ifelse(index_time == 1, 1, dp2$period_counts[1:(index_time - 1), sum(period_size)] + 1)

      # Extract the relevant row for current time period 'i' only once
      relevant_g_row <- dp2$crosstable_counts[i, ]

      # ------------------------------------------------
      # FOR CONTROL UNITS FIRST
      # ------------------------------------------------

      # Calculate the starting and ending index for control units (min_control_group to max_control_group)
      start_control <- start_time + sum(relevant_g_row[, col_names[1:min_idx - 1], with = FALSE], na.rm = TRUE)
      csize <- sum(relevant_g_row[, col_names[min_idx:max_idx], with = FALSE], na.rm = TRUE)
      end_control <- start_control + csize - 1

      # Impute control units with 0
      did_cohort_index[start_control:end_control] <- 0

      # ------------------------------------------------
      # FOR TREATED UNITS
      # ------------------------------------------------

      # Calculate correction_index based on the current g and the relevant row for the current time period
      index_cohort <- which(dp2$cohort_counts[, cohort] == dp2$reverse_mapping[group])
      correction_index <- ifelse(index_cohort == 1, 0, sum(relevant_g_row[, col_names[1:(index_cohort - 1)], with = FALSE], na.rm = TRUE))

      # Compute start and end points for treated units
      start_treat <- start_time + correction_index
      g_size <- relevant_g_row[, get(as.character(dp2$reverse_mapping[group]))]
      end_treat <- start_treat + g_size - 1

      # Impute treated units with 1
      did_cohort_index[start_treat:end_treat] <- 1
    }

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
#' @return Time period indicating the pre treatment period.
#' @noRd
run_DRDID <- function(cohort_data, covariates, dp2){

  if(dp2$panel){
    # --------------------------------------
    # Panel Data
    # --------------------------------------

    # still total number of units (not just included in G or C)
    n <- cohort_data[, .N]

    # pick up the indices for units that will be used to compute ATT(g,t)
    valid_obs <- which(cohort_data[, !is.na(D)])
    cohort_data <- cohort_data[valid_obs]
    # TODO; THIS HAS TO BE BETTER WRITTEN
    if(dp2$xformla != ~1){
      covariates <- covariates[valid_obs,]
    } else {
      covariates <- rep(1, length(valid_obs))
    }
    covariates <- as.matrix(covariates)

    # num obs. for computing ATT(g,t)
    n1 <- cohort_data[, .N]

    #-----------------------------------------------------------------------------
    # code for actually computing ATT(g,t)
    #-----------------------------------------------------------------------------


    if (inherits(dp2$est_method, "function")) {
      # user-specified function
      attgt <- dp2$est_method(y1=cohort_data[, y1],
                          y0=cohort_data[, y0],
                          D=cohort_data[, D],
                          covariates=covariates,
                          i.weights=cohort_data[, i.weights],
                          inffunc=TRUE)
    } else if (dp2$est_method == "ipw") {
      # inverse-probability weights
      attgt <- std_ipw_did_panel(y1=cohort_data[, y1],
                                        y0=cohort_data[, y0],
                                        D=cohort_data[, D],
                                        covariates=covariates,
                                        i.weights=cohort_data[, i.weights],
                                        boot=FALSE, inffunc=TRUE)
    } else if (dp2$est_method == "reg") {
      # regression
      attgt <- reg_did_panel(y1=cohort_data[, y1],
                                    y0=cohort_data[, y0],
                                    D=cohort_data[, D],
                                    covariates=covariates,
                                    i.weights=cohort_data[, i.weights],
                                    boot=FALSE, inffunc=TRUE)
    } else {
      # doubly robust, this is default
      attgt <- drdid_panel(y1=cohort_data[, y1],
                                  y0=cohort_data[, y0],
                                  D=cohort_data[, D],
                                  covariates=covariates,
                                  i.weights=cohort_data[, i.weights],
                                  boot=FALSE, inffunc=TRUE)
    }

    # adjust influence function to account for only using
    # subgroup to estimate att(g,t)
    inf_func_vector <- rep(0, n)
    inf_func_not_na <- (n/n1)*attgt$att.inf.func
    inf_func_vector[valid_obs] <- inf_func_not_na

  } else {
    # --------------------------------------
    # Repeated Cross-Section
    # --------------------------------------
    # still total number of units (not just included in G or C)
    n <- cohort_data[, .N]

    # pick up the indices for units that will be used to compute ATT(g,t)
    valid_obs <- which(cohort_data[, !is.na(D)])
    cohort_data <- cohort_data[valid_obs]

    # num obs. for computing ATT(g,t)
    n1 <- cohort_data[, .N]

    # TODO; THIS HAS TO BE BETTER WRITTEN
    if(dp2$xformla != ~1){
      covariates <- covariates[valid_obs,]
    } else {
      covariates <- covariates[valid_obs]
    }

    covariates <- as.matrix(covariates)

    #-----------------------------------------------------------------------------
    # code for actually computing ATT(g,t)
    #-----------------------------------------------------------------------------

    if (inherits(dp2$est_method, "function")) {
      # user-specified function
      attgt <- dp2$est_method(y=cohort_data[, y],
                          post=cohort_data[, post],
                          D=cohort_data[, D],
                          covariates=covariates,
                          i.weights=cohort_data[, i.weights],
                          inffunc=TRUE)
    } else if (dp2$est_method == "ipw") {
      # inverse-probability weights
      attgt <- std_ipw_did_rc(y=cohort_data[, y],
                             post=cohort_data[, post],
                             D=cohort_data[, D],
                             covariates=covariates,
                             i.weights=cohort_data[, i.weights],
                             boot=FALSE, inffunc=TRUE)
    } else if (dp2$est_method == "reg") {
      # regression
      attgt <- reg_did_rc(y=cohort_data[, y],
                         post=cohort_data[, post],
                         D=cohort_data[, D],
                         covariates=covariates,
                         i.weights=cohort_data[, i.weights],
                         boot=FALSE, inffunc=TRUE)
    } else {
      # doubly robust, this is default
      attgt <- drdid_rc(y=cohort_data[, y],
                       post=cohort_data[, post],
                       D=cohort_data[, D],
                       covariates=covariates,
                       i.weights=cohort_data[, i.weights],
                       boot=FALSE, inffunc=TRUE)
    }

    # n/n1 adjusts for estimating the
    # att_gt only using observations from groups
    # G and C
    # adjust influence function to account for only using
    # subgroup to estimate att(g,t)
    inf_func_vector <- rep(0, n)
    inf_func_not_na <- (n/n1)*attgt$att.inf.func
    inf_func_vector[valid_obs] <- inf_func_not_na

  }

  return(list(att = attgt$ATT, inf_func = inf_func_vector))

}


#' @title Run ATT estimation for a given group-time pair
#'
#' @description `run_att_gt_estimation` does the main work for computing
#'  multiperiod group-time average treatment effects
#' @param gt A numeric vector of the group and time period coming from a list of (g,t) cells
#' @param dp2 A DIDparams object v2.0
#'
#' @return a list with the gt cell and the results after performing estimation
#'
#' @keywords internal
run_att_gt_estimation <- function(gt, dp2){
  # get the current g,t values
  g <- gt[[1]]
  t <- gt[[2]]
  if(dp2$print_details){cat("\n", paste0("Evaluating (g,t) = (",g,",",t,")"))}
  pret <- get_pret(g, t, dp2$base_period, dp2$anticipation)

  # if we are in period (g-1) or base period out of bounds, normalize results to be equal to NULL
  # and break without computing anything
  if(t == pret | !pret %in% seq_along(dp2$time_periods)){
    # if(dp2$print_details){cat("\n Skipping (g,t) as base period is out of bounds or for varying base period.")}
    return(NULL)
  }

  # get units in treatment and control group
  did_cohort_index <- get_did_cohort_index(g, t, pret, dp2)
  # In case of not treatment and control group in the cohort, return NULL
  valid_did_cohort <- any(did_cohort_index == 1) & any(did_cohort_index == 0)
  if(!isTRUE(valid_did_cohort)){
    if(dp2$print_details){cat("\n Skipping (g,t) as no treatment group or control group found")}
    return(NULL)
  }

  # Get the matrix of covariates
  covariates <- dp2$covariates
  if(dp2$panel){
    cohort_data <- data.table(did_cohort_index, dp2$outcomes_tensor[[t]], dp2$outcomes_tensor[[pret]], dp2$weights_vector)
    names(cohort_data) <- c("D", "y1", "y0", "i.weights")
  } else {
    dp2$time_invariant_data[, post := fifelse(get(dp2$tname) == dp2$reverse_mapping[t], 1, 0)]
    cohort_data <- data.table(did_cohort_index, dp2$time_invariant_data[[dp2$yname]], dp2$time_invariant_data$post, dp2$time_invariant_data$weights)
    names(cohort_data) <- c("D", "y", "post", "i.weights")
  }


  # run estimation
  did_result <- tryCatch(run_DRDID(cohort_data, covariates, dp2),
                         error = function(e) {
                           warning("\n Error in computing internal 2x2 DiD for (g,t) = (",g,",",t,"): ", e$message)
                           return(NULL)
                         })
  return(did_result)

}

#' @title Compute Group-Time Average Treatment Effects
#'
#' @description `compute.att_gt` does the (g,t) cell and send it to estimation,
#' then do all the post process after estimation
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
  treated_groups <- dp2$treated_groups
  time_periods <- dp2$time_periods
  # standardize the times to indexes
  # Create a mapping for time_periods to standardized form
  time_mapping <- setNames(seq_along(time_periods), time_periods)
  # Create a reverse mapping from standardized values back to the original
  reverse_mapping <- setNames(as.numeric(names(time_mapping)), time_mapping)
  dp2$reverse_mapping <- reverse_mapping
  # Convert both time_periods and treated_group using the mapping
  time_periods <- time_mapping[as.character(time_periods)]
  treated_groups <- time_mapping[as.character(treated_groups)]
  # Get (g,t) cells to perform computations. This replace the nested for-loop.
  gt_cells <- expand.grid(g = treated_groups, t = time_periods, stringsAsFactors = FALSE)
  gt_cells <- gt_cells[order(gt_cells$g, gt_cells$t), ]
  total_gt_iterations <- nrow(gt_cells)


  # Running estimation using run_att_gt_estimation() function
  # Helper function for processing each (g,t) pair
  process_gt <- function(gt_cell, idx) {
    g <- gt_cell$g
    t <- gt_cell$t

    # Run estimation
    gt_result <- run_att_gt_estimation(c(g, t), dp2)

    # Compute post-treatment indicator
    post.treat <- as.integer(g <= t)

    if (is.null(gt_result) || is.null(gt_result$att)) {
      # Estimation failed or was skipped
      if(dp2$base_period == "universal"){
        inffunc_updates <- rep(0, n)
        gt_result <- list(att = 0, group = reverse_mapping[g], year = reverse_mapping[t], post = post.treat, inffunc_updates = inffunc_updates)
        return(gt_result)
      } else {
        return(NULL)
      }

    } else {
      att <- gt_result$att
      inf_func <- gt_result$inf_func

      # Handle NaN ATT
      if (is.nan(att)) {
        att <- 0
        inf_func <- rep(0, n)
      }

      # Save ATT and influence function
      inffunc_updates <- inf_func
      gt_result <- list(att = att, group = reverse_mapping[g], year = reverse_mapping[t], post = post.treat, inffunc_updates = inffunc_updates)
      return(gt_result)
    }
  }

  # run the estimation for each (g,t) pair with process_gt
  gt_results <- lapply(seq_len(total_gt_iterations), function(idx) {
    process_gt(gt_cells[idx, ], idx)
  })

  # Filter out NULL results
  gt_results <- Filter(Negate(is.null), gt_results)

  # Post processing: Apply the updates to the sparse matrix in one shot
  n_rows <- length(gt_results[[1]]$inffunc_updates)
  update_inffunc <- data.table(matrix(NA_real_, nrow = n_rows, ncol = length(gt_results)))
  for (i in seq_along(gt_results)) {
    update_inffunc[[i]] <- gt_results[[i]]$inffunc_updates
  }

  # Update the sparse matrix with the values collected
  inffunc <- as(Matrix::Matrix(as.matrix(update_inffunc), sparse = TRUE), "CsparseMatrix")

  return(list(attgt.list=gt_results, inffunc=inffunc))
}
