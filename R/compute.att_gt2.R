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
#' @export
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
#' @export
get_did_cohort_index <- function(group, time, pret, dp2){
  # return a vector of dimension id_size with 1, 0 or NA values
  treated_groups <- dp2$treated_groups
  time_periods <- dp2$time_periods
  # based on control_group option
  min_control_group_index <- ifelse(dp2$control_group == "notyettreated", max(time, pret) + dp2$anticipation + 1, Inf)
  min_control_group <- ifelse(min_control_group_index == Inf, Inf, time_periods[min_control_group_index])
  max_control_group <- Inf # always include the never treated units


  # select the DiD cohort
  did_cohort_index <- rep(NA, dp2$id_count)
  # getting the index to get units who will participate in the estimation for the (g,t) cell.
  start_control <- dp2$cohort_counts[cohort < min_control_group, sum(cohort_size)]+1
  end_control <- dp2$cohort_counts[cohort <= max_control_group, sum(cohort_size)]
  index <- which(dp2$cohort_counts[, cohort] == time_periods[group])
  start_treat <- ifelse(index == 1, 1, dp2$cohort_counts[1:(index-1), sum(cohort_size)]+1)
  end_treat <- dp2$cohort_counts[1:index, sum(cohort_size)]
  # set the cohort index; .C = 0 and .G = 1
  did_cohort_index[start_control:end_control] <- 0
  did_cohort_index[start_treat:end_treat] <- 1

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
#' @export
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

    # preparing vectors
    # Ypost <- cohort_data[, y1]
    # Ypre <- cohort_data[, y0]
    # w <- cohort_data[, i.weights]
    # G <- cohort_data[, D]

    if (inherits(dp2$est_method, "function")) {
      # user-specified function
      attgt <- est_method(y1=cohort_data[, y1],
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
    covariates <- covariates[valid_obs,]
    # add the intercept
    # Check if the ".intercept" name is already in the covariates matrix
    #if(".intercept" %in% names(covariates)){stop("did is trying to impute a new column .intercept, but this already exists. Please check your dataset")}
    #covariates[, .intercept := -1L]
    # coerce to a matrix
    covariates <- as.matrix(covariates)

    # num obs. for computing ATT(g,t)
    n1 <- cohort_data[, .N]

    # TODO; IMPLEMENT REPEATED CROSS-SECTION
    attgt <- NULL

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
#'
#' @export
run_att_gt_estimation <- function(gt, dp2){
  # get the current g,t values
  g <- gt[[1]]
  t <- gt[[2]]
  if(dp2$print_details){cat("\n Evaluating (g,t) = (",g,t,")")}
  pret <- get_pret(g, t, dp2$base_period, dp2$anticipation)

  # if we are in period (g-1) or base period out of bounds, normalize results to be equal to NULL
  # and break without computing anything
  if(t == pret | !pret %in% seq_along(dp2$time_periods)){
    if(dp2$print_details){cat("\n Skipping (g,t) = (",g,t,") as base period is out of bounds or for varying base period.")}
    return(NULL)
  }

  # get units in treatment and control group
  did_cohort_index <- get_did_cohort_index(g, t, pret, dp2)
  # In case of not treatment and control group in the cohort, return NULL
  valid_did_cohort <- any(did_cohort_index == 1) & any(did_cohort_index == 0)
  if(!isTRUE(valid_did_cohort)){
    if(dp2$print_details){cat("\n Skipping (g,t) = (",g,t,") as no treatment group or control group found")}
    return(NULL)
  }

  # Get the matrix of covariates
  covariates <- dp2$covariates
  cohort_data <- data.table(did_cohort_index, dp2$outcomes_tensor[[t]], dp2$outcomes_tensor[[pret]], dp2$weights_vector)
  names(cohort_data) <- c("D", "y1", "y0", "i.weights")

  # run estimation
  did_result <- tryCatch(run_DRDID(cohort_data, covariates, dp2),
                         error = function(e) {
                           warning("\n Error in computing internal 2x2 DiD for (g,t) = (",g,t,"):", e$message)
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
#'
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
