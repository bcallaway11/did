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
  treated_groups <- dp2$treated_groups
  time_periods <- dp2$time_periods
  # based on control_group option

  min_control_group <-  ifelse((dp2$control_group == "notyettreated"),
                               dp2$cohort_counts$cohort[which(dp2$cohort_counts$cohort > dp2$time_periods[max(time, pret) + tfac ]+ dp2$anticipation)][1],
                               Inf)
  max_control_group <- Inf # always include the never treated units as the maximum. We add a correction in case is needed afterwards.


  # select the DiD cohort
  ifelse(dp2$allow_unbalanced_panel, did_cohort_index <- rep(NA, dp2$time_invariant_data[, .N]), did_cohort_index <- rep(NA, dp2$id_count))

  if(dp2$panel){

    # adding some correction in control group to avoid weird behavior when the control group is not yet treated
    if(!max_control_group %in% dp2$cohort_counts$cohort){
      max_control_group <- tail(dp2$cohort_count$cohort,1)
    }

    # getting the index to get units who will participate in the estimation for the (g,t) cell.
    start_control <- dp2$cohort_counts[cohort < min_control_group, sum(cohort_size)]+1
    end_control <- dp2$cohort_counts[cohort <= max_control_group, sum(cohort_size)]
    index <- which(dp2$cohort_counts[, cohort] == dp2$treated_groups[group])
    start_treat <- ifelse(index == 1, 1, dp2$cohort_counts[1:(index-1), sum(cohort_size)]+1)
    end_treat <- dp2$cohort_counts[1:index, sum(cohort_size)]
    # set the cohort index; .C = 0 and .G = 1
    did_cohort_index[start_control:end_control] <- 0
    did_cohort_index[start_treat:end_treat] <- 1
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
    # if we are running unbalanced panel data, we get a temporary copy of cohort data to compute influence function
    if(dp2$allow_unbalanced_panel){
      cohort_data_init <- copy(cohort_data[, .(D, .rowid)])
    }
    # still total number of units (not just included in G or C)
    n <- cohort_data[, .N]

    # pick up the indices for units that will be used to compute ATT(g,t)
    valid_obs <- which(cohort_data[, !is.na(D)])
    cohort_data <- cohort_data[valid_obs]

    # num obs. for computing ATT(g,t)
    n1 <- cohort_data[, .N]

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
    if(dp2$allow_unbalanced_panel){
      # since this is technically a panel data but ran as RCS, we need to adjust the influence function
      # by aggregating influence value by .rowid (since several obs of one unit could be used to estimate ATT in each 2x2)
      cohort_data_init[, inf_func_long := 0]
      # Assign values from vec to the valid rows identified by valid_obs
      cohort_data_init[valid_obs, inf_func_long := (dp2$id_count/n1)*attgt$att.inf.func]
      inf_func_vector <- cohort_data_init[, .(inf_func_agg = sum(inf_func_long, na.rm = TRUE)), by = .rowid][ , inf_func_agg]

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
  tfac <- ifelse(dp2$base_period != "universal", 1, 0)
  # set pret
  # varying base period
  pret <- t


  # universal base period
  if (dp2$base_period == "universal") {
    # use same base period as for post-treatment periods
    pret <- tail(which( (dp2$time_periods + dp2$anticipation) < dp2$treated_groups[g]),1)
  }

  # check if in post-treatment period
  if ((dp2$treated_groups[g]<=dp2$time_periods[(t+tfac)])) {

    # update pre-period if in post-treatment period to
    # be  period (g-delta-1)
    pret <- tail(which( (dp2$time_periods+dp2$anticipation) < dp2$treated_groups[g]),1)

    # print a warning message if there are no pre-treatment period
    if (length(pret) == 0) {
      warning(paste0("There are no pre-treatment periods for the group first treated at ", dp2$treated_groups[g], "\nUnits from this group are dropped"))

      # if there are not pre-treatment periods, code will
      # jump out of this loop
      return(NULL)
    }
  }

  #-----------------------------------------------------------------------------
  # if we are in period (g-1), normalize results to be equal to 0
  # and break without computing anything
  if (dp2$base_period == "universal") {
    if (dp2$time_periods[pret] == dp2$time_periods[(t+tfac)]) {
      return(NULL)
    }
  }

  # get units in treatment and control group
  did_cohort_index <- get_did_cohort_index(group = g, time = t, tfac = tfac, pret = pret, dp2)
  # In case of not treatment and control group in the cohort, return NULL
  valid_did_cohort <- any(did_cohort_index == 1) & any(did_cohort_index == 0)
  if(!isTRUE(valid_did_cohort)){
    if(dp2$print_details){cat("\n Skipping (g,t) as no treatment group or control group found")}
    return(NULL)
  }



  if(dp2$panel){
    cohort_data <- data.table(did_cohort_index, dp2$outcomes_tensor[[t+tfac]], dp2$outcomes_tensor[[pret]], dp2$weights_vector)
    names(cohort_data) <- c("D", "y1", "y0", "i.weights")
    covariates <- dp2$covariates_tensor[[base::min(pret, t)]]
  } else {

    log_vec <- dp2$time_invariant_data[[ dp2$tname ]] == dp2$time_periods[t+tfac]
    # convert TRUE/FALSE to 1/0 in place (fastest)
    set(dp2$time_invariant_data, j = "post", value = as.integer(log_vec))
    cohort_data <- data.table(did_cohort_index, dp2$time_invariant_data[[dp2$yname]], dp2$time_invariant_data$post, dp2$time_invariant_data$weights, dp2$time_invariant_data$.rowid)
    names(cohort_data) <- c("D", "y", "post", "i.weights", ".rowid")
    covariates <- dp2$covariates_matrix
  }

  # run estimation
  did_result <- tryCatch(run_DRDID(cohort_data, covariates, dp2),
                         error = function(e) {
                           warning("\n Error in computing internal 2x2 DiD for (g,t) = (",dp2$treated_groups[g],",",dp2$time_periods[t],"): ", e$message)
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
  time_periods <- dp2$time_periods # tlist

  tlist.length <- ifelse(dp2$base_period != "universal", length(time_periods) - 1, length(time_periods))
  tfac <- ifelse(dp2$base_period != "universal", 1, 0)


  # in terms of indexes, not calendar times
  gt_cells <- expand.grid(g = 1:dp2$treated_groups_count, t = 1:tlist.length, stringsAsFactors = FALSE)
  gt_cells <- gt_cells[order(gt_cells$g, gt_cells$t), ]
  total_gt_iterations <- nrow(gt_cells)


  # Running estimation using run_att_gt_estimation() function
  # Helper function for processing each (g,t) pair
  process_gt <- function(gt_cell, idx) {
    g <- gt_cell$g
    t <- gt_cell$t

    # Run estimation
    gt_result <- run_att_gt_estimation(g, t, dp2)

    # Compute post-treatment indicator
    post.treat <- as.integer(g <= t)

    # Compute post-treatment indicator
    post.treat <- 1*(dp2$treated_groups[g] <= dp2$time_periods[t+tfac])

    if (is.null(gt_result) || is.null(gt_result$att)) {
      # Estimation failed or was skipped
      if(dp2$base_period == "universal"){
        inffunc_updates <- rep(0, n)
        gt_result <- list(att = 0, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat, inffunc_updates = inffunc_updates)
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
      gt_result <- list(att = att, group = dp2$treated_groups[g], year = dp2$time_periods[t+tfac], post = post.treat, inffunc_updates = inffunc_updates)
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
