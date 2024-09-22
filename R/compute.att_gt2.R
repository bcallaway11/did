get_pret <- function(g,t, base_period, anticipation){
  if(base_period == "universal"){
    pret <- g - 1 - anticipation
  } else {
    pret <- ifelse(t >= g, g - 1 - anticipation, t - 1)
  }

  return(pret)
}

# get the cohort for the current g,t values
get_did_cohort_index <- function(group, time, pret, dp){
  # return a vector of dimension id_size with 1, 0 or NA values
  treated_groups <- dp$treated_groups
  # based on control_group option
  min_control_group <- ifelse(dp$control_group == "notyettreated", max(time, pret) + dp$anticipation + 1, Inf)
  max_control_group <- ifelse(dp$control_group == "notyettreated", max(treated_groups), Inf)

  # select the DiD cohort
  did_cohort_index <- rep(NA, dp$id_count)
  # getting the index to get units who will participate in the estimation for the (g,t) cell.
  start_control <- dp$cohort_counts[cohort < min_control_group, sum(cohort_size)]+1
  end_control <- dp$cohort_counts[cohort <= max_control_group, sum(cohort_size)]
  index <- which(dp$cohort_counts[, cohort] == group)
  start_treat <- ifelse(index == 1, 1, dp$cohort_counts[1:(index-1), sum(cohort_size)]+1)
  end_treat <- dp$cohort_counts[1:index, sum(cohort_size)]
  # set the cohort index; .C = 0 and .G = 1
  did_cohort_index[start_control:end_control] <- 0
  did_cohort_index[start_treat:end_treat] <- 1

  return(did_cohort_index)
}

# Wrapper to run DRDID package
run_DRDID <- function(cohort_data, covariates, dp){

  if(dp$panel){
    # --------------------------------------
    # Panel Data
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

    #-----------------------------------------------------------------------------
    # code for actually computing ATT(g,t)
    #-----------------------------------------------------------------------------

    # preparing vectors
    Ypost <- cohort_data[, y1]
    Ypre <- cohort_data[, y0]
    w <- cohort_data[, i.weights]
    G <- cohort_data[, D]

    if (inherits(dp$est_method, "function")) {
      # user-specified function
      attgt <- est_method(y1=Ypost, y0=Ypre,
                          D=G,
                          covariates=covariates,
                          i.weights=w,
                          inffunc=TRUE)
    } else if (dp$est_method == "ipw") {
      # inverse-probability weights
      attgt <- DRDID::std_ipw_did_panel(Ypost, Ypre, G,
                                        covariates=covariates,
                                        i.weights=w,
                                        boot=FALSE, inffunc=TRUE)
    } else if (dp$est_method == "reg") {
      # regression
      attgt <- DRDID::reg_did_panel(Ypost, Ypre, G,
                                    covariates=covariates,
                                    i.weights=w,
                                    boot=FALSE, inffunc=TRUE)
    } else {
      # doubly robust, this is default
      attgt <- DRDID::drdid_panel(Ypost, Ypre, G,
                                  covariates=covariates,
                                  i.weights=w,
                                  boot=FALSE, inffunc=TRUE)
    }

    # adjust influence function to account for only using
    # subgroup to estimate att(g,t)
    attgt$att.inf.func <- (n/n1)*attgt$att.inf.func

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

  return(attgt = attgt)

}


#' @title Compute Group-Time Average Treatment Effects
#'
#' @description `run_att_gt_estimation` does the main work for computing
#'  multiperiod group-time average treatment effects
#' @param gt A numeric vector of the group and time period coming from a list of (g,t) cells
#' @param dp A DIDparams object
#'
#' @return a list with the gt cell and the results after performing estimation
#'
#' @keywords internal
#'
#' @export
run_att_gt_estimation <- function(gt, dp){
  # get the current g,t values
  g <- gt[1]
  t <- gt[2]
  pret <- get_pret(g, t, dp$base_period, dp$anticipation)

  # if we are in period (g-1) or base period out of bounds, normalize results to be equal to NULL
  # and break without computing anything
  if(t == pret | !pret %in% dp$time_periods){
    return(NULL)
  }

  # get treatment and control group
  did_cohort_index <- get_did_cohort_index(g, t, pret, dp)
  # In case of not treatment and control group in the cohort, return NULL
  if(all(is.na(did_cohort_index))){
    return(NULL)
  }

  # Get the matrix of covariates
  covariates <- dp$covariates
  cohort_data <- data.table(did_cohort_index, dp$outcomes_tensor[[t]], dp$outcomes_tensor[[pret]], dp$weights_vector)
  names(cohort_data) <- c("D", "y1", "y0", "i.weights")
  did_result <- run_DRDID(cohort_data, covariates, dp)
  return(did_result)

}

#' @title Compute Group-Time Average Treatment Effects
#'
#' @description `compute.att_gt` does the (g,t) cell and send it to estimation,
#' then do all the post process after estimation
#'
#' @param dp A DIDparams object
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
compute.att_gt <- function(dp) {

  # TODO; do we need the object "data" here if we already have the tensors??
  treated_groups <- dp$treated_groups
  time_periods <- dp$time_periods
  # Get (g,t) cells to perform computations. This replace the nested for-loop.
  gt_cells <- expand.grid(g = treated_groups, t = time_periods, stringsAsFactors = FALSE) |> transpose() |> as.list()
  # Running estimation using run_att_gt_estimation() function
  gt_results <- lapply(gt_cells, run_att_gt_estimation, dp)

  # Post processing creating attgt.list and inffunc

  return(list(attgt.list=attgt.list, inffunc=inffunc))
}
