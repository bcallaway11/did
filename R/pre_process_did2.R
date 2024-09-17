#' @title validate_args
#' @description A utility function to validate arguments passed to a att_gt()
#' @param data data.table used in function
#' @param args list of arguments to validate
#'
#' @return nothing, but throws an error if any of the arguments are not valid
#'
#' @export
#' @noRd
validate_args <- function(args, data){

  # ---------------------- Error Checking ----------------------

  # Flag for control group types
  if(!(args$control_group %in% c("nevertreated", "notyettreated"))){
    stop("control_group must be either 'nevertreated' or 'notyettreated'")
  }

  # Flag for yname
  if (!is.element(args$yname, base::colnames(data))) {
    stop("yname = ",yname,  " could not be found in the data provided.")
  }

  # Flag for tname
  if (!is.element(args$tname, base::colnames(data))) {
    stop("tname = ",tname,  " could not be found in the data provided.")
  }

  # check if times periods are numeric
  if (!all(sapply(data[, get(args$tname)], is.integer))) {
    stop("tname = ",tname,  " is not integer. Please convert it")
  }

  # Flag for gname
  if ( !is.element(gname, base::colnames(data))) {
    stop("gname = ",gname,  " could not be found in the data provided.")
  }

  # Check if gname is numeric
  if (!all(sapply(data[, get(args$gname)], is.numeric))) {
    stop("gname = ",gname,  " is not numeric. Please convert it")
  }

  # Flag for idname
  if(!is.null(args$idname)){
    # Check if idname is in the data
    if ( !is.element(args$idname, base::colnames(data))) {
      stop("idname = ",idname,  " could not be found in the data provided.")
    }

    #  check if idname is integer
    if (!all(sapply(data[, get(args$idname)], is.integer))) {
      stop("data[, idname] must be integer Please convert it.")
    }

    # Check if gname is unique by idname: irreversibility of the treatment
    checkTreatmentUniqueness(data, args$idname, args$gname)
  }


  # Check if any combination of idname and tname is duplicated
  n_id_year = anyDuplicated(data[, .(args$idname, args$tname)])
  # If any combination is duplicated, stop execution and throw an error
  if (n_id_year > 0) {
    stop("The value of idname must be unique (by tname)")
  }

  # Flag for weightsname
  if(!is.null(args$weightsname)){
    if (!is.element(args$weightsname, base::colnames(data))) {
      stop("weightsname = ", weightsname, " could not be found in the data provided.")
    }
  }


  # Flag for based period: not in c("universal", "varying"), stop
  if (!args$base_period %in% c("universal", "varying")) {
    stop("base_period must be either 'universal' or 'varying'.")
  }

  # Flags for cluster variable
  if (!is.null(args$clustervars)) {
    # dropping idname from cluster
    if (args$idname %in% args$clustervars) {
      args$clustervars <- setdiff(args$clustervars, args$idname)
    }

    # flag if cluster variable is not in dataset
    if (!is.element(args$cluster, base::colnames(data))) {
      stop("clustervars = ",cluster,  " could not be found in the data provided. Please check your arguments.")
    }

    # check if user is providing more than 2 cluster variables (different than idname)
    if (length(args$cluster) > 1) {
      stop("You can only provide 1 cluster variable additionally to the one provided in idname. Please check your arguments.")
    }

    # Check that cluster variables do not vary over time within each unit
    if (length(cluster) > 0) {
      # Efficiently check for time-varying cluster variables
      clust_tv <- data[, lapply(.SD, function(col) length(unique(col)) == 1), by = id, .SDcols = args$cluster]
      # If any cluster variable varies over time within any unit, stop execution
      if (!all(unlist(clust_tv[, -1, with = FALSE]))) {
        stop("did cannot handle time-varying cluster variables at the moment. Please check your cluster variable.")
      }
    }
  }

  # Check if anticipation is numeric
  if (!is.numeric(args$anticipation)) {
    stop("anticipation must be numeric. Please convert it.")
  }

  # Check if anticipation is positive
  if (args$anticipation < 0) {
    stop("anticipation must be positive. Please check your arguments.")
  }

  # Check if anticipation is integer
  if (!is.integer(args$anticipation)) {
    stop("anticipation must be an integer. Please check your arguments.")
  }

}

did_standarization <- function(data, args){
  # keep relevant columns in data
  cols_to_keep <-  c(args$idname, args$tname, args$gname, args$yname, args$weightsname, args$clustervars)

  model_frame <- model.frame(args$xformla, data = data, na.action = na.pass)
  # Subset the dataset to keep only the relevant columns
  data <- data[, ..cols_to_keep]

  # Column bind the model frame to the data
  data <- cbind(data, as.data.table(model_frame))

  # Check if any covariates were missing
  n_orig <- data[, .N]
  data <- data[complete.cases(data)]
  n_new <- data[, .N]
  n_diff <- n_orig - n_new
  if (n_diff != 0) {
    warning(paste0("dropped ", n_diff, " rows from original data due to missing data"))
  }

  # Set weights
  base::ifelse(is.null(args$weightsname), weights <- rep(1, n_new), weights <- data[[args$weightsname]])
  # enforcing normalization of weights. At this point we already drop any missing values in weights.
  weights <- weights/mean(weights)
  data$weights <- weights

  # get a list of dates from min to max
  tlist <- data[, sort(unique(get(args$tname)))]

  # Coerce control group identified with zero as Inf
  data[get(args$gname) == 0, (args$gname) := Inf]

  # Working out the dates
  # Identify groups with treatment time bigger than the maximum treatment time
  # calculate the maximum treatment time
  max_treatment_time <- max(tlist, na.rm = TRUE)
  data[, asif_never_treated := (get(args$gname) > max_treatment_time)]
  # replace NA values with FALSE in the logical vector
  data[is.na(asif_never_treated), asif_never_treated := FALSE]
  # set gname to 0 for those groups considered as never treated
  data[asif_never_treated == TRUE, (args$gname) := Inf]
  # remove the temporary column
  data[, asif_never_treated := NULL]

  # get list of treated groups (by time) from min to max
  glist <- sort(unique(data[[args$gname]]))

  # Check  if there is a never treated group in the data
  if (!(Inf %in% glist)) {
    if (args$control_group == "nevertreated") {
      stop("There is no available never-treated group")
    } else {
      # Drop all time periods with time periods >= latest treated
      max_treated_time <- max(glist[is.finite(glist)], na.rm = TRUE)
      latest_treated_time <- max_treated_time - args$anticipation
      data <- data[get(args$tname) < latest_treated_time]

      tlist <- sort(unique(data[[args$tname]]))
      glist <- sort(unique(data[[args$gname]]))

      # don't compute ATT(g,t) for groups that are only treated at end
      # and only play a role as a comparison group
      glist <- glist[glist < max_treated_time]
    }
  }

  # get only the first period
  first_period <- tlist[1]

  # drop groups treated in the first period or before; keep only treated groups
  glist <- glist[glist != Inf & glist > first_period + args$anticipation]

  # Check for groups treated in the first period and drop them
  # identify groups treated in the first period
  data[, treated_first_period := (get(gname) <= first_period)]
  data[is.na(treated_first_period), treated_first_period := FALSE]

  # count the number of units treated in the first period
  nfirstperiod <- uniqueN(data[treated_first_period == TRUE, get(args$idname)])

  # handle units treated in the first period
  if (nfirstperiod > 0) {
    warning(paste0("Dropped ", nfirstperiod, " units that were already treated in the first period."))
    data <- data[get(args$gname) %in% c(glist, Inf)]

    # update tlist and glist
    tlist <- data[, sort(unique(get(args$tname)))]
    glist <- data[, sort(unique(get(args$gname)))]

    # Drop groups treated in the first period or before
    first_period <- tlist[1]
    glist <- glist[glist != Inf & glist > first_period + args$anticipation]
  }
  # remove the temporary column
  data[, treated_first_period := NULL]

  # If user specifies repeated cross sections,
  # set that it really is repeated cross sections
  args$true_repeated_cross_sections <- FALSE
  if (!args$panel) {
    args$true_repeated_cross_sections <- TRUE
  }

  #-----------------------------------------------------------------------------
  # setup data in panel case
  #-----------------------------------------------------------------------------
  # Check if data is a balanced panel if panel = TRUE and allow_unbalanced_panel = TRUE
  bal_panel_test <- (args$panel)*(args$allow_unbalanced_panel)

  if (bal_panel_test) {
    # First, focus on complete cases and make a balanced dataset
    data_comp <- data[complete.cases(data)]

    # uniqueN for faster unique counts
    n_all <- uniqueN(data_comp[[args$idname]])

    # Make balanced panel
    data_bal <- BMisc::makeBalancedPanel(data_comp, args$idname, args$tname)
    n_bal <- uniqueN(data_bal[[args$idname]])

    # Determine if the panel is unbalanced
    args$allow_unbalanced_panel <- n_bal < n_all
    message(if (args$allow_unbalanced_panel)
      "You have an unbalanced panel. Proceeding as such."
      else
        "You have a balanced panel. Setting allow_unbalanced_panel = FALSE.")
  }

  if (args$panel) {
    # Check for unbalanced panel
    if (args$allow_unbalanced_panel) {
      # Flag for true repeated cross sections
      args$panel <- FALSE
      args$true_repeated_cross_sections <- FALSE
    } else {
      # Coerce balanced panel

      # Focus on complete cases
      keepers <- complete.cases(data)
      n <- uniqueN(data[[args$idname]])
      n_keep <- uniqueN(data[keepers, ][[args$idname]])
      if (n_keep < n) {
        warning(paste0("Dropped ", (n - n_keep), " observations that had missing data."))
        data <- data[keepers, ]
      }


      # Make balanced panel
      n_old <- uniqueN(data[[args$idname]])
      data <- as.data.table(BMisc::makeBalancedPanel(data, args$idname, args$tname)) # coerce to data.table again. This is ugly, find a better way to do it.
      n <- uniqueN(data[[args$idname]])

      if (n < n_old) {
        warning(paste0("Dropped ", n_old - n, " observations while converting to balanced panel."))
      }

      # If all data is dropped, stop execution
      if (nrow(data) == 0) {
        stop("All observations dropped while converting data to balanced panel. Consider setting `panel = FALSE` and/or revisiting 'idname'.")
      }

      n <- data[get(args$tname) == tlist[1], .N]

      # Ensure The value of gname must be the same across all periods for each particular individual.
      checkTreatmentUniqueness(data, args$idname, args$gname)
    }
  }

  #-----------------------------------------------------------------------------
  # setup data in repeated cross section
  #-----------------------------------------------------------------------------

  if (!args$panel) {
    # Focus on complete cases
    data <- data[complete.cases(data)]

    if (nrow(data) == 0) {
      stop("All observations dropped due to missing data problems.")
    }

    # n-row data.frame to hold the influence function
    if (args$true_repeated_cross_sections) {
      data[, .rowid := .I] # Create row index
      args$idname <- ".rowid"
    } else {
      # Set rowid to idname for repeated cross section/unbalanced
      data[, .rowid := get(args$idname)]
    }

    # Count unique number of cross section observations
    n <- uniqueN(data[[args$idname]])
  }

  # Check if groups is empty (usually a problem with the way people defined groups)
  if(length(glist)==0){
    stop("No valid groups. The variable in 'gname' should be expressed as the time a unit is first treated (0 if never-treated).")
  }

  # if there are only two time periods, then uniform confidence
  # bands are the same as pointwise confidence intervals
  if (length(tlist)==2) {
    args$cband <- FALSE
  }

  #-----------------------------------------------------------------------------
  # more error handling after we have balanced the panel

  # Check against very small groups
  # Calculate group sizes, dividing the count of each group by the length of tlist
  gsize <- data[, .N / length(tlist), by = get(args$gname)]

  # How many in each group before giving a warning
  reqsize <- length(BMisc::rhs.vars(xformla)) + 5

  # Filter groups smaller than reqsize
  gsize <- gsize[V1 < reqsize]

  # Warn if some groups are small
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[[args$gname]], collapse = ",")
    warning(paste0("Be aware that there are some small groups in your dataset.\n  Check groups: ", gpaste, "."))

    # Check if the never treated group is too small
    if (Inf %in% gsize[[gname]] & control_group == "nevertreated") {
      stop("Never treated group is too small, try setting control_group=\"notyettreated\"")
    }
  }
  #-----------------------------------------------------------------------------
  # Sort the data for easy access later on
  setorderv(data, c(args$tname, args$gname, args$idname), c(1,1,1))

  # Assign new args regarding number of time periods and groups
  # How many time periods
  args$time_periods_count <- length(tlist)
  # How many treated groups
  args$treated_groups_count <- length(glist)

  return(list(data = data, args = args))
}

get_did_partitions <- function(data, args){
  return(did_partitions)
}

augment_args <- function(did_partitions, args){
  return(args)
}

#' @title Process `did` Function Arguments
#'
#' @description Function to process arguments passed to the main methods in the
#'  `did` package as well as conducting some tests to make sure
#'  data is in proper format / try to throw helpful error messages.
#'
#' @inheritParams att_gt
#' @param call Function call to att_gt
#'
#' @return a [`DIDparams`] object
#'
#' @export
pre_process_did <- function(yname,
                            tname,
                            idname,
                            gname,
                            xformla = NULL,
                            data,
                            panel = TRUE,
                            allow_unbalanced_panel,
                            control_group = c("nevertreated","notyettreated"),
                            anticipation = 0,
                            weightsname = NULL,
                            alp = 0.05,
                            bstrap = FALSE,
                            cband = FALSE,
                            biters = 1000,
                            clustervars = NULL,
                            est_method = "dr",
                            base_period = "varying",
                            print_details = TRUE,
                            pl = FALSE,
                            cores = 1,
                            call = NULL) {


  # coerce data to data.table first
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  # gathering all the arguments except data
  args_names <- setdiff(names(formlas()), "data")
  args <- mget(arg_names, sys.frame(sys.nframe()))

  # run error checking on arguments
  validate_args(args, data)

  # put in blank xformla if no covariates
  if (is.null(args$xformla)) {
    args$xformla <- ~1
  }

  # Put the data in a standard format after some validation
  cleaned_did <- did_standarization(data, args)

  # Partition staggered did into a 2x2 DiD for faster implementation
  did_partitions <- get_did_partitions(cleaned_did$data, cleaned_did$args)

  # get the augmented arguments useful to perform estimation
  args <- augment_args(did_partitions, cleaned_did$args)

  # store parameters for passing around later

  dp <- DIDparams(did_partitions, args)
  dp <- DIDparams(yname=yname,
                  tname=tname,
                  idname=idname,
                  gname=gname,
                  xformla=xformla,
                  data=as.data.frame(data),
                  control_group=control_group,
                  anticipation=anticipation,
                  weightsname=weightsname,
                  alp=alp,
                  bstrap=bstrap,
                  biters=biters,
                  clustervars=clustervars,
                  cband=cband,
                  print_details=print_details,
                  pl=pl,
                  cores=cores,
                  est_method=est_method,
                  base_period=base_period,
                  panel=panel,
                  true_repeated_cross_sections=true_repeated_cross_sections,
                  n=n,
                  nG=nG,
                  nT=nT,
                  tlist=tlist,
                  glist=glist,
                  call=call)


}
