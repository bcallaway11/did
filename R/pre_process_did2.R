#' @title validate_args
#' @description A utility function to validate arguments passed to a att_gt()
#' @param data data.table used in function
#' @param args list of arguments to validate
#'
#' @return nothing, but throws an error if any of the arguments are not valid
#' @noRd
validate_args <- function(args, data){

  data_names <- names(data)

  # ---------------------- Error Checking ----------------------
  args$control_group <- args$control_group[1]
  # Flag for control group types
  control_group_message <- "control_group must be either 'nevertreated' or 'notyettreated'"
  dreamerr::check_set_arg(args$control_group, "match", .choices = c("nevertreated", "notyettreated"), .message = control_group_message, .up = 1)

  # Flag for tname, gname, yname
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  dreamerr::check_set_arg(args$tname, args$gname, args$yname, "match", .choices = data_names, .message = name_message, .up = 1)

  # Flag for clustervars and weightsname
  checkvar_message <- "__ARG__ must be NULL or a character scalar if a name of columns from the dataset."
  dreamerr::check_set_arg(args$weightsname, args$clustervars, "NULL | match", .choices = data_names, .message = checkvar_message, .up = 1)

  # check if times periods are numeric
  if(!data[, is.numeric(get(args$tname))]){stop("tname = ",args$tname,  " is not numeric. Please convert it")}

  # Check if gname is numeric
  if(!data[, is.numeric(get(args$gname))]){stop("gname = ",args$gname,  " is not numeric. Please convert it")}

  # Flag for idname
  if(!is.null(args$idname)){
    # Check if idname is in the data
    name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
    dreamerr::check_set_arg(args$idname, "match", .choices = data_names, .message = name_message, .up = 1)

    #  check if idname is numeric
    if(!data[, is.numeric(get(args$idname))]){stop("idname = ", args$idname,  " is not numeric. Please convert it")}

    # Check if gname is unique by idname: irreversibility of the treatment
    check_treatment_uniqueness <- data[, .(constant = all(get(args$gname)[1] == get(args$gname))), by = get(args$idname)][, all(constant, na.rm = TRUE)]
    if (!check_treatment_uniqueness) {
      stop("The value of gname (treatment variable) must be the same across all periods for each particular unit. The treatment must be irreversible.")
    }

    # Check if any combination of idname and tname is duplicated
    n_id_year = anyDuplicated(data[, .(args$idname, args$tname)])
    # If any combination is duplicated, stop execution and throw an error
    if (n_id_year > 0) {
      stop("The value of idname must be unique (by tname). Some units are observed more than once in a period.")
    }
  }

  # Flag for based period: not in c("universal", "varying"), stop
  args$base_period <- args$base_period[1]
  base_period_message <- "base_period must be either 'universal' or 'varying'."
  dreamerr::check_set_arg(args$base_period, "match", .choices = c("universal", "varying"), .message = base_period_message, .up = 1)

  # Flags for cluster variable
  if (!is.null(args$clustervars)) {
    # dropping idname from cluster
    if ((!is.null(args$idname)) && (args$idname %in% args$clustervars)){
      args$clustervars <- setdiff(args$clustervars, args$idname)
    }

    # check if user is providing more than 2 cluster variables (different than idname)
    if (length(args$clustervars) > 1) {
      stop("You can only provide 1 cluster variable additionally to the one provided in idname. Please check your arguments.")
    }

    # Check that cluster variables do not vary over time within each unit
    if (length(args$clustervars) > 0) {
      # Efficiently check for time-varying cluster variables
      clust_tv <- data[, lapply(.SD, function(col) length(unique(col)) == 1), by = get(args$idname), .SDcols = args$clustervars]
      # If any cluster variable varies over time within any unit, stop execution
      if (!all(unlist(clust_tv[, -1, with = FALSE]))) {
        stop("did cannot handle time-varying cluster variables at the moment. Please check your cluster variable.")
      }
    }
  }

  # Check if anticipation is numeric using
  if (!is.numeric(args$anticipation)) {
    stop("anticipation must be numeric. Please convert it.")
  }

  # Check if anticipation is positive
  if (args$anticipation < 0) {
    stop("anticipation must be positive. Please check your arguments.")
  }

}

#' @title DiD standarization
#' @description A utility function to coerce the data in standard format
#' @param data data.table used in function
#' @param args list of arguments to validate
#'
#' @return A List, containing order data and arguments
#' @noRd
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
  # if (!(Inf %in% glist)) {
  #   if (args$control_group == "nevertreated") {
  #     stop("There is no available never-treated group")
  #   } else {
  #     # Drop all time periods with time periods >= latest treated
  #     max_treated_time <- max(glist[is.finite(glist)], na.rm = TRUE)
  #     latest_treated_time <- max_treated_time - args$anticipation
  #     data <- data[get(args$tname) < latest_treated_time]
  #
  #     tlist <- sort(unique(data[[args$tname]]))
  #     glist <- sort(unique(data[[args$gname]]))
  #
  #     # don't compute ATT(g,t) for groups that are only treated at end
  #     # and only play a role as a comparison group
  #     glist <- glist[glist < max_treated_time]
  #   }
  # }

  if (!(Inf %in% glist)) {
    # Compute latest treated cohort once, and the cutoff time
    latest_g <- max(glist[is.finite(glist)], na.rm = TRUE)
    cutoff_t  <- latest_g - args$anticipation

    if (args$control_group == "nevertreated") {
      # Warn the user
      warning(
        "No never-treated group is available. ",
        "The last treated cohort is being coerced as 'never-treated' units."
      )

      # Drop all periods ≥ (latest_g - anticipation)
      data <- data[get(args$tname) < cutoff_t]

      # For any row where gname == latest_g, set gname := Inf
      data[, (args$gname) := as.numeric(get(args$gname))] # Convert the column to numeric so Inf can be stored
      data[get(args$gname) == latest_g, (args$gname) := Inf]
    } else {
      # 3. If not "nevertreated", we simply drop those periods and leave gnames alone
      data <- data[get(args$tname) < cutoff_t]
    }

    # Recompute tlist and glist from the filtered/modified data
    tlist <- sort(unique(data[[args$tname]]))
    glist <- sort(unique(data[[args$gname]]))

    # 5. If control_group != "nevertreated", drop the max cohort from glist
    if (args$control_group != "nevertreated") {
      glist <- glist[glist < latest_g]
    }
  }


  # get only the first period
  first_period <- tlist[1]

  # drop groups treated in the first period or before; keep only treated groups
  glist <- glist[glist != Inf & glist > first_period + args$anticipation]

  # Check for groups treated in the first period and drop them
  # identify groups treated in the first period
  data[, treated_first_period := (get(args$gname) <= first_period)]
  data[is.na(treated_first_period), treated_first_period := FALSE]

  # count the number of units treated in the first period
  nfirstperiod <- ifelse(args$panel, uniqueN(data[treated_first_period == TRUE, get(args$idname)]), nrow(data[treated_first_period == TRUE]))

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

    # Check if the panel is balanced
    is_balanced <- check_balance(data_comp, id_col = args$idname, time_col = args$tname)

    args$allow_unbalanced_panel <- !is_balanced
    message(if (!is_balanced)
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
      # n_old <- uniqueN(data[[args$idname]])
      # data <- BMisc::makeBalancedPanel(data, args$idname, args$tname, return_data.table = TRUE) # coerce to data.table again. This is ugly, find a better way to do it.
      # n <- uniqueN(data[[args$idname]])
      #
      # if (n < n_old) {
      #   warning(paste0("Dropped ", n_old - n, " observations while converting to balanced panel."))
      # }
      raw_time_size <- uniqueN(data[[args$tname]])
      unit_count <- data[, .(count = .N), by = get(args$idname)]
      if(any(unit_count[, count < raw_time_size])){
        mis_unit <- unit_count[count < raw_time_size]
        warning(nrow(mis_unit), " units are missing in some periods. Converting to balanced panel by dropping them.")
        data <- data[!get(args$idname) %in% mis_unit$get]
      }

      # If all data is dropped, stop execution
      if (nrow(data) == 0) {
        stop("All observations dropped while converting data to balanced panel. Consider setting `panel = FALSE` and/or revisiting 'idname'.")
      }

      n <- data[get(args$tname) == tlist[1], .N]

      # Check if the treatment is irreversible
      # Ensure The value of gname must be the same across all periods for each particular individual.
      check_treatment_uniqueness <- data[, .(constant = uniqueN(get(args$gname)) == 1), by = get(args$idname)][, all(constant)]
      if (!check_treatment_uniqueness) {
        stop("The value of gname (treatment variable) must be the same across all periods for each particular unit. The treatment must be irreversible.")
      }
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
  reqsize <- length(BMisc::rhs.vars(args$xformla)) + 5

  # Filter groups smaller than reqsize
  gsize <- gsize[V1 < reqsize]

  # Warn if some groups are small
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[,1], collapse = ",")
    warning(paste0("Be aware that there are some small groups in your dataset.\n  Check groups: ", gpaste, "."))

    # Check if the never treated group is too small
    if ((Inf %in% gsize[,1]) & (args$control_group == "nevertreated")) {
      stop("never treated group is too small, try setting control_group=\"notyettreated\"")
    }
  }

  #-----------------------------------------------------------------------------
  # Sort the data for easy access later on
  if (args$panel){
    setorderv(data, c(args$tname, args$gname, args$idname), c(1,1,1))
  } else {
    setorderv(data, c(args$tname, args$gname), c(1,1))
  }

  # Assign new args regarding number of time periods and groups
  # How many time periods
  args$time_periods_count <- length(tlist)
  # time periods
  args$time_periods <- tlist
  # How many treated groups
  args$treated_groups_count <- length(glist)
  # treated groups
  args$treated_groups <- glist
  # id size
  args$id_count <- n

  return(list(data = data, args = args))
}

#' @title DiD tensors
#' @description A utility function that split the data in a list of outcomes tensors and a list of arguments.
#' Tensor are objects of dimension id_count x 1 x time_periods_count and are used for faster filtering in the computation of the DiD estimator.
#' @param data data.table used in function
#' @param args list of arguments to validate
#'
#' @return A List, containing outcomes tensor, time invariant data, cohort counts, covariates matrix, cluster vector and weights vector.
#' @noRd
get_did_tensors <- function(data, args){

  # Getting the outcomes tensor: a vector a outcome variables per time period of dimension id_count x 1 x time_periods_count
  outcomes_tensor <- list()
  #if(!args$allow_unbalanced_panel){
  if(args$panel){
    for(time in seq_along(args$time_periods)){
      # iterating over index
      # creating buckets of outcomes of size from 1 to id_size, id_size+1 to 2*id_size, etc.
      start <- (time - 1) * args$id_count + 1
      end <- time * args$id_count
      outcomes_tensor[[time]] <- data[seq(start,end), get(args$yname)]
    }
  } else {
    # for(time in args$time_periods){
    #   outcome_vector_time <- rep(NA, args$id_count)  # Initialize vector with NAs
    #   # Create a new column with Y values if tname == current_time, otherwise NA
    #   data[, outcome_vector_time := fifelse(get(args$tname) == time, get(args$yname), NA_real_)]
    #   # Store the vector in the outcomes_tensor
    #   outcomes_tensor[[time]] <- data$outcome_vector_time
    #   # drop the column created
    #   data[, outcome_vector_time := NULL]
    # }
    outcomes_tensor <- NULL
  }

  # Getting the time invariant data
  time_invariant_cols <-  c(args$idname, args$gname, "weights", args$clustervars)
  #if(!args$allow_unbalanced_panel){
  if(args$panel){
    # We can do this filtering because the data is already sorted appropriately
    invariant_data <- data[1:args$id_count]
  } else {
    if (!args$allow_unbalanced_panel) {
      # get the first observation for each unit
      invariant_data <- data[data[, .I[1], by = get(args$idname)]$V1]
      # # order by idname
      # setorderv(invariant_data, c(args$idname), c(1))
    } else {
      invariant_data <- copy(data)
    }
  }

  # Get cohort counts
  cohorts <- c(args$treated_groups, Inf)
  cohort_counts <- invariant_data[, .(cohort_size = .N) , by = get(args$gname)]
  names(cohort_counts)[1] <- "cohort" # changing the name

  # Get period counts
  period_counts <- invariant_data[, .(period_size = .N), by = get(args$tname)]
  names(period_counts)[1] <- "period" # changing the name

  # Get crosstable counts
  crosstable_counts <- invariant_data[, .N, by = .(get(args$tname), get(args$gname))]   # Directly count occurrences by tname and gname
  names(crosstable_counts)[1] <- "period" # changing the name
  names(crosstable_counts)[2] <- "cohort" # changing the name
  crosstable_counts <- dcast(crosstable_counts, period ~ cohort, value.var = "N", fill = 0)  # Reshape


  # # Get covariates if any
  # if(args$xformla == ~1){
  #   covariates <- rep(1, args$id_count)
  # } else {
  #   covariates <- as.data.table(model.matrix(args$xformla, data = invariant_data, na.action = na.pass))
  # }

  # if (args$xformla == ~1) {
  #   covariates_tensor <- list(rep(1, args$id_count))
  # } else {
  #   if (args$panel){
  #     covariates_tensor <- vector("list", args$time_periods_count)
  #     for (tt in seq_along(args$time_periods)) {
  #       rng <- ((tt - 1) * args$id_count + 1):(tt * args$id_count)
  #       covariates_tensor[[tt]] <-
  #         model.matrix(args$xformla, data = data[rng], na.action = na.pass)
  #     }
  #   } else {
  #     covariates_tensor <- as.data.table(model.matrix(args$xformla, data = invariant_data, na.action = na.pass))
  #   }
  #
  # }
  if (args$xformla == ~1) {
    covariates_tensor <- replicate(args$time_periods_count,
                                   rep(1, args$id_count),
                                   simplify = FALSE)
    covariates_matrix <- matrix(1, nrow = nrow(invariant_data), ncol = 1)
  } else {
    if (args$panel) {
      covariates_tensor <- vector("list", args$time_periods_count)
      for (tt in seq_along(args$time_periods)) {
        rng <- ((tt - 1) * args$id_count + 1):(tt * args$id_count)
        covariates_tensor[[tt]] <-
          model.matrix(args$xformla, data = data[rng], na.action = na.pass)
      }
      covariates_matrix <- NULL        # not needed in balanced panel
    } else {                             # RCS / unbalanced
      covariates_tensor <- NULL        # slices not needed
      covariates_matrix <- model.matrix(args$xformla,
                                        data = invariant_data,
                                        na.action = na.pass)
    }
  }

  # Get the cluster variable only
  if(!is.null(args$clustervars)){
    cluster <- invariant_data[, .SD, .SDcols = args$clustervars] |> unlist()
  } else {
    cluster <- NA
  }

  # Get the weights only
  weights <- invariant_data[, .SD, .SDcols = "weights"] |> unlist()

  # Gather all the arguments to return
  return(list(outcomes_tensor = outcomes_tensor,
              data = data,
              time_invariant_data = invariant_data,
              cohort_counts = cohort_counts,
              period_counts = period_counts,
              crosstable_counts = crosstable_counts,
              # covariates = covariates,
              covariates_matrix  = covariates_matrix,
              covariates_tensor = covariates_tensor,
              cluster = cluster,
              weights = weights))
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
pre_process_did2 <- function(yname,
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
                            faster_mode=FALSE,
                            pl = FALSE,
                            cores = 1,
                            call = NULL) {


  # coerce data to data.table first
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  # gathering all the arguments except data
  args_names <- setdiff(names(formals()), "data")
  args <- mget(args_names, sys.frame(sys.nframe()))

  # pick a control_group by default
  args$control_group <- control_group[1]

  # run error checking on arguments
  validate_args(args, data)

  # put in blank xformla if no covariates
  if (is.null(args$xformla)) {
    args$xformla <- ~1
  } else {
    # extract variable names from the formula
    formula_vars <- all.vars(args$xformla)

    # identify variables in xformla not in data
    missing_vars <- setdiff(formula_vars, names(data))

    # error checking for missing variables in data
    if (length(missing_vars) > 0) {
      stop(paste("The following variables are not in data:", paste(missing_vars, collapse = ", ")), call. = FALSE)
    }
  }

  # Put the data in a standard format after some validation
  cleaned_did <- did_standarization(data, args)

  # Partition staggered DiD into a 2x2 DiD for faster implementation
  did_tensors <- get_did_tensors(cleaned_did$data, cleaned_did$args)

  # store parameters for passing around later
  dp <- DIDparams2(did_tensors, cleaned_did$args, call=call)

  return(dp)
}
