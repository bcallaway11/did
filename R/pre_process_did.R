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
                            fix_weights = NULL,
                            alp = 0.05,
                            bstrap = FALSE,
                            cband = FALSE,
                            biters = 1000,
                            clustervars = NULL,
                            est_method = "dr",
                            base_period = "varying",
                            print_details = TRUE,
                            faster_mode = FALSE,
                            pl = FALSE,
                            cores = 1,
                            call = NULL) {
  #-----------------------------------------------------------------------------
  # Data pre-processing and error checking
  #-----------------------------------------------------------------------------
  # set control group
  control_group <- control_group[1]
  if(!(control_group %in% c("nevertreated","notyettreated"))){
    stop("control_group must be either 'nevertreated' or 'notyettreated'")
  }
  base_period <- base_period[1]
  if (!(base_period %in% c("universal", "varying"))) {
    stop("base_period must be either 'universal' or 'varying'.")
  }
  # Check if anticipation is numeric and non-negative (same contract as the fast path)
  if (!is.numeric(anticipation)) {
    stop("anticipation must be numeric. Please convert it.")
  }
  if (anticipation < 0) {
    stop("anticipation must be non-negative. Please check your arguments.")
  }
  check_reserved_did_names(yname = yname, tname = tname, idname = idname,
                           gname = gname, xformla = xformla,
                           weightsname = weightsname,
                           clustervars = clustervars)
  # make sure dataset is a data.frame
  # this gets around RStudio's default of reading data as tibble
  if (!all( class(data) == "data.frame")) {
    data <- as.data.frame(data)
  }

  # validate that all required column names exist in the data
  required_cols <- c(yname, tname, idname, gname, weightsname, clustervars)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("The following column(s) are not found in the data: ",
         paste(missing_cols, collapse = ", "), ". ",
         "Please check the spelling of yname, tname, idname, gname, weightsname, and clustervars.")
  }

  # At most one cluster variable beyond idname is supported (clustering at the
  # unit level via idname is implicit). Enforced here -- mirroring the fast path
  # (pre_process_did2) and mboot(), with identical wording -- so every
  # faster_mode x bstrap combination rejects the input up front; the analytical
  # (bstrap = FALSE) path used to silently cluster on the first extra variable only.
  cv_check <- clustervars
  if (!is.null(cv_check) && !is.null(idname) && (idname %in% cv_check)) {
    cv_check <- setdiff(cv_check, idname)
  }
  if (length(cv_check) > 1) {
    stop("At most one cluster variable (beyond 'idname') is supported. Please reduce to one.")
  }

  # make sure time periods are numeric
  if (! (is.numeric(data[, tname])) ) stop("The time variable '", tname, "' must be numeric. Please convert it.")

  #  make sure gname is numeric
  if (! (is.numeric(data[, gname])) ) stop("The group variable '", gname, "' must be numeric. Please convert it.")

  #  make sure the outcome is numeric (logical 0/1 outcomes are also allowed)
  if (! (is.numeric(data[, yname]) || is.logical(data[, yname])) ) stop("The outcome variable '", yname, "' must be numeric. Please convert it.")

  # put in blank xformla if no covariates or check whether all variables are in data
  if (is.null(xformla)) {
    xformla <- ~1
  } else {
    # extract variable names from the formula
    formula_vars <- all.vars(xformla)

    # identify variables in xformla not in data
    missing_vars <- setdiff(formula_vars, names(data))

    # error checking for missing variables in data
    if (length(missing_vars) > 0) {
      stop(paste("The following variables are not in data:", paste(missing_vars, collapse = ", ")), call. = FALSE)
    }
  }

  # drop irrelevant columns from data. Keep the RAW covariate variables (all.vars of
  # xformla), NOT the evaluated model.frame, so model.matrix(xformla, .) can be
  # rebuilt downstream. This is what makes transform formulae work (e.g. ~I(X^2),
  # ~poly(X, 2), ~log(X)): the evaluated model.frame would store columns named by the
  # transformed expression -- losing the raw variable that model.matrix needs -- and,
  # for matrix-valued terms like poly(), a matrix column that later breaks the
  # data.table coercion in compute.att_gt(). For bare-variable formulae (e.g. ~X,
  # ~X1+X2) and factor covariates the kept columns are identical to before.
  xvars <- all.vars(xformla)
  keep_cols <- unique(c(idname, tname, yname, gname, weightsname, clustervars, xvars))
  data <- data[, keep_cols, drop = FALSE]

  # check if any covariates were missing
  n_orig <- nrow(data)
  # drop rows with any missing id / time / outcome / group / weight / cluster or any
  # missing RAW covariate value
  data <- data[complete.cases(data), ]
  # also drop rows whose EVALUATED design is non-finite (e.g. log of a non-positive
  # covariate), preserving the previous model.frame-based row dropping. We use
  # model.frame (NOT model.matrix) with na.action = na.pass: model.frame keeps EVERY
  # row -- including those where a term evaluates to NA/NaN -- so complete.cases()
  # flags them and the indicator stays aligned with `data`. (model.matrix would
  # instead silently drop the NaN rows, making the mask shorter than `data` and the
  # offending rows survive.) Inf-valued terms are kept, matching the prior behavior.
  # Safe to evaluate now that raw-covariate NAs have been removed (so poly()/ns()/...
  # will not error on NA input).
  if (length(xvars) > 0L && nrow(data) > 0L) {
    mf_check <- suppressWarnings(model.frame(xformla, data = data, na.action = na.pass))
    finite_rows <- complete.cases(mf_check)
    if (!all(finite_rows)) data <- data[finite_rows, ]
  }
  n_diff <- n_orig - nrow(data)
  if (n_diff != 0) {
    warning(paste0("dropped ", n_diff, " rows from original data due to missing data"))
  }

  # weights if null
  if (is.null(weightsname)) w <- rep(1, nrow(data)) else w <- data[,weightsname]

  # Validate user-supplied weights: negative weights flip signs and a non-positive
  # mean divides by ~0 during normalization, silently producing NA/NaN ATTs.
  if (!is.null(weightsname) && (any(w < 0, na.rm = TRUE) || isTRUE(mean(w, na.rm = TRUE) <= 0)))
    stop("The weights variable '", weightsname, "' must be non-negative with a positive mean.")

  if (".w" %in% colnames(data)) stop("Your data already contains a column named '.w', which is reserved for internal use by `did`. Please rename this column before calling att_gt().")
  data$.w <- w

  # Check for time-varying weights in panel data. Grouped max/min via data.table
  # (GForce-optimized) is much faster than tapply() with a per-unit closure and
  # produces the same per-unit ranges as diff(range(x)) on the raw weights.
  if (!is.null(weightsname) && panel) {
    dtw <- data.table(id = data[[idname]], .w = data[[weightsname]])
    w_rng <- dtw[, .(mx = max(.w), mn = min(.w)), by = "id"]
    if (any((w_rng$mx - w_rng$mn) > .Machine$double.eps^0.5, na.rm = TRUE)) {
      message(
        "Time-varying weights detected. For balanced panel data, the default ",
        "behavior uses the weight from the earlier of the two time periods in ",
        "each 2x2 comparison (the base period for post-treatment cells). ",
        "Use the 'fix_weights' argument to control this behavior. ",
        "See ?att_gt for details."
      )
    }
  }

  # Outcome variable will be denoted by y
  # data$.y <- data[, yname]

  # figure out the dates
  # list of dates from smallest to largest
  tlist <- sort(unique(data[,tname]))

  # Groups with treatment time bigger than max time period + anticipation are considered to be never treated.
  # We account for anticipation because units treated shortly after the last observed period
  # may already exhibit anticipatory effects during the observed time window, and should
  # not be used as controls.
  asif_never_treated <- (data[,gname] > max(tlist, na.rm = TRUE) + anticipation)
  asif_never_treated[is.na(asif_never_treated)] <- FALSE
  data[asif_never_treated, gname] <- 0

  # list of treated groups (by time) from smallest to largest
  glist <- sort(unique(data[,gname]))


  # Check if there is a never treated group

  # if ( length(glist[glist==0]) == 0) {
  #   if(control_group=="nevertreated"){
  #     stop("There is no available never-treated group")
  #   } else {
  #     # Drop all time periods with time periods >= latest treated
  #     data <- subset(data,(data[,tname] < (max(glist)-anticipation)))
  #     # Replace last treated time with zero
  #     # lines.gmax <- data[,gname]==max(glist, na.rm = TRUE)
  #     # data[lines.gmax,gname] <- 0
  #
  #     tlist <- sort(unique(data[,tname]))
  #     glist <- sort(unique(data[,gname]))
  #
  #     # don't comput ATT(g,t) for groups that are only treated at end
  #     # and only play a role as a comparison group
  #     glist <- glist[ glist < max(glist)]
  #   }
  # }

  if (!any(glist == 0)) {
    # Compute latest treated cohort once, and the cutoff time
    latest_g <- max(glist, na.rm = TRUE)
    cutoff_t  <- latest_g - anticipation

    if (control_group == "nevertreated") {
      # Warn the user
      warning(
        "No never-treated group is available. ",
        "The last treated cohort is being coerced as 'never-treated' units, and data from periods after that is being filtered out (no available comparison groups)."
      )

      # Drop all periods ≥ (latest_g - anticipation)
      data <- data[ data[[ tname ]] < cutoff_t, , drop = FALSE ]

      # For any row where gname == latest_g, set gname := 0
      lines.gmax <- data[, gname]==latest_g
      data[lines.gmax, gname] <- 0
    } else {
      # If "notyettreated", we simply drop those periods and leave gnames alone
      data <- data[ data[[ tname ]] < cutoff_t, , drop = FALSE ]
    }

    # Recompute tlist and glist from the filtered/modified data
    tlist <- sort(unique(data[,tname]))
    glist <- sort(unique(data[,gname]))


    # If control_group != "nevertreated", drop the max cohort from glist
    if (control_group != "nevertreated") {
      glist <- glist[glist < latest_g]
    }
  }

  # Only the treated groups
  glist <- glist[glist>0]

  # drop groups treated in the first period or before
  first.period <- tlist[1]
  glist <- glist[glist > first.period + anticipation]

  # check for groups treated in the first period (accounting for anticipation) and drop these
  # nfirstperiod <- length(unique(data[ !((data[,gname] > first.period) | (data[,gname]==0)), ] )[,idname])
  treated_first_period <- ( data[,gname] <= first.period + anticipation ) & ( !(data[,gname]==0) )
  treated_first_period[is.na(treated_first_period)] <- FALSE
  # if/else (not ifelse) so the data subset is built only for the relevant branch;
  # the result is identical to the previous ifelse(panel, ...) expression.
  nfirstperiod <- if (panel) length(unique(data[treated_first_period, idname])) else sum(treated_first_period)
  if ( nfirstperiod > 0 ) {
    warning(paste0("Dropped ", nfirstperiod, " units that were already treated in the first period",
                    if (anticipation > 0) paste0(" (accounting for anticipation = ", anticipation, ")") else "",
                    "."))
    data <- data[ data[,gname] %in% c(0,glist), ]
    # update tlist and glist
    tlist <- sort(unique(data[,tname]))
    glist <- sort(unique(data[,gname]))
    glist <- glist[glist>0]

    # drop groups treated in the first period or before
    first.period <- tlist[1]
    glist <- glist[glist > first.period + anticipation]

  }

  #  make sure id is numeric
  if (! is.null(idname)){
    #  make sure id is numeric
    if (! (is.numeric(data[, idname])) ) stop("The id variable '", idname, "' must be numeric. Please convert it.")

    # Validate treatment irreversibility and (id, period) uniqueness with a single
    # radix order over (id, t) plus O(n) adjacent-element scans, instead of the much
    # slower unique.data.frame()/anyDuplicated.data.frame() row-key machinery. After
    # ordering, each unit's rows are contiguous, so g varying within a unit shows up
    # as adjacent rows with equal id but differing g, and a duplicated (id, period)
    # pair as adjacent rows with equal id and equal t. All three columns are numeric
    # and NA-free at this point, so the comparisons are exact.
    idv <- data[, idname]
    nn <- length(idv)
    if (nn >= 2L) {
      tv <- data[, tname]
      o <- order(idv, tv, method = "radix")
      io <- idv[o]
      same_id <- io[-1L] == io[-nn]

      # Check that gname is time-invariant within each unit (treatment irreversibility).
      go <- data[, gname][o]
      if (any(same_id & (go[-1L] != go[-nn]))) {
        stop("The value of gname (treatment variable) must be the same across all periods for each particular unit. The treatment must be irreversible.")
      }

      # Check that (idname, tname) is unique: each unit observed at most once per
      # period. Mirrors the fast path (pre_process_did2.R) so both code paths reject
      # duplicated (id, period) rows identically -- without this guard the slow path
      # silently produced incorrect estimates on long-format data with duplicates.
      to <- tv[o]
      if (any(same_id & (to[-1L] == to[-nn]))) {
        stop("The value of idname must be unique (by tname). Some units are observed more than once in a period.")
      }

      # Check that cluster variables are time-invariant within each unit, mirroring
      # the fast path (validate_args) and mboot() -- with identical wording -- so
      # invalid clustering inputs are rejected up front regardless of bstrap.
      # Without this, the analytical (bstrap = FALSE) path fell back to i.i.d.
      # standard errors with a warning advising bstrap = TRUE, advice that then
      # fails in mboot() for the very same input. Reuses the (id, t) radix order
      # and adjacency scan above.
      for (cvar in setdiff(clustervars, idname)) {
        cvo <- data[o, cvar]
        if (any(same_id & (cvo[-1L] != cvo[-nn]))) {
          stop("Time-varying cluster variables are not supported. Please provide a time-invariant cluster variable.")
        }
      }
    }
  }



  # if user specifies repeated cross sections,
  # set that it really is repeated cross sections
  true_repeated_cross_sections <- FALSE
  if (!panel) {
    true_repeated_cross_sections <- TRUE
  }

  #-----------------------------------------------------------------------------
  # setup data in panel case
  #-----------------------------------------------------------------------------
  # Check if data is a balanced panel if panel = TRUE and allow_unbalanced_panel = TRUE
  bal_panel_test <- panel*allow_unbalanced_panel
  if (bal_panel_test) {
    # data is already complete-case filtered above and (id, period) uniqueness has
    # been validated, so the panel is balanced iff every unit appears in every
    # period, i.e. nrow(data) equals (number of units) x (number of periods). This
    # replaces a BMisc::makeBalancedPanel() round-trip that built a balanced copy
    # of the data just to compare unit counts and then threw it away.
    allow_unbalanced_panel <-
      nrow(data) != length(unique(data[[idname]])) * as.numeric(length(unique(data[[tname]])))
  }



  if (panel) {

    # check for unbalanced panel
    if (allow_unbalanced_panel) {

      # Flag for true repeated cross sections
      panel <- FALSE
      true_repeated_cross_sections <- FALSE

    } else {

      # this is the case where we coerce balanced panel

      # check for complete cases (rows with missing data were already dropped
      # above, so anyNA() short-circuits the redundant complete.cases() pass and
      # full-table subset copies in the common case)
      if (anyNA(data)) {
        keepers <- complete.cases(data)
        n <- length(unique(data[,idname]))
        n.keep <- length(unique(data[keepers,idname]))
        if (!all(keepers)) {
          warning(paste0("Dropped ", (n-n.keep), " observations that had missing data."))
          data <- data[keepers,]
        }
      }

      # make it a balanced data set: keep only units observed in every period.
      # A unit is fully observed iff its row count equals the number of distinct
      # periods ((id, period) uniqueness was validated above), so this
      # count-and-filter keeps exactly the same rows as BMisc::makeBalancedPanel()
      # without the per-group .SD materialization; any row-order difference is
      # normalized by the (id, period) sort below.
      uid <- unique(data[[idname]])
      n.old <- length(uid)
      cnt <- tabulate(match(data[[idname]], uid))
      keep_bal <- data[[idname]] %in% uid[cnt == length(unique(data[[tname]]))]
      if (!all(keep_bal)) data <- data[keep_bal, , drop = FALSE]
      n <- length(unique(data[,idname]))
      if (n < n.old) {
        # n.old - n counts dropped UNITS (unique ids), not unit-time rows; word the
        # warning accordingly, with the same text as the fast path (pre_process_did2)
        warning(n.old - n, " units are missing in some periods. Converting to balanced panel by dropping them.")
      }

      # If drop all data, you do not have a panel.
      if (nrow(data)==0) {
        stop("All observations dropped while converting data to balanced panel. Consider setting `panel = FALSE` and/or revisiting 'idname'.")
      }

      n <- nrow(data[ data[,tname]==tlist[1], ])

      # slow, repeated check here...
      ## # check that first.treat doesn't change across periods for particular individuals
      ## if (!all(sapply( split(data, data[,idname]), function(df) {
      ##   length(unique(df[,gname]))==1
      ## }))) {
      ##   stop("The value of gname must be the same across all periods for each particular individual.")
      ## }

    }
  }

  #-----------------------------------------------------------------------------
  # code for setting up repeated cross sections (and unbalanced panel)
  #-----------------------------------------------------------------------------
  if (!panel) {

    # check for complete cases (rows with missing data were already dropped
    # above, so anyNA() short-circuits the redundant complete.cases() pass and
    # full-table subset copies in the common case)
    if (anyNA(data)) {
      keepers <- complete.cases(data)
      if (!all(keepers)) {
        warning(paste0("Dropped ", sum(!keepers), " observations that had missing data."))
        data <- data[keepers,]
      }
    }

    # If drop all data, you do not have a panel.
    if (nrow(data)==0) {
      stop("All observations were dropped due to missing data. Check your outcome, group, time, and covariate variables for missing values.")
    }

    # n-row data.frame to hold the influence function
    if (true_repeated_cross_sections) {
      data$.rowid <- seq(1:nrow(data))
      idname <- ".rowid"
    } else {
      # set rowid to idname for repeated cross section/unbalanced
      data$.rowid <- data[, idname]
    }

    # n is unique number of cross section observations
    # this is different for repeated cross sections and unbalanced panel
    n <- length(unique(data[,idname]))
  }

  ## # Update tlist and glist because of data handling
  ## # figure out the dates
  ## # list of dates from smallest to largest
  ## tlist <- unique(data[,tname])[order(unique(data[,tname]))]
  ## # list of treated groups (by time) from smallest to largest
  ## glist <- unique(data[,gname])[order(unique(data[,gname]))]

  ## # Only the treated groups
  ## glist <- glist[glist>0]

  ## # drop groups treated in the first period or before
  ## first.period <- tlist[1]
  ## glist <- glist[glist > first.period + anticipation]

  # Check if groups is empty (usually a problem with the way people defined groups)
  if(length(glist)==0){
    stop("No valid groups. The variable in 'gname' should be expressed as the time a unit is first treated (0 if never-treated).")
  }

  # if there are only two time periods, then uniform confidence
  # bands are the same as pointwise confidence intervals
  if (length(tlist)==2) {
    cband <- FALSE
  }

  #-----------------------------------------------------------------------------
  # more error handling after we have balanced the panel

  # check against very small groups. tabulate(match(.)) yields the same per-group
  # row counts as the previous aggregate() call without invoking an R closure per
  # group; sorted gvals reproduces aggregate()'s ascending group order, and the
  # Group.1/x column names are kept for the subset()/paste logic below.
  gvals <- sort(unique(data[[gname]]))
  gcnt <- tabulate(match(data[[gname]], gvals), nbins = length(gvals))
  gsize <- data.frame(Group.1 = gvals, x = gcnt / length(tlist))

  # how many in each group before give warning
  # 5 is just a buffer, could pick something else, but seems to work fine
  reqsize <- length(BMisc::rhs.vars(xformla)) + 5

  # which groups to warn about
  gsize <- subset(gsize, x < reqsize) # x is name of column from aggregate

  # warn if some groups are small
  if (nrow(gsize) > 0) {
    gpaste <-  paste(gsize[,1], collapse=",")
    warning(paste0("Some groups in your dataset have very few observations, which may cause estimation problems.\n  Check groups: ", gpaste, "."))

    if ( (0 %in% gsize[,1]) & (control_group == "nevertreated") ) {
      stop("The never-treated group is too small to serve as a reliable control. Try setting `control_group = 'notyettreated'` to include not-yet-treated units as controls.")
    }
  }
  #----------------------------------------------------------------------------

  # How many time periods
  nT <- length(tlist)
  # How many treated groups
  nG <- length(glist)

  # order dataset wrt idname and tname, in place. The sort keys are tie-free
  # ((id, period) uniqueness was validated above; .rowid is 1:n for repeated cross
  # sections), so the permutation is unique and the row order matches the previous
  # order()-based subset without a full-table copy. `data` is function-local (it
  # was subset/copied above), so the by-reference conversion cannot touch the
  # user's input; setDF() returns a plain data.frame for downstream consumers
  # (DIDparams, compute.att_gt).
  setDT(data)
  setorderv(data, c(idname, tname))
  setDF(data)

  # store parameters for passing around later
  dp <- DIDparams(yname=yname,
                  tname=tname,
                  idname=idname,
                  gname=gname,
                  xformla=xformla,
                  data=data,
                  control_group=control_group,
                  anticipation=anticipation,
                  weightsname=weightsname,
                  fix_weights=fix_weights,
                  alp=alp,
                  bstrap=bstrap,
                  biters=biters,
                  clustervars=clustervars,
                  cband=cband,
                  print_details=print_details,
                  faster_mode=faster_mode,
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
