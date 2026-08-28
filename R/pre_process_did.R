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
  if (missing(control_group)) control_group <- "nevertreated"
  validate_choice_scalar(
    control_group,
    "control_group",
    c("nevertreated", "notyettreated"),
    "control_group must be either 'nevertreated' or 'notyettreated'"
  )
  validate_choice_scalar(
    base_period,
    "base_period",
    c("universal", "varying"),
    "base_period must be either 'universal' or 'varying'."
  )
  validate_logical_scalar(panel, "panel")
  validate_logical_scalar(allow_unbalanced_panel, "allow_unbalanced_panel")
  validate_logical_scalar(bstrap, "bstrap")
  validate_logical_scalar(cband, "cband")
  validate_logical_scalar(faster_mode, "faster_mode")
  validate_logical_scalar(print_details, "print_details")
  validate_logical_scalar(pl, "pl")
  validate_positive_whole_number(cores, "cores")
  validate_anticipation(anticipation)
  validate_alp(alp)
  if (bstrap) validate_positive_whole_number(biters, "biters")
  validate_xformla(xformla)
  # make sure dataset is a data.frame
  # this gets around RStudio's default of reading data as tibble
  if (!all( class(data) == "data.frame")) {
    data <- as.data.frame(data)
  }

  data_names <- names(data)
  validate_column_name(yname, "yname", data_names)
  validate_column_name(tname, "tname", data_names)
  validate_column_name(gname, "gname", data_names)
  validate_column_name(idname, "idname", data_names, allow_null = !panel)
  validate_column_name(weightsname, "weightsname", data_names, allow_null = TRUE)
  validate_column_names(clustervars, "clustervars", data_names, allow_null = TRUE)

  check_reserved_did_names(yname = yname, tname = tname, idname = idname,
                           gname = gname, xformla = xformla,
                           weightsname = weightsname,
                           clustervars = clustervars)

  # validate that all required column names exist in the data
  required_cols <- c(yname, tname, idname, gname, weightsname, clustervars)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("The following column(s) are not found in the data: ",
         paste(missing_cols, collapse = ", "), ". ",
         "Please check the spelling of yname, tname, idname, gname, weightsname, and clustervars.")
  }

  # Strip idname from clustervars before validation: users may naturally pass
  # clustervars = c(idname, extra_var), and idname clustering is implicit/redundant.
  # Mirrors the fast path (pre_process_did2), so the DIDparams object -- and hence
  # mboot() -- carries the same stripped value in both code paths.
  if (!is.null(clustervars) && !is.null(idname) && (idname %in% clustervars)) {
    clustervars <- setdiff(clustervars, idname)
    if (length(clustervars) == 0L) clustervars <- NULL
  }

  # At most one cluster variable beyond idname is supported (clustering at the
  # unit level via idname is implicit). Enforced here -- mirroring the fast path
  # (pre_process_did2) and mboot(), with identical wording -- so every
  # faster_mode x bstrap combination rejects the input up front; the analytical
  # (bstrap = FALSE) path used to silently cluster on the first extra variable only.
  if (length(clustervars) > 1) {
    stop("At most one cluster variable (beyond 'idname') is supported. Please reduce to one.")
  }

  # make sure time periods are numeric
  if (! (is.numeric(data[, tname])) ) stop("The time variable '", tname, "' must be numeric. Please convert it.")

  #  make sure gname is numeric
  if (! (is.numeric(data[, gname])) ) stop("The group variable '", gname, "' must be numeric. Please convert it.")
  # store gname as double (the fast path does the same) so cohort labels in messages
  # are identical across modes for integer input (att_gt()'s group vector was already
  # double on both paths; the exported pre_process_did() now returns double codes too)
  if (is.integer(data[, gname])) data[, gname] <- as.numeric(data[, gname])

  # gname must be 0 (never-treated) or a positive treatment-timing value.
  # Negative codes are not supported: 0 is reserved for never-treated, so a
  # non-positive time scale is ambiguous. Reject them up front so both code
  # paths behave identically (the fast path previously accepted negative codes
  # silently while this slow path later errored with "No valid groups").
  if (any(data[, gname] < 0, na.rm = TRUE)) {
    stop("The group variable '", gname, "' must be 0 (never-treated) or a ",
         "positive treatment-timing value; negative values are not supported. ",
         "If your time periods are non-positive, shift them so the earliest ",
         "period is >= 1.")
  }

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
  # drop rows with any missing or non-finite id / time / outcome / weight /
  # cluster or RAW covariate value. gname is excluded from the finite check
  # because Inf is a valid never-treated code there (see complete_finite_cases);
  # missing/NaN gname is still dropped via complete.cases().
  data <- data[complete_finite_cases(data, finite_exclude = gname), ]
  # also drop rows whose EVALUATED design is missing/non-finite (e.g. log of a non-positive
  # covariate), preserving the previous model.frame-based row dropping. We use
  # model.frame (NOT model.matrix) with na.action = na.pass: model.frame keeps EVERY
  # row -- including those where a term evaluates to NA/NaN -- so complete.cases()
  # flags them and the indicator stays aligned with `data`. (model.matrix would
  # instead silently drop the NaN rows, making the mask shorter than `data` and the
  # offending rows survive.)
  # Safe to evaluate now that raw-covariate NAs have been removed (so poly()/ns()/...
  # will not error on NA input).
  if (length(xvars) > 0L && nrow(data) > 0L) {
    mf_check <- suppressWarnings(model.frame(xformla, data = data, na.action = na.pass))
    finite_rows <- complete_finite_cases(mf_check)
    if (!all(finite_rows)) data <- data[finite_rows, ]
  }
  n_diff <- n_orig - nrow(data)
  if (n_diff != 0) {
    warning(paste0("dropped ", n_diff, " rows from original data due to missing or non-finite data"))
  }
  # Nothing left to estimate on: stop here, on both paths, with the missing-data
  # message (the cohort logic below would otherwise fabricate a "last treated
  # cohort" of -Inf from an empty cohort list).
  if (nrow(data) == 0) {
    if (n_orig == 0) stop("The data has no rows.")
    stop("All observations were dropped due to missing data. Check your outcome, group, time, weights, cluster, and covariate variables for missing or non-finite values.")
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
      # wording is deliberately path-neutral: whether the balanced-panel or the
      # unbalanced (per-observation) path runs is only known further below, once
      # the panel has been checked for balance
      message(
        "Time-varying weights detected. For balanced panel data, the default ",
        "behavior uses the weight from the earlier of the two time periods in ",
        "each 2x2 comparison (the base period for post-treatment cells); for ",
        "unbalanced panel data, each observation carries its own ",
        "period-specific weight. ",
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

  # Fallback applied when no never-treated group is available -- in the raw data here,
  # or in the balanced sample (see the re-check after the balanced-panel coercion
  # below). Periods from the last treated cohort's treatment date (net of
  # anticipation) on are dropped, and that cohort serves only as the comparison group:
  # coerced to never-treated (0) under "nevertreated", left in the data but excluded
  # from glist by the caller under "notyettreated". Returns the filtered data; the
  # caller recomputes tlist/glist and issues any warning.
  no_never_treated_fallback <- function(data, latest_g) {
    cutoff_t <- latest_g - anticipation
    # Drop all periods >= (latest_g - anticipation)
    data <- data[ data[[ tname ]] < cutoff_t, , drop = FALSE ]
    # Nothing left: every period is at or after the cutoff. Reachable with a single
    # cohort treated in the first period (net of anticipation), with gname coded on a
    # different scale than tname (e.g. every unit coded 1 in a gname meant as a 0/1
    # treatment indicator, against calendar years), or in
    # degenerate re-check inputs; stop here with the cause instead of letting the
    # balancing code blame idname/panel. Same text as the fast path (pre_process_did2).
    if (nrow(data) == 0) {
      stop("No valid groups: there is no never-treated group in the estimation sample, and every observed period is at or after the treatment date (net of anticipation) of the last treated cohort (first treated at ", fmt_g(latest_g), "), so no period remains once that cohort is used as the comparison group. Check that 'gname' is coded as the period of first treatment on the same scale as 'tname' (0 for never-treated), and that at least one cohort is observed before it is treated.")
    }
    if (control_group == "nevertreated") {
      # For any row where gname == latest_g, set gname := 0
      lines.gmax <- data[, gname]==latest_g
      data[lines.gmax, gname] <- 0
    }
    data
  }

  # Sentence announcing the fallback (byte-identical in the fast path, pre_process_did2:
  # tests pin "filtered out" and cross-mode identity). The "notyettreated" case is
  # deliberately silent on the raw data: the fallback is the routine design there, so
  # it is announced only when the balanced-panel coercion is what removed the
  # never-treated group (see the re-check below).
  no_nt_text <- function(latest_g) {
    if (control_group == "nevertreated") {
      paste0("The last treated cohort (first treated at ", fmt_g(latest_g), ") is being coerced as 'never-treated' units, and data from periods at or after its treatment date (net of anticipation) is being filtered out (no available comparison groups).")
    } else {
      paste0("The last treated cohort (first treated at ", fmt_g(latest_g), ") is used only as a not-yet-treated comparison group (no ATT(g,t) is computed for it), and data from periods at or after its treatment date (net of anticipation) is being filtered out (no available comparison groups).")
    }
  }

  # Latest treated cohort the fallback was applied to (NA while a never-treated group
  # is available); consulted by the re-check after the balanced-panel coercion below.
  no_nt_latest_g <- NA_real_
  # Treated cohorts removed entirely by the balanced-panel coercion (see below), and
  # whether that left no treated cohort at all.
  gone_cohorts <- numeric(0)
  emptied_by_balancing <- FALSE
  # What the balanced-panel coercion removed (noun phrases), for the "No valid groups"
  # error: at top level R prints the error before the deferred warnings, so the error
  # itself has to name the loss and the remedy.
  bal_lost <- character(0)

  if (!any(glist == 0)) {
    # Compute latest treated cohort once
    latest_g <- max(glist, na.rm = TRUE)
    no_nt_latest_g <- latest_g
    if (control_group == "nevertreated") {
      warning("No never-treated group is available. ", no_nt_text(latest_g))
    }
    data <- no_never_treated_fallback(data, latest_g)

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
    # Drop ONLY the first-period-treated units, by row identity. The previous
    # `data[gname %in% c(0, glist)]` dropped by cohort membership in glist, which --
    # when there is no never-treated group -- also deleted the latest cohort that was
    # deliberately removed from glist above (the `glist[glist < latest_g]` trim) so it
    # could serve as a not-yet-treated control. That silently deleted a valid
    # comparison cohort and corrupted ATT(g,t) for the other groups; the
    # treated_first_period mask removes exactly the already-treated units, nothing else.
    data <- data[ !treated_first_period, , drop = FALSE ]
    # update tlist and glist
    tlist <- sort(unique(data[,tname]))
    glist <- sort(unique(data[,gname]))
    glist <- glist[glist>0]

    # drop groups treated in the first period or before
    first.period <- tlist[1]
    glist <- glist[glist > first.period + anticipation]

    # The latest cohort stays in the data as a not-yet-treated control but, when there
    # is still no never-treated group, must remain excluded from glist (it gets no ATT
    # of its own) -- mirroring the exclusion above for the nfirstperiod == 0 case.
    if (control_group != "nevertreated" && !any(data[,gname] == 0)) {
      glist <- glist[glist < latest_g]
    }

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

      # With user-level panel = FALSE the data are genuine repeated cross sections:
      # every observation is its own sampling unit, so a supplied idname must not
      # repeat in ANY way (within or across periods). `panel` is still the user's
      # argument here -- the allow_unbalanced_panel flip happens further below -- so
      # unbalanced panels (panel = TRUE) are unaffected. Checked before the
      # irreversibility and (idname, tname) scans, mirroring the fast path
      # (pre_process_did2) with identical wording, so repeated cross sections get
      # this message instead of the panel-flavored ones.
      if (!panel && any(same_id)) {
        stop("The value of idname must be unique when panel = FALSE. Repeated cross sections treat each observation as a distinct sampling unit, but some values of '", idname, "' appear in more than one row. If the same units are observed in multiple periods, use panel = TRUE (with allow_unbalanced_panel = TRUE if the panel is incomplete). If the data are genuine repeated cross sections, give each observation its own unique value of '", idname, "' (or omit idname).")
      }

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
    # tell the user which branch was taken, with the same wording as the fast path
    # (pre_process_did2), so the silent reset is visible in both code paths
    message(if (allow_unbalanced_panel)
      "You have an unbalanced panel. Proceeding as such."
      else
        "You have a balanced panel. Setting allow_unbalanced_panel = FALSE.")
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
      # never-treated units going into the coercion (for the re-check's message below)
      n_nt_before <- length(unique(data[[idname]][data[[gname]] == 0]))

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
        stop("Converting to a balanced panel removed every unit: no unit is observed in every period once rows with missing data are removed. Set allow_unbalanced_panel = TRUE to keep units that are not observed in every period, or address the missing data.")
      }

      n <- nrow(data[ data[,tname]==tlist[1], ])

      # Balancing drops units WHOLE, so it can remove every unit of a TREATED cohort
      # (one covariate missing in one period for all of them is enough). glist was
      # fixed before balancing; left stale, this path returns NA for every cell of
      # that cohort with no explanation and the fast path stops with an internal
      # error. Reconcile glist with the balanced sample and say which cohorts were
      # lost. Same logic and text as the fast path (pre_process_did2).
      gone_cohorts <- setdiff(glist, unique(data[, gname]))
      if (length(gone_cohorts) > 0) {
        warning("Converting to a balanced panel removed every unit of the cohort(s) first treated at ", paste(fmt_g(gone_cohorts), collapse = ", "), " (not observed in every period once rows with missing data are removed); no ATT(g,t) is computed for them. To keep these units, set allow_unbalanced_panel = TRUE or address the missing data.")
        glist <- setdiff(glist, gone_cohorts)
        bal_lost <- c(bal_lost, paste0("every unit of the cohort(s) first treated at ", paste(fmt_g(gone_cohorts), collapse = ", ")))
      }
      emptied_by_balancing <- length(gone_cohorts) > 0 && length(glist) == 0

      # Likewise it can remove every never-treated unit, or the whole latest cohort
      # that the fallback above had just turned into the comparison group. The
      # availability check above ran on the pre-balancing data, so nothing noticed:
      # no control units remained and every ATT(g,t) came back NA, reported here as a
      # misleading "overlap condition violated" (a logit of G on X with G all ones).
      # Re-run the same rule on the balanced sample and announce it, saying what was
      # lost and how to keep it.
      gvec_bal <- data[, gname]
      if (!any(gvec_bal == 0) &&
          (is.na(no_nt_latest_g) || max(gvec_bal) < no_nt_latest_g)) {
        # n_nt_before counts the never-treated units left after the missing-data row
        # drop above (units missing in every period are already gone by then)
        lost_what <- if (is.na(no_nt_latest_g)) {
          paste0("all never-treated units (", n_nt_before, " remained after dropping rows with missing data)")
        } else {
          paste0("every unit of the comparison cohort (first treated at ", fmt_g(no_nt_latest_g), ")")
        }
        lost_desc <- paste0(lost_what, if (is.na(no_nt_latest_g)) " were removed" else " was removed")
        bal_lost <- c(bal_lost, lost_what)
        latest_g <- max(gvec_bal)
        no_nt_latest_g <- latest_g
        warning(
          "No never-treated group is available after converting to a balanced panel: ",
          lost_desc, " because none of them is observed in every period. ",
          no_nt_text(latest_g),
          " To keep those units instead, set allow_unbalanced_panel = TRUE (cells with no comparison units in their base or current period will then be NA) or address the missing data."
        )
        # (the helper stops with an informative message if no period remains)
        data <- no_never_treated_fallback(data, latest_g)
        tlist <- sort(unique(data[,tname]))
        glist <- sort(unique(data[,gname]))
        glist <- glist[glist > 0 & glist > tlist[1] + anticipation]
        if (control_group != "nevertreated") {
          glist <- glist[glist < latest_g]
        }
        n <- nrow(data[ data[,tname]==tlist[1], ])
      }

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
    if (emptied_by_balancing) {
      stop("No valid groups: converting to a balanced panel removed every unit of every treated cohort (first treated at ", paste(fmt_g(gone_cohorts), collapse = ", "), "). Set allow_unbalanced_panel = TRUE to keep units that are not observed in every period, or address the missing data.",
           if (!is.na(no_nt_latest_g)) paste0(" (There is no never-treated group, so the last treated cohort, first treated at ", fmt_g(no_nt_latest_g), ", serves only as the comparison group and is not counted.)") else "")
    }
    if (!is.na(no_nt_latest_g)) {
      # the fallback consumed the last cohort that could have had an ATT(g,t); balancing
      # may have removed the others (gone_cohorts), in which case say so here too (at
      # top level R prints the error before the deferred warnings)
      stop("No valid groups: there is no never-treated group in the estimation sample, so the last treated cohort (first treated at ", fmt_g(no_nt_latest_g), ") serves only as the comparison group, and no other treated cohort remains to compute ATT(g,t) for (cohorts already treated in the first period, net of anticipation, are dropped and do not count",
           if (length(bal_lost) > 0) paste0("; converting to a balanced panel removed ", paste(bal_lost, collapse = " and "), " -- set allow_unbalanced_panel = TRUE to keep them, or address the missing data") else "",
           "). At least two treated cohorts observed before treatment, or a never-treated group, are required.")
    }
    stop("No valid groups. The variable in 'gname' should be expressed as the time a unit is first treated (0 if never-treated).")
  }

  # A single remaining period cannot support any 2x2 comparison (reachable only through
  # the fallback above); stop identically on both paths instead of failing downstream.
  if (length(tlist) < 2) {
    stop("Only one time period remains after dropping the periods from the treatment date of the last treated cohort (net of anticipation) onward; at least two are required.")
  }

  # if there are only two time periods, then uniform confidence
  # bands are the same as pointwise confidence intervals
  if (length(tlist)==2) {
    # only announce the override when the user actually asked for a band
    if (cband) {
      message("Only two time periods are available; uniform confidence bands coincide with pointwise confidence intervals. Setting cband = FALSE.")
    }
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
  reqsize <- length(BMisc::rhs_vars(xformla)) + 5

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
