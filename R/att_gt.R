#' @title Group-Time Average Treatment Effects
#'
#' @description `att_gt` computes average treatment effects in DID
#'  setups where there are more than two periods of data and allowing for
#'  treatment to occur at different points in time and allowing for
#'  treatment effect heterogeneity and dynamics.
#'  See Callaway and Sant'Anna (2021) for a detailed description.
#'
#' @param yname The name of the outcome variable
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param idname The individual (cross-sectional unit) id name
#' @param gname The name of the variable in `data` that
#'  contains the first period when a particular observation is treated.
#'  This should be a positive number for all observations in treated groups.
#'  It defines which "group" a unit belongs to.  It should be 0 for units
#'  in the untreated group.
#' @param weightsname The name of the column containing the sampling weights.
#'  If not set, all observations have same weight. When weights are
#'  time-invariant (constant within each unit across periods), all
#'  \code{fix_weights} options produce identical results and no special
#'  handling is needed.
#'
#'  When weights vary across time (e.g., time-varying population sizes),
#'  the default behavior differs by panel type:
#'  \describe{
#'    \item{Balanced panel}{Each 2x2 DiD comparison uses the weight from the
#'      earlier of the two time periods involved. For post-treatment cells,
#'      this is the base period (g-1). For pre-treatment cells with
#'      \code{base_period="varying"}, this is the pre-treatment period itself.
#'      The panel DRDID estimators are used.}
#'    \item{Repeated cross sections and unbalanced panels}{Both periods'
#'      per-observation weights are passed directly to the RC DRDID estimators,
#'      so each observation carries its own period-specific weight.}
#'  }
#'  Use the \code{fix_weights} argument to override the default behavior.
#' @param fix_weights Controls how time-varying sampling weights are resolved.
#'  Only relevant when weights vary across time; with time-invariant weights,
#'  all options produce identical results. Options:
#'  \describe{
#'    \item{\code{NULL} (default)}{For balanced panel: uses the weight from
#'      the earlier of the two time periods in each 2x2 comparison. For
#'      post-treatment cells, this is the base period (g-1). For
#'      pre-treatment cells, this depends on the \code{base_period} setting.
#'      For RC/unbalanced panel: uses per-observation weights from both
#'      periods.}
#'    \item{\code{"varying"}}{Uses per-observation, period-specific weights
#'      for all panel types. For balanced panel data, this switches to the
#'      repeated cross-section DRDID estimators so that pre-period and
#'      post-period observations each carry their own weight. Covariates
#'      are held fixed at their pre-period values (same as the default
#'      panel estimator). This is the most flexible option for weights but
#'      sacrifices the efficiency of the panel estimator. For RC/unbalanced
#'      panel, this is identical to the default. Not supported with custom
#'      \code{est_method} functions.}
#'    \item{\code{"base_period"}}{Fixes weights at the base period (g-1) for
#'      all (g,t) cells within a group, for both pre-treatment and
#'      post-treatment comparisons. Ensures all cells within a group use the
#'      same weights. For unbalanced panels, units not observed in the base
#'      period are dropped with a warning. Not supported for repeated cross
#'      sections (\code{panel = FALSE}).}
#'    \item{\code{"first_period"}}{Fixes weights at the first time period in
#'      the dataset for all (g,t) cells. For unbalanced panels, units not
#'      observed in the first period are dropped with a warning. Not supported
#'      for repeated cross sections (\code{panel = FALSE}).}
#'  }
#' @param alp the significance level, default is 0.05
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier bootstrap.  Default is `TRUE` (in addition, cband
#'  is also by default `TRUE` indicating that uniform confidence bands
#'  will be returned).  If `bstrap=FALSE`, analytical standard errors are
#'  reported; these are cluster-robust when `clustervars` is supplied.
#' @param biters The number of bootstrap iterations to use.  The default is 1000,
#'  and this is only applicable if `bstrap=TRUE`.
#' @param clustervars A vector of variables names to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level. Clustered standard errors are available with the multiplier bootstrap
#'  (`bstrap=TRUE`) or analytically (`bstrap=FALSE`).
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability `1-alp`.  In order to compute uniform confidence
#'  bands, `bstrap` must also be set to `TRUE`.  The default is
#' `TRUE`.
#' @param print_details Whether or not to show details/progress of computations.
#'   Default is `FALSE`.
#' @param pl Whether or not to use parallel processing
#' @param cores The number of cores to use for parallel processing
#' @param est_method the method to compute group-time average treatment effects.  The default is "dr" which uses the doubly robust
#' approach in the `DRDID` package.  Other built-in methods
#' include "ipw" for inverse probability weighting and "reg" for
#' first step regression estimators.  The user can also pass their
#' own function for estimating group time average treatment
#' effects.  The required signature depends on the data structure:
#'
#' **Panel data** (`panel=TRUE`): `f(y1, y0, D, covariates,
#' i.weights, inffunc, ...)` where `y1` is an `n x 1` vector of
#' post-treatment outcomes, `y0` is an `n x 1` vector of
#' pre-treatment outcomes, `D` is a binary vector indicating
#' treatment group membership, `covariates` is an `n x k` matrix,
#' `i.weights` is a vector of sampling weights, and `inffunc` is a
#' logical requesting influence-function computation.
#'
#' **Repeated cross sections / unbalanced panel** (`panel=FALSE`):
#' `f(y, post, D, covariates, i.weights, inffunc, ...)` where `y` is
#' the outcome vector (length `n`), `post` is a binary indicator for
#' the post-treatment period, `D` is a binary treatment indicator,
#' `covariates` is an `n x k` matrix, `i.weights` is a vector of
#' sampling weights, and `inffunc` is a logical.
#'
#' In both cases the function should return a list that includes
#' `ATT` (the estimated group-time average treatment effect) and
#' `att.inf.func` (an `n x 1` influence function — one entry per
#' observation passed into the estimator).
#' The function can return other things as well, but these are
#' the only two that are required. `est_method` is only used
#' if covariates are included.
#' @param xformla A formula for the covariates to include in the
#'  model.  It should be of the form `~ X1 + X2`.  Default
#'  is NULL which is equivalent to `xformla=~1`.  This is
#'  used to create a matrix of covariates which is then passed
#'  to the 2x2 DID estimator chosen in `est_method`.
#'
#'  For time-varying covariates: (1) With balanced panel data,
#'  in each 2x2 comparison, the covariates
#'  are taken to be the value of the covariates in the earlier time
#'  period, and all of the underlying computations involve changes in Y
#'  as a function of those values of covariates.  (2) With repeated cross
#'  sections data and unbalanced panel data, the covariates are taken
#'  from each time period and computations involve Y_post conditional
#'  on X_post minus Y_pre conditional on X_pre.  A byproduct of this
#'  is that, with balanced panel data and in the presence of
#'  time-varying covariates, it is possible to get different numerical
#'  results according to whether or not `allow_unbalanced_panel=TRUE` or
#'  `FALSE`.
#' @param panel Whether or not the data is a panel dataset.
#'  The panel dataset should be provided in long format -- that
#'  is, where each row corresponds to a unit observed at a
#'  particular point in time.  The default is TRUE.  When
#'  using a panel dataset, the variable `idname` must
#'  be set.  When `panel=FALSE`, the data is treated
#'  as repeated cross sections.
#' @param allow_unbalanced_panel Whether or not function should
#'  "balance" the panel with respect to time and id.  The default
#'  value is `FALSE` which means that [att_gt()] will drop
#'  all units where data is not observed in all periods.
#'  The advantage of this is that the computations are faster
#'  (sometimes substantially).
#' @param control_group Which units to use as the control group.
#'  The default is "nevertreated" which sets the control group
#'  to be the group of units that never participate in the
#'  treatment.  This group does not change across groups or
#'  time periods.  The other option is to set
#'  `group="notyettreated"`.  In this case, the control group
#'  is set to the group of units that have not yet participated
#'  in the treatment in that time period.  This includes all
#'  never treated units, but it includes additional units that
#'  eventually participate in the treatment, but have not
#'  participated yet.
#' @param anticipation The number of time periods before participating
#'  in the treatment where units can anticipate participating in the
#'  treatment and therefore it can affect their untreated potential outcomes
#' @param faster_mode This option enables a faster version of `did`, optimizing
#' computation time for large datasets by improving data management within the package.
#' The default is set to `FALSE`. While the difference is minimal for small datasets,
#' it is recommended for use with large datasets.
#' @param base_period Whether to use a "varying" base period or a
#'  "universal" base period.  Either choice results in the same
#'  post-treatment estimates of ATT(g,t)'s.  In pre-treatment
#'  periods, using a varying base period amounts to computing a
#'  pseudo-ATT in each treatment period by comparing the change
#'  in outcomes for a particular group relative to its comparison
#'  group in the pre-treatment periods (i.e., in pre-treatment
#'  periods this setting computes changes from period t-1 to period
#'  t, but repeatedly changes the value of t)
#'
#'  A universal base period fixes the base period to always be
#'  (g-anticipation-1).  This does not compute
#'  pseudo-ATT(g,t)'s in pre-treatment periods, but rather
#'  reports average changes in outcomes from period t to
#'  (g-anticipation-1) for a particular group relative to its comparison
#'  group.  This is analogous to what is often reported in event
#'  study regressions.
#'
#'  Using a varying base period results in an estimate of
#'  ATT(g,t) being reported in the period immediately before
#'  treatment.  Using a universal base period normalizes the
#'  estimate in the period right before treatment (or earlier when
#'  the user allows for anticipation) to be equal to 0, but one
#'  extra estimate in an earlier period.
#'
#' @param ... Additional arguments to be passed to a custom `est_method`
#'  function. These are ignored when using built-in estimation methods
#'  (`"dr"`, `"ipw"`, `"reg"`).
#'
#' @references Callaway, Brantly and Pedro H.C. Sant'Anna.  \"Difference-in-Differences with Multiple Time Periods.\" Journal of Econometrics, Vol. 225, No. 2, pp. 200-230, 2021. \doi{10.1016/j.jeconom.2020.12.001}, <https://arxiv.org/abs/1803.09015>
#'
#' @return an [`MP`] object containing all the results for group-time average
#'  treatment effects
#'
#' @details # Examples:
#'
#' **Basic [att_gt()] call:**
#' ```{r, comment = "#>", collapse = TRUE}
#' # Example data
#' data(mpdta)
#' set.seed(09152024)
#' out1 <- att_gt(yname="lemp",
#'                tname="year",
#'                idname="countyreal",
#'                gname="first.treat",
#'                xformla=NULL,
#'                data=mpdta)
#' summary(out1)
#' ```
#'
#' **Using covariates:**
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' out2 <- att_gt(yname="lemp",
#'                tname="year",
#'                idname="countyreal",
#'                gname="first.treat",
#'                xformla=~lpop,
#'                data=mpdta)
#' summary(out2)
#' ```
#'
#' **Specify comparison units:**
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' out3 <- att_gt(yname="lemp",
#'                tname="year",
#'                idname="countyreal",
#'                gname="first.treat",
#'                xformla=~lpop,
#'                control_group = "notyettreated",
#'                data=mpdta)
#' summary(out3)
#' ```
#'
#' @export


att_gt <- function(yname,
                   tname,
                   idname = NULL,
                   gname,
                   xformla = NULL,
                   data,
                   panel = TRUE,
                   allow_unbalanced_panel = FALSE,
                   control_group = c("nevertreated", "notyettreated"),
                   anticipation = 0,
                   weightsname = NULL,
                   fix_weights = NULL,
                   alp = 0.05,
                   bstrap = TRUE,
                   cband = TRUE,
                   biters = 1000,
                   clustervars = NULL,
                   est_method = "dr",
                   base_period = "varying",
                   faster_mode = TRUE,
                   print_details = FALSE,
                   pl = FALSE,
                   cores = 1,
                   ...) {
  # Capture extra arguments for custom est_method
  extra_args <- list(...)

  # Warn if extra arguments passed with built-in est_method
  if (length(extra_args) > 0 && !inherits(est_method, "function")) {
    warning("Extra arguments (", paste(names(extra_args), collapse = ", "),
            ") are ignored when using built-in est_method = \"", est_method,
            "\". Extra arguments are only passed to custom est_method functions.")
  }

  # Validate fix_weights
  if (!is.null(fix_weights)) {
    if (!is.character(fix_weights) || length(fix_weights) != 1 ||
        !(fix_weights %in% c("varying", "base_period", "first_period"))) {
      stop("fix_weights must be NULL or one of \"varying\", \"base_period\", or \"first_period\".")
    }
    if (!panel && fix_weights %in% c("base_period", "first_period")) {
      stop("fix_weights = \"", fix_weights, "\" is not supported for repeated cross sections ",
           "(panel = FALSE) because units are not tracked across periods. ",
           "Use fix_weights = \"varying\" or NULL instead.")
    }
    if (fix_weights == "varying" && panel && inherits(est_method, "function")) {
      stop("fix_weights = \"varying\" is not currently supported with custom est_method functions ",
           "on panel data. The \"varying\" option uses repeated cross-section estimators internally, ",
           "which require a different function signature (y, post, D) than the panel signature ",
           "(y1, y0, D). Use fix_weights = NULL, \"base_period\", or \"first_period\" instead.")
    }
  }

  # Validate est_method
  if (!inherits(est_method, "function")) {
    if (!is.character(est_method) || length(est_method) != 1) {
      stop("est_method must be a character string (\"dr\", \"ipw\", or \"reg\") or a custom function. ",
           "Received an object of class '", class(est_method)[1], "'.")
    }
    if (!(est_method %in% c("dr", "ipw", "reg"))) {
      stop("est_method must be one of \"dr\", \"ipw\", or \"reg\". Received: \"", est_method, "\".")
    }
  }

  # Warn users about anticipation and never-treated units
  if (anticipation > 0) {
    message("Note: anticipation = ", anticipation, ". Never-treated units (with group status 0 or Inf) ",
            "are assumed to never anticipate treatment. Anticipation only applies to eventually-treated units.")
  }

  # Check if user wants to run faster mode:
  if (faster_mode) {
    # this is a DIDparams2 object
    dp <- pre_process_did2(
      yname = yname,
      tname = tname,
      idname = idname,
      gname = gname,
      xformla = xformla,
      data = data,
      panel = panel,
      allow_unbalanced_panel = allow_unbalanced_panel,
      control_group = control_group,
      anticipation = anticipation,
      weightsname = weightsname,
      fix_weights = fix_weights,
      alp = alp,
      bstrap = bstrap,
      cband = cband,
      biters = biters,
      clustervars = clustervars,
      est_method = est_method,
      base_period = base_period,
      print_details = print_details,
      faster_mode = faster_mode,
      pl = pl,
      cores = cores,
      call = match.call()
    )

    # attach extra args for custom est_method
    dp$extra_args <- extra_args

    #-----------------------------------------------------------------------------
    # Compute all ATT(g,t)
    #-----------------------------------------------------------------------------
    results <- compute.att_gt2(dp)
  } else {
    # this is a DIDparams object
    dp <- pre_process_did(
      yname = yname,
      tname = tname,
      idname = idname,
      gname = gname,
      xformla = xformla,
      data = data,
      panel = panel,
      allow_unbalanced_panel = allow_unbalanced_panel,
      control_group = control_group,
      anticipation = anticipation,
      weightsname = weightsname,
      fix_weights = fix_weights,
      alp = alp,
      bstrap = bstrap,
      cband = cband,
      biters = biters,
      clustervars = clustervars,
      est_method = est_method,
      base_period = base_period,
      print_details = print_details,
      pl = pl,
      cores = cores,
      call = match.call()
    )

    # attach extra args for custom est_method
    dp$extra_args <- extra_args

    #-----------------------------------------------------------------------------
    # Compute all ATT(g,t)
    #-----------------------------------------------------------------------------
    results <- compute.att_gt(dp)
  }

  # extract ATT(g,t) and influence functions
  attgt.list <- results$attgt.list
  inffunc <- results$inffunc

  # process results
  # attgt.results <- process_attgt(attgt.list)
  tryCatch(
    {
      # Attempt to run this line for process results
      attgt.results <- process_attgt(attgt.list)
    },
    error = function(e) {
      # Handle the error
      if (faster_mode) {
        # If faster_mode is TRUE, send this stop message
        stop("An unexpected error occurred, normally associated with a singular matrix due to not enough control units. Try changing faster_mode=FALSE.")
      } else {
        # If faster_mode is FALSE, send this stop message
        stop("An unexpected error occurred, normally associated with a singular matrix due to not enough control units.")
      }
    }
  )
  group <- attgt.results$group
  att <- attgt.results$att
  tt <- attgt.results$tt



  # analytical standard errors
  # estimate variance. The i.i.d. form is clustered at the unit level; with a coarser cluster variable
  # the variance is formed from cluster sums instead (the cluster_analytic branch below).
  n <- ifelse(faster_mode, dp$id_count, dp$n)
  # Analytical variance of the ATT(g,t)'s. With a cluster variable (beyond idname) and no bootstrap, use
  # the cluster-robust form -- the cluster sums of the influence function -- mirroring the cluster-sum
  # aggregation of the multiplier bootstrap (Callaway & Sant'Anna 2021, Remark 10). It is kept on the same
  # 1/n scaling as the i.i.d. V so that se = sqrt(diag(V)/n) and the Wald statistic n * b' V^{-1} b stay
  # valid for either form.
  # cluster variables beyond the unit id; clustering on the unit id alone is just the unit-level (i.i.d.)
  # SE. Use the *internal* unit id dp$idname -- the user's idname for panels, and ".rowid" for repeated
  # cross-sections / unbalanced panels (where the user may omit idname), mirroring how mboot() identifies units.
  unit_id <- dp$idname
  extra_clustervars <- clustervars[!(clustervars %in% c(unit_id, idname, ""))]
  # Per-unit cluster identifiers aligned with the rows of inffunc. faster_mode supplies dp$cluster_vector,
  # but for unbalanced panels (internally panel = FALSE, so the time-invariant data keeps every
  # observation) that vector is observation-length rather than unit-length and no longer aligns with the
  # unit-level influence function. Derive the vector from the data exactly as mboot() does -- one row per
  # cross-sectional unit, keyed on the internal unit id -- and store it back whenever it is missing OR does
  # not align with inffunc. This keeps the analytical clustered SE and aggte() working for faster_mode TRUE
  # and FALSE and for balanced and unbalanced panels; repeated cross sections (where the per-observation
  # vector already has one row per unit and aligns) fall through unchanged.
  if (length(extra_clustervars) > 0 && !is.null(dp$data) && !is.null(unit_id) &&
      (is.null(dp$cluster_vector) || length(dp$cluster_vector) != nrow(inffunc))) {
    cdat <- as.data.frame(dp$data)
    if (all(c(unit_id, extra_clustervars[1L]) %in% names(cdat))) {
      rebuilt_cv <- as.vector(unique(cdat[, c(unit_id, extra_clustervars[1L])])[, 2L])
      # only adopt the rebuilt vector if it lines up with the influence-function rows
      if (length(rebuilt_cv) == nrow(inffunc)) dp$cluster_vector <- rebuilt_cv
    }
  }
  cluster_vec <- dp$cluster_vector
  cluster_analytic <- (length(extra_clustervars) > 0) && !bstrap &&
    !is.null(cluster_vec) && length(cluster_vec) == nrow(inffunc)
  if (cluster_analytic) {
    Sc <- rowsum(as.matrix(inffunc), cluster_vec)   # n_clusters x k cluster sums
    V  <- Matrix::Matrix(crossprod(Sc) / n)
  } else {
    V <- Matrix::t(inffunc) %*% inffunc / (n)
  }
  se <- sqrt(Matrix::diag(V) / n)

  # Zero standard error replaced by NA
  se[se <= sqrt(.Machine$double.eps) * 10] <- NA

  # If a cluster variable beyond idname is requested without the bootstrap but the cluster-robust variance
  # could not be formed (e.g. cluster identifiers unavailable), fall back to i.i.d. and let the user know.
  if ((length(extra_clustervars) > 0) & !bstrap & !cluster_analytic) {
    warning("Clustered standard errors could not be computed analytically; the reported standard errors do NOT account for clustering. Set bstrap = TRUE for the cluster-robust multiplier bootstrap.")
  }

  # Identify entries of main diagonal V that are zero or NA
  zero_na_sd_entry <- unique(which(is.na(se)))

  # bootstrap variance matrix
  if (bstrap) {
    bout <- mboot(inffunc, DIDparams = dp, pl = pl, cores = cores)
    bres <- bout$bres

    if (length(zero_na_sd_entry) > 0) {
      se[-zero_na_sd_entry] <- bout$se[-zero_na_sd_entry]
    } else {
      se <- bout$se
    }
  }
  # Zero standard error replaced by NA
  se[se <= sqrt(.Machine$double.eps) * 10] <- NA


  #-----------------------------------------------------------------------------
  # compute Wald pre-test
  #-----------------------------------------------------------------------------

  # Determine whether the variance matrix V is a valid basis for the Wald pre-test. The i.i.d. form
  # V = t(inffunc) %*% inffunc / n is valid for balanced panels, unbalanced panels (influence functions
  # are aggregated to the unit level), and repeated cross-sections (CLT applies per independent
  # observation). When a cluster variable beyond idname is specified and V is built from cluster sums
  # (cluster_analytic, above), V is cluster-robust and the Wald test remains valid. It is only skipped
  # when an extra cluster variable is present but V is the i.i.d. form (e.g. bstrap = TRUE), since then
  # V ignores between-cluster correlation; in that case we rely on the bootstrap confidence bands.

  wald_invalid <- NULL
  if (length(extra_clustervars) > 0 && !cluster_analytic) {
    wald_invalid <- paste0(
      "The Wald pre-test is not reported when clustering beyond the unit level ",
      "(clustervars = '", paste(extra_clustervars, collapse = "', '"), "') ",
      "because the analytical variance matrix does not account for between-cluster ",
      "correlation. Use the bootstrap confidence intervals to assess pre-trends."
    )
  }

  if (!is.null(wald_invalid)) {
    message(wald_invalid)
    W <- NULL
    Wpval <- NULL
  }

  if (is.null(wald_invalid)) {
    # select which periods are pre-treatment
    pre <- which(group > tt)

    # Drop group-periods that have variance equal to zero (singularity problems)
    if (length(zero_na_sd_entry) > 0) {
      pre <- pre[!(pre %in% zero_na_sd_entry)]
    }
    # pseudo-atts in pre-treatment periods
    preatt <- as.matrix(att[pre])

    # covariance matrix of pre-treatment atts
    preV <- as.matrix(V[pre, pre])

    # check if there are actually any pre-treatment periods
    if (length(preV) == 0) {
      msg <- paste0(
        "No pre-treatment periods available for the Wald pre-test of parallel trends. ",
        "This can happen when all groups are first treated early in the panel ",
        "(e.g., in the second time period) so that no pre-treatment ATT(g,t) estimates exist."
      )
      if (anticipation > 0) {
        msg <- paste0(msg, " Note: anticipation=", anticipation, " further reduces the number of available pre-treatment periods.")
      }
      warning(msg)
      W <- NULL
      Wpval <- NULL
    } else if (sum(is.na(preV))) {
      warning("Not returning pre-test Wald statistic due to NA pre-treatment values")
      W <- NULL
      Wpval <- NULL
    } else if (rcond(preV) <= .Machine$double.eps) {
      # singular covariance matrix for pre-treatment periods
      warning("Not returning pre-test Wald statistic due to singular covariance matrix")
      W <- NULL
      Wpval <- NULL
    } else {
      # everything is working...
      W <- n * t(preatt) %*% solve(preV) %*% preatt
      q <- length(pre) # number of restrictions
      Wpval <- round(1 - pchisq(W, q), 5)
    }
  }


  #-----------------------------------------------------------------------------
  # compute confidence intervals / bands
  #-----------------------------------------------------------------------------

  # critical value from N(0,1), for pointwise
  cval <- qnorm(1 - alp / 2)

  # in order to get uniform confidence bands
  # HAVE to use the bootstrap
  if (bstrap) {
    if (cband) {
      # for uniform confidence band
      # compute new critical value
      # see paper for details
      bSigma <- apply(
        bres, 2,
        function(b) {
          (quantile(b, .75, type = 1, na.rm = T) -
            quantile(b, .25, type = 1, na.rm = T)) / (qnorm(.75) - qnorm(.25))
        }
      )

      bSigma[bSigma <= sqrt(.Machine$double.eps) * 10] <- NA

      # sup-t confidence band
      bT <- apply(bres, 1, function(b) max(abs(b / bSigma), na.rm = TRUE))
      cval <- quantile(bT, 1 - alp, type = 1, na.rm = T)
      if (cval >= 7) {
        warning("Simultaneous critical value is arguably `too large' to be reliable. This usually happens when the number of observations per group is small and/or there is not much variation in outcomes.")
      }
    }
  }


  # Return this list
  return(MP(group = group, t = tt, att = att, V_analytical = V, se = se, c = cval, inffunc = inffunc, n = n, W = W, Wpval = Wpval, alp = alp, DIDparams = dp))
}
