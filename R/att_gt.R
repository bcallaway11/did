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
#'  If not set, all observations have same weight.
#' @param alp the significance level, default is 0.05
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier bootstrap.  If standard errors are clustered, then one
#'  must set `bstrap=TRUE`. Default is `TRUE` (in addition, cband
#'  is also by default `TRUE` indicating that uniform confidence bands
#'  will be returned.  If bstrap is `FALSE`, then analytical
#'  standard errors are reported.
#' @param biters The number of bootstrap iterations to use.  The default is 1000,
#'  and this is only applicable if `bstrap=TRUE`.
#' @param clustervars A vector of variables names to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level. By default, we cluster at individual level (when `bstrap=TRUE`).
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability `1-alp`.  In order to compute uniform confidence
#'  bands, `bstrap` must also be set to `TRUE`.  The default is
#' `TRUE`.
#' @param print_details Whether or not to show details/progress of computations.
#'   Default is `FALSE`.
#' @param pl Whether or not to use parallel processing
#'  (not implemented yet)
#' @param cores The number of cores to use for parallel processing
#'  (not implemented yet)
#' @param est_method the method to compute group-time average treatment effects.  The default is "dr" which uses the doubly robust
#' approach in the `DRDID` package.  Other built-in methods
#' include "ipw" for inverse probability weighting and "reg" for
#' first step regression estimators.  The user can also pass their
#' own function for estimating group time average treatment
#' effects.  This should be a function
#' `f(Y1,Y0,treat,covariates)` where `Y1` is an
#' `n` x `1` vector of outcomes in the post-treatment
#' outcomes, `Y0` is an `n` x `1` vector of
#' pre-treatment outcomes, `treat` is a vector indicating
#' whether or not an individual participates in the treatment,
#' and `covariates` is an `n` x `k` matrix of
#' covariates.  The function should return a list that includes
#' `ATT` (an estimated average treatment effect), and
#' `inf.func` (an `n` x `1` influence function).
#' The function can return other things as well, but these are
#' the only two that are required. `est_method` is only used
#' if covariates are included.
#' @param xformla A formula for the covariates to include in the
#'  model.  It should be of the form `~ X1 + X2`.  Default
#'  is NULL which is equivalent to `xformla=~1`.  This is
#'  used to create a matrix of covariates which is then passed
#'  to the 2x2 DID estimator chosen in `est_method`.
#' @param panel Whether or not the data is a panel dataset.
#'  The panel dataset should be provided in long format -- that
#'  is, where each row corresponds to a unit observed at a
#'  particular point in time.  The default is TRUE.  When
#'  is using a panel dataset, the variable `idname` must
#'  be set.  When `panel=FALSE`, the data is treated
#'  as repeated cross sections.
#' @param allow_unbalanced_panel Whether or not function should
#'  "balance" the panel with respect to time and id.  The default
#'  values if `FALSE` which means that [att_gt()] will drop
#'  all units where data is not observed in all periods.
#'  The advantage of this is that the computations are faster
#'  (sometimes substantially).
#' @param control_group Which units to use the control group.
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
#'
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
                   idname=NULL,
                   gname,
                   xformla=NULL,
                   data,
                   panel=TRUE,
                   allow_unbalanced_panel=FALSE,
                   control_group=c("nevertreated","notyettreated"),
                   anticipation=0,
                   weightsname=NULL,
                   alp=0.05,
                   bstrap=TRUE,
                   cband=TRUE,
                   biters=1000,
                   clustervars=NULL,
                   est_method="dr",
                   base_period="varying",
                   print_details=FALSE,
                   pl=FALSE,
                   cores=1) {

  # this is a DIDparams object
  dp <- pre_process_did(yname=yname,
                        tname=tname,
                        idname=idname,
                        gname=gname,
                        xformla=xformla,
                        data=data,
                        panel=panel,
                        allow_unbalanced_panel=allow_unbalanced_panel,
                        control_group=control_group,
                        anticipation=anticipation,
                        weightsname=weightsname,
                        alp=alp,
                        bstrap=bstrap,
                        cband=cband,
                        biters=biters,
                        clustervars=clustervars,
                        est_method=est_method,
                        base_period=base_period,
                        print_details=print_details,
                        pl=pl,
                        cores=cores,
                        call=match.call()
  )

  #-----------------------------------------------------------------------------
  # Compute all ATT(g,t)
  #-----------------------------------------------------------------------------
  results <- compute.att_gt(dp)


  # extract ATT(g,t) and influence functions
  attgt.list <- results$attgt.list
  inffunc <- results$inffunc

  # process results
  attgt.results <- process_attgt(attgt.list)
  group <- attgt.results$group
  att <- attgt.results$att
  tt <- attgt.results$tt



  # analytical standard errors
  # estimate variance
  # this is analogous to cluster robust standard errors that
  # are clustered at the unit level

  # note to self: this def. won't work with unbalanced panel,
  # same with clustered standard errors
  # but it is always ignored b/c bstrap has to be true in that case
  n <- dp$n
  V <- Matrix::t(inffunc)%*%inffunc/n
  se <- sqrt(Matrix::diag(V)/n)

  # Zero standard error replaced by NA
  se[se <= sqrt(.Machine$double.eps)*10] <- NA

  # if clustering along another dimension...we require using the
  # bootstrap (in principle, could come up with analytical standard
  # errors here though)
  if ( (length(clustervars) > 0) & !bstrap) {
    warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  }

  # Identify entries of main diagonal V that are zero or NA
  zero_na_sd_entry <- unique(which(is.na(se))) 

  # bootstrap variance matrix
  if (bstrap) {

    bout <- mboot(inffunc, DIDparams=dp)
    bres <- bout$bres
    
    if(length(zero_na_sd_entry)>0) {
      se[-zero_na_sd_entry] <- bout$se[-zero_na_sd_entry]
    } else {
      se <- bout$se
    }
  }
  # Zero standard error replaced by NA
  se[se <= sqrt(.Machine$double.eps)*10] <- NA


  #-----------------------------------------------------------------------------
  # compute Wald pre-test
  #-----------------------------------------------------------------------------

  # select which periods are pre-treatment
  pre <- which(group > tt)

  # Drop group-periods that have variance equal to zero (singularity problems)
  if(length(zero_na_sd_entry)>0){
    pre <- pre[!(pre %in% zero_na_sd_entry)]
  }
  # pseudo-atts in pre-treatment periods
  preatt <- as.matrix(att[pre])

  # covariance matrix of pre-treatment atts
  preV <- as.matrix(V[pre,pre])

  # check if there are actually any pre-treatment periods
  if (length(preV) == 0) {
    message("No pre-treatment periods to test")
    W  <- NULL
    Wpval <- NULL
  } else if(sum(is.na(preV))) {
    warning("Not returning pre-test Wald statistic due to NA pre-treatment values")
    W <- NULL
    Wpval <- NULL
  } else if (rcond(preV) <= .Machine$double.eps) {
    # singluar covariance matrix for pre-treatment periods
    warning("Not returning pre-test Wald statistic due to singular covariance matrix")
    W <- NULL
    Wpval <- NULL
  } else {
    # everything is working...
    W <- n*t(preatt)%*%solve(preV)%*%preatt
    q <- length(pre) # number of restrictions
    Wpval <- round(1-pchisq(W,q),5)
  }


  #-----------------------------------------------------------------------------
  # compute confidence intervals / bands
  #-----------------------------------------------------------------------------

  # critical value from N(0,1), for pointwise
  cval <- qnorm(1-alp/2)

  # in order to get uniform confidence bands
  # HAVE to use the bootstrap
  if (bstrap){
    if (cband) {
      # for uniform confidence band
      # compute new critical value
      # see paper for details
      bSigma <- apply(bres, 2,
                      function(b) (quantile(b, .75, type=1, na.rm = T) -
                                     quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))

      bSigma[bSigma <= sqrt(.Machine$double.eps)*10] <- NA

      # sup-t confidence band
      bT <- apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = TRUE))
      cval <- quantile(bT, 1-alp, type=1, na.rm = T)
      if(cval >= 7){
        warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
      }
    }
  }


  # Return this list
  return(MP(group=group, t=tt, att=att, V_analytical=V, se=se, c=cval, inffunc=inffunc, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp))

}
