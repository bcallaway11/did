#' @title Group-Time Average Treatment Effects
#'
#' @description \code{att_gt} computes average treatment effects in DID
#'  setups where there are more than two periods of data and allowing for
#'  treatment to occur at different points in time and allowing for
#'  treatment effect heterogeneity and dynamics.
#'  See Callaway and Sant'Anna (2019) for a detailed description.
#'
#' @param yname The name of the outcome variable
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param idname The individual (cross-sectional unit) id name
#' @param first.treat.name The name of the variable in \code{data} that
#'  contains the first period when a particular observation is treated.
#'  This should be a positive number for all observations in treated groups.
#'  It defines which "group" a unit belongs to.  It should be 0 for units
#'  in the untreated group.
#' @param weightsname The name of the column containing the sampling weights.
#'  If not set, all observations have same weight.
#' @param alp the significance level, default is 0.05
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier boostrap.  If standard errors are clustered, then one
#'  must set \code{bstrap=TRUE}. Default is \code{TRUE} (in addition, cband
#'  is also by default \code{TRUE} indicating that uniform confidence bands
#'  will be returned.  If bstrap is \code{FALSE}, then analytical
#'  standard errors are reported.
#' @param biters The number of boostrap iterations to use.  The default is 1000,
#'  and this is only applicable if \code{bstrap=TRUE}.
#' @param clustervars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level.
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability \code{1-alp}.  In order to compute uniform confidence
#'  bands, \code{bstrap} must also be set to \code{TRUE}.  The default is
#' \code{TRUE}.
#' @param printdetails Boolean for showing detailed results or not
#' @param pl Boolean for whether or not to use parallel processing
#'  (not implemented yet)
#' @param cores The number of cores to use for parallel processing
#'  (not implemented yet)
#' @param estMethod the method to compute group-time average treatment effects.  The default is "dr" which uses the doubly robust
#' approach in the \code{DRDID} package.  Other built-in methods
#' include "ipw" for inverse probability weighting and "reg" for
#' first step regression estimators.  The user can also pass their
#' own function for estimating group time average treatment
#' effects.  This should be a function
#' \code{f(Y1,Y0,treat,covariates)} where \code{Y1} is an
#' \code{n} x \code{1} vector of outcomes in the post-treatment
#' outcomes, \code{Y0} is an \code{n} x \code{1} vector of
#' pre-treatment outcomes, \code{treat} is a vector indicating
#' whether or not an individual participates in the treatment,
#' and \code{covariates} is an \code{n} x \code{k} matrix of
#' covariates.  The function should return a list that includes
#' \code{ATT} (an estimated average treatment effect), and
#' \code{inf.func} (an \code{n} x \code{1} influence function).
#' The function can return other things as well, but these are
#' the only two that are required. \code{estMethod} is only used
#' if covariates are included.
#' @param xformla A formula for the covariates to include in the
#'  model.  It should be of the form \code{~ X1 + X2}.  Default
#'  is NULL which is equivalent to \code{xformla=~1}.  This is
#'  used to create a matrix of covariates which is then passed
#'  to the 2x2 DID estimator chosen in \code{estMethod}.
#' @param panel Whether or not the data is a panel dataset.
#'  The panel dataset should be provided in long format -- that
#'  is, where each row corresponds to a unit observed at a
#'  particular point in time.  The default is TRUE.  When
#'  is using a panel dataset, the variable \code{idname} must
#'  be set.  When \code{panel=FALSE}, the data is treated
#'  as repeated cross sections.
#' @param control.group Which units to use the control group.
#'  The default is "nevertreated" which sets the control group
#'  to be the group of units that never participate in the
#'  treatment.  This group does not change across groups or
#'  time periods.  The other option is to set
#'  \code{group="notyettreated"}.  In this case, the control group
#'  is set to the group of units that have not yet participated
#'  in the treatment in that time period.  This includes all
#'  never treated units, but it includes additional units that
#'  eventually participate in the treatment, but have not
#'  participated yet.
#' @references Callaway, Brantly and Sant'Anna, Pedro H. C.. "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2019).
#'
#' @examples
#' data(mpdta)
#'
#' # with covariates
#' out1 <- att_gt(yname="lemp",
#'                tname="year",
#'                idname="countyreal",
#'                first.treat.name="first.treat",
#'                xformla=~lpop,
#'                data=mpdta,
#'                printdetails=FALSE)
#' summary(out1)
#'
#' # without covariates
#' out2 <- att_gt(yname="lemp",
#'                tname="year",
#'                idname="countyreal",
#'                first.treat.name="first.treat",
#'                xformla=NULL,
#'                data=mpdta,
#'                printdetails=FALSE)
#' summary(out2)
#'
#' @return an \code{\link{MP}} object containing all the results for group-time average
#'  treatment effects
#'
#' @export
att_gt <- function(yname,
                   tname,
                   idname=NULL,
                   first.treat.name,
                   xformla=NULL,
                   data,
                   panel=TRUE,
                   control.group=c("nevertreated","notyettreated"),
                   weightsname=NULL,
                   alp=0.05,
                   bstrap=TRUE,
                   cband=TRUE,
                   biters=1000,
                   clustervars=NULL,
                   estMethod="dr",
                   printdetails=TRUE,
                   pl=FALSE,
                   cores=1) {

  # this is a DIDparams object
  dp <- pre_process_did(yname=yname,
                        tname=tname,
                        idname=idname,
                        first.treat.name=first.treat.name,
                        xformla=xformla,
                        data=data,
                        panel=panel,
                        control.group=control.group,
                        weightsname=weightsname,
                        alp=alp,
                        bstrap=bstrap,
                        cband=cband,
                        biters=biters,
                        clustervars=clustervars,
                        estMethod=estMethod,
                        printdetails=printdetails,
                        pl=pl,
                        cores=cores
                        )

  #-----------------------------------------------------------------------------
  # Compute all ATT(g,t)
  #-----------------------------------------------------------------------------
  results <- compute.att_gt(dp)


  # extract ATT(g,t) and influence functions
  attgt.list <- results$attgt.list
  inffunc <- results$inffunc

  # process results
  attgt.results <- process_attgt(results)
  group <- attgt.results$group
  att <- attgt.results$att
  tt <- attgt.results$tt
  inffunc1 <- attgt.results$inf.func



  # estimate variance
  # this is analogous to cluster robust standard errors that
  # are clustered at the unit level
  n <- dp$n
  V <- t(inffunc1)%*%inffunc1/n

  # if clustering along another dimension...we require using the
  # bootstrap (in principle, could come up with analytical standard
  # errors here though)
  if ( (length(clustervars) > 0) & !bstrap) {
    warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  }


  # bootstrap variance matrix
  if (bstrap) {

    bout <- mboot(inffunc1, DIDparams=dp)
    bres <- bout$bres
    V[!is.na(V)] <- bout$V
  }



  #-----------------------------------------------------------------------------
  # compute Wald pre-test
  #-----------------------------------------------------------------------------

  # select which periods are pre-treatment
  pre <- which(group > tt)

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
    } else if (det(preV) == 0) {
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

  # in order to get uniform confidencs bands
  # HAVE to use the bootstrap
  if (bstrap){
    if (cband) {
      # for uniform confidence band
      # compute new critical value
      # see paper for details
      bSigma <- apply(bres, 2,
                      function(b) (quantile(b, .75, type=1, na.rm = T) -
                                     quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
      # sup-t confidence band
      bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
      cval <- quantile(bT, 1-alp, type=1, na.rm = T)
    }
  }

  # Return this list
  return(MP(group=group, t=tt, att=att, V=V, c=cval, inffunc=inffunc1, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp))

}
