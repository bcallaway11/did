#' @title reset.sim
#' @description a function to create a "reasonable" set of parameters
#'  to create simulated panel data that obeys a parallel trends assumption.
#'  In particular, it provides parameters where the the effect of participating
#'  in the treatment is equal to one in all post-treatment time periods.
#' 
#'  After calling this function, the user can change particular values of the
#'  parameters in order to generate dynamics, heterogeneous effects across
#'  groups, etc.
#'
#' @param time.periods The number of time periods to include
#' @param n The total number of observations
#' @param ipw If TRUE, sets parameters so that DGP is
#'  compatible with recovering ATT(g,t)'s using IPW (i.e.,
#'  where logit that just includes a linear term in X works).  If
#'  FALSE, sets parameters that will be incompatible with IPW.
#'  Either way, these parameters can be specified by the user
#'  if so desired.
#' @param reg If TRUE, sets parameters so that DGP is compatible
#'  with recovering ATT(g,t)'s using regressions on untreated
#'  untreated potential outcomes.  If FALSE, sets parameters that
#'  will be incompatible with using regressions (i.e., regressions
#'  that include only linear term in X).  Either way, these
#'  parameters can be specified by the user if so desired.
#'
#' @return list of simulation parameters
#' 
#' @export
reset.sim <- function(time.periods=4, n=5000, ipw=TRUE, reg=TRUE) {
  #-----------------------------------------------------------------------------
  # set parameters
  #-----------------------------------------------------------------------------
  # coefficient on X 
  bett <- seq(1:time.periods)
  # time fixed effect
  thet <- seq(1:time.periods)
  # time fixed effect
  theu <- thet # changing this creates violations of parallel trends
  # covariate effect
  betu <- bett # changing this creates violations of conditional parallel trends
  #-----------------------------------------------------------------------------
  # parameters for treated potential outcomes
  #-----------------------------------------------------------------------------
  te.bet.ind <- rep(1,time.periods) # no selective treatment timing
  te.bet.X <- bett #no heterogeneous effects by X
  te.t <- thet # no calendar time effects
  te.e <- rep(0,time.periods) # no dynamic effects
  te <- 1 # overall basic effect
  # parameters in generalized propensity score
  # don't make them too big otherwise can get divide by 0
  gamG <- c(0,1:time.periods)/(2*time.periods)

  # return list of parameters
  list(time.periods=time.periods,
       bett=bett,
       thet=thet,
       theu=theu,
       betu=betu,
       te.bet.ind=te.bet.ind,
       te.bet.X=te.bet.X,
       te.t=te.t,
       te.e=te.e,
       te=te,
       n=n,
       gamG=gamG,
       ipw=ipw,
       reg=reg
       )
}

#' @title build_sim_dataset
#'
#' @description A function for building simulated data
#'
#' @param sp_list A list of simulation parameters.  See `reset.sim` to generate
#'  some default values for parameters
#' @param panel whether to construct panel data (the default) or repeated
#'  cross sections data
#'
#' @return a data.frame with the following columns
#'   \itemize{
#'     \item G observations group
#'     \item X value of covariate
#'     \item id observation's id
#'     \item cluster observation's cluster (by construction there is no within-cluster correlation)
#'     \item period time period for current observation
#'     \item Y outcome
#'     \item treat whether or not this unit is ever treated
#'   }
#'
#' @export
build_sim_dataset <- function(sp_list, panel=TRUE) {
  #-----------------------------------------------------------------------------
  # build dataset
  #-----------------------------------------------------------------------------
  time.periods <- sp_list$time.periods
  nt <- sp_list$nt
  bett <- sp_list$bett
  thet=sp_list$thet
  nu <- sp_list$nu
  theu <- sp_list$theu
  betu <- sp_list$betu
  te.bet.ind <- sp_list$te.bet.ind
  te.bet.X <- sp_list$te.bet.X
  te.t <- sp_list$te.t
  te.e <- sp_list$te.e
  te <- sp_list$te
  n <- sp_list$n
  gamG <- sp_list$gamG
  ipw <- sp_list$ipw
  reg <- sp_list$reg

  X <- rnorm(n)

  if (ipw) {
    pr <- exp(outer(X,gamG)) / apply( exp(outer(X,gamG)), 1, sum)
  } else {
    pr <- exp(outer((pnorm(X)+0.5)^2,gamG)) / apply( exp(outer((pnorm(X)+0.5)^2,gamG)), 1, sum)
  }

  G <- apply(pr, 1, function(pvec) sample(seq(0,time.periods), size=1, prob=pvec))

  Gt <- G[G>0]
  nt <- length(Gt)

  if (reg) {
    Xmodel <- X
  } else {
    Xmodel <- X^2
  }

  Xt <- Xmodel[G>0]

  # draw individual fixed effect
  Ct <- rnorm(nt, mean=G)

  # generate untreated potential outcomes in each time period
  Ynames <- paste0("Y",1:time.periods)
  #Ynames <- paste0(1:time.periods)
  Y0tmat <- sapply(1:time.periods, function(t) {
    thet[t] + Ct + Xt*bett[t] + rnorm(nt)
  })
  Y0tdf <- as.data.frame(Y0tmat)
  
  # generate treated potential outcomes
  Y1tdf <- sapply(1:time.periods, function(t) {
    te.t[t] + te.bet.ind[Gt]*Ct + Xt*te.bet.X[t] + (Gt <= t)*te.e[sapply(1:nt, function(i) max(t-Gt[i]+1,1))] + te + rnorm(nt) # hack for the dynamic effects but ok
  })

  # generate observed data
  Ytdf <- sapply(1:time.periods, function(t) {
    (Gt<=t)*Y1tdf[,t] + (Gt>t)*Y0tdf[,t]
  })
  colnames(Ytdf) <- Ynames

  # store observed data for treated group
  dft <- cbind.data.frame(G=Gt,X=X[G>0],Ytdf)

  # untreated units

  # draw untreated covariate
  nu <- sum(G==0)
  Xu <- Xmodel[G==0]

  # draw untreated fixed effect
  Cu <- rnorm(nu, mean=0)


  # generate untreated potential outcomes
  Y0umat <- sapply(1:time.periods, function(t) {
    theu[t] + Cu + rnorm(nu) + Xu*betu[t]
  })
  Y0udf <- as.data.frame(Y0umat)
  colnames(Y0udf) <- Ynames

  # store dataset of observed outcomes for untreated units
  dfu <- cbind.data.frame(G=0,X=X[G==0],Y0udf)

  # store overall dataset
  df <- rbind.data.frame(dft, dfu)

  # generate id variable
  df$id <- 1:nrow(df)
  # generate clusters (there's no actual within-cluster correlation)
  df$cluster <- sample(1:50, size=nrow(df), replace=TRUE)

  # convert data from wide to long format
  ddf <- tidyr::pivot_longer(df,
                             cols=tidyr::starts_with("Y"),
                             names_to="period",
                             names_prefix="Y",
                             values_to="Y")
  ddf$period <- as.numeric(ddf$period)
  ddf$treat <- 1*(ddf$G > 0)
  ddf <- ddf[order(ddf$id, ddf$period),] # reorder data

  if (!panel) { # repeated cross sections
    n <- nt+nu
    Time <- sample(1:time.periods, size=n, replace=TRUE, prob=rep(1/time.periods, time.periods))
    right.row <- sapply( unique(ddf$id), function(i) {
      which(ddf$id==i & ddf$period==Time[i])
    })
    ddf <- ddf[right.row,]
  }

  ddf <- subset(ddf, G != 1)
  ddf
}


#' @title sim
#' @description An internal function that builds simulated data, computes
#'  ATT(g,t)'s and some aggregations.  It is useful for testing the inference
#'  procedures in the `did` function.
#'
#' @inheritParams reset.sim
#' @inheritParams build_sim_dataset
#'
#' @param ret which type of results to return.  The options are `Wpval` (returns
#'  1 if the p-value from a Wald test that all pre-treatment ATT(g,t)'s are equal
#'  is less than .05),
#'  `cband`  (returns 1 if a uniform confidence band covers 0 for groups and times),
#'  `simple` (returns 1 if, using the simple treatment effect aggregation results
#'  in rejecting that this aggregated treatment effect parameter is equal to 0),
#'  `dynamic` (returns 1 if the uniform confidence band from the dynamic treatment
#'  effect aggregation covers 0 in all pre- and post-treatment periods).  The default
#'  value is NULL, and in this case the function will just return the results from
#'  the call to `att_gt`.
#' @param bstrap whether or not to use the bootstrap to conduct inference (default is TRUE)
#' @param cband whether or not to compute uniform confidence bands in the call to `att_gt`
#'  (the default is TRUE)
#' @param control_group Whether to use the "nevertreated" comparison group (the default)
#'  or the "notyettreated" as the comparison group
#' @param xformla Formula for covariates in `att_gt` (default is `~X`)
#' @param est_method Which estimation method to use in `att_gt` (default is "dr")
#' @param clustervars Any additional variables which should be clustered on
#' @param panel whether to simulate panel data (the default) or otherwise repeated
#'  cross sections data
#' 
#' @return When `ret=NULL`, returns the results of the call to `att_gt`, otherwise returns
#'  1 if the specified test rejects or 0 if not.
#' 
#' @export
sim <- function(sp_list,
                ret=NULL,
                bstrap=TRUE,
                cband=TRUE,
                control_group="nevertreated",
                xformla=~X,
                est_method="dr",
                clustervars=NULL,
                panel=TRUE) {

  ddf <- build_sim_dataset(sp_list=sp_list,
                           panel=panel)

  time.periods <- sp_list$time.periods
  te.e <- sp_list$te.e
  te <- sp_list$te
  
  # get results
  res <- att_gt(yname="Y", xformla=xformla, data=ddf, tname="period", idname="id",
                gname="G", 
                bstrap=bstrap, cband=cband, control_group=control_group,
                est_method=est_method,
                clustervars=clustervars,
                panel=panel)


  if (is.null(ret)) {
    return(res)
  } else if (ret=="Wpval") {
    rej <- 1*(res$Wpval < .05)
    return(rej)
  } else if (ret=="cband") {
    cu <- res$att + res$c*res$se 
    cl <- res$att - res$c*res$se
    covers0 <- 1*(all( (cu > 0) & (cl < 0)))
    return(covers0)
  } else if (ret=="simple") {
    agg <- aggte(res)
    rej <- 1*( abs(agg$overall.att / agg$overall.se) > qnorm(.975) )
    return(rej)
  } else if (ret=="dynamic") {
    truth <- c(rep(0,(time.periods-2)),te+te.e[1:(time.periods-1)])
    agg <- aggte(res, type="dynamic")
    cu <- agg$att.egt + agg$crit.val.egt * agg$se.egt
    cl <- agg$att.egt - agg$crit.val.egt * agg$se.egt
    coverstruth <- 1*(all( (cu > truth) & (cl < truth)))
    return(coverstruth)
  } else if (ret=="notyettreated") {
    agg <- aggte(res)
    rej <- 1*( abs(agg$overall.att / agg$overall.se) > qnorm(.975) )
    return(rej)
  } else {
    return(res)
  }

}

## pretest_sim <- function(ret=NULL, bstrap=FALSE, cband=FALSE,
##                         control.group="nevertreated", panel=TRUE, xformla=~X, cores=1) {

##   ddf <- build_ipw_dataset(panel=panel)

##   # get results
##   res <- conditional_did_pretest(yname="Y", xformla=xformla, data=ddf,
##                                  tname="period", idname="id",
##                                  first.treat.name="G", estMethod="ipw",
##                                  printdetails=FALSE,
##                                  bstrap=bstrap, cband=cband,
##                                  control.group=control.group,
##                                  panel=panel, 
##                                  pl=TRUE, cores=cores)

##   res$CvMpval
## }
