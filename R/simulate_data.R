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
#' @param p The probability of being treated
#'
#' @return list of simulation parameters
#' 
#' @keywords internal
#' @export
reset.sim <- function(time.periods=4, n=5000, p=0.5) {
  #-----------------------------------------------------------------------------
  # set parameters
  #-----------------------------------------------------------------------------
  # number of time periods
  # number of treated units
  n <- 5000
  p <- 0.5
  nt <- rbinom(1, n, p)
  # coefficient on X 
  bett <- seq(1:time.periods)
  # time fixed effect
  thet <- seq(1:time.periods)
  # number of untreated units
  nu <- n - nt
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
  #-----------------------------------------------------------------------------
  # extra parameters for an ipw sim (otherwise, these are not used)
  #-----------------------------------------------------------------------------
  # this ignores p from earlier and just
  # these are parameters in generalized propensity score
  # don't make them too big otherwise can get divide by 0
  gamG <- c(0,1:time.periods)/(2*time.periods)

  # return list of parameters
  list(time.periods=time.periods,
       nt=nt,
       bett=bett,
       thet=thet,
       nu=nu,
       theu=theu,
       betu=betu,
       te.bet.ind=te.bet.ind,
       te.bet.X=te.bet.X,
       te.t=te.t,
       te.e=te.e,
       te=te,
       n=n,
       gamG=gamG)
}


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
  
  # randomly assign treated units to groups
  # random groups
  G <- sample(1:time.periods, size=nt, replace=TRUE)
  # fixed groups
  # G <- unlist(lapply(1:T, function(g) rep(g, nt/T)))
  
  # draw a single covariate
  Xt <- rnorm(nt)

  # draw individual fixed effect
  Ct <- rnorm(nt, mean=G)

  # generate untreated potential outcomes in each time period
  Ynames <- paste0("Y",1:time.periods)
  Ynames <- paste0(1:time.periods)
  Y0tmat <- sapply(1:time.periods, function(t) {
    thet[t] + Ct + Xt*bett[t] + rnorm(nt)
  })
  Y0tdf <- as.data.frame(Y0tmat)
  
  # generate treated potential outcomes
  Y1tdf <- sapply(1:time.periods, function(t) {
    te.t[t] + te.bet.ind[G]*Ct + Xt*te.bet.X[t] + (G <= t)*te.e[sapply(1:nt, function(i) max(t-G[i]+1,1))] + te + rnorm(nt) # hack for the dynamic effects but ok
  })

  # generate observed data
  Ytdf <- sapply(1:time.periods, function(t) {
    (G<=t)*Y1tdf[,t] + (G>t)*Y0tdf[,t]
  })
  colnames(Ytdf) <- Ynames

  # store observed data for treated group
  dft <- cbind.data.frame(G,X=Xt,Ytdf)

  # untreated units

  # draw untreated covariate
  Xu <- rnorm(nu, mean=1)

  # draw untreated fixed effect
  Cu <- rnorm(nu)


  # generate untreated potential outcomes
  Y0umat <- sapply(1:time.periods, function(t) {
    theu[t] + Cu + rnorm(nu) + Xu*betu[t]
  })
  Y0udf <- as.data.frame(Y0umat)
  colnames(Y0udf) <- Ynames

  # store dataset of observed outcomes for untreated units
  dfu <- cbind.data.frame(G=0,X=Xu,Y0udf)

  # store overall dataset
  df <- rbind.data.frame(dft, dfu)

  # generate id variable
  df$id <- 1:nrow(df)

  # convert data from wide to long format
  ddf <- tidyr::gather(df, period, Y, -G, -X, -id)
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



build_ipw_dataset <- function(sp_list, panel=TRUE) {

  #-----------------------------------------------------------------------------
  # build dataset
  # build things up from the generalized propensity score
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
  
  X <- rnorm(n)

  pr <- exp(outer(X,gamG)) / apply( exp(outer(X,gamG)), 1, sum)
  GG <- apply(pr, 1, function(pvec) sample(seq(0,time.periods), size=1, prob=pvec))

  G <- GG[GG>0]
  nt <- length(G)
  Xt <- X[GG>0]

  # draw individual fixed effect
  Ct <- rnorm(nt, mean=G)

  # generate untreated potential outcomes in each time period
  Ynames <- paste0("Y",1:time.periods)
  Ynames <- paste0(1:time.periods)
  Y0tmat <- sapply(1:time.periods, function(t) {
    thet[t] + Ct + Xt*bett[t] + rnorm(nt)
  })
  Y0tdf <- as.data.frame(Y0tmat)
  
  # generate treated potential outcomes
  Y1tdf <- sapply(1:time.periods, function(t) {
    te.t[t] + te.bet.ind[G]*Ct + Xt*te.bet.X[t] + (G <= t)*te.e[sapply(1:nt, function(i) max(t-G[i]+1,1))] + te + rnorm(nt) # hack for the dynamic effects but ok
  })

  # generate observed data
  Ytdf <- sapply(1:time.periods, function(t) {
    (G<=t)*Y1tdf[,t] + (G>t)*Y0tdf[,t]
  })
  colnames(Ytdf) <- Ynames

  # store observed data for treated group
  dft <- cbind.data.frame(G,X=Xt,Ytdf)

  # untreated units
  nu <- n-nt
  
  # draw untreated covariate
  Xu <- X[GG==0]

  # draw untreated fixed effect
  Cu <- rnorm(nu)


  # generate untreated potential outcomes
  Y0umat <- sapply(1:time.periods, function(t) {
    theu[t] + Cu + rnorm(nu) + Xu*betu[t]
  })
  Y0udf <- as.data.frame(Y0umat)
  colnames(Y0udf) <- Ynames

  # store dataset of observed outcomes for untreated units
  dfu <- cbind.data.frame(G=0,X=Xu,Y0udf)

  # store overall dataset
  df <- rbind.data.frame(dft, dfu)

  # generate id variable
  df$id <- 1:nrow(df)

  # convert data from wide to long format
  ddf <- tidyr::gather(df, period, Y, -G, -X, -id)
  ddf$period <- as.numeric(ddf$period)
  ddf$treat <- 1*(ddf$G > 0)
  ddf <- ddf[order(ddf$id, ddf$period),] # reorder data

  if (!panel) { # repeated cross sections
    Time <- sample(1:time.periods, size=n, replace=TRUE, prob=rep(1/time.periods, time.periods))
    right.row <- sapply( unique(ddf$id), function(i) {
      which(ddf$id==i & ddf$period==Time[i])
    })
    ddf <- ddf[right.row,]
  }
  
  ddf <- subset(ddf, G != 1)
  
  ddf
}


sim <- function(ret=NULL, bstrap=FALSE, cband=FALSE, control.group="nevertreated",
                panel=TRUE) {

  ddf <- build_sim_dataset(panel)
  
  # get results
  res <- att_gt(yname="Y", xformla=~X, data=ddf, tname="period", idname="id",
                gname="G", est_method="reg", 
                bstrap=bstrap, cband=cband, control_group=control.group,
                panel=panel)


  if (is.null(ret)) {
    return(res)
  } else if (ret=="Wpval") {
    rej <- 1*(res$Wpval < .05)
    return(rej)
  } else if (ret=="cband") {
    cu <- res$att + res$c * sqrt(diag(res$V))/sqrt(res$n)
    cl <- res$att - res$c * sqrt(diag(res$V))/sqrt(res$n)
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
