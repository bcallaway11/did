reset.sim <- function() {
  #-----------------------------------------------------------------------------
  # set parameters
  #-----------------------------------------------------------------------------
  # number of time periods
  # time.periods <<- 4
  # number of treated units
  nt <<- 4000
  # coefficient on X 
  bett <<- seq(1:time.periods)
  # time fixed effect
  thet <<- seq(1:time.periods)
  # number of untreated units
  nu <<- 4000
  # time fixed effect
  theu <<- thet # changing this creates violations of parallel trends
  # covariate effect
  betu <<- bett # changing this creates violations of conditional parallel trends
  #-----------------------------------------------------------------------------
  # parameters for treated potential outcomes
  #-----------------------------------------------------------------------------
  te.bet.ind <<- rep(1,time.periods) # no selective treatment timing
  te.bet.X <<- bett #no heterogeneous effects by X
  te.t <<- thet # no calendar time effects
  te.e <<- rep(0,time.periods) # no dynamic effects
  te <<- 1 # overall basic effect
  #-----------------------------------------------------------------------------
  # extra parameters for an ipw sim (otherwise, these are not used)
  #-----------------------------------------------------------------------------
  n <<- 4000
  # these are parameters in generalized propensity score
  # don't make them too big otherwise can get divide by 0
  gamG <<- c(0,1:time.periods)/(2*time.periods)
}


build_sim_dataset <- function(panel=TRUE) {
  #-----------------------------------------------------------------------------
  # build dataset
  #-----------------------------------------------------------------------------
  
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
  library(tidyr)
  ddf <- gather(df, period, Y, -G, -X, -id)
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



build_ipw_dataset <- function(panel=TRUE) {

  #-----------------------------------------------------------------------------
  # build dataset
  # build things up from the generalized propensity score
  #-----------------------------------------------------------------------------
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
  library(tidyr)
  ddf <- gather(df, period, Y, -G, -X, -id)
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

pretest_sim <- function(ret=NULL, bstrap=FALSE, cband=FALSE,
                        control.group="nevertreated", panel=TRUE, xformla=~X, cores=1) {

  ddf <- build_ipw_dataset(panel=panel)

  # get results
  res <- conditional_did_pretest(yname="Y", xformla=xformla, data=ddf,
                                 tname="period", idname="id",
                                 first.treat.name="G", estMethod="ipw",
                                 printdetails=FALSE,
                                 bstrap=bstrap, cband=cband,
                                 control.group=control.group,
                                 panel=panel, 
                                 pl=TRUE, cores=cores)

  res$CvMpval
}
