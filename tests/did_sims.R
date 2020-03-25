#-----------------------------------------------------------------------------

# load did code
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
# remotes::install_github("bcallaway11/did", ref="pedro")
remotes::install_github("bcallaway11/DRDID")
remotes::install_github("bcallaway11/BMisc")

reset.sim <- function() {
  #-----------------------------------------------------------------------------
  # set parameters
  #-----------------------------------------------------------------------------
  # number of time periods
  time.periods <<- 5
  # number of treated units
  nt <<- 1000
  # coefficient on X 
  bett <<- seq(1:time.periods)
  # time fixed effect
  thet <<- seq(1:time.periods)
  # number of untreated units
  nu <<- 1000
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
}

reset.sim()

sim <- function(ret=NULL, bstrap=FALSE, cband=FALSE, control.group="nevertreated") {
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

  ddf <- subset(ddf, G != 1)

  # get results
  res <- att_gt(yname="Y", xformla=~X, data=ddf, tname="period", idname="id",
                first.treat.name="G", estMethod="reg", printdetails=FALSE,
                bstrap=bstrap, cband=cband, control.group=control.group)


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

#-----------------------------------------------------------------------------
# check if Wald pre-tests are working
#-----------------------------------------------------------------------------
reset.sim()
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="Wpval"), cl=3)
mean( bout )
# expect to reject about 5% of the time


#-----------------------------------------------------------------------------
# check if uniform confidence bands are working
#-----------------------------------------------------------------------------
reset.sim()
te <- 0
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="cband", bstrap=TRUE, cband=TRUE), cl=3)
mean( bout )
# expect to cover all about 95% of the time


#-----------------------------------------------------------------------------
# check if simple att is working
#-----------------------------------------------------------------------------
reset.sim()
te <- 0
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="simple"), cl=3)
mean(bout)
# expect to reject that simple att = 0 about 5% of the time


#-----------------------------------------------------------------------------
# check if dynamic effects are working
#-----------------------------------------------------------------------------
reset.sim()
te <- 0
te.e <- 1:time.periods
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="dynamic"), cl=3)
mean(bout)
# expect to cover the truth about 95% of the time

#-----------------------------------------------------------------------------
# check not yet treated as control
#-----------------------------------------------------------------------------
reset.sim()
te <- 0
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="notyettreated", control.group="notyettreated"))
mean( bout )
# should reject about 5% of the time




## # check if dynamic effects are working
## te.e <- time.periods:1 / (time.periods/4)
## te.bet.ind <- time.periods:1 / (time.periods/2)

res <- sim()
ggdid(res, ylim=c(-2,14))
ate <- aggte(res, type="calendar")
ggdid(ate)

## p <- ggdid(ate)
## ggsave("~/Downloads/event_study4_balanced5.pdf", p, width=12)

library(ggplot2)
library(gridExtra)
## ggdid(res, ylim=c(0,2))
