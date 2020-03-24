#-----------------------------------------------------------------------------

# load did code
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
# remotes::install_github("bcallaway11/did", ref="pedro")
remotes::install_github("bcallaway11/DRDID")
remotes::install_github("bcallaway11/BMisc")

#-----------------------------------------------------------------------------
# set parameters
#-----------------------------------------------------------------------------
# number of time periods
T <- 5
# number of treated units
nt <- 1000
# coefficient on X 
bett <- seq(1:T)
# time fixed effect
thet <- seq(1:T)
# treatment effect
# (could make this vary across groups / times / length of exposure / be random
te  <- 1
# number of untreated units
nu <- 1000
# time fixed effect
theu <- thet # changing this creates violations of parallel trends
# covariate effect
betu <- bett # changing this creates violations of conditional parallel trends



sim <- function(ret=NULL, bstrap=FALSE, cband=FALSE) {
  #-----------------------------------------------------------------------------
  # build dataset
  #-----------------------------------------------------------------------------
  
  # randomly assign treated units to groups
  # random groups
  G <- sample(1:T, size=nt, replace=TRUE)
  # fixed groups
  # G <- unlist(lapply(1:T, function(g) rep(g, nt/T)))
  
  # draw a single covariate
  Xt <- rnorm(nt)

  # draw individual fixed effect
  Ct <- rnorm(nt, mean=G)

  # generate untreated potential outcomes in each time period
  Ynames <- paste0("Y",1:T)
  Ynames <- paste0(1:T)
  Y0tmat <- sapply(1:T, function(t) {
    thet[t] + Ct + rnorm(nt) + Xt*bett[t]
  })
  Y0tdf <- as.data.frame(Y0tmat)

  # generate treated potential outcomes
  Y1tdf <- Y0tdf + te

  # generate observed data
  Ytdf <- sapply(1:T, function(t) {
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
  Y0umat <- sapply(1:T, function(t) {
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
                bstrap=bstrap, cband=cband)


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
    rej <- 1*( abs(res$aggte$simple.att / res$aggte$simple.se) > qnorm(.975) )
    return(rej)
  } else {
    return(res)
  }

}


# check if Wald pre-tests are working
biters <- 1000
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="Wpval"), cl=3)
mean( bout )

# check if uniform confidence bands are working
te <- 0
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="cband", bstrap=TRUE, cband=TRUE), cl=3)
mean( bout )

# check if simple att is working
te <- 0
biters <- 1000
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="simple"), cl=3)
mean(bout)



res <- sim()
ggdid(res)
ate <- aggte(res, type="dynamic", balance.e=3)
ggdid(ate)


library(ggplot2)
library(gridExtra)
ggdid(res, ylim=c(0,2))
