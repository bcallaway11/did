#-----------------------------------------------------------------------------

# load did code
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
# remotes::install_github("bcallaway11/did", ref="pedro")
remotes::install_github("bcallaway11/DRDID")

#-----------------------------------------------------------------------------
# build treated groups
#-----------------------------------------------------------------------------

# number of time periods
T <- 10

# number of treated units
nt <- 1000

# randomly assign treated units to groups
G <- sample(1:T, size=nt, replace=TRUE)

# draw a single covariate
Xt <- rnorm(nt)

# draw individual fixed effect
Ct <- rnorm(nt, mean=G)

# coefficient on X 
bett <- seq(1:T)

# time fixed effect
thet <- seq(1:T)

# generate untreated potential outcomes in each time period
Ynames <- paste0("Y",1:T)
Ynames <- paste0(1:T)
Y0tmat <- sapply(1:T, function(t) {
  thet[t] + Ct + rnorm(nt) + Xt*bett[t]
})
Y0tdf <- as.data.frame(Y0tmat)

# treatment effect
# (could make this vary across groups / times / length of exposure
te  <- 1

# generate treated potential outcomes
Y1tdf <- Y0tdf + te

# generate observed data
Ytdf <- sapply(1:T, function(t) {
  (G<=t)*Y1tdf[,t] + (G>t)*Y0tdf[,t]
})
colnames(Ytdf) <- Ynames

# store observed data for treated group
dft <- cbind.data.frame(G,X=Xt,Ytdf)

#-----------------------------------------------------------------------------
# build untreated group
#-----------------------------------------------------------------------------

# number of untreated units
nu <- 1000

# draw untreated covariate
Xu <- rnorm(nu, mean=1)

# draw untreated fixed effect
Cu <- rnorm(nu)

# time fixed effect
theu <- thet # changing this creates violations of parallel trends

# covariate effect
betu <- bett # changing this creates violations of conditional parallel trends

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
ddf <- gather(df, period, Y, as.character(1:10))
ddf$period <- as.numeric(ddf$period)
ddf$treat <- 1*(ddf$G > 0)
ddf <- subset(ddf, G != 1) ## not dropping this group causes code to fail
ddf <- ddf[order(ddf$id, ddf$period),]

# results
##res <- mp.spatt(Y ~ treat, xformla=~X, data=ddf, tname="period",
##                first.treat.name="G")

res <- att_gt(outcome="Y", data=ddf, tname="period", idname="id",
              first.treat.name="G")

library(ggplot2)
library(gridExtra)
ggdid(res)
