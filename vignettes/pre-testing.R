## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)

## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
#fldr <- "~/Dropbox/did/R/"
#sapply(paste0(fldr,list.files(fldr)), source)
library(did)
library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)

## ----echo=FALSE---------------------------------------------------------------
source("setup_sims.R")
time.periods <- 4
reset.sim()
bett <- betu <- rep(0, time.periods)
te <- 0
set.seed(1814)

## -----------------------------------------------------------------------------
# generate dataset with 4 time periods
time.periods <- 4

# generate dynamic effects
te.e <- time.periods:1

# generate data set with these parameters
# (main thing: it generates a dataset that satisfies
# parallel trends in all periods...including pre-treatment)
data <- build_sim_dataset()

head(data)

## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# modify the dataset a bit so that we can run an event study
#-----------------------------------------------------------------------------

# generate leads and lags of the treatment
Dtl <- sapply(-(time.periods-1):(time.periods-2), function(l) {
    dtl <- 1*( (data$period == data$G + l) & (data$G > 0) )
    dtl
})
Dtl <- as.data.frame(Dtl)
cnames1 <- paste0("Dtmin",(time.periods-1):1)
colnames(Dtl) <- c(cnames1, paste0("Dt",0:(time.periods-2)))
data <- cbind.data.frame(data, Dtl)
row.names(data) <- NULL

head(data)

#-----------------------------------------------------------------------------
# run the event study regression
#-----------------------------------------------------------------------------

# load plm package
library(plm)

# run event study regression
# normalize effect to be 0 in pre-treatment period
es <- plm(Y ~ Dtmin3 + Dtmin2 + Dt0 + Dt1 + Dt2, data=data, model="within", effect="twoways", index=c("id","period"))

summary(es)

#-----------------------------------------------------------------------------
# make an event study plot
#-----------------------------------------------------------------------------

# some housekeeping for making the plot
# add 0 at event time -1
coefs1 <- coef(es)
ses1 <- sqrt(diag(summary(es)$vcov))
idx.pre <- 1:(time.periods-2)
idx.post <- (time.periods-1):length(coefs1)
coefs <- c(coefs1[idx.pre], 0, coefs1[idx.post])
ses <- c(ses1[idx.pre], 0, ses1[idx.post])
exposure <- -(time.periods-1):(time.periods-2)

cmat <- data.frame(coefs=coefs, ses=ses, exposure=exposure)

library(ggplot2)
library(gridExtra)

ggplot(data=cmat, mapping=aes(y=coefs, x=exposure)) +
  geom_line(linetype="dashed") +
  geom_point() + 
  geom_errorbar(aes(ymin=(coefs-1.96*ses), ymax=(coefs+1.96*ses)), width=0.2) +
  ylim(c(-2,5)) +
  theme_bw()

