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
te <- 0
set.seed(1814)

## -----------------------------------------------------------------------------
# generate dataset with 4 time periods
time.periods <- 4

# add dynamic effects
te.e <- 1:time.periods

# generate data set with these parameters
data <- build_sim_dataset()

nrow(data)
head(data)

## -----------------------------------------------------------------------------
# estimate group-time average treatment effects using att_gt method
example.attgt <- att_gt(yname="Y",
                        tname="period",
                        idname="id",
                        gname="G",
                        xformla=~X,
                        data=data
                        )

# summarize the results
summary(example.attgt)

## -----------------------------------------------------------------------------
# plot the results
ggdid(example.attgt)

## -----------------------------------------------------------------------------
agg.simple <- aggte(example.attgt, type="simple")
summary(agg.simple)

## -----------------------------------------------------------------------------
agg.es <- aggte(example.attgt, type="dynamic")
summary(agg.es)
ggdid(agg.es)

## -----------------------------------------------------------------------------
agg.gs <- aggte(example.attgt, type="group")
summary(agg.gs)
ggdid(agg.gs)

## -----------------------------------------------------------------------------
agg.ct <- aggte(example.attgt, type="calendar")
summary(agg.ct)
ggdid(agg.ct)

## ----eval=FALSE---------------------------------------------------------------
#  example.attgt.altcontrol <- att_gt(yname="Y",
#                                     tname="period",
#                                     idname="id",
#                                     gname="G",
#                                     xformla=~X,
#                                     data=data,
#                                     control_group="notyettreated"			
#                                     )
#  summary(example.attgt.altcontrol)

## ----eval=FALSE---------------------------------------------------------------
#  example.attgt.reg <- att_gt(yname="Y",
#                              tname="period",
#                              idname="id",
#                              gname="G",
#                              xformla=~X,
#                              data=data,
#                              est_method="reg"
#                              )
#  summary(example.attgt.reg)

## -----------------------------------------------------------------------------
library(did)
data(mpdta)

## -----------------------------------------------------------------------------
head(mpdta)

## -----------------------------------------------------------------------------
# estimate group-time average treatment effects without covariates
mw.attgt <- att_gt(yname="lemp",
                   gname="first.treat",
                   idname="countyreal",
                   tname="year",
                   xformla=~1,
                   data=mpdta,
                   )

# summarize the results
summary(mw.attgt)

# plot the results
# set ylim so that all plots have the same scale along y-axis
ggdid(mw.attgt, ylim=c(-.3,.3))

## -----------------------------------------------------------------------------
# aggregate the group-time average treatment effects
mw.dyn <- aggte(mw.attgt, type="dynamic")
summary(mw.dyn)
ggdid(mw.dyn, ylim=c(-.3,.3))

## -----------------------------------------------------------------------------
mw.dyn.balance <- aggte(mw.attgt, type="dynamic", balance_e=1)
summary(mw.dyn.balance)
ggdid(mw.dyn.balance, ylim=c(-.3,.3))

## ----eval=FALSE---------------------------------------------------------------
#  mw.attgt.X <- att_gt(yname="lemp",
#                     gname="first.treat",
#                     idname="countyreal",
#                     tname="year",
#                     xformla=~lpop,
#                     data=mpdta,
#                     )

