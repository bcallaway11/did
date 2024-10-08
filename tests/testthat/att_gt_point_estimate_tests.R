## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
# source("setup_sims.R")
# fldr <- "~/Dropbox/did/R/"
#fldr <- "/Users/santanph/Dropbox/Co-authored Projects/did/R/"
# sapply(paste0(fldr,list.files(fldr)), source)
library(DRDID)
library(BMisc)
# library(ggplot2)
# library(ggpubr)


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, ipw model is incorectly specified here
#-----------------------------------------------------------------------------
set.seed(09142024)
time.periods <- 4
did::reset.sim()
data <- did::build_sim_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr")
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg")
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw")
res

#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, reg model is incorectly specified here
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_ipw_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr")
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg")
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw")
res


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test if 2 period case works (possible to introduce bugs that affect this
# case only)
# Expected results: warning about no pre-treatment periods to test
#-----------------------------------------------------------------------------
time.periods <- 2
did::reset.sim()
data <- did::build_sim_dataset()

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw")
res

summary(aggte(res, type="simple"))
summary(aggte(res, type="group"))
summary(aggte(res, type="dynamic"))
summary(aggte(res, type="calendar"))


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test no covariates case
# Expected Result: te=1, p-value for pre-test uniformly distributed,
#  identical results for different estimation methods
#-----------------------------------------------------------------------------
time.periods <- 4
did::reset.sim()
bett <- betu <- rep(0,time.periods)
data <- did::build_sim_dataset()

res <- att_gt(yname="Y", xformla=~1, data=data, tname="period", idname="id",
              gname="G", est_method="dr")
res

res <- att_gt(yname="Y", xformla=~1, data=data, tname="period", idname="id",
              gname="G", est_method="reg")
res


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test repeated cross sections, regression sims
# Expected result: te=1, p-value for pre-test uniformly distributed
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_sim_dataset(panel=FALSE)

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", panel=FALSE)
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg", panel=FALSE)
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw", panel=FALSE)
res

#-----------------------------------------------------------------------------
# test repeated cross sections, ipw sims
# Expected result: te=1, p-value for pre-test uniformly distributed
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_ipw_dataset(panel=FALSE)

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", panel=FALSE)
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg", panel=FALSE)
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw", panel=FALSE)
res

#-----------------------------------------------------------------------------
# test repeated cross sections, test aggregations
# Expected result: te=length of exposure, p-value for pre-test uniformly distributed
#-----------------------------------------------------------------------------
did::reset.sim()
te.e <- 1:time.periods
data <- did::build_sim_dataset(panel=FALSE)

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", panel=FALSE)
res

summary(aggte(res))
summary(aggte(res, type="dynamic"))
summary(aggte(res, type="group"))
summary(aggte(res, type="calendar"))


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# these are same test cases as for panel data
# but estimate using allow_unbalanced_panel = TRUE
# but setting an id which gives a way to incorporate
# unbalanced panel
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test uniform[0,1]
#-----------------------------------------------------------------------------
time.periods <- 4
did::reset.sim()
data <- did::build_sim_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", allow_unbalanced_panel=TRUE)
res


#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, reg model is incorectly specified here
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_ipw_dataset()

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw", allow_unbalanced_panel=TRUE)
res

#-----------------------------------------------------------------------------
# try it with an actual unbalanced panel
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, reg model is incorectly specified here
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_ipw_dataset()
data <- data[sample(1:nrow(data),  size=floor(.9*nrow(data))),]

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", panel=TRUE, allow_unbalanced_panel=TRUE)
res

#-----------------------------------------------------------------------------
# version that should error
# have to have an idname if you use an unbalanced panel
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_sim_dataset()
data <- data[sample(1:nrow(data),  size=floor(.9*nrow(data))),]

tryCatch(res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname=NULL,
                       gname="G", est_method="reg", panel=TRUE, allow_unbalanced_panel=TRUE),
         error=function(cond) {
           message("expected error:")
           message(cond)
           message("\n")
           return(NA)
         })



## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test not yet treated as control
# Expected result: te=1, p-value for pre-test U[0,1]
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_ipw_dataset(panel=FALSE)

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="notyettreated",
              gname="G", est_method="dr", panel=FALSE)
res


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test not yet treated as control in case w/o never treated group
# Expected result: te=1, p-value for pre-test U[0,1]
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_sim_dataset()
data <- subset(data, G > 0) # drop nevertreated

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="notyettreated",
              gname="G", est_method="dr", panel=FALSE)
res


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test nevertreated as control in case w/o never treated group
# Expected result: te=1, p-value for pre-test U[0,1], error on no never treated
#  units
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_sim_dataset()
data <- subset(data, G > 0) # drop nevertreated

# dr
tryCatch(res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="dr", panel=FALSE),
         error=function(cond) {
           message("expected error:")
           message(cond)
           message("\n")
           return(NA)
         }
         )


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# *test dynamic effects*
# expected result: te=length of exposure
#-----------------------------------------------------------------------------
did::reset.sim()
te <- 0
te.e <- 1:time.periods
data <- did::build_sim_dataset()

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="reg", panel=FALSE)
res
summary(aggte(res, type="dynamic"))


#-----------------------------------------------------------------------------
# test group treatment timing
# Expected result: te constant within group / varies across groups
#-----------------------------------------------------------------------------
did::reset.sim()
te <- 0
te.bet.ind <- 1:time.periods
data <- did::build_ipw_dataset(panel=FALSE)

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="ipw", panel=FALSE)
res
summary(aggte(res, type="group"))


#-----------------------------------------------------------------------------
# test calendar time effects
# expected result: te=time
#-----------------------------------------------------------------------------
did::reset.sim()
te <- 0
te.t <- thet + 1:time.periods
data <- did::build_sim_dataset(panel=FALSE)

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="dr", panel=FALSE)
res
summary(aggte(res, type="calendar"))

#-----------------------------------------------------------------------------
# test balancing with respect to length of exposure
# expected result: balancing fixes treatment effect dynamics
#-----------------------------------------------------------------------------
did::reset.sim()
te <- 0
te.e <- 1:time.periods
te.bet.ind <- 1:time.periods
data <- did::build_sim_dataset()

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="dr", panel=FALSE)
res
summary(aggte(res, type="dynamic"))
summary(aggte(res, type="dynamic", balance_e=1))


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test that att_gt and aggte work with unequally spaced periods
# expected result: te=length of exposure
#-----------------------------------------------------------------------------
time.periods <- 8
did::reset.sim()
te <- 0
te.e <- 1:time.periods
data <- did::build_sim_dataset()
keep.periods <- c(1,2,5,7)
data <- subset(data, G %in% c(0, keep.periods))
data <- subset(data, period %in% keep.periods)

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="reg", panel=FALSE)
res
summary(aggte(res, type="dynamic"))
summary(aggte(res, type="group"))
summary(aggte(res, type="calendar"))


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test that att_gt and aggte work with unequally spaced groups
# expected result: te=length of exposure
#-----------------------------------------------------------------------------
time.periods <- 5
did::reset.sim()
te <- 0
te.e <- 1:time.periods
data <- did::build_sim_dataset()
keep.groups <- c(3,5)
data <- subset(data, G %in% c(0, keep.groups))

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="notyettreated",
              gname="G", est_method="reg", panel=FALSE)
res
summary(aggte(res, type="dynamic", balance_e=0))
summary(aggte(res, type="group"))
summary(aggte(res, type="calendar"))


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test that att_gt works when some units are treated in first period
# expected result: te=length of exposure, code runs with warning message about
#  dropped units
#-----------------------------------------------------------------------------
time.periods <- 4
did::reset.sim()
te <- 1
data <- did::build_sim_dataset()
data <- subset(data, period >= 2)

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="reg", panel=FALSE)
res


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# *test dynamic effects*
# expected result: te=length of exposure
#-----------------------------------------------------------------------------
did::reset.sim()
te <- 0
te.e <- 1:time.periods
data <- did::build_sim_dataset()

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control_group="nevertreated",
              gname="G", est_method="dr", panel=FALSE)
res
summary(aggte(res, type="dynamic", min_e=-1, max_e=1))


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# *test dynamic effects*
# expected result: te=length of exposure - 1 (w/ one period -1 anticipation)
#-----------------------------------------------------------------------------
time.periods <- 5
did::reset.sim()
te <- 0
te.e <- -1:(time.periods-2)
data <- did::build_sim_dataset()
data$G <- data$G + 1 # add anticipation

#-----------------------------------------------------------------------------
# will get results wrong due to anticipation effect for g=6
# which shows up in the comparison group
# this will only affect g=5 with 1 period anticipation
#-----------------------------------------------------------------------------
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              idname="id",
              control_group="nevertreated",
              gname="G", est_method="dr",
              print_details=TRUE,
              anticipation=1
              )
res
summary(aggte(res, type="dynamic"))

#-----------------------------------------------------------------------------
# drop last time period and results should be correct
#-----------------------------------------------------------------------------
data <- subset(data, period < time.periods)
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              idname="id",
              control_group="nevertreated",
              gname="G", est_method="dr",
              print_details=TRUE,
              anticipation=1
              )
res
summary(aggte(res, type="dynamic"))

#-----------------------------------------------------------------------------
# incorrectly ignore anticipation, dynamic effects are incorrect due to ignoring
# anticipation
#-----------------------------------------------------------------------------
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              idname="id",
              control_group="nevertreated",
              gname="G", est_method="dr",
              print_details=TRUE,
              anticipation=0
              )
res
summary(aggte(res, type="dynamic"))


## -----------------------------------------------------------------------------
time.periods <- 4
did::reset.sim()
data <- did::build_sim_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", alp=0.01)
res


## -----------------------------------------------------------------------------
time.periods <- 4
did::reset.sim()
data <- did::build_sim_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", alp=0.01, cband=FALSE)
res


## -----------------------------------------------------------------------------
## some groups later than last treated period
## plus missing groups
time.periods <- 7
did::reset.sim()
data <- did::build_sim_dataset()
data <- subset(data, period <= 4)
missingG_ids <- sample(unique(data$id), size=10)
data[data$id %in% missingG_ids,"G"] <- NA

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", cband=FALSE)
res


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# incorrectly specified id
#-----------------------------------------------------------------------------
time.periods <- 4
did::reset.sim()
data <- did::build_sim_dataset()

# dr
tryCatch(res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="brant",
                       gname="G", est_method="dr"),
         error=function(cond) {
           message("expected error:")
           message(cond)
           message("\n")
           return(NA)
         }
         )


#-----------------------------------------------------------------------------
# incorrectly specified time period
#-----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# custom estimation method
# Expected results: te=1, pre-test p-value uniformly distributed, code runs
#-----------------------------------------------------------------------------
did::reset.sim()
data <- did::build_sim_dataset(panel=TRUE)

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method=DRDID::drdid_imp_panel, panel=TRUE)
res

