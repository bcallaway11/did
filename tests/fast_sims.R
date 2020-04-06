source("setup_sims.R")
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
remotes::install_github("pedrohcgs/DRDID")
library(BMisc)
library(ggplot2)
library(gridExtra)


#-----------------------------------------------------------------------------
# These are simulations that will run fast
# and serve as basic checks that we have not
# introduced any bugs into the code
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, ipw model is incorectly specified here
#-----------------------------------------------------------------------------
reset.sim()
data <- build_sim_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="dr", printdetails=FALSE)
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="reg", printdetails=FALSE)
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="ipw", printdetails=FALSE)
res

#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, reg model is incorectly specified here
#-----------------------------------------------------------------------------
reset.sim()
data <- build_ipw_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="dr", printdetails=FALSE)
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="reg", printdetails=FALSE)
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="ipw", printdetails=FALSE)
res



#-----------------------------------------------------------------------------
# test if 2 period case works (possible to introduce bugs that affect this
# case only)
# Expected results: warning about no pre-treatment periods to test
#-----------------------------------------------------------------------------
reset.sim()
time.periods <- 2
data <- build_sim_dataset()

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="ipw", printdetails=FALSE)
res


#-----------------------------------------------------------------------------
# test no covariates case
# Expected Result: te=1, p-value for pre-test uniformly distributed
#-----------------------------------------------------------------------------
reset.sim()
bett <- betu <- rep(0,time.periods)
data <- build_sim_dataset()

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="reg", printdetails=FALSE)
res


#-----------------------------------------------------------------------------
# test repeated cross sections, regression sims
# Expected result: te=1, p-value for pre-test uniformly distributed
#-----------------------------------------------------------------------------
reset.sim()
data <- build_sim_dataset(panel=FALSE)

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="dr", printdetails=FALSE, panel=FALSE)
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="reg", printdetails=FALSE, panel=FALSE)
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="ipw", printdetails=FALSE, panel=FALSE)
res

#-----------------------------------------------------------------------------
# test repeated cross sections, ipw sims
# Expected result: te=1, p-value for pre-test uniformly distributed
#-----------------------------------------------------------------------------
reset.sim()
data <- build_ipw_dataset(panel=FALSE)

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="dr", printdetails=FALSE, panel=FALSE)
res

# reg
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="reg", printdetails=FALSE, panel=FALSE)
res

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod="ipw", printdetails=FALSE, panel=FALSE)
res

#-----------------------------------------------------------------------------
# test not yet treated as control
# Expected result: te=1, p-value for pre-test U[0,1]
#-----------------------------------------------------------------------------
reset.sim()
data <- build_ipw_dataset(panel=FALSE)

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
              control.group="nevertreated",
              first.treat.name="G", estMethod="dr", printdetails=FALSE, panel=FALSE)
res



#-----------------------------------------------------------------------------
# incorrectly specified id
#-----------------------------------------------------------------------------
reset.sim()
data <- build_sim_dataset()

# dr
res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="brant",
              first.treat.name="G", estMethod="dr", printdetails=FALSE)
res


#-----------------------------------------------------------------------------
# incorrectly specified time period
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# custom estimation method
# Expected results: te=1, pre-test p-value uniformly distributed, code runs
#-----------------------------------------------------------------------------
reset.sim()
data <- build_sim_dataset(panel=TRUE)

res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              first.treat.name="G", estMethod=DRDID::drdid_imp_panel, printdetails=FALSE, panel=TRUE)
res
