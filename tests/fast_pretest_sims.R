source("setup_sims.R")
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
remotes::install_github("pedrohcgs/DRDID")
library(BMisc)
library(ggplot2)
library(gridExtra)


#-----------------------------------------------------------------------------
# pre-test simulations
# ...not actually that fast, but not nearly as slow as the other ones!
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# panel with covariates, 2 periods
# Expected result: Throw helpful error cause no periods to test
#-----------------------------------------------------------------------------
reset.sim()
n <- 1000
time.periods <- 2
gamG <- c(0,1:time.periods)/(2*time.periods)
pretest_sim(ret=NULL, bstrap=FALSE, cband=FALSE,
            control.group="nevertreated", panel=TRUE, xformla=~X, cores=3) 


#-----------------------------------------------------------------------------
# panel with covariates, 3 periods
# Note: this case needs to be separately tested from more periods to test
# Expected result: p-value uniformly distributed between 0 and 1
#-----------------------------------------------------------------------------
reset.sim()
n <- 1000
time.periods <- 3
gamG <- c(0,1:time.periods)/(2*time.periods)
pretest_sim(control.group="nevertreated", panel=TRUE, xformla=~X, cores=3) 

#-----------------------------------------------------------------------------
# panel with covariate, 4 periods
# Expected result: p-value uniformly distributed between 0 and 1
#-----------------------------------------------------------------------------
reset.sim()
n <- 1000
time.periods <- 4
gamG <- c(0,1:time.periods)/(2*time.periods)
pretest_sim(ret=NULL, bstrap=FALSE, cband=FALSE,
            control.group="nevertreated", panel=TRUE, xformla=~X, cores=3) 



#-----------------------------------------------------------------------------
# panel w/o covariates, 4 periods
# Expected result: p-value uniformly between 0 and 1
#-----------------------------------------------------------------------------
reset.sim()
n <- 1000
time.periods <- 4
gamG <- c(0,1:time.periods)/(2*time.periods)
bett <- betu <- rep(0,time.periods)
pretest_sim(ret=NULL, 
            control.group="nevertreated", panel=TRUE, cores=3) 


#-----------------------------------------------------------------------------
# repeated cross sections with covariates
# Expected result: p-value uniformly between 0 and 1
#-----------------------------------------------------------------------------
reset.sim()
n <- 1000
time.periods <- 4
gamG <- c(0,1:time.periods)/(2*time.periods)
pretest_sim(xformla=~X, cores=3, panel=FALSE)


#-----------------------------------------------------------------------------
# panel with covariate, 4 periods, parallel trends violated
# Expected result: rej parallel trends
#-----------------------------------------------------------------------------
reset.sim()
n <- 1000
time.periods <- 4
gamG <- c(0,1:time.periods)/(2*time.periods)
betu <- rep(0,time.periods)
pretest_sim(ret=NULL, 
            control.group="nevertreated", panel=TRUE, xformla=~X, cores=3) 



#-----------------------------------------------------------------------------
# panel with covariate, 4 periods, not-yet-treated as control
# Expected result: p-value uniformly between 0 and 1
#-----------------------------------------------------------------------------
reset.sim()
n <- 1000
time.periods <- 4
gamG <- c(0,1:time.periods)/(2*time.periods)
pretest_sim(ret=NULL, 
            control.group="notyettreated", panel=TRUE, xformla=~X, cores=3) 
