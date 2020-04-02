#-----------------------------------------------------------------------------

# load did code
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
# remotes::install_github("bcallaway11/did", ref="pedro")
remotes::install_github("bcallaway11/DRDID")
remotes::install_github("bcallaway11/BMisc")



reset.sim()
sim(panel=FALSE)



#-----------------------------------------------------------------------------
# check if Wald pre-tests are working
#-----------------------------------------------------------------------------
reset.sim()
te <- 0
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="Wpval", panel=FALSE), cl=3)
mean( bout )
# expect to reject about 5% of the time


#-----------------------------------------------------------------------------
# check if uniform confidence bands are working
#-----------------------------------------------------------------------------
reset.sim()
te <- 0
biters <- 100
bout <- pbapply::pbsapply(1:biters, function(b) sim(ret="cband", bstrap=TRUE, cband=TRUE, panel=FALSE), cl=3)
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

res <- sim(panel=FALSE, cband=TRUE, bstrap=TRUE)
res
ggdid(res)#, ylim=c(-2,14))
ate <- aggte(res, type="dynamic")
ggdid(ate)

## p <- ggdid(ate)
## ggsave("~/Downloads/event_study4_balanced5.pdf", p, width=12)

library(ggplot2)
library(gridExtra)
## ggdid(res, ylim=c(0,2))


#-----------------------------------------------------------------------------
# test if 2 period case works (possible to introduce bugs that affect this
# case only)
#-----------------------------------------------------------------------------
reset.sim()
time.periods <- 2
sim()
