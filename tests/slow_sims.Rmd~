
#-----------------------------------------------------------------------------
# check if Wald pre-tests are working
#-----------------------------------------------------------------------------
reset.sim()
te <- 1
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



