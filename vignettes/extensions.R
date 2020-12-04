## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)

## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
library(did)
# Source the currently version of the did package (based on our Dropbox)
fldr <- here::here("R/")
sapply(paste0(fldr,list.files(fldr)), source)
# Source simulation designs
source(here::here("vignettes/setup_sims.R"))

## ----echo=FALSE---------------------------------------------------------------
time.periods <- 5
reset.sim()
te <- 0
te.e <- c(-1,rep(1,time.periods-1))
bett <- betu <- rep(0,time.periods)
#set.seed(1814)

# generate data set with these parameters
dta <- build_sim_dataset()
dta$G <- ifelse(dta$G==0, 0, dta$G+1)
dta$G <- ifelse(dta$G==6, 0, dta$G)

## -----------------------------------------------------------------------------
nrow(dta)
head(dta)

## ---- fig.width=12,fig.height=8, fig.align="center"---------------------------
# estimate group-time average treatment effects using att_gt method
# (and ignoring pre-treatment "dip")
attgt.ignoredip <- att_gt(yname = "Y",
                          tname = "period",
                          idname = "id",
                          gname = "G",
                          xformla = ~1,
                          data = dta,
                          )

# summarize the results
summary(attgt.ignoredip)

# make dynamic effects plot
p <- ggdid(aggte(attgt.ignoredip, "dynamic"))

# add actual treatment effects to the plot
truth <- cbind.data.frame(e = seq(-3,2), att.e = c(0,0,-1,1,1,1))
p <- p + geom_line(data = truth, aes(x = e, y = att.e), inherit.aes = FALSE, color = "blue")
p <- p + geom_point(data = truth, aes(x = e, y = att.e), inherit.aes = FALSE, color = "blue")
p

## -----------------------------------------------------------------------------
compute.attgt <- function(dta) {
  # pick up all groups
  groups <- unique(dta$G)

  # pick up all time periods
  time.periods <- unique(dta$period)

  # sort the groups and drop the untreated group
  groups <- sort(groups)[-1]

  # sort the time periods and drop the first two
  # (can't compute treatment effects for these two
  # periods with one-period anticipation -- w/o anticipation
  # we would just drop one period here)
  time.periods <- sort(time.periods)[-c(1,2)]

  # drop last time period (because we don't know if
  # these units are affected by anticipation or not
  # and we are being very careful)
  # (if you were worried about more than one anticipation
  # period here, would need to drop more time periods
  # from the end)
  time.periods <- time.periods[-length(time.periods)]

  # list to store all group-time average treatment effects
  # that we calculate
  attgt.list <- list()
  counter <- 1

  # loop over all groups
  for (g in groups) {

    # get the correct "base" period for this group
    # (subtract 2 to avoid anticipation)
    main.base.period <- g - 2
    
    # loop over all time periods
    for (tp in time.periods) {

      #----------------------------------------------------
      # if it's a pre-treatment time period (used for the
      # pre-test, we need to adjust the base period)

      # group not treated yet
      if (tp < g) {
        # move two periods before
        base.period <- tp - 2
      } else {
        # this is a post-treatment period
        base.period <- main.base.period
      }
      #----------------------------------------------------



      #----------------------------------------------------
      # now, all we need to do is collect the right subset
      # of the data and estimate a 2x2 DiD

      # get group g and untreated group
      this.data <- subset(dta, G==g | G==0)

      # get current period and base period data
      this.data <- subset(this.data, period==tp | period==base.period)

      # set up to compute 2x2 estimator
      Ypost <- subset(this.data, period==tp)$Y
      Ypre <- subset(this.data, period==base.period)$Y

      # dummy variable being in group g
      G <- 1*(subset(this.data, period==tp)$G == g)

      # compute 2x2 estimator using DRDID package
      # (in this unconditional case, it would be straightforward
      # to calculate the 2x2 att just using averages, but we
      # like the DRDID package as it will work for more complicated
      # cases as well)
      attgt <- DRDID::reg_did_panel(Ypost, Ypre, G, covariates=NULL)$ATT

      # save results
      attgt.list[[counter]] <- list(att=attgt, group=g, time.period=tp)
      counter <- counter+1
      #----------------------------------------------------

    }
  }

  #-----------------------------------------------------------------------------
  # aggregate into dynamic effects
  
  # turn results into a data.frame
  attgt.results <- do.call("rbind.data.frame", attgt.list)

  # add event time to the results
  attgt.results$e <- attgt.results$time.period - attgt.results$group
  
  # calculate relative sizes of each group
  # (will be used as weights)
  n.group <- sapply(groups, function(gg) nrow(subset(dta, G==gg)))
  # merge in group sizes
  ngroup.mat <- cbind(groups, n.group)
  attgt.results <- merge(attgt.results, ngroup.mat, by.x = "group", by.y = "groups")

  # event times to calculate dynamic effects
  eseq <- unique(attgt.results$e) 
  eseq <- sort(eseq)

  # calculate average effects by event time
  att.e <- c()
  counter <- 1
  for (this.e in eseq) {
    # get subset of results at this event time
    res.e <- subset(attgt.results, e==this.e)

    # calculate weights by group size
    res.e$weight <- res.e$n.group / sum(res.e$n.group)

    # calculate dynamic effect as weighted average
    att.e[counter] <- sum(res.e$att * res.e$weight)

    # on to the next one
    counter <- counter+1
  }

  # store dynamic effects results
  dyn.results <- cbind.data.frame(e = eseq, att.e = att.e)

  # return group-time average treatment effects and dynamic effects
  return(list(attgt.results=attgt.results[,c("group","att","time.period")],
              dyn.results=dyn.results))
}

## ----antests, cache=TRUE------------------------------------------------------
anticipation.results <- compute.attgt(dta)
anticipation.results

## ----bootstrap, cache=TRUE----------------------------------------------------
# the number of bootstrap iterations
biters <- 100

# list to store bootstrap results
boot.res <- list()

# loop for each nonparametric bootstrap iteration
for (b in 1:biters) {
  # draw a bootstrap sample; here, we'll call an outside function
  bdata <- BMisc::blockBootSample(dta, "id")

  # call our function for estimating dynamic effects on the
  # bootstrapped data
  boot.res[[b]] <- compute.attgt(bdata)$dyn.results$att.e
}

# use the bootstrapped results to compute standard errors
boot.res <- t(simplify2array(boot.res))
boot.se <- apply(boot.res, 2, sd)

# add the standard errors to the main results
dyn.results <- anticipation.results$dyn.results
dyn.results$att.se <- boot.se

## -----------------------------------------------------------------------------
p <- ggplot(data=dyn.results, aes(x=e, y=att.e)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=(att.e-1.96*att.se), ymax=(att.e+1.96*att.se)), width=0.1) +
  theme_bw()

p <- p + geom_line(data=truth, aes(x=e, y=att.e), inherit.aes=FALSE, color="blue")
p <- p + geom_point(data=truth, aes(x=e, y=att.e), inherit.aes=FALSE, color="blue")
p


