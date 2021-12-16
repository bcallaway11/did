library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, ipw model is incorectly specified here
#-----------------------------------------------------------------------------
test_that("att_gt works w/o dynamics, time effects, or group effects", {
  sp <- reset.sim()
  sp$ipw <- FALSE
  data <- build_sim_dataset(sp)

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr")
  # reg
  res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg")


  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_reg$att[1], 1, tol=.5)
})


test_that("att_gt works using ipw", {
  sp <- reset.sim()
  sp$reg <- FALSE
  data <- build_sim_dataset(sp)

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr")


  # ipw
  res_ipw <- att_gt(yname="Y", xformla=~1, data=data, tname="period", idname="id",
                gname="G", est_method="ipw")

  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_ipw$att[1], 1, tol=.5)
})

test_that("two period case", {
  sp <- reset.sim(time.periods=2)
  sp$ipw <- FALSE
  sp$n <- 10000
  data <- build_sim_dataset(sp)
  
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg")
  res
  
  agg_simple <- aggte(res, type="simple")
  agg_group <- aggte(res, type="group")
  agg_dynamic <- aggte(res, type="dynamic")
  agg_calendar <- aggte(res, type="calendar")

  expect_equal(agg_simple$overall.att, 1, tol=.5)
  expect_equal(agg_group$overall.att, 1, tol=.5)
  expect_equal(agg_dynamic$overall.att, 1, tol=.5)
  expect_equal(agg_calendar$overall.att, 1, tol=.5)
})

test_that("no covariates case", {
  time.periods <- 4
  sp <- reset.sim(time.periods=time.periods)
  
  # no effect of covariates
  sp$bett <- sp$betu <- rep(0,time.periods)
  data <- build_sim_dataset(sp)

  res_dr <- att_gt(yname="Y", xformla=~1, data=data, tname="period", idname="id",
                gname="G", est_method="dr")

  res_reg <- att_gt(yname="Y", xformla=~1, data=data, tname="period", idname="id",
                gname="G", est_method="reg")

  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_reg$att[1], 1, tol=.5)
})

test_that("repeated cross section", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp, panel=FALSE)

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=FALSE)

  # reg
  res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", panel=FALSE)

  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_reg$att[1], 1, tol=.5)
})


test_that("ipw repeated cross sections", {
  sp <- reset.sim()
  sp$reg <- FALSE
  sp$n <- 20000 # these are noisy
  data <- build_sim_dataset(sp, panel=FALSE)

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="dr", panel=FALSE)
  

  #ipw
  res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw", panel=FALSE)

  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_ipw$att[1], 1, tol=.5)
})


test_that("repeated cross sections dynamic effects", {
  time.periods <- 4
  sp <- reset.sim(time.periods=time.periods)
  sp$te.e <- 1:time.periods
  data <- build_sim_dataset(sp, panel=FALSE)

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=FALSE)
  

  agg_dynamic <- aggte(res_dr, type="dynamic")
  agg_idx <- agg_dynamic$egt==2

  expect_equal(agg_dynamic$att.egt[agg_idx], 3, tol=.5)
})

test_that("unbalanced panel", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  # drop second row to create unbalanced panel
  data <- data[-2,]

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="dr", allow_unbalanced_panel=TRUE)

  expect_equal(res_dr$att[1], 1, tol=.5)

  expect_warning(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="dr", allow_unbalanced_panel=FALSE))

  # ipw version
  sp <- reset.sim()
  sp$reg <- FALSE
  data <- build_sim_dataset(sp)
  data <- data[-2,]

  res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="ipw", allow_unbalanced_panel=TRUE)

  expect_equal(res_dr$att[1], 1, tol=.5)

  # unbalanced paenl without providing id, should error
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data <- data[sample(1:nrow(data),  size=floor(.9*nrow(data))),]

  expect_error(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname=NULL,
                      gname="G", est_method="reg", panel=TRUE, allow_unbalanced_panel=TRUE))
})

test_that("not yet treated comparison group", {
  sp <- reset.sim()
  sp$reg <- FALSE
  data <- build_sim_dataset(sp, panel=FALSE)

  # dr
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="notyettreated", 
                gname="G", est_method="dr", panel=FALSE)

  expect_equal(res$att[1], 1, tol=.5)

  # no never treated group
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data <- subset(data, G > 0) # drop nevertreated

  # dr
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="notyettreated",
                gname="G", est_method="dr", panel=FALSE)
  expect_equal(res$att[1], 1, tol=.5)


  # try to use never treated group as comparison group, should error
  expect_error(att_gt(yname="Y", xformla=~X, data=data, tname="period",
                      control_group="nevertreated",
                      gname="G", est_method="dr", panel=FALSE))

})

test_that("aggregations", {

  # dynamic effects
  time.periods <- 4
  sp <- reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- build_sim_dataset(sp)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="nevertreated",
                gname="G", est_method="reg", panel=FALSE)
  
  agg_dynamic <- aggte(res, type="dynamic")
  agg_idx <- agg_dynamic$egt==2

  expect_equal(agg_dynamic$att.egt[agg_idx], 2, tol=.5)


  # group effects

  time.periods <- 4
  sp <- reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.bet.ind <- 1:time.periods
  sp$reg <- FALSE
  data <- build_sim_dataset(sp, panel=FALSE)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="notyettreated",
                gname="G", est_method="ipw", panel=FALSE)

  agg_group <- aggte(res, type="group")

  expect_equal(agg_group$att.egt[2], 2*2, tol=.5)


  # calendar time effects
  time.periods <- 4
  sp <- reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.t <- sp$thet + 1:time.periods
  data <- build_sim_dataset(sp, panel=FALSE)
  
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="nevertreated",
                gname="G", est_method="dr", panel=FALSE)
  
  agg_calendar <- aggte(res, type="calendar")
  expect_equal(agg_calendar$att.egt[2], 2, tol=.5)


  # balancing with respect to event time

  sp <- reset.sim()
  sp$te <- 0
  sp$te.e <- 1:time.periods
  sp$te.bet.ind <- 1:time.periods
  data <- build_sim_dataset(sp)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="nevertreated",
                gname="G", est_method="dr", panel=FALSE)

  
  agg_dynamic <- aggte(res, type="dynamic")
  agg_dynamic_balance <- aggte(res, type="dynamic", balance_e=1)

  ad_idx <- which(agg_dynamic$egt == 1)
  adb_idx <- which(agg_dynamic_balance$egt == 1)

  expect_equal(agg_dynamic_balance$att.egt[adb_idx] - agg_dynamic_balance$att.egt[adb_idx-1], 1, tol=.5)
})

test_that("unequally spaced groups", {
  time.periods <- 8
  sp <- reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- build_sim_dataset(sp)
  keep.periods <- c(1,2,5,7)
  data <- subset(data, G %in% c(0, keep.periods))
  data <- subset(data, period %in% keep.periods)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="nevertreated",
                gname="G", est_method="reg", panel=FALSE)

  agg_dynamic <- aggte(res, type="dynamic")
  agg_idx <- agg_dynamic$egt==2

  expect_equal(agg_dynamic$att.egt[agg_idx], 3, tol=.5)

  agg_dynamic_balance <- aggte(res, type="dynamic", balance_e=0)
  agg_idx2 <- which(agg_dynamic_balance$egt==0)
  expect_equal(agg_dynamic_balance$att.egt[agg_idx2], 1, tol=.5)
})

test_that("some units treated in first period", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data <- subset(data, period >= 2)

  expect_warning(att_gt(yname="Y", xformla=~X, data=data, tname="period",
                        control_group="nevertreated",
                        gname="G", est_method="reg", panel=FALSE))
})

test_that("min and max length of exposures", {
  sp <- reset.sim()
  time.periods <- 4
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- build_sim_dataset(sp)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="nevertreated",
                gname="G", est_method="dr", panel=FALSE)
  
  agg_dynamic <- aggte(res, type="dynamic", min_e=-1, max_e=1)
  agg_idx <- which(agg_dynamic$egt == 1)
  expect_equal(agg_dynamic$att.egt[agg_idx], 2, tol=.5)
})


test_that("anticipation", {
  time.periods <- 5
  sp <- reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- -1:(time.periods-2)
  data <- build_sim_dataset(sp)
  data$G <- ifelse(data$G==0, 0, data$G + 1) # add anticipation
  data <- subset(data, G <= time.periods) # drop last period (due to way data is constructed)
  # this will have an anticipation effect=-1, no effect at exposure,
  # and treatment effects increasing by one in subsequent periods

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                idname="id",
                control_group="nevertreated",
                gname="G", est_method="dr",
                anticipation=1
                )

  agg_dynamic <- aggte(res, type="dynamic")
  agg_idx <- which(agg_dynamic$egt==2)
  expect_equal(agg_dynamic$att.egt[agg_idx], 2, tol=.5)

  # incorrectly ignore anticipation
  # causes over-stating treatment effects
  # due to using incorrect base-period
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                idname="id",
                control_group="nevertreated",
                gname="G", est_method="dr",
                anticipation=0
                )

  agg_dynamic <- aggte(res, type="dynamic")
  agg_idx <- which(agg_dynamic$egt==2)
  expect_equal(agg_dynamic$att.egt[agg_idx], 3, tol=.5)

})

test_that("significance level and uniform confidence bands", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)

  # 5% significance level
  set.seed(1234)
  res05 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                  gname="G", est_method="dr", alp=0.05)
  # 1% significance level
  set.seed(1234)
  res01 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                  gname="G", est_method="dr", alp=0.01)
  # 5% pointwise
  set.seed(1234)
  res_pw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                  gname="G", est_method="dr", alp=0.05, cband=FALSE)

  expect_lte(res05$att[1] + res05$c*res05$se[1],
             res01$att[1] + res01$c*res01$se[1])
  expect_gte(res05$att[1] + res05$c*res05$se[1],
             res_pw$att[1] + res_pw$c*res_pw$se[1])

})

test_that("malformed data", {

  # some groups later than last treated period
  # plus missing groups
  time.periods <- 7
  sp <- reset.sim(time.periods=time.periods)
  data <- build_sim_dataset(sp)
  data <- subset(data, period <= 4)
  missingG_ids <- sample(unique(data$id), size=10)
  data[data$id %in% missingG_ids,"G"] <- NA

  expect_warning(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="dr"))
  

  #-----------------------------------------------------------------------------
  # incorrectly specified id
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  
  expect_error(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="brant",
                      gname="G", est_method="dr"))
})

test_that("varying or universal base period", {

  time.periods <- 8
  sp <- reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- build_sim_dataset(sp)
  data <- subset(data, (G<=5) | G==0 )
  # add pre-treatment effects 
  data$G <- ifelse(data$G==0, 0, data$G+3)
  
  
  # dr
  res_varying <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr", base_period="varying")

  agg_dynamic_varying <- aggte(res_varying, type="dynamic")

  res_universal <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                          gname="G", est_method="dr", base_period="universal")

  agg_dynamic_universal <- aggte(res_universal, type="dynamic")

  agg_idx <- which(agg_dynamic_varying$egt == -3)

  expect_equal(agg_dynamic_varying$att.egt[agg_idx], 1, tol=.5)
  expect_equal(agg_dynamic_universal$att.egt[agg_idx], -2, tol=.5)
})


test_that("small groups", {
  # code should still compute in this case (as comparison
  # group is large, but should give a warning about small groups)

  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  # keep only one observation from group 2
  G2_keep_id <- unique(subset(data, G==2)$id)[1]
  data <- subset(data, (G != 2) | (id == G2_keep_id))

  # dr
  expect_warning(res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr"), "small groups")
  # reg
  expect_warning(res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg"), "small groups")
  # ipw
  expect_warning(res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw"), "small groups")

  # estimates will be imprecise for group 2
  idx <- which(res_dr$group == 3 & res_dr$t==3)
  expect_equal(res_dr$att[idx], 1, tol=.5)
  expect_equal(res_reg$att[idx], 1, tol=.5)
})


test_that("small comparison group", {
  # code doesn't run here if use never treated comparison group
  # but should run for all groups except the last one when
  # the not-yet-treated comparison group
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  # keep only one observation from untreated group
  G0_keep_id <- unique(subset(data, G==0)$id)[1]
  data <- subset(data, (G != 0) | (id == G0_keep_id))

  #-----------------------------------------------------------------------------
  # never treated comparison group
  #-----------------------------------------------------------------------------
  # dr
  expect_error(expect_warning(res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", 
              gname="G", est_method="dr"), "small groups"), "never treated group is too small")
  # reg
  expect_error(expect_warning(res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", 
              gname="G", est_method="reg"), "small groups"), "never treated group is too small")
  # ipw
  expect_error(expect_warning(res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", 
              gname="G", est_method="ipw"), "small groups"), "never treated group is too small")

  #-----------------------------------------------------------------------------
  # not-yet-treated comparison group
  #-----------------------------------------------------------------------------
  # dr
  expect_warning(res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated", 
              gname="G", est_method="dr"), "small group")
  # reg
  expect_warning(res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
              gname="G", est_method="reg"), "Not enough control units")
  # ipw
  expect_warning(res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
              gname="G", est_method="ipw"), "overlap condition violated")

  # code should still work for some (g,t)'s
  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_reg$att[1], 1, tol=.5)
  expect_equal(res_ipw$att[1], 1, tol=.5)

  #-----------------------------------------------------------------------------
  # aggregations
  #-----------------------------------------------------------------------------

  agg_dyn <- aggte(res_dr, type="dynamic", na.rm=TRUE)
  agg_group <- aggte(res_reg, type="group", na.rm=TRUE)

  expect_equal(agg_dyn$att.egt[3], 1, tol=.5)
  expect_equal(agg_group$att.egt[1], 1, tol=.5)

  # make sure that standard errors are computed too
  expect_false(is.na(agg_dyn$se.egt[3]))
  expect_false(is.na(agg_group$se.egt[1]))

  skip_if(TRUE, message="known bug about calendar time aggregations, not high priority to fix")
  agg_cal <- aggte(res_ipw, type="calendar", na.rm=TRUE)
  expect_equal(agg_cal$att.egt[1], 1, tol=.5)
  expect_false(is.na(agg_cal$se.egt[1]))
})

test_that("custom estimation method", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method=DRDID::drdid_imp_panel, panel=TRUE)
  expect_equal(res$att[1], 1, tol=.5)  
})



test_that("sampling weights", {
  # the idea here is that we can re-weight and should
  # get the same thing as if we subset
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data2 <- data
  keepids <- sample(unique(data$id), length(unique(data$id)))
  data$w <- 1*(data$id %in% keepids) # weights shouldn't have to have mean/sum 1
  data2 <- subset(data, id %in% keepids)

  res_weights <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                        gname="G", est_method="reg", weightsname="w")

  res_subset <- att_gt(yname="Y", xformla=~X, data=data2, tname="period", idname="id", control_group="notyettreated",
                        gname="G", est_method="reg")

  # test for same att's
  expect_equal(res_weights$att[1], res_subset$att[1])
  # test for same standard errors
  expect_equal(res_weights$se[1], res_subset$se[1], tol=.02)

})

test_that("clustered standard errors", {
  # check that we can compute when clustered standard errors are supplied
  # either as numeric or as factor
  sp <- reset.sim()
  data <- build_sim_dataset(sp)

  data$cluster <- as.numeric(data$cluster)
  res_numeric <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", clustervars="cluster")

  data$cluster <- as.factor(data$cluster)
  res_factor <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", 
                        gname="G", est_method="dr", clustervars="cluster")

  # test for same att's
  expect_equal(res_factor$att[1], res_numeric$att[1])
  # test for same standard errors
  expect_equal(res_factor$se[1], res_numeric$se[1], tol=.02)

  #-----------------------------------------------------------------------------
  # clustered standard errors with unbalanced panel
  data <- data[-3,] # drop one observation
  res_ub <- att_gt(yname="Y",
              tname="period",
              idname="id",
              gname="G",
              xformla=~X,
              data=data,
              panel=TRUE,
              allow_unbalanced_panel=TRUE,
              clustervars="cluster")
  expect_equal(res_ub$att[1], 1, tol=.5)

  #-----------------------------------------------------------------------------
  # also, check that we error when clustering variable varies within unit
  # over time
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data$cluster <- as.numeric(data$cluster)
  data[1,]$cluster <- data[1,]$cluster+1

  expect_error(res_vc <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                        gname="G", est_method="dr", clustervars="cluster"), "handle time-varying cluster variables")

  #-----------------------------------------------------------------------------
  # clustered standard errors with repeated cross sections data
  data <- build_sim_dataset(sp, panel=FALSE)
  res_rc <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                   gname="G", est_method="dr", clustervars="cluster", panel=FALSE)
  expect_equal(res_rc$att[1], 1, tol=.5)
})
