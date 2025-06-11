library(DRDID)
library(BMisc)
library(tidyverse)
library(data.table)
#library(ggplot2)
#library(ggpubr)


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, ipw model is incorrectly specified here
#-----------------------------------------------------------------------------
test_that("att_gt works w/o dynamics, time effects, or group effects", {
  set.seed(09142024)
  sp <- did::reset.sim()
  sp$ipw <- FALSE
  data <- did::build_sim_dataset(sp)

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
  set.seed(09142024)
  sp <- did::reset.sim()
  sp$reg <- FALSE
  data <- did::build_sim_dataset(sp)

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
  set.seed(09142024)
  sp <- did::reset.sim(time.periods=2)
  sp$ipw <- FALSE
  sp$n <- 10000
  data <- did::build_sim_dataset(sp)

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
  set.seed(09142024)
  time.periods <- 4
  sp <- did::reset.sim(time.periods=time.periods)

  # no effect of covariates
  sp$bett <- sp$betu <- rep(0,time.periods)
  data <- did::build_sim_dataset(sp)

  res_dr <- att_gt(yname="Y", xformla=~1, data=data, tname="period", idname="id",
                gname="G", est_method="dr")

  res_reg <- att_gt(yname="Y", xformla=~1, data=data, tname="period", idname="id",
                gname="G", est_method="reg")

  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_reg$att[1], 1, tol=.5)
})

test_that("repeated cross section", {
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp, panel=FALSE)

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
  set.seed(09142024)
  sp <- did::reset.sim()
  sp$reg <- FALSE
  sp$n <- 20000 # these are noisy
  data <- did::build_sim_dataset(sp, panel=FALSE)

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
  set.seed(09142024)
  time.periods <- 4
  sp <- did::reset.sim(time.periods=time.periods)
  sp$te.e <- 1:time.periods
  data <- did::build_sim_dataset(sp, panel=FALSE)

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=FALSE)


  agg_dynamic <- aggte(res_dr, type="dynamic")
  agg_idx <- agg_dynamic$egt==2

  expect_equal(agg_dynamic$att.egt[agg_idx], 3, tol=.5)
})

test_that("unbalanced panel", {
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  # drop second row to create unbalanced panel
  data <- data[-2,]

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="dr", allow_unbalanced_panel=TRUE)

  expect_equal(res_dr$att[1], 1, tol=.5)

  expect_warning(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="dr", allow_unbalanced_panel=FALSE))

  # ipw version
  set.seed(09142024)
  sp <- did::reset.sim()
  sp$reg <- FALSE
  data <- did::build_sim_dataset(sp)
  data <- data[-2,]

  res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="ipw", allow_unbalanced_panel=TRUE)

  expect_equal(res_dr$att[1], 1, tol=.5)

  # unbalanced paenl without providing id, should error
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  data <- data[sample(1:nrow(data),  size=floor(.9*nrow(data))),]

  expect_error(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname=NULL,
                      gname="G", est_method="reg", panel=TRUE, allow_unbalanced_panel=TRUE))
})

test_that("not yet treated comparison group", {
  set.seed(09142024)
  sp <- did::reset.sim()
  sp$reg <- FALSE
  data <- did::build_sim_dataset(sp, panel=FALSE)

  # dr
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="notyettreated",
                gname="G", est_method="dr", panel=FALSE)

  expect_equal(res$att[1], 1, tol=.5)

  # no never treated group
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  data <- subset(data, G > 0) # drop nevertreated

  # dr
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="notyettreated",
                gname="G", est_method="dr", panel=FALSE)
  expect_equal(res$att[1], 1, tol=.5)


  # try to use never treated group as comparison group, should warn
  expect_warning(nonev_orig <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                      control_group="nevertreated",
                      gname="G", est_method="dr", panel=FALSE, faster_mode = FALSE))

  # try to use never treated group as comparison group with faster mode, should warn
  expect_warning(nonev_faster <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                        control_group="nevertreated",
                        gname="G", est_method="dr", panel=FALSE, faster_mode = TRUE))

  # make use both methods give same ATT(g,t) with no never treated group
  expect_equal(nonev_orig$att, nonev_faster$att)

})

test_that("aggregations", {
  set.seed(09142024)
  # dynamic effects
  time.periods <- 4
  sp <- did::reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- did::build_sim_dataset(sp)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="nevertreated",
                gname="G", est_method="reg", panel=FALSE)

  agg_dynamic <- aggte(res, type="dynamic")
  agg_idx <- agg_dynamic$egt==2

  expect_equal(agg_dynamic$att.egt[agg_idx], 2, tol=.5)


  # group effects
  set.seed(09142024)
  time.periods <- 4
  sp <- did::reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.bet.ind <- 1:time.periods
  sp$reg <- FALSE
  data <- did::build_sim_dataset(sp, panel=FALSE)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="notyettreated",
                gname="G", est_method="ipw", panel=FALSE)

  agg_group <- aggte(res, type="group")

  expect_equal(agg_group$att.egt[2], 2*2, tol=.5)


  # calendar time effects
  set.seed(09142024)
  time.periods <- 4
  sp <- did::reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.t <- sp$thet + 1:time.periods
  data <- did::build_sim_dataset(sp, panel=FALSE)

  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                control_group="nevertreated",
                gname="G", est_method="dr", panel=FALSE)

  agg_calendar <- aggte(res, type="calendar")
  expect_equal(agg_calendar$att.egt[2], 2, tol=.5)


  # balancing with respect to event time
  set.seed(09142024)
  sp <- did::reset.sim()
  sp$te <- 0
  sp$te.e <- 1:time.periods
  sp$te.bet.ind <- 1:time.periods
  data <- did::build_sim_dataset(sp)

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
  set.seed(09142024)
  time.periods <- 8
  sp <- did::reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- did::build_sim_dataset(sp)
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
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  data <- subset(data, period >= 2)

  expect_warning(att_gt(yname="Y", xformla=~X, data=data, tname="period",
                        control_group="nevertreated",
                        gname="G", est_method="reg", panel=FALSE))
})

test_that("min and max length of exposures", {
  set.seed(09142024)
  sp <- did::reset.sim()
  time.periods <- 4
  sp$te <- 0
  sp$te.e <- 1:time.periods
  sp$bett <- sp$betu <- rep(0,time.periods)
  data <- did::build_sim_dataset(sp)

  res <- att_gt(yname="Y", xformla=~1, data=data, tname="period",
                idname="id",
                control_group="nevertreated",
                gname="G", est_method="reg", panel=TRUE,
                base_period="varying")

  agg_dynamic <- aggte(res, type="dynamic", min_e=-1, max_e=1)
  agg_idx <- which(agg_dynamic$egt == 1)
  expect_equal(agg_dynamic$att.egt[agg_idx], 2, tol=.5)
})


test_that("anticipation", {
  set.seed(09142024)
  time.periods <- 5
  sp <- did::reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- -1:(time.periods-2)
  data <- did::build_sim_dataset(sp)
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


  # test for previous bug when using anticipation and
  # a notyettreated comparison group
  data <- subset(data, G != 0)
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period",
                idname="id",
                control_group="notyettreated",
                gname="G", est_method="dr",
                anticipation=1
                )

  agg_dynamic <- aggte(res, type="dynamic")
  agg_idx <- which(agg_dynamic$egt==0)
  expect_equal(agg_dynamic$att.egt[agg_idx], 0, tol=.3)
})

test_that("significance level and uniform confidence bands", {
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)

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
  set.seed(09142024)
  # some groups later than last treated period
  # plus missing groups
  time.periods <- 7
  sp <- did::reset.sim(time.periods=time.periods)
  data <- did::build_sim_dataset(sp)
  data <- subset(data, period <= 4)
  missingG_ids <- sample(unique(data$id), size=10)
  data[data$id %in% missingG_ids,"G"] <- NA

  expect_warning(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method="dr"))


  #-----------------------------------------------------------------------------
  # incorrectly specified id
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)

  expect_error(att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="brant",
                      gname="G", est_method="dr"))
})

test_that("varying or universal base period", {
  set.seed(09142024)
  time.periods <- 8
  sp <- did::reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- did::build_sim_dataset(sp)
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
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
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
  set.seed(09142024)
  # code doesn't run here if use never treated comparison group
  # but should run for all groups except the last one when
  # the not-yet-treated comparison group
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
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
  # not-yet-treated comparison group with faster_mode = TRUE
  #-----------------------------------------------------------------------------

  # dr
  expect_error(expect_warning(res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                                               gname="G", est_method="dr", faster_mode = TRUE), "matrix is singular"), "faster_mode=FALSE")
  # reg
  expect_error(expect_warning(res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                                                gname="G", est_method="reg", faster_mode = TRUE), "matrix is singular"), "faster_mode=FALSE")
  # ipw
  expect_warning(res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                                   gname="G", est_method="ipw", faster_mode = TRUE), "small groups")

  #-----------------------------------------------------------------------------
  # not-yet-treated comparison group
  #-----------------------------------------------------------------------------
  # dr
  expect_warning(res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
              gname="G", est_method="dr", faster_mode = FALSE), "small group")
  # reg
  expect_warning(res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
              gname="G", est_method="reg", faster_mode = FALSE), "Not enough control units")
  # ipw
  expect_warning(res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
              gname="G", est_method="ipw", faster_mode = FALSE), "overlap condition violated")

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
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  res <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                gname="G", est_method=DRDID::drdid_imp_panel, panel=TRUE)
  expect_equal(res$att[1], 1, tol=.5)
})



test_that("sampling weights", {
  set.seed(09142024)
  # the idea here is that we can re-weight and should
  # get the same thing as if we subset
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
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
  set.seed(09142024)
  # check that we can compute when clustered standard errors are supplied
  # either as numeric or as factor
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)

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
              faster_mode = FALSE,
              clustervars="cluster")
  expect_equal(res_ub$att[1], 1, tol=.5)

  # clustered standard errors with unbalanced panel with faster mode
  res_ub_faster <- att_gt(yname="Y",tname="period",idname="id",gname="G",
              xformla=~X,data=data,panel=TRUE, allow_unbalanced_panel=TRUE,
              faster_mode = TRUE, clustervars="cluster")
  expect_equal(res_ub_faster$att[1], 1, tol=.5)

  #-----------------------------------------------------------------------------
  # also, check that we error when clustering variable varies within unit
  # over time
  set.seed(09142024)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  data$cluster <- as.numeric(data$cluster)
  data[1,]$cluster <- data[1,]$cluster+1

  expect_error(res_vc <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                        gname="G", est_method="dr", clustervars="cluster"), "handle time-varying cluster variables")

  #-----------------------------------------------------------------------------
  # clustered standard errors with repeated cross sections data
  data <- did::build_sim_dataset(sp, panel=FALSE)
  res_rc <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", control_group="notyettreated",
                   gname="G", est_method="dr", clustervars="cluster", panel=FALSE)
  expect_equal(res_rc$att[1], 1, tol=.5)
})

test_that("faster mode enabled for panel data", {
  data <- did::mpdta
  out <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "universal",
                control_group = "nevertreated", est_method = "dr", faster_mode = FALSE)
  out2 <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "universal",
                control_group = "nevertreated", est_method = "dr", faster_mode = TRUE)

  # check if results are equal.
  expect_equal(out$att, out2$att)
  expect_equal(out$se, as.numeric(out2$se))
  # --------------------------------------------------------------------------------------------------------
  out <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "varying",
                control_group = "nevertreated", est_method = "dr", faster_mode = FALSE)
  out2 <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                 xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "varying",
                 control_group = "nevertreated", est_method = "dr", faster_mode = TRUE)

  # check if results are equal.
  expect_equal(out$att, out2$att)
  expect_equal(out$se, as.numeric(out2$se))

  # --------------------------------------------------------------------------------------------------------
  out <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "varying",
                control_group = "notyettreated", est_method = "dr", faster_mode = FALSE)
  out2 <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                 xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "varying",
                 control_group = "notyettreated", est_method = "dr", faster_mode = TRUE)

  # check if results are equal.
  expect_equal(out$att, out2$att)
  expect_equal(out$se, as.numeric(out2$se))

  # --------------------------------------------------------------------------------------------------------
  out <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "universal",
                control_group = "notyettreated", est_method = "dr", faster_mode = FALSE)
  out2 <- att_gt(yname = "lemp", gname = "first.treat", idname = "countyreal", tname = "year",
                 xformla = ~1, data = data, bstrap = FALSE, cband = FALSE, base_period = "universal",
                 control_group = "notyettreated", est_method = "dr", faster_mode = TRUE)

  # check if results are equal.
  expect_equal(out$att, out2$att)
  expect_equal(out$se, as.numeric(out2$se))

})

test_that("faster model enabled for repeated cross sectional data", {

  data_rcs <- as.data.table(did::build_sim_dataset(reset.sim(time.periods=4, n=1000), panel=FALSE))
  data_rcs$period <- as.integer(data_rcs$period)
  data_rcs[G == 0, G := Inf]

  # ----------------------------------------------------------------------------------------------------------------------------------
  out_rcs <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                    bstrap = FALSE, cband = FALSE, base_period = "universal", control_group = "nevertreated",
                    est_method = "dr", faster_mode = FALSE)

  out_rcs2 <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                      bstrap = FALSE, cband = FALSE, base_period = "universal", control_group = "nevertreated",
                      est_method = "dr", faster_mode = TRUE )

  # check if results are equal.
  expect_equal(out_rcs$att, out_rcs2$att)
  expect_equal(out_rcs$se, as.numeric(out_rcs2$se))

  # ----------------------------------------------------------------------------------------------------------------------------------
  out_rcs <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                    bstrap = FALSE, cband = FALSE, base_period = "varying", control_group = "nevertreated",
                    est_method = "dr", faster_mode = FALSE)

  out_rcs2 <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                      bstrap = FALSE, cband = FALSE, base_period = "varying", control_group = "nevertreated",
                      est_method = "dr", faster_mode = TRUE )

  # check if results are equal.
  expect_equal(out_rcs$att, out_rcs2$att)
  expect_equal(out_rcs$se, as.numeric(out_rcs2$se))

  # ----------------------------------------------------------------------------------------------------------------------------------
  out_rcs <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                    bstrap = FALSE, cband = FALSE, base_period = "varying", control_group = "notyettreated",
                    est_method = "dr", faster_mode = FALSE)

  out_rcs2 <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                      bstrap = FALSE, cband = FALSE, base_period = "varying", control_group = "notyettreated",
                      est_method = "dr", faster_mode = TRUE )

  # check if results are equal.
  expect_equal(out_rcs$att, out_rcs2$att)
  expect_equal(out_rcs$se, as.numeric(out_rcs2$se))

  # ----------------------------------------------------------------------------------------------------------------------------------
  out_rcs <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                    bstrap = FALSE, cband = FALSE, base_period = "universal", control_group = "notyettreated",
                    est_method = "dr", faster_mode = FALSE)

  out_rcs2 <- att_gt(yname = "Y", gname = "G", tname = "period", xformla = ~1, data = data_rcs, panel = FALSE,
                      bstrap = FALSE, cband = FALSE, base_period = "universal", control_group = "notyettreated",
                      est_method = "dr", faster_mode = TRUE )

  # check if results are equal.
  expect_equal(out_rcs$att, out_rcs2$att)
  expect_equal(out_rcs$se, as.numeric(out_rcs2$se))

})

test_that("faster model enabled for unbalanced panel data and time-varying covariates", {

  set.seed(05202025)
  # balanced panel dimensions
  N  <- 1000        # individuals
  TT <- 6           # time periods (t = 1, …, 6)
  # create id–time grid
  dt <- CJ(id = 1:N, t = 1:TT)   # CJ:=Cartesian join == balanced panel
  # assign cohort (first time treated) or never-treated
  # 40% treated in period 4, 30% in 5, 10% in 6, rest never treated
  cohort_vals  <- c(0, 4, 5, 6)                     # 0 = never treated
  cohort_probs <- c(1 - .40 - .30 - .10, .40, .30, .10)

  dt[, g := sample(
    cohort_vals,
    size = 1,          # one draw per id → time-invariant
    prob = cohort_probs),
    by = id]  # group by id

  # treatment indicator: D = 1 if t ≥ g and g > 0
  dt[, D := as.integer(t >= g & g > 0)]

  # four time-varying covariates
  dt[, `:=`(
    x1 = rnorm(.N),
    x2 = runif(.N, 0, 10) + 0.1 * t,
    x3 = rnorm(N)[id] + 0.2 * t + rnorm(.N),
    x4 = sin(2 * pi * t / TT) + rnorm(.N)
  )]

  # outcome model: baseline + covariate effects + treatment effect
  # note: just for debugging purposes, not a well-defined DGP for theoretical results!!
  beta  <- c(1.5, -0.8,  0.4, 2.0)          # coefficients on x1–x4
  tau   <- 3                                # true treatment effect
  alpha_i <- rnorm(N)[dt$id]                # id fixed effect
  gamma_t <- seq(-1, 1, length.out = TT)[dt$t]  # common time trend

  dt[, y := alpha_i + gamma_t +
       beta[1]*x1 + beta[2]*x2 + beta[3]*x3 + beta[4]*x4 +
       tau * D + rnorm(.N, sd = 2)]

  # outcome model: baseline + covariate effects + treatment effect
  # note: just for debugging purposes, not a well-defined DGP for theoretical results!!
  beta  <- c(1.5, -0.8,  0.4, 2.0)          # coefficients on x1–x4
  tau   <- 3                                # true treatment effect
  alpha_i <- rnorm(N)[dt$id]                # id fixed effect
  gamma_t <- seq(-1, 1, length.out = TT)[dt$t]  # common time trend

  dt[, y := alpha_i + gamma_t +
       beta[1]*x1 + beta[2]*x2 + beta[3]*x3 + beta[4]*x4 +
       tau * D + rnorm(.N, sd = 2)]

  #################################################################
  ##  Make the panel UNBALANCED
  ##  drop 25 % of rows at random, but leave ≥1 row per id
  #################################################################
  drop_rate <- 0.25
  dt_unbal <- dt[, .SD[runif(.N) > drop_rate], by = id]     # keep rows with prob 0.75
  dt_unbal  <- dt_unbal[, if (.N > 0) .SD, by = id]         # (safety) ids with ≥1 obs

  # ----------------------------------------------------------------------------------------------------------------------------------
  # let's run did with covariates using faster_mode=FALSE
  att_slower <- att_gt(
    yname   = "y",
    tname   = "t",
    idname  = "id",
    gname   = "g",
    xformla = ~ x1 + x2 + x3 + x4,   # four time-varying covariates
    data    = dt,
    base_period = "universal", # "universal" or "varying"
    panel   = TRUE,
    faster_mode = FALSE,
    bstrap = FALSE,
  )

  # now let's run did with covariates using faster_mode=TRUE
  att_faster <- att_gt(
    yname   = "y",
    tname   = "t",
    idname  = "id",
    gname   = "g",
    xformla = ~ x1 + x2 + x3 + x4,   # four time-varying covariates
    data    = dt,
    panel   = TRUE,
    base_period = "universal", # "universal" or "varying"
    faster_mode = TRUE,
    bstrap = FALSE,
  )

  # check if results are equal.
  expect_equal(att_slower$att, att_faster$att)
  expect_equal(att_slower$se, as.numeric(att_faster$se))

  # get event study estimates
  out1 = att_slower %>%
    aggte(type = "dynamic",  cband = FALSE, bstrap = FALSE)

  out2 = att_faster %>%
    aggte(type = "dynamic", cband = FALSE, bstrap = FALSE)

  # check if results are equal.
  expect_equal(out1$att.egt, out2$att.egt)
  expect_equal(out1$se.egt, as.numeric(out2$se.egt))

  # ----------------------------------------------------------------------------------------------------------------------------------
  # running the same but with unbalanced panel data
  att_slow <- att_gt(
    yname   = "y",
    tname   = "t",
    idname  = "id",
    gname   = "g",
    xformla = ~ x1 + x2 + x3 + x4,
    data    = dt_unbal,
    panel   = TRUE,
    allow_unbalanced_panel = TRUE,
    faster  = FALSE,
    bstrap = FALSE
  )

  expect_message(att_fast <- att_gt(
    yname   = "y",
    tname   = "t",
    idname  = "id",
    gname   = "g",
    xformla = ~ x1 + x2 + x3 + x4,
    data    = dt_unbal,
    panel   = TRUE,
    allow_unbalanced_panel = TRUE,
    faster  = TRUE,
    bstrap = FALSE
  ), "unbalanced panel")

  # check if results are equal.
  expect_equal(att_slow$att, att_fast$att)
  expect_equal(att_slow$se, as.numeric(att_fast$se))


  # get event study estimates
  out1 = att_slow %>%
    aggte(type = "dynamic",  cband = FALSE, bstrap = FALSE)


  out2 = att_fast %>%
    aggte(type = "dynamic", cband = FALSE, bstrap = FALSE)

  # check if results are equal.
  expect_equal(out1$att.egt, out2$att.egt)
  expect_equal(out1$se.egt, out2$se.egt, tol=.0005)

})


test_that("faster_mode = TRUE matches baseline on filtered sim dataset when there are not subsequent cohort and time periods", {
  # simulate full panel
  set.seed(09142024)
  sp <- reset.sim()
  dt <- build_sim_dataset(sp)

  # filter down to two periods
  # here we know build_sim_dataset() has periods 1:4, so pick 2 & 4
  # cohorts -> [3,4], periods -> [2,4]
  dt2 <- dt[dt$period %in% c(2, 4), ]

  # run att_gt with both modes (no errors)
  expect_warning({
    res_slow <- att_gt(
      yname         = "Y",
      tname         = "period",
      idname        = "id",
      gname         = "G",
      data          = dt2,
      panel         = TRUE,
      control_group = "nevertreated",
      xformla       = NULL,
      est_method    = "dr",
      base_period   = "universal",
      faster_mode   = FALSE
    )
    res_fast <- att_gt(
      yname         = "Y",
      tname         = "period",
      idname        = "id",
      gname         = "G",
      data          = dt2,
      panel         = TRUE,
      control_group = "nevertreated",
      xformla       = NULL,
      est_method    = "dr",
      base_period   = "universal",
      faster_mode   = TRUE
    )
  }, "Dropped 999 units that were already treated in the first period")


  # they should have the same length and (within tol) the same values
  expect_length(res_slow$att, 4)
  expect_length(res_fast$att, 4)
  expect_equal(res_slow$att, res_fast$att, tolerance = 1e-8)
})
