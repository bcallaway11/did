library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)


#-----------------------------------------------------------------------------
#
# These are tests (primarily coming from github) related to bugs
# encountered by users.
#
#-----------------------------------------------------------------------------

test_that("having column named t1 causes code to crash", {
  data(mpdta, package="did")
  out <- att_gt(yname = "lemp",
                gname = "first.treat",
                idname = "countyreal",
                tname = "year",
                xformla = ~1,
                data = mpdta,
                est_method = "reg",
                control_group="notyettreated"
                )
  mpdta$t1 <- 1

  out <- att_gt(yname = "lemp",
                gname = "first.treat",
                idname = "countyreal",
                tname = "year",
                xformla = ~1,
                data = mpdta,
                est_method = "reg",
                control_group="notyettreated"
                )
  expect_false(is.null(out), "code crashed due to stange variable names")
})

test_that("missing covariates", {
  # should warn about missing data but otherwise run
  data(mpdta, package="did")

  mpdta[1, "lpop"] <- NA
  expect_warning(out <- att_gt(yname = "lemp",
                gname = "first.treat",
                idname = "countyreal",
                tname = "year",
                xformla = ~lpop,
                data = mpdta,
                est_method = "reg",
                control_group="notyettreated"
                ))

})


test_that("repeated cross sections small groups with covariates", {
  # from https://github.com/bcallaway11/did/issues/64
  sp <- reset.sim(time.periods=3)
  data <- build_sim_dataset(sp, panel=FALSE)
  data$X2 <- rnorm(nrow(data))
  data$X3 <- rnorm(nrow(data))
  dropids <- unique(subset(data, G==2)$id)
  # interestingly, this seems to return all NA's (which "works")
  # if you only keep one or two observations
  dropids <- dropids[-c(1,2,3)] # keep three observations from group 2
  data <- subset(data, !(id %in% dropids))

  expect_warning(res_dr <- att_gt(yname="Y", xformla=~X+X2+X3, data=data, tname="period", idname="id",
                   gname="G", est_method=DRDID::drdid_rc1, panel=FALSE), "there are some small groups")

  expect_true(is.numeric(res_dr$se[1]))

  skip_if(TRUE, message="known bug, code crashes in this case, fix is probably in DRDID package")
  res_dr <- att_gt(yname="Y", xformla=~X+X2+X3, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=FALSE)

  expect_true(is.numeric(res_dr$se[1]))
})


test_that("fewer time periods than groups", {
  # from https://github.com/bcallaway11/did/issues/56
  # not sure if this is actually a bug,
  # can easily circumvent all of these issues by
  # manually recoding the groups
  time.periods <- 6
  sp <- reset.sim(time.periods=time.periods)
  sp$te <- 0
  sp$te.e <- 1:time.periods
  data <- build_sim_dataset(sp)
  data <- subset(data, !(period %in% c(2,5)))

  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr")
  res_idx <- which(res_dr$group==2 & res_dr$t==3)
  expect_equal(res_dr$att[res_idx], 2, tol=.5)

  dyn_agg <- aggte(res_dr, type="dynamic")
  dyn_idx <- which(dyn_agg$egt==3)
  expect_equal(dyn_agg$att.egt[dyn_idx], 4, tol=.5)
  expect_false(any(is.na(dyn_agg$att.egt)))

  skip_on_cran()
  group_agg <- aggte(res_dr, type="group")
  group_idx <- which(group_agg$egt==3)
  #this seems to fail for groups that are exactly equal to
  #missing time periods -- that seems like a bug!
  expect_equal(group_agg$att.egt[group_idx], mean(c(1,2,4)), tol=.5)
  expect_false(any(is.na(group_agg$att.egt)))

  # calendar aggregation does not compute at all in this case
  # low priority to fix this
  cal_agg <- aggte(res_dr, type="calendar")
  # improve this test if ever get this working
  expect_false(is.na(cal_agg$att.egt[1]))
})
