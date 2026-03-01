library(DRDID)
library(BMisc)
# library(ggplot2)
# library(ggpubr)


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
  expect_false(is.null(out), "code crashed due to strange variable names")
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
  sp <- did::reset.sim(time.periods=3)
  data <- build_sim_dataset(sp, panel=FALSE)
  data$X2 <- rnorm(nrow(data))
  data$X3 <- rnorm(nrow(data))
  dropids <- unique(subset(data, G==2)$id)
  # interestingly, this seems to return all NA's (which "works")
  # if you only keep one or two observations
  dropids <- dropids[-c(1,2,3)] # keep three observations from group 2
  data <- subset(data, !(id %in% dropids))

  expect_warning(res_dr <- att_gt(yname="Y", xformla=~X+X2+X3, data=data, tname="period", idname="id",
                   gname="G", est_method=DRDID::drdid_rc1, panel=FALSE), "very few observations")

  expect_true(is.numeric(res_dr$se[1]))

  skip("known bug, code crashes in this case, fix is probably in DRDID package")
})


test_that("fewer time periods than groups", {
  # from https://github.com/bcallaway11/did/issues/56
  # not sure if this is actually a bug,
  # can easily circumvent all of these issues by
  # manually recoding the groups
  time.periods <- 6
  sp <- did::reset.sim(time.periods=time.periods)
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

  group_agg <- aggte(res_dr, type="group")
  group_idx <- which(group_agg$egt==3)
  expect_equal(group_agg$att.egt[group_idx], mean(c(1,2,4)), tol=.5)
  expect_false(any(is.na(group_agg$att.egt)))

  cal_agg <- aggte(res_dr, type="calendar")
  expect_false(is.na(cal_agg$att.egt[1]))
})


test_that("0 pre-treatment estimates when outcomes are 0", {
  # from https://github.com/bcallaway11/did/issues/126
  sp <- did::reset.sim(time.periods=10)
  data <- build_sim_dataset(sp)
  data <- subset(data, G != 0) # drop never treated
  data <- subset(data, G > 6)
  data <- subset(data, period > 5)
  data$Y[(data$period < data$G)] <- 0 # set pre-treatment = 0
  res <- att_gt(yname="Y",
                tname="period",
                idname="id",
                gname="G",
                data=data,
                control_group = "notyettreated",
                base_period="universal")
  res_idx <- which(res$group==9 & res$t==7)
  expect_equal(res$att[res_idx],0)

  res <- att_gt(yname="Y",
                tname="period",
                idname="id",
                gname="G",
                data=data,
                control_group = "notyettreated",
                base_period="varying")
  res_idx <- which(res$group==9 & res$t==7)
  expect_equal(res$att[res_idx],0)
})

test_that("variables not in dataset", {
  sp <- did::reset.sim(time.periods=3)
  data <- build_sim_dataset(sp)

  X2  <- factor(data$cluster)

  expect_error(att_gt(yname="Y", xformla=~X2, data=data, tname="period", idname="id", control_group="notyettreated",
                      gname="G", est_method="dr", clustervars="cluster"), " variables are not in data")

})

test_that("groups treated after max(t) but within anticipation window are not coerced to never-treated", {
  # Bug: with anticipation > 0, groups treated shortly after the last observed

  # period may have anticipatory effects during the observed window.
  # These should NOT be coerced to never-treated (and used as controls).
  #
  # Example: tlist = [1,2,3,4,5], gname = 6, anticipation = 2
  # Anticipation effects start at period 6 - 2 = 4 (within observed data).
  # If coerced to never-treated, these contaminated units serve as controls.

  set.seed(20250228)
  n <- 600
  dt <- data.frame(
    id = rep(1:n, each = 5),
    time = rep(1:5, n)
  )
  # Groups: 0 (never-treated), 4 (treated at 4), 6 (treated at 6, after max(t)=5)
  group_vals <- c(rep(0, 200), rep(4, 200), rep(6, 200))
  dt$group <- group_vals[dt$id]
  dt$y <- rnorm(nrow(dt))

  # With anticipation = 0: group 6 is beyond max(t)=5, should be coerced to never-treated
  suppressWarnings({
    res_ant0 <- att_gt(yname = "y", gname = "group", idname = "id", tname = "time",
                       data = dt, anticipation = 0, bstrap = FALSE,
                       faster_mode = FALSE, print_details = FALSE,
                       control_group = "nevertreated")
  })

  # With anticipation = 2: group 6 anticipates at period 4, should NOT be coerced
  # to never-treated. Instead, group 6 should remain as a treated group.
  suppressWarnings({
    res_ant2 <- att_gt(yname = "y", gname = "group", idname = "id", tname = "time",
                       data = dt, anticipation = 2, bstrap = FALSE,
                       faster_mode = FALSE, print_details = FALSE,
                       control_group = "nevertreated")
  })

  # With anticipation = 0, group 6 is coerced to never-treated (only group 4 in results)
  expect_true(all(res_ant0$group == 4))

  # With anticipation = 2, group 6 should NOT be in the control group.
  # It may or may not appear in the results as a treated group (depending on
  # whether there are enough pre-treatment periods), but the key check is
  # that the number of control units differs between the two runs.
  # Group 6 units (200 units) should not be part of the never-treated pool.

  # Run the same test with faster_mode = TRUE for consistency
  suppressWarnings({
    res_ant0_fast <- att_gt(yname = "y", gname = "group", idname = "id", tname = "time",
                            data = dt, anticipation = 0, bstrap = FALSE,
                            faster_mode = TRUE, print_details = FALSE,
                            control_group = "nevertreated")
  })
  suppressWarnings({
    res_ant2_fast <- att_gt(yname = "y", gname = "group", idname = "id", tname = "time",
                            data = dt, anticipation = 2, bstrap = FALSE,
                            faster_mode = TRUE, print_details = FALSE,
                            control_group = "nevertreated")
  })

  # Both paths should agree
  expect_equal(res_ant0$att, res_ant0_fast$att)
  expect_equal(res_ant2$att, res_ant2_fast$att)

  # With anticipation = 2, group 6 should appear as a treated group in results
  # (it has enough pre-treatment periods: periods 1,2,3 are all before 6-2=4)
  expect_true(6 %in% res_ant2$group)
  expect_true(6 %in% res_ant2_fast$group)
})

