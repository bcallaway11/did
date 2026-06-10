# =============================================================================
# Coverage for output/summary methods and exported utilities that were
# previously unexercised: summary.AGGTEobj / print.AGGTEobj (all aggregation
# types + control-group / est-method / cband branches), summary.MP.TEST,
# trimmer(), and the get_wide_data() input guards. These are behaviour tests
# of existing code paths (no new behaviour).
# =============================================================================

test_that("summary.AGGTEobj / print.AGGTEobj cover every aggregation type", {
  set.seed(20260609)
  sp <- did::reset.sim(time.periods = 4)
  data <- did::build_sim_dataset(sp)
  mp <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
           gname = "G", bstrap = TRUE, biters = 100, cband = TRUE,
           control_group = "nevertreated")))
  for (ty in c("simple", "group", "dynamic", "calendar")) {
    agg <- suppressWarnings(aggte(mp, type = ty))
    expect_s3_class(agg, "AGGTEobj")
    expect_output(print(agg))      # print.AGGTEobj delegates to summary.AGGTEobj
    expect_output(summary(agg))    # exercises overall block + per-type egt table
  }
})

test_that("summary.AGGTEobj covers pointwise (non-cband) and not-yet-treated / est_method branches", {
  set.seed(20260609)
  sp <- did::reset.sim(time.periods = 4)
  data <- did::build_sim_dataset(sp)
  # bstrap = FALSE -> the "Pointwise" cband_text1b branch; notyettreated control;
  # cycle est_method to hit the dr / ipw / reg text branches.
  for (em in c("dr", "ipw", "reg")) {
    mp <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", bstrap = FALSE, est_method = em,
             control_group = "notyettreated")))
    agg <- suppressWarnings(aggte(mp, type = "dynamic"))
    expect_output(summary(agg))
  }
})

test_that("summary.MP.TEST prints the conditional pre-test results", {
  skip_on_cran()
  sp <- did::reset.sim(time.periods = 3, n = 200)
  data <- did::build_sim_dataset(sp)
  cdp <- suppressWarnings(suppressMessages(conditional_did_pretest(
    yname = "Y", tname = "period", idname = "id", gname = "G",
    xformla = ~X, data = data)))
  expect_s3_class(cdp, "MP.TEST")
  expect_output(summary(cdp), "Cramer von Mises")
})

test_that("trimmer() runs for both control-group definitions", {
  set.seed(20260609)
  sp <- did::reset.sim(time.periods = 3)
  data <- did::build_sim_dataset(sp)
  g <- 2  # an early treatment group present in the sim
  # trimmer uses print() for any flagged units; swallow it. It returns the
  # flagged ids or NULL when none violate the support condition.
  out_ny <- capture.output(res_ny <- trimmer(g = g, tname = "period", idname = "id",
    gname = "G", xformla = ~X, data = data, control_group = "notyettreated"))
  out_nt <- capture.output(res_nt <- trimmer(g = g, tname = "period", idname = "id",
    gname = "G", xformla = ~X, data = data, control_group = "nevertreated"))
  expect_true(is.null(res_ny) || is.numeric(res_ny) || is.data.frame(res_ny) || is.atomic(res_ny))
  expect_true(is.null(res_nt) || is.numeric(res_nt) || is.data.frame(res_nt) || is.atomic(res_nt))
})

test_that("get_wide_data() guards reject non-data.table and non-2-period input", {
  dt2 <- data.table::data.table(id = c(1, 1, 2, 2), period = c(1, 2, 1, 2), Y = rnorm(4))
  dt3 <- data.table::data.table(id = rep(1:2, each = 3), period = rep(1:3, 2), Y = rnorm(6))
  expect_error(get_wide_data(as.data.frame(dt2), "Y", "id", "period"),
               "data must be a data.table")
  expect_error(get_wide_data(dt3, "Y", "id", "period"),
               "data must have exactly 2 periods")
})
