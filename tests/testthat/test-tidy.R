
#-----------------------------------------------------------------------------
#
# Tests for tidy(), glance(), and nobs() S3 methods
#
#-----------------------------------------------------------------------------

data(mpdta, package = "did")

set.seed(1234)
mp <- att_gt(
  yname   = "lemp",
  gname   = "first.treat",
  idname  = "countyreal",
  tname   = "year",
  xformla = ~1,
  data    = mpdta,
  bstrap  = FALSE
)

agg_dyn      <- aggte(mp, type = "dynamic",  bstrap = FALSE, cband = FALSE)
agg_group    <- aggte(mp, type = "group",    bstrap = FALSE, cband = FALSE)
agg_calendar <- aggte(mp, type = "calendar", bstrap = FALSE, cband = FALSE)
agg_simple   <- aggte(mp, type = "simple",   bstrap = FALSE, cband = FALSE)

# ---------------------------------------------------------------------------
# tidy.MP
# ---------------------------------------------------------------------------

test_that("tidy.MP returns expected columns", {
  td <- broom::tidy(mp)
  expect_s3_class(td, "data.frame")
  expect_true(all(c("term", "group", "time", "estimate", "std.error",
                    "statistic", "p.value",
                    "conf.low", "conf.high",
                    "point.conf.low", "point.conf.high") %in% names(td)))
})

test_that("tidy.MP statistic equals estimate / std.error", {
  td <- broom::tidy(mp)
  expect_equal(td$statistic, td$estimate / td$std.error)
})

test_that("tidy.MP p.value matches normal approximation", {
  td <- broom::tidy(mp)
  expected <- 2 * (1 - pnorm(abs(td$statistic)))
  expect_equal(td$p.value, expected)
})

test_that("tidy.MP p.value is between 0 and 1", {
  td <- broom::tidy(mp)
  expect_true(all(td$p.value >= 0 & td$p.value <= 1, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# tidy.AGGTEobj
# ---------------------------------------------------------------------------

test_that("tidy.AGGTEobj (dynamic) returns expected columns", {
  td <- broom::tidy(agg_dyn)
  expect_true(all(c("type", "term", "event.time", "estimate", "std.error",
                    "statistic", "p.value",
                    "conf.low", "conf.high",
                    "point.conf.low", "point.conf.high") %in% names(td)))
})

test_that("tidy.AGGTEobj (group) returns expected columns", {
  td <- broom::tidy(agg_group)
  expect_true(all(c("type", "term", "group", "estimate", "std.error",
                    "statistic", "p.value",
                    "conf.low", "conf.high",
                    "point.conf.low", "point.conf.high") %in% names(td)))
})

test_that("tidy.AGGTEobj (calendar) returns expected columns", {
  td <- broom::tidy(agg_calendar)
  expect_true(all(c("type", "term", "estimate", "std.error",
                    "statistic", "p.value",
                    "conf.low", "conf.high",
                    "point.conf.low", "point.conf.high") %in% names(td)))
})

test_that("tidy.AGGTEobj (simple) returns expected columns", {
  td <- broom::tidy(agg_simple)
  expect_true(all(c("type", "term", "estimate", "std.error",
                    "statistic", "p.value",
                    "conf.low", "conf.high",
                    "point.conf.low", "point.conf.high") %in% names(td)))
})

test_that("tidy.AGGTEobj statistic and p.value are consistent", {
  for (agg in list(agg_dyn, agg_group, agg_calendar, agg_simple)) {
    td <- broom::tidy(agg)
    expect_equal(td$statistic, td$estimate / td$std.error)
    expect_equal(td$p.value, 2 * (1 - pnorm(abs(td$statistic))))
  }
})

# ---------------------------------------------------------------------------
# nobs
# ---------------------------------------------------------------------------

test_that("nobs.MP returns number of unique units", {
  n <- stats::nobs(mp)
  expect_type(n, "integer")
  expect_equal(n, length(unique(mpdta$countyreal)))
})

test_that("nobs.AGGTEobj returns number of unique units", {
  for (agg in list(agg_dyn, agg_group, agg_calendar, agg_simple)) {
    n <- stats::nobs(agg)
    expect_equal(n, length(unique(mpdta$countyreal)))
  }
})

test_that("nobs.MP and nobs.AGGTEobj agree", {
  expect_equal(stats::nobs(mp), stats::nobs(agg_dyn))
})
