# Anticipation periods (g - anticipation <= t < g) use the base period
# g - anticipation - 1 under either base_period and are not used in the pre-test.
# Under base_period = "varying" they used t - 1 (so with anticipation >= 2 they were
# changes between two anticipation periods), and they entered the pre-test, which
# then rejected parallel trends because of the anticipation the user had declared.

# true parallel trends; cohort g has anticipation effect a_eff in the delta periods
# before g and treatment effect tau from g on
make_antic <- function(seed, n_per, periods = 1:6, g = 4, delta = 2, a_eff = 1, tau = 2) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(2 * n_per), time = periods)
  d$g <- ifelse(d$id <= n_per, 0, g)
  d$y <- rnorm(nrow(d)) + 0.5 * d$time + (d$id %% 7) / 10 +
    ifelse(d$g == g & d$time >= g - delta & d$time < g, a_eff, 0) + ifelse(d$g == g & d$time >= g, tau, 0)
  d
}
fit_a <- function(d, ...) suppressWarnings(suppressMessages(
  att_gt(yname = "y", tname = "time", idname = "id", gname = "g", data = d, bstrap = FALSE, cband = FALSE, ...)))

test_that("anticipation periods use the base period g - anticipation - 1 under both base periods", {
  d <- make_antic(seed = 1, n_per = 2000)
  for (delta in c(1, 2)) for (fm in c(TRUE, FALSE)) {
    v <- fit_a(d, anticipation = delta, base_period = "varying", faster_mode = fm)
    u <- fit_a(d, anticipation = delta, base_period = "universal", faster_mode = fm)
    wv <- which(v$t >= v$group - delta); wu <- which(u$t >= u$group - delta)
    expect_gte(length(wv), 3L)
    expect_identical(v$t[wv], u$t[wu])
    expect_equal(v$att[wv], u$att[wu], tolerance = 1e-10)
    expect_equal(v$se[wv], u$se[wu], tolerance = 1e-10)
  }
  # anticipation effect 1 in periods 2-3, treatment effect 2 from period 4
  v <- fit_a(d, anticipation = 2)
  expect_lt(abs(v$att[v$t == 2] - 1), 0.2)
  expect_lt(abs(v$att[v$t == 3] - 1), 0.2)   # was ~0: the change between two anticipation periods
  for (tt in 4:6) expect_lt(abs(v$att[v$t == tt] - 2), 0.2)
})

test_that("the pre-test uses only cells with t < g - anticipation", {
  # periods 1:7, cohort 5: anticipation 2 leaves the single pre-treatment cell (5,2)
  d <- make_antic(seed = 2, n_per = 500, periods = 1:7, g = 5, delta = 2)
  for (fm in c(TRUE, FALSE)) {
    f2 <- fit_a(d, anticipation = 2, faster_mode = fm)
    i <- which(f2$group == 5 & f2$t == 2)
    expect_length(as.numeric(f2$W), 1L)
    expect_equal(as.numeric(f2$W), (f2$att[i] / f2$se[i])^2, tolerance = 1e-8)
    expect_identical(as.numeric(f2$Wpval), round(1 - pchisq(as.numeric(f2$W), df = 1), 5))
    f0 <- fit_a(d, anticipation = 0, faster_mode = fm)
    expect_length(as.numeric(f0$W), 1L)
    expect_identical(as.numeric(f0$Wpval), round(1 - pchisq(as.numeric(f0$W), df = 3), 5))
  }
  # universal base: the normalized base cell (att 0, se NA) is not a pre-treatment estimate
  w <- capture_warnings(suppressMessages(
    fu <- att_gt(yname = "y", tname = "time", idname = "id", gname = "g", data = make_antic(3, 200, periods = 1:6, g = 4),
                 anticipation = 2, base_period = "universal", bstrap = FALSE, cband = FALSE)))
  expect_null(fu$W)
  expect_true(any(grepl("No pre-treatment periods available", w)))
  expect_false(any(grepl("missing or zero variance", w)))
})

test_that("declared anticipation no longer makes the pre-test reject true parallel trends", {
  # periods 1:8, cohort 6, anticipation in periods 4-5: the pre-test uses (6,2), (6,3)
  p_declared <- vapply(1:20, function(s) fit_a(make_antic(s, 300, periods = 1:8, g = 6), anticipation = 2)$Wpval, numeric(1))
  expect_lte(mean(p_declared < 0.05), 0.3)   # was 20/20
  # the same effects, not declared, are pre-trend violations and are still detected
  p_ignored <- vapply(1:20, function(s) fit_a(make_antic(s, 300, periods = 1:8, g = 6), anticipation = 0)$Wpval, numeric(1))
  expect_gte(mean(p_ignored < 0.05), 0.7)
})

test_that("the pretest-only cell filter follows the same rule", {
  fit <- fit_a(make_antic(4, 100, periods = 1:7, g = 5), anticipation = 2, faster_mode = FALSE)
  dp <- fit$DIDparams
  dp$pretreatment_cells_only <- TRUE
  cells <- suppressWarnings(did:::compute.att_gt(dp))$attgt.list
  expect_gte(length(cells), 1L)
  expect_true(all(vapply(cells, function(x) x$year < x$group - 2, logical(1))))
})

test_that("event-study aggregation reports the anticipation effects; the overall ATT averages t >= g only", {
  d <- make_antic(seed = 4, n_per = 2000)
  for (fm in c(TRUE, FALSE)) {
    fit <- fit_a(d, anticipation = 2, faster_mode = fm)
    dyn <- suppressWarnings(aggte(fit, type = "dynamic", cband = FALSE))
    expect_true(all(c(-2, -1, 0, 1, 2) %in% dyn$egt))
    expect_lt(abs(dyn$att.egt[dyn$egt == -1] - 1), 0.2)
    expect_lt(abs(dyn$overall.att - 2), 0.2)
    expect_lt(abs(suppressWarnings(aggte(fit, type = "simple", cband = FALSE))$overall.att - 2), 0.2)
  }
})

test_that("anticipation = 0 keeps the plain definitions (pseudo-ATTs with base t - 1; pre-test over t < g)", {
  d <- make_antic(seed = 5, n_per = 400, delta = 0, a_eff = 0)
  for (fm in c(TRUE, FALSE)) {
    v <- fit_a(d, anticipation = 0, base_period = "varying", faster_mode = fm)
    u <- fit_a(d, anticipation = 0, base_period = "universal", faster_mode = fm)
    expect_false(isTRUE(all.equal(v$att[v$t == 2], u$att[u$t == 2])))
    expect_equal(v$att[v$t >= 4], u$att[u$t >= 4], tolerance = 1e-10)
    expect_length(as.numeric(v$W), 1L)
    expect_identical(as.numeric(v$Wpval), round(1 - pchisq(as.numeric(v$W), df = sum(v$t < v$group)), 5))
  }
})
