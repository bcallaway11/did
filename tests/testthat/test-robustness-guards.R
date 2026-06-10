# =============================================================================
# Robustness guards: malformed inputs must be rejected with a clear, identical
# message in BOTH the fast (faster_mode = TRUE) and slow (faster_mode = FALSE)
# paths. Regression tests for two confirmed silent-failure bugs:
#   - duplicate (id, period) rows (slow path used to return wrong numbers)
#   - negative / zero-mean weights (used to silently produce NA)
# =============================================================================

test_that("duplicate (id, period) rows are rejected identically by both modes", {
  set.seed(1)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  ddup <- rbind(d, d[d$id == d$id[1] & d$period == 2, ])  # one duplicated (id, period)
  msg <- "must be unique \\(by tname\\)"
  expect_error(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = ddup,
    tname = "period", idname = "id", gname = "G", faster_mode = TRUE)), msg)
  expect_error(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = ddup,
    tname = "period", idname = "id", gname = "G", faster_mode = FALSE)), msg)
  # The guard is intentionally unconditional on panel: with idname supplied it
  # also rejects duplicated (id, period) rows on repeated cross sections and
  # unbalanced panels (the slow path used to silently accept them there).
  for (fm in c(TRUE, FALSE)) {
    expect_error(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = ddup,
      tname = "period", idname = "id", gname = "G", panel = FALSE,
      faster_mode = fm)), msg)
    expect_error(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = ddup,
      tname = "period", idname = "id", gname = "G", allow_unbalanced_panel = TRUE,
      faster_mode = fm)), msg)
  }
})

test_that("negative / zero-mean weights are rejected identically by both modes", {
  set.seed(1)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  d$wgt <- ifelse(d$id %% 5 == 0, -1, 1)   # some negative weights
  msg <- "must be non-negative with a positive mean"
  expect_error(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
    tname = "period", idname = "id", gname = "G", weightsname = "wgt", faster_mode = TRUE)), msg)
  expect_error(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
    tname = "period", idname = "id", gname = "G", weightsname = "wgt", faster_mode = FALSE)), msg)

  d$wzero <- 0   # all-zero -> non-positive mean
  expect_error(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
    tname = "period", idname = "id", gname = "G", weightsname = "wzero")), msg)
})

test_that("valid positive weights are still accepted in both modes", {
  set.seed(2)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  d$wgt <- runif(nrow(d), 0.5, 2)
  for (fm in c(TRUE, FALSE)) {
    r <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
      tname = "period", idname = "id", gname = "G", weightsname = "wgt",
      faster_mode = fm, bstrap = FALSE)))
    expect_s3_class(r, "MP")
    expect_true(any(is.finite(r$att)))
  }
})

test_that("fast RC path returns NA instead of crashing when a cell has no treated post observations", {
  set.seed(43)
  ids <- 1:90
  d <- expand.grid(id = ids, period = 1:4)
  d$G <- ifelse(d$id <= 40, 0, ifelse(d$id <= 65, 2, 3))
  d <- d[!(d$G == 2 & d$period >= 3), ]
  d$X <- rnorm(nrow(d))
  d$Y <- 0.5 * d$X + 0.2 * d$period + (d$G > 0 & d$period >= d$G) +
    rnorm(nrow(d), 0, 0.2)

  expect_warning(
    res <- suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
      tname = "period", idname = "id", gname = "G", panel = FALSE,
      faster_mode = TRUE, est_method = "dr", bstrap = FALSE)),
    "No units in group 2 in time period 3")
  expect_s3_class(res, "MP")
  expect_true(any(is.na(res$att[res$group == 2])))
})

test_that("never-treated small-group stop fires in both modes when multiple groups are small", {
  set.seed(44)
  d <- expand.grid(id = 1:6, period = 1:3)
  d$G <- ifelse(d$id <= 3, 0, 2)
  d$Y <- rnorm(nrow(d))
  for (j in 1:6) d[[paste0("X", j)]] <- rnorm(nrow(d))
  xf <- ~X1 + X2 + X3 + X4 + X5 + X6
  msg <- "never-treated group is too small"

  expect_error(suppressWarnings(att_gt(yname = "Y", xformla = xf, data = d,
    tname = "period", idname = "id", gname = "G", faster_mode = FALSE,
    bstrap = FALSE)), msg)
  expect_error(suppressWarnings(att_gt(yname = "Y", xformla = xf, data = d,
    tname = "period", idname = "id", gname = "G", faster_mode = TRUE,
    bstrap = FALSE)), msg)
})

test_that("slow unbalanced influence aggregation handles fractional numeric ids", {
  set.seed(45)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  d <- d[-sample(seq_len(nrow(d)), floor(0.05 * nrow(d))), ]
  d$id <- d$id + 0.123456789012345

  slow <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
    data = d, tname = "period", idname = "id", gname = "G",
    allow_unbalanced_panel = TRUE, faster_mode = FALSE, bstrap = FALSE)))
  fast <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
    data = d, tname = "period", idname = "id", gname = "G",
    allow_unbalanced_panel = TRUE, faster_mode = TRUE, bstrap = FALSE)))

  expect_equal(slow$att, fast$att, tolerance = 1e-10)
  expect_equal(slow$se, fast$se, tolerance = 1e-8)
})

test_that("fast path preserves user columns named weights", {
  set.seed(46)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)

  d$weights <- d$Y
  slow_y <- suppressWarnings(suppressMessages(att_gt(yname = "weights", xformla = ~X,
    data = d, tname = "period", idname = "id", gname = "G",
    faster_mode = FALSE, bstrap = FALSE)))
  fast_y <- suppressWarnings(suppressMessages(att_gt(yname = "weights", xformla = ~X,
    data = d, tname = "period", idname = "id", gname = "G",
    faster_mode = TRUE, bstrap = FALSE)))
  expect_equal(slow_y$att, fast_y$att, tolerance = 1e-10)
  expect_false(all(abs(fast_y$att) < 1e-12))

  d$weights <- d$X
  slow_x <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~weights,
    data = d, tname = "period", idname = "id", gname = "G",
    faster_mode = FALSE, bstrap = FALSE)))
  fast_x <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~weights,
    data = d, tname = "period", idname = "id", gname = "G",
    faster_mode = TRUE, bstrap = FALSE)))
  expect_equal(slow_x$att, fast_x$att, tolerance = 1e-10)
  expect_equal(slow_x$se, fast_x$se, tolerance = 1e-8)
})

test_that("slow RC path NA-cells a throwing preliminary logit instead of aborting att_gt", {
  # Regression test: a -Inf covariate (log(0), reachable via transform-formula
  # support) makes overlap_logit_fit() throw. The slow RC branch used to run the
  # overlap/rcond guards OUTSIDE the per-cell tryCatch, hard-aborting the whole
  # att_gt() call while the fast path and the slow panel path degraded to NA
  # cells with a warning. Both modes must now fail identically, cell by cell.
  set.seed(20260609)
  sp <- did::reset.sim(time.periods = 4, n = 400)
  d <- did::build_sim_dataset(sp)
  d$Xpos <- exp(d$X)
  d$Xpos[d$id == unique(d$id)[1]] <- 0  # log(0) = -Inf for one unit

  w_slow <- capture_warnings(
    slow <- suppressMessages(att_gt(yname = "Y", xformla = ~log(Xpos), data = d,
      tname = "period", idname = "id", gname = "G", panel = FALSE,
      est_method = "dr", faster_mode = FALSE, bstrap = FALSE)))
  w_fast <- capture_warnings(
    fast <- suppressMessages(att_gt(yname = "Y", xformla = ~log(Xpos), data = d,
      tname = "period", idname = "id", gname = "G", panel = FALSE,
      est_method = "dr", faster_mode = TRUE, bstrap = FALSE)))

  expect_true(any(grepl("Error computing internal 2x2 DiD", w_slow)))
  expect_true(any(grepl("Error computing internal 2x2 DiD", w_fast)))
  expect_true(any(is.na(slow$att)))      # affected cells degrade to NA
  expect_true(any(is.finite(slow$att)))  # healthy cells still estimated
  expect_equal(is.na(slow$att), is.na(fast$att))
  expect_equal(slow$att, fast$att, tolerance = 1e-10)
  expect_equal(slow$se, fast$se, tolerance = 1e-10)
})

test_that("overlap_logit_fit matches fastglm::fastglm exactly, including the 0.999 decision", {
  # overlap_logit_fit() (fastglmPure with method = 0L, tol = 1e-8, maxit = 100L)
  # replaced fastglm::fastglm(x, y, family = binomial()) at every overlap-check
  # call site with a docstring claim of bit-identical fitted values. Lock that
  # contract -- and the max(fitted) >= 0.999 decision that drives the
  # att = NA-with-warning behavior -- on the degenerate designs where IRLS is
  # most fragile: near separation, complete separation, rank deficiency, and an
  # all-zero column (the sparse-factor-level case).
  set.seed(20260610)
  n <- 200
  x <- rnorm(n)
  designs <- list(
    near_sep = list(X = cbind(1, x), D = rbinom(n, 1, plogis(8 * x))),
    comp_sep = list(X = cbind(1, x), D = as.numeric(x > 0)),
    rank_def = list(X = cbind(1, x, x), D = rbinom(n, 1, plogis(2 * x))),
    zero_col = list(X = cbind(1, x, 0), D = rbinom(n, 1, plogis(2 * x)))
  )
  for (nm in names(designs)) {
    X <- designs[[nm]]$X; D <- designs[[nm]]$D
    # expected "fitted probabilities numerically 0 or 1" noise under separation
    f1 <- suppressWarnings(fastglm::fastglm(X, D, family = binomial()))
    f2 <- suppressWarnings(did:::overlap_logit_fit(X, D))
    expect_identical(f1$fitted.values, f2$fitted.values)
    expect_identical(max(f1$fitted.values) >= 0.999, max(f2$fitted.values) >= 0.999)
  }
})

test_that("slow panel path converts estimator NaN cells to NA like fast mode", {
  set.seed(47)
  d <- expand.grid(id = 1:90, period = 1:4)
  d$G <- ifelse(d$id <= 30, 0, ifelse(d$id <= 60, 2, 3))
  unit_x <- rnorm(90)
  d$X <- unit_x[d$id]
  d$Y <- 0.5 * d$X + 0.2 * d$period + (d$G > 0 & d$period >= d$G) +
    rnorm(nrow(d), 0, 0.2)
  d$w <- ifelse(d$G == 3, 0, 1)

  slow <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
    data = d, tname = "period", idname = "id", gname = "G",
    weightsname = "w", faster_mode = FALSE, bstrap = FALSE)))
  fast <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
    data = d, tname = "period", idname = "id", gname = "G",
    weightsname = "w", faster_mode = TRUE, bstrap = FALSE)))

  expect_true(any(is.na(slow$att[slow$group == 3])))
  expect_false(any(is.nan(slow$att)))
  expect_equal(is.na(slow$att), is.na(fast$att))
})
