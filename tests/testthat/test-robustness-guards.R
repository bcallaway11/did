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
