# =============================================================================
# Tests for glance.MP and glance.AGGTEobj S3 methods
# =============================================================================

# Shared setup: generate MP and AGGTEobj results
set.seed(20260401)
sp <- reset.sim()
data_gl <- build_sim_dataset(sp)

mp_slow <- suppressWarnings(suppressMessages(
  att_gt(yname = "Y", xformla = ~X, data = data_gl, tname = "period",
         idname = "id", gname = "G", est_method = "dr",
         bstrap = FALSE, faster_mode = FALSE)
))

mp_fast <- suppressWarnings(suppressMessages(
  att_gt(yname = "Y", xformla = ~X, data = data_gl, tname = "period",
         idname = "id", gname = "G", est_method = "dr",
         bstrap = FALSE, faster_mode = TRUE)
))

agg_types <- c("simple", "dynamic", "group", "calendar")
agg_slow <- lapply(setNames(agg_types, agg_types), function(tp) {
  suppressWarnings(aggte(mp_slow, type = tp))
})
agg_fast <- lapply(setNames(agg_types, agg_types), function(tp) {
  suppressWarnings(aggte(mp_fast, type = tp))
})

# =============================================================================
# glance.MP tests
# =============================================================================

test_that("glance.MP returns a one-row data.frame", {
  gl <- glance(mp_slow)
  expect_equal(nrow(gl), 1)
  expect_s3_class(gl, "data.frame")
})

test_that("glance.MP has expected columns", {
  gl <- glance(mp_slow)
  expected_cols <- c("nobs", "ngroup", "ntime", "control.group", "est.method")
  expect_true(all(expected_cols %in% names(gl)))
})

test_that("glance.MP values are reasonable", {
  gl <- glance(mp_slow)
  expect_true(gl$nobs > 0)
  expect_true(gl$ngroup > 0)
  expect_true(gl$ntime > 0)
  expect_equal(gl$control.group, "nevertreated")
  expect_equal(gl$est.method, "dr")
})

test_that("glance.MP nobs matches nobs.MP", {
  gl <- glance(mp_slow)
  expect_equal(gl$nobs, nobs(mp_slow))
})

test_that("glance.MP works with faster_mode = TRUE", {
  gl <- glance(mp_fast)
  expect_equal(nrow(gl), 1)
  expect_true(gl$nobs > 0)
  expect_true(gl$ngroup > 0)
  expect_true(gl$ntime > 0)
})

test_that("glance.MP agrees between faster_mode settings", {
  gl_slow <- glance(mp_slow)
  gl_fast <- glance(mp_fast)
  expect_equal(gl_slow$nobs, gl_fast$nobs)
  expect_equal(gl_slow$ngroup, gl_fast$ngroup)
  expect_equal(gl_slow$ntime, gl_fast$ntime)
  expect_equal(gl_slow$control.group, gl_fast$control.group)
  expect_equal(gl_slow$est.method, gl_fast$est.method)
})

# =============================================================================
# glance.AGGTEobj tests
# =============================================================================

test_that("glance.AGGTEobj returns a one-row data.frame for all 4 types", {
  for (tp in agg_types) {
    gl <- glance(agg_slow[[tp]])
    expect_equal(nrow(gl), 1, label = paste("nrow for", tp))
    expect_s3_class(gl, "data.frame")
  }
})

test_that("glance.AGGTEobj has expected columns", {
  gl <- glance(agg_slow[["dynamic"]])
  expected_cols <- c("type", "nobs", "ngroup", "ntime", "control.group", "est.method")
  expect_true(all(expected_cols %in% names(gl)))
})

test_that("glance.AGGTEobj type column matches requested type", {
  for (tp in agg_types) {
    gl <- glance(agg_slow[[tp]])
    expect_equal(gl$type, tp, label = paste("type for", tp))
  }
})

test_that("glance.AGGTEobj values are not NULL or NA", {
  for (tp in agg_types) {
    gl <- glance(agg_slow[[tp]])
    expect_false(any(sapply(gl, is.null)), label = paste("no NULLs for", tp))
    expect_false(any(is.na(gl)), label = paste("no NAs for", tp))
  }
})

test_that("glance.AGGTEobj works with faster_mode = TRUE", {
  for (tp in agg_types) {
    gl <- glance(agg_fast[[tp]])
    expect_equal(nrow(gl), 1, label = paste("nrow for", tp, "fast"))
    expect_false(any(sapply(gl, is.null)), label = paste("no NULLs for", tp, "fast"))
  }
})

test_that("glance.MP and glance.AGGTEobj agree on nobs", {
  gl_mp <- glance(mp_slow)
  for (tp in agg_types) {
    gl_agg <- glance(agg_slow[[tp]])
    expect_equal(gl_mp$nobs, gl_agg$nobs, label = paste("nobs for", tp))
  }
})
