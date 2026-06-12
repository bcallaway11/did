library(testthat)

# ============================================================
# Shared helpers
# ============================================================

make_panel_cov <- function(n = 60, n_periods = 6, seed = 42) {
  set.seed(seed)
  ids     <- rep(1:n, each = n_periods)
  times   <- rep(1:n_periods, times = n)
  cohorts <- rep(c(3, 5, Inf), each = n / 3)
  g_vec   <- cohorts[ids]
  x1_unit <- rnorm(n)
  x2_unit <- rnorm(n)
  x1      <- x1_unit[ids]   # time-invariant
  x2      <- x2_unit[ids]
  y       <- 0.5 * times + 0.3 * x1 +
             as.numeric(times >= g_vec) * (1 + 0.2 * x1_unit[ids]) +
             rnorm(n * n_periods, sd = 0.5)
  data.frame(id = ids, t = times, y = y,
             g = g_vec, x1 = x1, x2 = x2)
}

# ============================================================
# Deprecated covariates argument
# ============================================================

test_that("covariates= argument errors with redirect message", {
  df <- make_panel_cov()
  expect_error(
    edid(df, "y", "id", "t", "g", covariates = c("x1")),
    regexp = "replaced by.*xformla"
  )
})

# ============================================================
# xformla type check
# ============================================================

test_that("non-formula xformla errors immediately", {
  df <- make_panel_cov()
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = "x1"),
    regexp = "one-sided formula"
  )
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = 1L),
    regexp = "one-sided formula"
  )
})

# ============================================================
# Missing variable in xformla
# ============================================================

test_that("xformla with missing column errors informatively", {
  df <- make_panel_cov()
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = ~ nonexistent_var),
    regexp = "nonexistent_var"
  )
})

# ============================================================
# NA in covariates: any period
# ============================================================

test_that("NA in covariate column (any row) is rejected", {
  df <- make_panel_cov()

  # NA in period 1 for unit 1
  df_na1 <- df
  df_na1$x1[1] <- NA_real_
  expect_error(
    edid(df_na1, "y", "id", "t", "g", xformla = ~ x1),
    regexp = "NA values"
  )

  # NA in a later period (period 4) for unit 5
  df_na2 <- df
  df_na2$x1[df_na2$id == 5 & df_na2$t == 4] <- NA_real_
  expect_error(
    edid(df_na2, "y", "id", "t", "g", xformla = ~ x1),
    regexp = "NA values"
  )
})

# ============================================================
# Time-varying covariate rejection
# ============================================================

test_that("time-varying covariate is rejected with informative error", {
  df <- make_panel_cov()

  # Introduce variation for unit 1 in period 2 (keep period 1 value the same)
  period1_val <- df$x1[df$id == 1 & df$t == 1]
  df$x1[df$id == 1 & df$t == 2] <- period1_val + 1.0

  expect_error(
    edid(df, "y", "id", "t", "g", xformla = ~ x1),
    regexp = "time-varying"
  )
})

# ============================================================
# ~1 formula: no error, routes to no-cov path
# ============================================================

test_that("xformla = ~1 runs without error and returns edid_fit", {
  df  <- make_panel_cov()
  fit <- edid(df, "y", "id", "t", "g", xformla = ~1)
  expect_s3_class(fit, "edid_fit")
})

# ============================================================
# Empty formula expansion (only intercept): routes to no-cov path
# ============================================================

test_that("xformla with no variables silently routes to no-covariate path", {
  # ~1 + 0 has no variables (all.vars() = character(0)), so
  # it is treated the same as xformla = ~1: no covariate path invoked.
  df  <- make_panel_cov()
  fit_nocov <- edid(df, "y", "id", "t", "g")
  fit_10    <- edid(df, "y", "id", "t", "g", xformla = ~1 + 0)
  # Both should give the same result (no-cov path)
  expect_equal(fit_nocov$overall$att, fit_10$overall$att, tolerance = 1e-10)
})

# ============================================================
# model.matrix() expansion errors are caught
# ============================================================

test_that("model.matrix() failure in xformla is caught with informative error", {
  df       <- make_panel_cov()
  # Make x1 character (model.matrix would fail or coerce unexpectedly)
  df$x1chr <- as.character(df$x1)
  # character cols fail the numeric-or-factor check
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = ~ x1chr),
    regexp = "numeric or factor"
  )
})

# ============================================================
# Curse-of-dimensionality warning: continuous columns only
# ============================================================

# Small panel with one continuous covariate and 7 binary dummies (the
# Dobkin/HRS-style design where the d >= 5 warning previously counted the
# dummies as continuous), plus 4 extra continuous covariates for the
# genuinely high-dimensional case.
make_panel_curse <- function(n = 60L, n_periods = 3L, seed = 7L) {
  set.seed(seed)
  ids   <- rep(seq_len(n), each = n_periods)
  times <- rep(seq_len(n_periods), times = n)
  g_vec <- rep(c(3, Inf), each = n / 2L)[ids]
  xc    <- replicate(5L, rnorm(n))               # continuous: xc1..xc5
  xd    <- replicate(7L, rbinom(n, 1L, 0.4))     # binary dummies: xd1..xd7
  df <- data.frame(id = ids, t = times,
                   y = 0.4 * times + 0.3 * xc[ids, 1L] +
                     as.numeric(times >= g_vec) + rnorm(n * n_periods, sd = 0.5),
                   g = g_vec)
  for (j in 1:5) df[[paste0("xc", j)]] <- xc[ids, j]
  for (j in 1:7) df[[paste0("xd", j)]] <- xd[ids, j]
  df
}

test_that("d >= 5 curse-of-dimensionality warning counts only continuous covariates", {
  df <- make_panel_curse()

  # 1 continuous + 7 binary dummies (8 columns after expansion): the kernel
  # rate is driven by d = 1 continuous dimension -> NO curse warning.
  w_dummies <- capture_warnings(
    edid(df, "y", "id", "t", "g",
         xformla = ~ xc1 + xd1 + xd2 + xd3 + xd4 + xd5 + xd6 + xd7,
         cband = FALSE, misspec_robust = FALSE)
  )
  expect_false(any(grepl("curse of dimensionality", w_dummies)))

  # 5 genuinely continuous covariates -> the warning fires and reports d = 5.
  w_cont <- capture_warnings(
    edid(df, "y", "id", "t", "g",
         xformla = ~ xc1 + xc2 + xc3 + xc4 + xc5,
         cband = FALSE, misspec_robust = FALSE)
  )
  expect_true(any(grepl("5 continuous covariates", w_cont)))
  expect_true(any(grepl("curse of dimensionality", w_cont)))

  # Mixed: 5 continuous + dummies still warns (dummies don't mask the curse).
  w_mixed <- capture_warnings(
    edid(df, "y", "id", "t", "g",
         xformla = ~ xc1 + xc2 + xc3 + xc4 + xc5 + xd1 + xd2,
         cband = FALSE, misspec_robust = FALSE)
  )
  expect_true(any(grepl("5 continuous covariates", w_mixed)))

  # Constant-weight schemes are exempt regardless of dimension.
  w_avg <- capture_warnings(
    edid(df, "y", "id", "t", "g",
         xformla = ~ xc1 + xc2 + xc3 + xc4 + xc5,
         weight_scheme = "averaged", cband = FALSE, misspec_robust = FALSE)
  )
  expect_false(any(grepl("curse of dimensionality", w_avg)))
})
