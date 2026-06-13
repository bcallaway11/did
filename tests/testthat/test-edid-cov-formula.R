library(testthat)

# ============================================================
# Shared helpers
# ============================================================

make_panel_2cov <- function(n = 90, n_periods = 6, seed = 1) {
  set.seed(seed)
  ids     <- rep(1:n, each = n_periods)
  times   <- rep(1:n_periods, times = n)
  g_unit  <- rep(c(3, 5, Inf), each = n / 3)
  g_vec   <- g_unit[ids]
  x1u     <- rnorm(n)
  x2u     <- rnorm(n)
  x1      <- rep(x1u, each = n_periods)
  x2      <- rep(x2u, each = n_periods)
  y       <- 0.5 * times + x1 + x2^2 +
             as.numeric(times >= g_vec) * 1.0 +
             rnorm(n * n_periods, sd = 0.5)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1, x2 = x2)
}

# ============================================================
# Time-varying covariate: must fail regardless of which period varies
# ============================================================

test_that("time-varying covariate (post-period change) fails before estimation", {
  df <- make_panel_2cov(seed = 10)
  # Only post-period change for one unit — still must be rejected
  df$x1[df$id == 2 & df$t == 4] <- df$x1[df$id == 2 & df$t == 4] + 1.5
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = ~ x1),
    regexp = "time-varying"
  )
})

test_that("time-varying covariate (early-period change) fails before estimation", {
  df <- make_panel_2cov(seed = 11)
  # Change in period 2 for one unit
  df$x1[df$id == 3 & df$t == 2] <- df$x1[df$id == 3 & df$t == 1] + 0.5
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = ~ x1),
    regexp = "time-varying"
  )
})

# ============================================================
# NA covariate: both period-1 and late-period NAs rejected
# ============================================================

test_that("NA in covariate period 1 for one unit is rejected", {
  df <- make_panel_2cov(seed = 20)
  df$x1[df$id == 1 & df$t == 1] <- NA_real_
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = ~ x1),
    regexp = "NA"
  )
})

test_that("NA in covariate late period for one unit is rejected", {
  df <- make_panel_2cov(seed = 21)
  df$x2[df$id == 5 & df$t == 5] <- NA_real_
  expect_error(
    edid(df, "y", "id", "t", "g", xformla = ~ x2),
    regexp = "NA"
  )
})

# ============================================================
# formula semantics: model.matrix() expansion
# ============================================================

test_that("covariate_matrix uses model.matrix() and handles I() correctly", {
  df <- make_panel_2cov(seed = 30)
  p0 <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1)
  p1 <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + I(x1^2))

  n_units <- length(unique(df$id))
  expect_equal(nrow(p0$covariate_matrix), n_units)
  expect_equal(ncol(p0$covariate_matrix), 1L)   # just x1

  expect_equal(nrow(p1$covariate_matrix), n_units)
  expect_equal(ncol(p1$covariate_matrix), 2L)   # x1, I(x1^2)

  # I(x1^2) column equals x1^2
  x1_unit  <- df$x1[df$t == df$t[1]][match(sort(unique(df$id)), df$id[df$t == df$t[1]])]
  expect_equal(p1$covariate_matrix[, 2L], x1_unit^2, tolerance = 1e-10)
})

test_that("interaction in xformla produces expected number of columns", {
  df <- make_panel_2cov(seed = 31)
  p  <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 * x2)
  n_units <- length(unique(df$id))
  expect_equal(nrow(p$covariate_matrix), n_units)
  # x1, x2, x1:x2 = 3 columns
  expect_equal(ncol(p$covariate_matrix), 3L)
})

# ============================================================
# factor covariate: supported via model.matrix
# ============================================================

test_that("factor covariate is accepted and covariate_matrix has dummy columns", {
  df       <- make_panel_2cov(seed = 40)
  fac_unit <- factor(rep(c("A", "B", "C"), each = nrow(df) / 6 / 3 + 1))[1:(nrow(df) / 6)]
  df$fac   <- rep(fac_unit, each = 6)

  p <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + fac)
  # model.matrix with factor "fac" (3 levels) and x1 produces 3 columns
  # (x1, facB, facC) with treatment contrast; intercept is removed
  n_units <- length(unique(df$id))
  expect_equal(nrow(p$covariate_matrix), n_units)
  expect_true(ncol(p$covariate_matrix) >= 2L)  # at least x1 + one dummy
})

# ============================================================
# Formula with x1:x2 vs x1 * x2 parity (interaction)
# ============================================================

test_that("xformla=~x1+x2+x1:x2 and xformla=~x1*x2 produce identical covariate matrices", {
  df <- make_panel_2cov(seed = 50)
  p1 <- prepare_edid_panel(df, "y", "id", "t", "g",
                                  xformla = ~ x1 + x2 + x1:x2)
  p2 <- prepare_edid_panel(df, "y", "id", "t", "g",
                                  xformla = ~ x1 * x2)
  expect_equal(p1$covariate_matrix, p2$covariate_matrix, tolerance = 1e-10)
})

# ============================================================
# Reproducibility via seed
# ============================================================

test_that("two calls with same seed produce identical results on covariate path", {
  df   <- make_panel_2cov(seed = 60)
  fit1 <- edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 42L)
  fit2 <- edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 42L)

  expect_equal(fit1$overall$att, fit2$overall$att, tolerance = 1e-12)
  expect_equal(fit1$overall$se,  fit2$overall$se,  tolerance = 1e-12)
  expect_equal(fit1$att_gt$att,  fit2$att_gt$att,  tolerance = 1e-12)
})

test_that("different seeds produce different fold assignments (probabilistically)", {
  # With n=90 units and K=5 folds, different seeds should almost always give
  # different fold vectors.
  set.seed(999)
  folds1 <- build_crossfit_folds_edid(90L, 5L, seed = 1L)
  folds2 <- build_crossfit_folds_edid(90L, 5L, seed = 2L)
  expect_false(identical(folds1, folds2))
})

test_that("same seed always produces same fold assignments", {
  folds1 <- build_crossfit_folds_edid(90L, 5L, seed = 77L)
  folds2 <- build_crossfit_folds_edid(90L, 5L, seed = 77L)
  expect_identical(folds1, folds2)
})
