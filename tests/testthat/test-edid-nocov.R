library(testthat)

# ============================================================
# 5.1 compute_omega_star_nocov_edid(): return type and dimensions
# ============================================================
test_that("compute_omega_star_nocov_edid() returns H x H symmetric numeric matrix", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 1)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  H <- nrow(pairs)
  omega <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  expect_true(is.matrix(omega))
  expect_equal(dim(omega), c(H, H))
  expect_true(isSymmetric(omega, tol = 1e-10))
})

test_that("compute_omega_star_nocov_edid() is positive semi-definite (eigenvalues >= 0)", {
  df      <- make_panel_1cohort(n_treat = 40, n_never = 40, n_periods = 5, seed = 2)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  eigs    <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))  # PSD up to numerical noise
})

# ============================================================
# 5.2 compute_omega_star_nocov_edid(): PT-Post returns 1x1 matrix
# ============================================================
test_that("compute_omega_star_nocov_edid() returns 1x1 matrix under PT-Post", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 3)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "post")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "post")
  expect_equal(dim(omega), c(1L, 1L))
  # The 1x1 Omega* must be positive (variance is non-negative)
  expect_true(omega[1, 1] >= 0)
})

# ============================================================
# 5.3 compute_efficient_weights_edid(): properties
# ============================================================
test_that("compute_efficient_weights_edid() weights sum to 1", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 4)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  w       <- compute_efficient_weights_edid(omega)
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_equal(length(w), nrow(pairs))
})

test_that("compute_efficient_weights_edid() returns w=1 for a single pair (H=1)", {
  # Construct a 1x1 omega_star
  omega_1x1 <- matrix(0.5)
  w <- compute_efficient_weights_edid(omega_1x1)
  expect_equal(w, 1.0, tolerance = 1e-12)
})

test_that("compute_efficient_weights_edid() returns uniform weights when Omega* is all zeros", {
  H <- 4L
  omega_zero <- matrix(0, H, H)
  w <- compute_efficient_weights_edid(omega_zero)
  expect_equal(w, rep(1/H, H), tolerance = 1e-12)
})

test_that("compute_efficient_weights_edid() uses pseudoinverse fallback for singular Omega*", {
  # Singular 2x2 omega (rank 1)
  v <- c(1, 2)
  omega_sing <- outer(v, v) * 0.1
  # Should not error; should return weights summing to 1
  w <- compute_efficient_weights_edid(omega_sing)
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_equal(length(w), 2L)
})

test_that("compute_efficient_weights_edid() returns numeric vector with no NA or NaN", {
  df      <- make_panel_1cohort(n_treat = 25, n_never = 25, n_periods = 5, seed = 5)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  w       <- compute_efficient_weights_edid(omega)
  expect_true(all(is.finite(w)))
})

# ============================================================
# 5.4 compute_generated_outcomes_nocov_edid(): shape and finiteness
# ============================================================
test_that("compute_generated_outcomes_nocov_edid() returns length-H finite vector", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 6)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  y_hat   <- compute_generated_outcomes_nocov_edid(3L, 4L, pairs, panel, "all")
  expect_equal(length(y_hat), nrow(pairs))
  expect_true(all(is.finite(y_hat)))
})

test_that("compute_generated_outcomes_nocov_edid() PT-Post returns length-1 vector", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 7)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "post")
  y_hat   <- compute_generated_outcomes_nocov_edid(3L, 4L, pairs, panel, "post")
  expect_equal(length(y_hat), 1L)
})

test_that("compute_generated_outcomes_nocov_edid() ATT=2 panel: generated outcome close to 2 for post period", {
  # Panel with known ATT=2 (from make_panel_1cohort default)
  df    <- make_panel_1cohort(n_treat = 200, n_never = 200, n_periods = 5, seed = 99)
  panel <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "post")
  y_hat <- compute_generated_outcomes_nocov_edid(3L, 3L, pairs, panel, "post")
  # With 200 treated, 200 never-treated, and true ATT=2, should be within 0.5 of 2
  expect_equal(y_hat[1], 2, tolerance = 0.5)
})

# ============================================================
# 5.5 compute_eif_nocov_edid(): shape, finiteness, zero-mean
# ============================================================
test_that("compute_eif_nocov_edid() returns length-n finite vector", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 8)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  n       <- panel$n
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  w       <- compute_efficient_weights_edid(omega)
  y_hat   <- compute_generated_outcomes_nocov_edid(3L, 4L, pairs, panel, "all")
  att_gt  <- sum(w * y_hat)
  eif     <- compute_eif_nocov_edid(3L, 4L, pairs, w, panel, att_gt, "all")
  expect_equal(length(eif), n)
  expect_true(all(is.finite(eif)))
})

test_that("compute_eif_nocov_edid() has zero mean (up to numerical precision)", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 9)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  w       <- compute_efficient_weights_edid(omega)
  y_hat   <- compute_generated_outcomes_nocov_edid(3L, 4L, pairs, panel, "all")
  att_gt  <- sum(w * y_hat)
  eif     <- compute_eif_nocov_edid(3L, 4L, pairs, w, panel, att_gt, "all")
  expect_equal(mean(eif), 0, tolerance = 1e-10)
})

test_that("compute_eif_nocov_edid() PT-Post: EIF has zero mean", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 10)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "post")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "post")
  w       <- compute_efficient_weights_edid(omega)
  y_hat   <- compute_generated_outcomes_nocov_edid(3L, 4L, pairs, panel, "post")
  att_gt  <- sum(w * y_hat)
  eif     <- compute_eif_nocov_edid(3L, 4L, pairs, w, panel, att_gt, "post")
  expect_equal(mean(eif), 0, tolerance = 1e-10)
})

test_that("compute_eif_nocov_edid() sum of squared EIF is positive (non-degenerate)", {
  df      <- make_panel_1cohort(n_treat = 30, n_never = 30, n_periods = 5, seed = 11)
  panel   <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  pairs   <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                        panel$period_1, "all")
  omega   <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  w       <- compute_efficient_weights_edid(omega)
  y_hat   <- compute_generated_outcomes_nocov_edid(3L, 4L, pairs, panel, "all")
  att_gt  <- sum(w * y_hat)
  eif     <- compute_eif_nocov_edid(3L, 4L, pairs, w, panel, att_gt, "all")
  expect_true(sum(eif^2) > 0)
})
