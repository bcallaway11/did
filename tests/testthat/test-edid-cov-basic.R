library(testthat)

# ============================================================
# Shared helpers
# ============================================================

make_panel_cov <- function(n = 90, n_periods = 6, seed = 1) {
  set.seed(seed)
  ids     <- rep(1:n, each = n_periods)
  times   <- rep(1:n_periods, times = n)
  g_unit  <- rep(c(3, 5, Inf), each = n / 3)
  g_vec   <- g_unit[ids]
  x1      <- rep(rnorm(n), each = n_periods)
  y       <- 0.5 * times + 0.2 * x1 +
             as.numeric(times >= g_vec) * 1.0 +
             rnorm(n * n_periods, sd = 0.5)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1)
}

# ============================================================
# Regression: xformla=NULL and xformla=~1 must be identical
# ============================================================

test_that("xformla=NULL and xformla=~1 produce bit-for-bit identical results", {
  df     <- make_panel_cov(seed = 10)
  fit0   <- edid(df, "y", "id", "t", "g")
  fit1   <- edid(df, "y", "id", "t", "g", xformla = ~1)

  expect_equal(fit0$overall$overall.att,    fit1$overall$overall.att,    tolerance = 1e-12)
  expect_equal(fit0$overall$overall.se,     fit1$overall$overall.se,     tolerance = 1e-12)
  expect_equal(fit0$att_gt$att,     fit1$att_gt$att,     tolerance = 1e-12)
  expect_equal(fit0$att_gt$se,      fit1$att_gt$se,      tolerance = 1e-12)
})

# Confirm the ~1 path truly skips covariate estimation (covariate_matrix is NULL)
test_that("xformla=~1 routes to no-covariate path (covariate_matrix is NULL)", {
  df     <- make_panel_cov(seed = 11)
  panel0 <- did:::prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~1)
  expect_null(panel0$covariate_matrix)
})

test_that("xformla=NULL routes to no-covariate path (covariate_matrix is NULL)", {
  df     <- make_panel_cov(seed = 12)
  panel0 <- did:::prepare_edid_panel(df, "y", "id", "t", "g", xformla = NULL)
  expect_null(panel0$covariate_matrix)
})

# ============================================================
# Output structure: covariate path returns same class/slots
# ============================================================

test_that("covariate path returns edid_fit with all required slots", {
  df  <- make_panel_cov(seed = 20)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 1L)

  expect_s3_class(fit, "edid_fit")
  expect_true(is.data.frame(fit$att_gt))
  expect_true(all(c("group", "time", "att", "se", "ci_lower", "ci_upper") %in%
                    names(fit$att_gt)))
  expect_true(!is.null(fit$overall))
  expect_true(is.numeric(fit$overall$overall.att))
  expect_true(is.finite(fit$overall$overall.att))
  expect_true(is.numeric(fit$overall$overall.se))
  expect_true(fit$overall$overall.se > 0)
})

# ============================================================
# Covariate path: ATTs are finite for all post-treatment cells
# ============================================================

test_that("covariate path: all post-treatment ATTs are finite", {
  df  <- make_panel_cov(seed = 30)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 1L)
  post_cells <- fit$att_gt[!fit$att_gt$is_pre, ]
  expect_true(nrow(post_cells) > 0L)
  expect_true(all(is.finite(post_cells$att)))
  expect_true(all(is.finite(post_cells$se)))
  expect_true(all(post_cells$se > 0))
})

# ============================================================
# No-covariate regression: no change after adding irrelevant xformla
# NOTE: we do NOT expect them to be equal, only that both are valid
# ============================================================

test_that("covariate path runs without error on 2D covariate formula", {
  df  <- make_panel_cov(seed = 40)
  df$x2 <- rep(rnorm(nrow(df) / 6), each = 6)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2, seed = 1L)
  expect_s3_class(fit, "edid_fit")
  expect_true(is.finite(fit$overall$overall.att))
})

# ============================================================
# Transformed formula: I(x1^2) supported via model.matrix()
# ============================================================

test_that("xformla with I(x1^2) runs and differs from ~x1 on nonlinear DGP", {
  set.seed(55)
  n <- 90; T <- 6
  ids   <- rep(1:n, each = T)
  times <- rep(1:T, times = n)
  x1u   <- rnorm(n)
  g_u   <- rep(c(3, 5, Inf), each = n / 3)
  x1    <- rep(x1u, each = T)
  g     <- g_u[ids]
  # nonlinear effect of x1
  y <- 0.5 * times + x1^2 + as.numeric(times >= g) + rnorm(n * T, sd = 0.5)
  df <- data.frame(id = ids, t = times, y = y, g = g, x1 = x1)

  # This small-n nonlinear DGP deliberately stresses overlap, so the quadratic
  # fit can trip the (legitimate) extreme-propensity-ratio diagnostic. Whether it
  # crosses the threshold is seed/platform-dependent, so suppress rather than
  # assert it -- the test's purpose is that both fits run and differ.
  fit_lin  <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1,            seed = 1L))
  fit_quad <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + I(x1^2), seed = 1L))

  # Both should run; they should not be identical (different model matrix)
  expect_s3_class(fit_lin,  "edid_fit")
  expect_s3_class(fit_quad, "edid_fit")
  expect_false(isTRUE(all.equal(fit_lin$overall$overall.att, fit_quad$overall$overall.att,
                                 tolerance = 1e-6)))
})

# ============================================================
# Factor covariate: supported via model.matrix()
# ============================================================

test_that("factor covariate is accepted and produces finite results", {
  df  <- make_panel_cov(seed = 60)
  df$fac <- as.factor(rep(c("A", "B", "C"), length.out = nrow(df)))
  # First: make factor time-invariant within unit
  fac_unit <- df$fac[df$t == 1]
  df$fac   <- fac_unit[df$id]

  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1 + fac, seed = 1L)
  expect_s3_class(fit, "edid_fit")
  expect_true(is.finite(fit$overall$overall.att))
})
