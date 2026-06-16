library(testthat)

# ============================================================
# Phase 2 gate for the weighted-covariate rollout: the three
# Omega*(X) conditional-covariance builders gained an obs-weight
# path (weighted Nadaraya-Watson / WLS conditional moment + weighted
# pooling over the marginal X).
#
# Omega*(X) is built from the conditional covariance of OUTCOMES given
# X (kernel / sieve smoothing) scaled by 1/p prefactors; it does not use
# prop_ratios / cond_means in its value, so these call the builders
# directly with NULL nuisances (and inv_propensities = NULL -> the
# pi-based prefactor fallback, which is itself weight-consistent).
#
# Two invariants:
#   (A) weighted branch with uw == 1 reproduces the unweighted build
#       BYTE-IDENTICALLY (the constant weight column normalizes to all
#       ones, so the weighted path must collapse to the legacy path).
#       Covers all 3 smoothers x {averaged, per-unit} builds.
#   (B) a genuinely dispersed mean-1 weight MOVES Omega (weights are not
#       inert) -- the smoother + pooling actually consume the weights.
#
# The unweighted path (uw == NULL) byte-identity across all smoothers is
# guarded separately by quality_reports/wcov/invariant_harness.R.
# ============================================================

make_omega_panel <- function(n = 300, seed = 21) {
  set.seed(seed)
  Tn  <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P   <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  wu   <- exp(0.8 * rnorm(n)); wu <- wu / mean(wu)         # dispersed, mean-1
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    # Never-treated MUST stay Inf: these tests call prepare_edid_panel / the Omega
    # builders DIRECTLY, which expect the internal never-treated convention (Inf).
    # (edid() recodes a user-facing g = 0 to Inf at its boundary; a direct panel build
    # does not, so coding never-treated as 0 here would make 0 a phantom treatment group
    # with pi_inf = 0 -> 1/pi_inf = Inf -> Inf * 0 = NaN in the Eq.(3.12) prefactors.)
    data.frame(id = 1:n, t = tt, g = gc,
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

# Build a (g, t, pairs) cell whose Omega is non-degenerate (H >= 1).
cell_pairs <- function(panel, g, t) {
  enumerate_valid_pairs_edid(
    target_g = g, treatment_groups = panel$treatment_groups,
    time_periods = as.numeric(names(panel$period_to_col)),
    period_1 = panel$period_1, pt_assumption = "all",
    anticipation = if (is.null(panel$anticipation)) 0L else panel$anticipation)
}

# Run one builder with uw = NULL vs uw = rep(1, n); arrays/matrices compared sans attributes.
expect_w1_identical <- function(builder, panel, g, t, pairs, pointwise, ...) {
  o0 <- builder(panel, g, t, pairs, NULL, NULL, NULL, return_pointwise = pointwise, ...)
  p1 <- panel; p1$unit_weights <- rep(1, panel$n)
  o1 <- builder(p1, g, t, pairs, NULL, NULL, NULL, return_pointwise = pointwise, ...)
  expect_equal(as.numeric(o1), as.numeric(o0), tolerance = 1e-12,
               info = sprintf("pointwise=%s", pointwise))
  invisible(list(o0 = o0, o1 = o1))
}

df    <- make_omega_panel()
panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1)
g_use <- panel$treatment_groups[1L]                 # = 2
t_use <- max(as.numeric(names(panel$period_to_col)))# = 4 (post-treatment)
prs   <- cell_pairs(panel, g_use, t_use)

test_that("the test cell has a non-empty moment set", {
  expect_gt(nrow(prs), 0L)
})

test_that("kernel-fast Omega*(X): weighted branch (w==1) == unweighted (averaged + per-unit)", {
  expect_w1_identical(compute_omega_star_kernel_fast_edid, panel, g_use, t_use, prs, FALSE)
  expect_w1_identical(compute_omega_star_kernel_fast_edid, panel, g_use, t_use, prs, TRUE)
})

test_that("kernel-orig Omega*(X): weighted branch (w==1) == unweighted (averaged + per-unit)", {
  expect_w1_identical(compute_omega_star_cov_edid, panel, g_use, t_use, prs, FALSE)
  expect_w1_identical(compute_omega_star_cov_edid, panel, g_use, t_use, prs, TRUE)
})

test_that("sieve Omega*(X): weighted branch (w==1) == unweighted (averaged + per-unit)", {
  expect_w1_identical(compute_omega_star_sieve_edid, panel, g_use, t_use, prs, FALSE, bs_df = 4L)
  expect_w1_identical(compute_omega_star_sieve_edid, panel, g_use, t_use, prs, TRUE,  bs_df = 4L)
})

test_that("a dispersed mean-1 weight MOVES Omega*(X) (weights are not inert)", {
  panel_w <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1, weightsname = "w")
  for (builder in list(compute_omega_star_kernel_fast_edid,
                       compute_omega_star_cov_edid,
                       compute_omega_star_sieve_edid)) {
    o_unw <- builder(panel,   g_use, t_use, prs, NULL, NULL, NULL, return_pointwise = FALSE)
    o_w   <- builder(panel_w, g_use, t_use, prs, NULL, NULL, NULL, return_pointwise = FALSE)
    expect_false(isTRUE(all.equal(as.numeric(o_w), as.numeric(o_unw), tolerance = 1e-6)))
  }
})

test_that("constant weight column normalizes to all-ones -> Omega == unweighted (end-to-end panel build)", {
  df_c <- df; df_c$w <- 3                                 # constant weight column
  panel_c <- prepare_edid_panel(df_c, "y", "id", "t", "g", xformla = ~ x1, weightsname = "w")
  expect_equal(panel_c$unit_weights, rep(1, panel_c$n), tolerance = 1e-12)
  o_unw <- compute_omega_star_kernel_fast_edid(panel,   g_use, t_use, prs, NULL, NULL, NULL, return_pointwise = FALSE)
  o_c   <- compute_omega_star_kernel_fast_edid(panel_c, g_use, t_use, prs, NULL, NULL, NULL, return_pointwise = FALSE)
  expect_equal(as.numeric(o_c), as.numeric(o_unw), tolerance = 1e-12)
})
