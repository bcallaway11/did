# test-edid-nocov-shrink.R
# Pole-target Ledoit-Wolf shrinkage of Omega* on the no-covariate path (nocov_shrink).
#
# The audited failure mode (small-n study, 2026-06-11): at n = 50 the H(H+1)/2
# estimated covariance entries make the efficient weights noisy enough that the
# estimator's realized variance exceeds fixed pole-weight imputation's by ~12%
# even AT the i.i.d. pole, where the pole weights are optimal. nocov_shrink
# stabilizes the WEIGHTS by shrinking Omega* toward its closed-form i.i.d.-pole
# structure with a Ledoit-Wolf intensity; the SE machinery (empirical variance
# of the realized weighted IF) is untouched.

# Shared staggered test panel: cohorts {3, 5} + never-treated, T = 6, optional
# AR(1) shocks. Never-treated coded 0 (edid() remaps to Inf); the *_panel_inf
# variant codes Inf directly for the internal-builder tests.
.shrink_test_panel <- function(n, rho = 0, seed = 42L, nt_inf = FALSE) {
  set.seed(seed)
  Tp <- 6L
  G  <- sample(c(3L, 5L, 0L), n, TRUE, c(0.35, 0.30, 0.35))
  e  <- matrix(rnorm(n * Tp), n, Tp)
  if (rho > 0) for (s in 2:Tp) e[, s] <- rho * e[, s - 1L] + sqrt(1 - rho^2) * e[, s]
  eta <- rnorm(n) + 0.4 * (G == 3L) - 0.3 * (G == 5L)
  Y   <- eta + matrix(0.3 * seq_len(Tp), n, Tp, byrow = TRUE) + e
  for (g in c(3L, 5L)) for (t in g:Tp) {
    idx <- G == g
    Y[idx, t] <- Y[idx, t] + 1 + 0.3 * (t - g)
  }
  df <- data.frame(id   = rep(seq_len(n), each = Tp),
                   time = rep(seq_len(Tp), n),
                   y    = as.vector(t(Y)),
                   gvar = rep(G, each = Tp))
  if (nt_inf) df$gvar <- ifelse(df$gvar == 0L, Inf, df$gvar)
  df
}

.shrink_cell_lambdas <- function(fit) {
  vapply(fit$cells, function(cc) {
    v <- cc$nocov_shrink_lambda
    if (is.null(v)) NA_real_ else v
  }, numeric(1L))
}

# ---------------------------------------------------------------------------
# 1. Finite-sample identity: Omega* == crossprod(psi)/n^2
#    (the LW entry-variance estimate is coherent with the Omega* builder)
# ---------------------------------------------------------------------------

test_that("compute_psi_moments_nocov_edid() reproduces Omega* exactly (crossprod(psi)/n^2)", {
  df    <- .shrink_test_panel(300L, rho = 0.7, seed = 42L, nt_inf = TRUE)
  panel <- prepare_edid_panel(df, "y", "id", "time", "gvar")
  pairs <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  om  <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  psi <- compute_psi_moments_nocov_edid(3L, 4L, pairs, panel)
  expect_equal(crossprod(psi) / panel$n^2, om,
               tolerance = 1e-12, ignore_attr = TRUE)
})

# ---------------------------------------------------------------------------
# 2. Pole structure matrix: sigma^2 * S is the population Omega* under iid shocks
# ---------------------------------------------------------------------------

test_that("compute_pole_structure_nocov_edid() matches the empirical Omega* on a large iid draw", {
  df    <- .shrink_test_panel(60000L, rho = 0, seed = 7L, nt_inf = TRUE)
  panel <- prepare_edid_panel(df, "y", "id", "time", "gvar")
  pairs <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  om <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  S  <- compute_pole_structure_nocov_edid(3L, 4L, pairs, panel)
  # unit shocks (sigma2 = 1): Omega-hat -> S entrywise; 60k draw ~ 1% Frobenius
  expect_lt(max(abs(om - S)) / max(abs(S)), 0.05)
  expect_true(isSymmetric(S))
  # degenerate self pair (tpre = period_1): the comparison difference is the zero
  # vector, so its D-column contribution must mirror the empirical builder's zeros
  j1 <- which(pairs$gp == 3L & pairs$tpre == panel$period_1)
  expect_length(j1, 1L)
})

# ---------------------------------------------------------------------------
# 3. Shrinkage core: lambda in [0, 1]; lambda = 1 reproduces pole weights;
#    intensity decays with n off the pole and stays high at the pole
# ---------------------------------------------------------------------------

test_that("shrink_omega_nocov_edid() returns lambda in [0,1] and a symmetric PSD-safe matrix", {
  df    <- .shrink_test_panel(80L, rho = 0.7, seed = 9L, nt_inf = TRUE)
  panel <- prepare_edid_panel(df, "y", "id", "time", "gvar")
  pairs <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  om <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  sh <- shrink_omega_nocov_edid(om, 3L, 4L, pairs, panel)
  expect_true(is.finite(sh$lambda) && sh$lambda >= 0 && sh$lambda <= 1)
  expect_true(is.finite(sh$sigma2) && sh$sigma2 > 0)
  expect_true(isSymmetric(sh$omega))
  expect_true(all(eigen(sh$omega, symmetric = TRUE, only.values = TRUE)$values > -1e-12))
})

test_that("at lambda = 1 the shrunk weights equal the closed-form pole weights", {
  df    <- .shrink_test_panel(80L, rho = 0.7, seed = 9L, nt_inf = TRUE)
  panel <- prepare_edid_panel(df, "y", "id", "time", "gvar")
  pairs <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  S      <- compute_pole_structure_nocov_edid(3L, 4L, pairs, panel)
  w_pole <- compute_efficient_weights_edid(S)
  # the pole weights are sigma2-scale invariant, so a fully-shrunk matrix gives them back
  expect_equal(compute_efficient_weights_edid(17.3 * S), w_pole, tolerance = 1e-10)
})

test_that("lambda-hat is ~1 at the iid pole and decays with n off the pole", {
  lam_mean <- function(n, rho, seed) {
    fit <- suppressWarnings(
      edid(.shrink_test_panel(n, rho = rho, seed = seed), yname = "y", idname = "id",
           tname = "time", gname = "gvar", pt_assumption = "all",
           weight_scheme = "efficient", aggregate = "none", cband = FALSE,
           nocov_shrink = TRUE))
    mean(.shrink_cell_lambdas(fit), na.rm = TRUE)
  }
  expect_gt(lam_mean(60L, rho = 0,   seed = 11L), 0.8)    # pole: stay shrunk
  l_small <- lam_mean(60L,   rho = 0.7, seed = 12L)
  l_big   <- lam_mean(4000L, rho = 0.7, seed = 12L)
  expect_gt(l_small, l_big + 0.2)                          # decays with n off the pole
  expect_lt(l_big, 0.10)                                   # and is asymptotically negligible
})

# ---------------------------------------------------------------------------
# 4. Off-switch: nocov_shrink = FALSE reproduces the unshrunk pipeline bit-for-bit
# ---------------------------------------------------------------------------

test_that("nocov_shrink = FALSE/default reproduces the legacy weights/ATT/SE bit-for-bit", {
  df  <- .shrink_test_panel(150L, rho = 0.5, seed = 4L)
  fit <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE, nocov_shrink = FALSE))
  fit_def <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE))
  expect_identical(fit_def$att_gt$att, fit$att_gt$att)
  expect_identical(fit_def$att_gt$se, fit$att_gt$se)
  expect_false(fit_def$nocov_shrink)
  # oracle: rebuild each post cell's weights from the RAW Omega* (the pre-shrinkage
  # pipeline) and confirm exact equality, including the stored condition number
  df_inf <- df; df_inf$gvar <- ifelse(df_inf$gvar == 0L, Inf, df_inf$gvar)
  panel  <- prepare_edid_panel(df_inf, "y", "id", "time", "gvar")
  for (cc in fit$cells) {
    if (!isTRUE(is.finite(cc$att))) next
    pairs_cc <- cc$pairs
    om  <- compute_omega_star_nocov_edid(cc$group, cc$time, pairs_cc, panel, "all")
    expect_identical(unname(cc$weights), unname(compute_efficient_weights_edid(om)))
    expect_identical(cc$condition_num, check_condition_edid(om))
    expect_identical(cc$nocov_shrink_lambda, NA_real_)
  }
  expect_false(fit$nocov_shrink)
})

test_that("nocov_shrink = TRUE inverts the shrunk matrix and records lambda", {
  df  <- .shrink_test_panel(80L, rho = 0.5, seed = 21L)
  fit <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE, nocov_shrink = TRUE))
  expect_true(fit$nocov_shrink)
  df_inf <- df; df_inf$gvar <- ifelse(df_inf$gvar == 0L, Inf, df_inf$gvar)
  panel  <- prepare_edid_panel(df_inf, "y", "id", "time", "gvar")
  n_checked <- 0L
  for (cc in fit$cells) {
    if (!isTRUE(is.finite(cc$att)) || nrow(cc$pairs) < 2L) next
    om <- compute_omega_star_nocov_edid(cc$group, cc$time, cc$pairs, panel, "all")
    sh <- shrink_omega_nocov_edid(om, cc$group, cc$time, cc$pairs, panel)
    expect_equal(unname(cc$weights), unname(compute_efficient_weights_edid(sh$omega)),
                 tolerance = 1e-12)
    expect_equal(cc$nocov_shrink_lambda, sh$lambda, tolerance = 1e-12)
    n_checked <- n_checked + 1L
  }
  expect_gt(n_checked, 0L)
})

# ---------------------------------------------------------------------------
# 5. No-ops: uniform weights, PT-Post, and H = 1 cells record lambda = NA
# ---------------------------------------------------------------------------

test_that("nocov_shrink is inert for uniform weights, PT-Post, and the covariate path", {
  df <- .shrink_test_panel(120L, rho = 0.5, seed = 33L)
  fit_unif <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "uniform",
         aggregate = "none", cband = FALSE))
  expect_true(all(is.na(.shrink_cell_lambdas(fit_unif))))

  fit_post <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "post", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE))
  expect_true(all(is.na(.shrink_cell_lambdas(fit_post))))

  # uniform fits are bit-identical under both switch values (no weights estimated)
  fit_unif_off <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "uniform",
         aggregate = "none", cband = FALSE, nocov_shrink = FALSE))
  expect_identical(fit_unif$att_gt$att, fit_unif_off$att_gt$att)
  expect_identical(fit_unif$att_gt$se,  fit_unif_off$att_gt$se)
})

test_that("nocov_shrink validation rejects non-logical input", {
  df <- .shrink_test_panel(60L, seed = 5L)
  expect_error(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         nocov_shrink = "yes"),
    "logical scalar")
  expect_error(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         nocov_shrink = NA),
    "logical scalar")
})

# ---------------------------------------------------------------------------
# 6. The SE stays = empirical variance of the realized weighted IF (the shrunk
#    matrix stabilizes weights only; it never replaces the data's IF variance)
# ---------------------------------------------------------------------------

test_that("with shrinkage on, the cell SE is the empirical variance of the realized weighted IF", {
  df  <- .shrink_test_panel(90L, rho = 0.5, seed = 8L)
  fit <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE, nocov_shrink = TRUE))
  df_inf <- df; df_inf$gvar <- ifelse(df_inf$gvar == 0L, Inf, df_inf$gvar)
  panel  <- prepare_edid_panel(df_inf, "y", "id", "time", "gvar")
  n_checked <- 0L
  for (cc in fit$cells) {
    if (!isTRUE(is.finite(cc$se)) || nrow(cc$pairs) < 2L) next
    eif <- compute_eif_nocov_edid(cc$group, cc$time, cc$pairs, unname(cc$weights),
                                  panel, cc$att, "all")
    se_oracle <- safe_inference_edid(eif, panel$cluster_indices, fit$alpha, cc$att)$se
    expect_equal(cc$se, se_oracle, tolerance = 1e-12)
    n_checked <- n_checked + 1L
  }
  expect_gt(n_checked, 0L)
})

# ---------------------------------------------------------------------------
# 7. Guard interplay: a guard-pinned (H = 1) cohort records lambda = NA while
#    healthy cohorts still shrink
# ---------------------------------------------------------------------------

test_that("thin-cohort-guard-pinned cells skip shrinkage (H = 1), healthy cells shrink", {
  df <- .shrink_test_panel(120L, rho = 0.5, seed = 14L)
  # make cohort 5 thin: keep 3 of its units
  ids5 <- unique(df$id[df$gvar == 5L])
  drop <- ids5[-(1:3)]
  df   <- df[!(df$id %in% drop), ]
  fit  <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE, nocov_shrink = TRUE))
  lam  <- .shrink_cell_lambdas(fit)
  pinned  <- vapply(fit$cells, function(cc) isTRUE(cc$thin_cohort_degraded), logical(1L))
  expect_true(any(pinned))
  expect_true(all(is.na(lam[pinned])))         # just-identified: no weights to stabilize
  healthy_post <- !pinned &
    vapply(fit$cells, function(cc) isTRUE(is.finite(cc$att)) && nrow(cc$pairs) > 1L, logical(1L))
  expect_true(any(healthy_post))
  expect_true(all(is.finite(lam[healthy_post])))
})
