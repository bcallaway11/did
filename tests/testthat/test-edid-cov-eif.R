library(testthat)

# ============================================================
# Tests for the covariate-path EIF and generated outcomes.
# These tests verify the FORMULA, not just zero-mean property.
# ============================================================

make_simple_panel <- function(n = 120, seed = 1) {
  set.seed(seed)
  T      <- 5
  ids    <- rep(1:n, each = T)
  times  <- rep(1:T, times = n)
  g_unit <- rep(c(3, Inf), each = n / 2)
  g_vec  <- g_unit[ids]
  x1u    <- rnorm(n)
  x1     <- rep(x1u, each = T)
  y      <- 0.5 * times + x1 +
            as.numeric(times >= g_vec) * 1.0 +
            rnorm(n * T, sd = 0.4)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1)
}

# ============================================================
# EIF has (approximately) zero mean
# ============================================================

test_that("covariate-path EIF has near-zero mean for each cell", {
  df  <- make_simple_panel(n = 120, seed = 10)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1,
               store_eif = TRUE, seed = 1L)
  # eif_matrix is n x n_cells
  eif_mat <- fit$eif
  if (!is.null(eif_mat)) {
    col_means <- colMeans(eif_mat, na.rm = TRUE)
    expect_true(all(abs(col_means) < 1e-10),
                info = paste("EIF column means:", paste(round(col_means, 8), collapse = ", ")))
  }
})

# ============================================================
# SE = sqrt(sum(eif^2)/n^2): verify via manual calculation
# ============================================================

test_that("reported SE matches manual EIF plug-in formula for valid-inference cells", {
  # Note: cells where inference_valid = FALSE (SE below eps threshold, e.g. exact
  # pre-treatment zeros) will have reported SE = NA even though sqrt(sum(eif^2)/n^2)
  # gives a finite (possibly 0) value. We compare only cells with finite reported SE.
  df  <- make_simple_panel(n = 120, seed = 20)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1,
               store_eif = TRUE, aggregate = "none", seed = 1L)
  eif_mat <- fit$eif
  if (!is.null(eif_mat)) {
    n            <- fit$n
    manual_ses   <- sqrt(colSums(eif_mat^2) / n^2)
    reported_ses <- fit$att_gt$se
    # Compare only where both are finite and reported SE > 0
    valid        <- is.finite(reported_ses) & reported_ses > 0
    if (sum(valid) > 0L) {
      expect_equal(manual_ses[valid], reported_ses[valid], tolerance = 1e-8,
                   info = "SE from EIF^2/n^2 must match reported SE for valid cells")
    }
  }
})

# ============================================================
# EIF formula correctness: direct vs compute_eif_cov_edid()
# ============================================================
# For a SINGLE pair (one column in gen_out_mat), the EIF is:
#   EIF_i = w * phi_i - att_gt
# where w = 1 (only one pair, weights sum to 1).
# This is the reference for checking the formula exactly.

test_that("compute_eif_cov_edid formula: weighted_phi minus att_gt, then centered", {
  # Build a minimal scenario with exactly 1 pair
  df <- make_simple_panel(n = 120, seed = 30)
  devtools::load_all(quiet = TRUE)

  panel_obj <- did:::prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1)
  g <- 3; t <- 3
  pairs <- did:::enumerate_valid_pairs_edid(g, panel_obj$treatment_groups,
                                             panel_obj$time_periods,
                                             panel_obj$period_1, "post",
                                             panel_obj$anticipation)
  fold_id <- did:::build_crossfit_folds_edid(panel_obj$n, 5L, seed = 1L)

  # Build nuisances
  pairs_for_nuis <- pairs
  pairs_for_nuis$gp[is.finite(pairs_for_nuis$gp) & pairs_for_nuis$gp == g] <- Inf
  prop_r  <- did:::estimate_all_propensity_ratios(panel_obj, g, pairs_for_nuis,
                                                   4L, 5L, fold_id)
  cond_m  <- did:::estimate_all_conditional_means(panel_obj, pairs_for_nuis, t,
                                                   4L, 5L, fold_id)

  gen_out <- did:::compute_generated_outcomes_cov_edid(panel_obj, g, t, pairs,
                                                        prop_r, cond_m, "post")
  H       <- ncol(gen_out)
  omega   <- did:::compute_omega_star_cov_edid(panel_obj, g, t, pairs, prop_r, cond_m)
  weights <- did:::compute_efficient_weights_edid(omega)
  att_gt  <- sum(weights * colMeans(gen_out, na.rm = TRUE))

  # Compute EIF via the function
  eif_fn <- did:::compute_eif_cov_edid(panel_obj, gen_out, weights, att_gt, g)

  # Correct EIF is the ratio-estimator influence function (fixes the prior over-coverage):
  #   EIF_i = w' phi_i - (G_{g,i} / pi_g) * att_gt    (mean-zero in sample; no de-meaning).
  # The old constant centering (w' phi - att_gt) omitted the influence of the estimated
  # treated-cohort share pi_hat_g and inflated the variance by att^2 (1/pi_g - 1).
  Gg      <- as.numeric(panel_obj$cohort_masks[[as.character(g)]])
  pi_g    <- panel_obj$cohort_fractions[[as.character(g)]]
  eif_ref <- drop(gen_out %*% weights) - (Gg / pi_g) * att_gt

  expect_equal(eif_fn, eif_ref, tolerance = 1e-10,
               info = "compute_eif_cov_edid must implement w^T phi - (G_g/pi_g) att_gt (ratio-estimator IF)")
})

# ============================================================
# Wrong-sign EIF diagnostic: must be caught by SE comparison
# ============================================================
# This is a "canary" test: if the EIF sign/scale were wrong, the implied
# SE would differ from the correct value.  We verify that a deliberately
# wrong EIF produces a detectably different SE.

test_that("wrong-sign EIF produces materially different SE", {
  df  <- make_simple_panel(n = 120, seed = 40)
  panel_obj <- did:::prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1)
  g <- 3; t <- 3
  pairs <- did:::enumerate_valid_pairs_edid(g, panel_obj$treatment_groups,
                                             panel_obj$time_periods,
                                             panel_obj$period_1, "post",
                                             panel_obj$anticipation)
  fold_id <- did:::build_crossfit_folds_edid(panel_obj$n, 5L, seed = 1L)
  pairs_nuis <- pairs
  pairs_nuis$gp[is.finite(pairs_nuis$gp) & pairs_nuis$gp == g] <- Inf
  prop_r <- did:::estimate_all_propensity_ratios(panel_obj, g, pairs_nuis,
                                                  4L, 5L, fold_id)
  cond_m <- did:::estimate_all_conditional_means(panel_obj, pairs_nuis, t,
                                                  4L, 5L, fold_id)
  gen_out <- did:::compute_generated_outcomes_cov_edid(panel_obj, g, t, pairs,
                                                        prop_r, cond_m, "post")
  omega   <- did:::compute_omega_star_cov_edid(panel_obj, g, t, pairs, prop_r, cond_m)
  weights <- did:::compute_efficient_weights_edid(omega)
  att_gt  <- sum(weights * colMeans(gen_out, na.rm = TRUE))
  n       <- panel_obj$n

  # Correct EIF
  eif_correct <- did:::compute_eif_cov_edid(panel_obj, gen_out, weights, att_gt, g)
  se_correct  <- sqrt(sum(eif_correct^2) / n^2)

  # Wrong EIF: adds (Ig/pi_g) * att_gt extra term (old bug)
  pi_g  <- panel_obj$cohort_fractions[[as.character(g)]]
  Ig    <- as.numeric(panel_obj$cohort_masks[[as.character(g)]])
  eif_wrong <- drop(gen_out %*% weights) + (Ig / pi_g) * att_gt
  eif_wrong <- eif_wrong - mean(eif_wrong)
  se_wrong  <- sqrt(sum(eif_wrong^2) / n^2)

  # The wrong SE should differ from the correct one when ATT != 0
  if (abs(att_gt) > 0.05) {
    expect_false(isTRUE(all.equal(se_correct, se_wrong, tolerance = 1e-4)),
                 info = "Wrong EIF (old bug) must produce detectably different SE")
  }
})

# ============================================================
# Generated outcomes: self-comparison pair has correct structure
# ============================================================

test_that("generated outcome for self-comparison pair: zero for non-g/non-inf units", {
  # For a self-comparison pair (gp=g remapped to Inf):
  # phi_i != 0 only for G_g or G_inf units.
  # For units in other cohorts, both Ig=0 and I_inf=0, so phi_i = 0.
  df <- make_simple_panel(n = 120, seed = 50)
  # Add a third cohort so we have "other" units
  set.seed(50)
  n_units  <- 120; T <- 5
  ids      <- rep(1:n_units, each = T)
  times    <- rep(1:T, times = n_units)
  g_unit   <- rep(c(3, 5, Inf), each = n_units / 3)
  g_vec    <- g_unit[ids]
  x1u      <- rnorm(n_units)
  x1       <- rep(x1u, each = T)
  y        <- 0.5 * times + x1 + as.numeric(times >= g_vec) + rnorm(n_units * T, sd = 0.4)
  df2      <- data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1)

  panel_obj <- did:::prepare_edid_panel(df2, "y", "id", "t", "g", xformla = ~ x1)
  g    <- 3; t <- 3
  # Use PT-Post pair: single self-comparison pair (gp=g, tpre=g-1)
  pairs <- did:::enumerate_valid_pairs_edid(g, panel_obj$treatment_groups,
                                             panel_obj$time_periods,
                                             panel_obj$period_1, "post",
                                             panel_obj$anticipation)
  # Self-comparison pair: gp == g; the code remaps to Inf
  fold_id    <- did:::build_crossfit_folds_edid(panel_obj$n, 5L, seed = 1L)
  pairs_nuis <- pairs
  pairs_nuis$gp[is.finite(pairs_nuis$gp) & pairs_nuis$gp == g] <- Inf
  prop_r <- did:::estimate_all_propensity_ratios(panel_obj, g, pairs_nuis,
                                                  4L, 5L, fold_id)
  cond_m <- did:::estimate_all_conditional_means(panel_obj, pairs_nuis, t,
                                                  4L, 5L, fold_id)

  gen_out <- did:::compute_generated_outcomes_cov_edid(panel_obj, g, t, pairs,
                                                        prop_r, cond_m, "post")
  # Cohort 5 units are neither G=g nor G=Inf, so their phi should be ~0
  mask_g5   <- (panel_obj$unit_cohorts == 5)
  phi_col1  <- gen_out[, 1L]
  phi_g5    <- phi_col1[mask_g5]
  expect_true(all(abs(phi_g5) < 1e-10),
              info = paste("phi for cohort-5 units in self-pair:", round(phi_g5, 4)))
})

# ============================================================
# Generated outcomes: E[phi] ≈ ATT for post-treatment cells
# ============================================================

test_that("mean generated outcome (weighted) approximates ATT for post-treatment cell", {
  # With a linear DGP where the true ATT is known (approx 1.0),
  # mean(gen_out %*% w) should be close to 1 for large enough n.
  set.seed(100)
  n_units <- 200; T <- 5
  ids   <- rep(1:n_units, each = T)
  times <- rep(1:T, times = n_units)
  g_unit <- rep(c(3, Inf), each = n_units / 2)
  g_vec  <- g_unit[ids]
  x1u    <- rnorm(n_units)
  x1     <- rep(x1u, each = T)
  y      <- 0.5 * times + 0.5 * x1 + as.numeric(times >= g_vec) + rnorm(n_units * T, sd = 0.5)
  df     <- data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1)

  fit  <- edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 1L)
  # ATT for cell (g=3, t=3) should be around 1.0
  cell_att <- fit$att_gt$att[fit$att_gt$group == 3 & fit$att_gt$time == 3]
  if (length(cell_att) == 1L) {
    expect_true(abs(cell_att - 1.0) < 0.5,
                info = paste("ATT(3,3) =", round(cell_att, 3), "; expected ~1.0"))
  }
})

# ============================================================
# Regression: covariate-path Omega* self-pair term5 must condition on G=g (Eq 3.12),
# NOT on G=Inf. The bug (self-pairs remapped to Inf) corrupted Omega* and the efficient
# weights (even negative weights) under cohort-specific pre-period heteroskedasticity with
# >=3 pre-periods. Fix: term5 conditions on the true cohort label g'_j (=g for self-pairs),
# matching the no-covariate path. On near-constant X the cov-path must match the no-cov path.
# ============================================================
test_that("cov-path Omega* self-pair term5 conditions on G=g (matches no-cov; no negative weights)", {
  set.seed(101); n <- 6000L; TP <- 5L; g_cohort <- 4L
  x1 <- rnorm(n, sd = 0.02)                                  # near-constant -> cond cov = uncond
  G  <- ifelse(runif(n) < 0.5, Inf, g_cohort); mu <- rnorm(n)
  sdc <- ifelse(is.finite(G), 1.6, 0.5); rhoc <- ifelse(is.finite(G), 0.7, 0.2)  # cohort heterosk.
  eps <- matrix(0, n, TP); eps[, 1] <- rnorm(n, sd = sdc)
  for (k in 2:TP) eps[, k] <- rhoc * eps[, k - 1] + sqrt(1 - rhoc^2) * rnorm(n, sd = sdc)
  rows <- lapply(1:TP, function(k) {
    tr <- as.numeric(is.finite(G) & k >= G)
    data.frame(id = 1:n, t = k, y = mu + 0.3 * k + 1.0 * tr + eps[, k], g = G, x1 = x1)
  })
  df <- do.call(rbind, rows)
  panel <- did:::prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1)
  g <- g_cohort; t <- 4L
  pairs <- did:::enumerate_valid_pairs_edid(g, panel$treatment_groups, panel$time_periods,
                                            panel$period_1, "all", panel$anticipation)
  expect_gte(nrow(pairs), 3L)                                # >=3 pre-periods -> self-pairs at t' != 1
  fold_id <- did:::build_crossfit_folds_edid(panel$n, 5L, seed = 1L)
  pn <- pairs; pn$gp[is.finite(pn$gp) & pn$gp == g] <- Inf   # nuisances use G_inf comparison (paper text)
  prop_r <- did:::estimate_all_propensity_ratios(panel, g, pn, 4L, 5L, fold_id)
  cond_m <- did:::estimate_all_conditional_means(panel, pn, t, 4L, 5L, fold_id)
  Om_cov   <- did:::compute_omega_star_cov_edid(panel, g, t, pairs, prop_r, cond_m)
  Om_nocov <- did:::compute_omega_star_nocov_edid(g, t, pairs, panel, "all")
  w_cov   <- did:::compute_efficient_weights_edid(Om_cov)
  w_nocov <- did:::compute_efficient_weights_edid(Om_nocov)
  expect_true(all(w_cov > -1e-6),
              info = paste("cov-path efficient weights must be non-negative; got",
                           paste(round(w_cov, 3), collapse = ", ")))
  expect_equal(w_cov, w_nocov, tolerance = 0.05,
               info = "cov-path and no-cov Omega* must give matching efficient weights")
})
