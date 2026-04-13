library(testthat)

# ============================================================
# 8.1 generate_multiplier_weights_edid(): Rademacher
# ============================================================
test_that("generate_multiplier_weights_edid() Rademacher: values are in {-1, +1}", {
  set.seed(42)
  W <- generate_multiplier_weights_edid(n = 50L, n_bootstrap = 200L,
                                         type = "rademacher", seed = 42L)
  expect_equal(dim(W), c(50L, 200L))
  unique_vals <- unique(as.vector(W))
  expect_true(all(unique_vals %in% c(-1, 1)))
})

test_that("generate_multiplier_weights_edid() Rademacher: roughly 50/50 split", {
  set.seed(1)
  W <- generate_multiplier_weights_edid(n = 1000L, n_bootstrap = 1L,
                                         type = "rademacher", seed = 1L)
  prop_pos <- mean(W == 1)
  expect_true(prop_pos > 0.4 && prop_pos < 0.6)
})

# ============================================================
# 8.2 generate_multiplier_weights_edid(): Mammen
# ============================================================
test_that("generate_multiplier_weights_edid() Mammen: values are in the two Mammen values", {
  set.seed(42)
  W <- generate_multiplier_weights_edid(n = 100L, n_bootstrap = 100L,
                                         type = "mammen", seed = 42L)
  mammen_lo <- -(sqrt(5) - 1) / 2
  mammen_hi <- (sqrt(5) + 1) / 2
  unique_vals <- unique(as.vector(W))
  expect_true(all(abs(unique_vals - mammen_lo) < 1e-10 |
                    abs(unique_vals - mammen_hi) < 1e-10))
})

# ============================================================
# 8.3 generate_multiplier_weights_edid(): Webb
# ============================================================
test_that("generate_multiplier_weights_edid() Webb: values are among 6 Webb values", {
  set.seed(42)
  W <- generate_multiplier_weights_edid(n = 200L, n_bootstrap = 100L,
                                         type = "webb", seed = 42L)
  webb_vals <- c(-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2))
  unique_vals <- unique(as.vector(W))
  for (v in unique_vals) {
    expect_true(any(abs(v - webb_vals) < 1e-10))
  }
})

# ============================================================
# 8.4 generate_multiplier_weights_edid(): cluster expansion
# ============================================================
test_that("generate_multiplier_weights_edid() cluster: units in same cluster get identical weight", {
  set.seed(42)
  G <- 5L
  n <- 20L  # 4 units per cluster
  cluster_idx <- rep(1:G, each = 4L)
  W <- generate_multiplier_weights_edid(n = n, n_bootstrap = 50L,
                                         type = "rademacher",
                                         cluster_indices = cluster_idx,
                                         seed = 42L)
  expect_equal(dim(W), c(n, 50L))
  # All units in cluster 1 must have the same weight for each draw
  for (b in seq_len(50)) {
    expect_true(length(unique(W[cluster_idx == 1L, b])) == 1L)
  }
})

# ============================================================
# 8.5 generate_multiplier_weights_edid(): reproducibility with seed
# ============================================================
test_that("generate_multiplier_weights_edid() is reproducible with seed", {
  W1 <- generate_multiplier_weights_edid(50L, 100L, "rademacher", seed = 777L)
  W2 <- generate_multiplier_weights_edid(50L, 100L, "rademacher", seed = 777L)
  expect_identical(W1, W2)
})

# ============================================================
# 8.6 compute_bootstrap_stats_edid(): output structure
# ============================================================
test_that("compute_bootstrap_stats_edid() returns named list with correct fields", {
  set.seed(1)
  att_hat    <- 1.5
  boot_draws <- rnorm(200, mean = att_hat, sd = 0.3)
  res <- compute_bootstrap_stats_edid(boot_draws, att_hat, alpha = 0.05)
  expect_named(res, c("se_boot", "ci_lower", "ci_upper", "p_value_boot"),
               ignore.order = TRUE)
})

test_that("compute_bootstrap_stats_edid() SE is positive", {
  set.seed(2)
  att_hat <- 0.8
  boot_draws <- rnorm(500, att_hat, 0.5)
  res <- compute_bootstrap_stats_edid(boot_draws, att_hat)
  expect_true(res$se_boot > 0)
})

test_that("compute_bootstrap_stats_edid() CI lower < upper", {
  set.seed(3)
  att_hat <- 2.0
  boot_draws <- rnorm(200, att_hat, 0.4)
  res <- compute_bootstrap_stats_edid(boot_draws, att_hat)
  expect_true(res$ci_lower < res$ci_upper)
})

test_that("compute_bootstrap_stats_edid() p_value in [0, 1]", {
  set.seed(4)
  att_hat <- 0.0
  boot_draws <- rnorm(300, 0, 0.5)
  res <- compute_bootstrap_stats_edid(boot_draws, att_hat)
  expect_true(res$p_value_boot >= 0 && res$p_value_boot <= 1)
})

# ============================================================
# 8.7 run_multiplier_bootstrap_edid(): basic smoke test
# ============================================================
test_that("run_multiplier_bootstrap_edid() runs without error and returns bootstrap draws", {
  df    <- make_panel_1cohort(seed = 77)
  panel <- prepare_edid_panel(df, yname = "outcome", idname = "unit",
                               tname = "time", gname = "first_treat")
  fit   <- fit_edid_cells(panel, pt_assumption = "all", alpha = 0.05,
                           store_eif = TRUE, covariates = NULL)
  boot_res <- run_multiplier_bootstrap_edid(
    cells           = fit$cells,
    eif_matrix      = fit$eif_matrix,
    cell_index      = fit$cell_index,
    panel_obj       = panel,
    n_bootstrap     = 100L,
    bootstrap_weights = "rademacher",
    seed            = 42L,
    aggregate       = "overall",
    alpha           = 0.05
  )
  expect_true(!is.null(boot_res$overall_b))
  expect_equal(length(boot_res$overall_b), 100L)
  expect_true(all(is.finite(boot_res$overall_b)))
})
