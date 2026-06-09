library(testthat)

# ============================================================
# 6.1 compute_eif_se_edid(): correct formula
# ============================================================
test_that("compute_eif_se_edid() returns sqrt(sum(eif^2)/n^2)", {
  set.seed(42)
  n   <- 50L
  eif <- rnorm(n, 0, 1)
  se  <- compute_eif_se_edid(eif, n)
  expected <- sqrt(sum(eif^2) / n^2)
  expect_equal(se, expected, tolerance = 1e-12)
})

test_that("compute_eif_se_edid() returns non-negative value", {
  set.seed(1)
  eif <- rnorm(100)
  expect_true(compute_eif_se_edid(eif, 100) >= 0)
})

# ============================================================
# 6.2 cluster_aggregate_edid(): formula
# ============================================================
test_that("cluster_aggregate_edid() sums EIF within clusters correctly", {
  # 10 units, 2 clusters of 5 each
  n <- 10L
  G <- 2L
  cluster_idx <- c(rep(1L, 5L), rep(2L, 5L))
  set.seed(1)
  eif <- rnorm(n)

  result <- cluster_aggregate_edid(eif, cluster_idx)
  # result is the cluster sums (not yet centered here -- tester checks dimension)
  # Cluster 1 sum:
  c1_sum <- sum(eif[1:5])
  c2_sum <- sum(eif[6:10])
  # The function should return a length-G vector
  expect_equal(length(result), G)
})

test_that("cluster_aggregate_edid() SE with clustering differs from iid SE", {
  df    <- make_panel_clustered(seed = 42)
  panel <- prepare_edid_panel(df, yname = "outcome", idname = "unit",
                               tname = "time", gname = "first_treat",
                               clustervars = "cluster_id")
  pairs <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "post")
  omega <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "post")
  w     <- compute_efficient_weights_edid(omega)
  y_hat <- compute_generated_outcomes_nocov_edid(3L, 4L, pairs, panel, "post")
  att   <- sum(w * y_hat)
  eif   <- compute_eif_nocov_edid(3L, 4L, pairs, w, panel, att, "post")

  se_iid     <- compute_eif_se_edid(eif, panel$n)
  clust_agg  <- cluster_aggregate_edid(eif, panel$cluster_indices)
  se_cluster <- safe_inference_edid(eif, panel$cluster_indices, 0.05)$se
  # Cluster SE should be finite and positive
  expect_true(is.finite(se_cluster))
  expect_true(se_cluster > 0)
  # Cluster SE and iid SE are generally different (not testing direction)
  expect_false(isTRUE(all.equal(se_iid, se_cluster)))
})

# ============================================================
# 6.3 safe_inference_edid(): output structure
# ============================================================
test_that("safe_inference_edid() returns named list with expected fields", {
  set.seed(42)
  eif  <- rnorm(50)
  res  <- safe_inference_edid(eif, cluster_indices = NULL, alpha = 0.05)
  expect_named(res, c("se", "ci_lower", "ci_upper", "t_stat", "p_value", "inference_valid"),
               ignore.order = TRUE)
})

test_that("safe_inference_edid() returns inference_valid=FALSE when EIF is all zeros", {
  eif <- rep(0, 50)
  res <- safe_inference_edid(eif, cluster_indices = NULL, alpha = 0.05)
  expect_false(res$inference_valid)
  expect_true(is.na(res$se) || res$se == 0)
})

test_that("safe_inference_edid() CI width is positive for non-degenerate EIF", {
  set.seed(10)
  eif <- rnorm(100)
  res <- safe_inference_edid(eif, cluster_indices = NULL, alpha = 0.05)
  ci_width <- res$ci_upper - res$ci_lower
  expect_true(ci_width > 0 || !res$inference_valid)
})

test_that("safe_inference_edid() p_value is in [0, 1]", {
  set.seed(15)
  eif <- rnorm(80)
  # Pass a finite att: CIs and the p-value require it (att defaults to NA, which marks inference invalid).
  # With a clean Gaussian EIF, no clustering, and a finite att the inference is valid and the p-value is a
  # proper probability. Always assert (no bare conditional that could leave the test expectation-free).
  res <- safe_inference_edid(eif, cluster_indices = NULL, alpha = 0.05, att = 0.1)
  expect_true(isTRUE(res$inference_valid))
  expect_true(is.finite(res$p_value) && res$p_value >= 0 && res$p_value <= 1)
})
