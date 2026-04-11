library(testthat)

# Helper: run fit_edid_cells on one-cohort panel
fit_one_cohort <- function(seed = 42, n_treat = 30, n_never = 30, n_periods = 5) {
  df    <- make_panel_1cohort(n_treat = n_treat, n_never = n_never,
                               n_periods = n_periods, seed = seed)
  panel <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  fit   <- fit_edid_cells(panel, pt_assumption = "all", alpha = 0.05,
                           store_eif = TRUE, covariates = NULL)
  list(fit = fit, panel = panel)
}

# ============================================================
# 7.1 aggregate_overall_edid(): basic properties
# ============================================================
test_that("aggregate_overall_edid() returns finite ATT", {
  res   <- fit_one_cohort()
  agg   <- aggregate_overall_edid(res$fit$cells, res$fit$eif_matrix,
                                   res$fit$cell_index, res$panel, alpha = 0.05)
  expect_true(is.finite(agg$att))
})

test_that("aggregate_overall_edid() ATT is within [-10, 10] for ATT=2 DGP", {
  res <- fit_one_cohort(seed = 42, n_treat = 50, n_never = 50)
  agg <- aggregate_overall_edid(res$fit$cells, res$fit$eif_matrix,
                                 res$fit$cell_index, res$panel, alpha = 0.05)
  expect_true(agg$att > -10 && agg$att < 10)
})

test_that("aggregate_overall_edid() SE is positive and finite", {
  res <- fit_one_cohort()
  agg <- aggregate_overall_edid(res$fit$cells, res$fit$eif_matrix,
                                 res$fit$cell_index, res$panel, alpha = 0.05)
  expect_true(is.finite(agg$se))
  expect_true(agg$se > 0)
})

test_that("aggregate_overall_edid() EIF aggregated vector has zero mean", {
  res <- fit_one_cohort()
  agg <- aggregate_overall_edid(res$fit$cells, res$fit$eif_matrix,
                                 res$fit$cell_index, res$panel, alpha = 0.05)
  expect_equal(mean(agg$eif_agg), 0, tolerance = 1e-8)
})

# ============================================================
# 7.2 aggregate_event_study_edid(): structure
# ============================================================
test_that("aggregate_event_study_edid() returns a list with entries for each relative time", {
  res <- fit_one_cohort()
  es  <- aggregate_event_study_edid(res$fit$cells, res$fit$eif_matrix,
                                     res$fit$cell_index, res$panel, alpha = 0.05)
  # Should have entries for each unique e = t - g in the cells
  expect_true(is.list(es))
  expect_true(length(es) >= 1)
  for (entry in es) {
    expect_true(is.finite(entry$att))
    expect_true(is.finite(entry$se))
  }
})

# ============================================================
# 7.3 aggregate_group_edid(): equal-time weights
# ============================================================
test_that("aggregate_group_edid() returns one entry per treated cohort", {
  df    <- make_panel_2cohort()
  panel <- prepare_edid_panel(df, "outcome", "unit", "time", "first_treat")
  fit   <- fit_edid_cells(panel, pt_assumption = "all", alpha = 0.05,
                           store_eif = TRUE, covariates = NULL)
  grp <- aggregate_group_edid(fit$cells, fit$eif_matrix, fit$cell_index, panel, alpha = 0.05)
  expect_equal(length(grp), length(panel$treatment_groups))
  for (entry in grp) {
    expect_true(is.finite(entry$att))
  }
})

# ============================================================
# 7.4 WIF correction: EIF of aggregate has correct mean
# ============================================================
test_that("compute_wif_contribution_edid() produces a length-n numeric vector", {
  res <- fit_one_cohort()
  wif <- compute_wif_contribution_edid(
    weight_fn  = function(cells, cell_idx, po) {
      # simple uniform weight function for test
      post_cells <- Filter(function(c) !c$is_pre && !is.na(c$att), cells)
      n_post <- length(post_cells)
      setNames(rep(1/n_post, n_post), seq_len(n_post))
    },
    cells      = res$fit$cells,
    eif_matrix = res$fit$eif_matrix,
    cell_index = res$fit$cell_index,
    panel_obj  = res$panel
  )
  expect_equal(length(wif), res$panel$n)
  expect_true(is.numeric(wif))
})
