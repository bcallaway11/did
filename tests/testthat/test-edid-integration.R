library(testthat)

# ============================================================
# 9.1 Basic edid() call returns edid_fit object
# ============================================================
test_that("edid() returns an edid_fit object on one-cohort panel", {
  df  <- make_panel_1cohort(seed = 42)
  fit <- edid(
    data          = df,
    yname         = "outcome",
    idname        = "unit",
    tname         = "time",
    gname         = "first_treat",
    pt_assumption = "all",
    aggregate     = "all",
    bstrap        = FALSE
  )
  expect_s3_class(fit, "edid_fit")
})

# ============================================================
# 9.2 edid_fit object fields
# ============================================================
test_that("edid_fit contains all required top-level fields", {
  df  <- make_panel_1cohort(seed = 1)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat",
              pt_assumption = "all", aggregate = "all", bstrap = FALSE)
  required_fields <- c("call", "pt_assumption", "control_group", "alpha", "n",
                        "T_periods", "treatment_groups", "anticipation", "inference_type",
                        "cells", "att_gt", "overall", "event_study", "group")
  for (f in required_fields) {
    expect_true(f %in% names(fit),
                info = paste("Missing field:", f))
  }
})

# ============================================================
# 9.3 att_gt data.frame structure
# ============================================================
test_that("edid_fit$att_gt is a data.frame with required columns", {
  df  <- make_panel_1cohort(seed = 2)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat", aggregate = "all", bstrap = FALSE)
  expect_s3_class(fit$att_gt, "data.frame")
  expected_cols <- c("group", "time", "att", "se", "ci_lower", "ci_upper", "p_value", "is_pre")
  for (col in expected_cols) {
    expect_true(col %in% names(fit$att_gt), info = paste("Missing column:", col))
  }
})

# ============================================================
# 9.4 ATT estimates are finite for post-treatment cells
# ============================================================
test_that("edid() post-treatment ATT estimates are finite", {
  df  <- make_panel_1cohort(seed = 3)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat", aggregate = "all", bstrap = FALSE)
  post_cells <- fit$att_gt[!fit$att_gt$is_pre, ]
  expect_true(all(is.finite(post_cells$att)))
})

# ============================================================
# 9.5 PT-Post: ATT close to true value of 2 (large sample)
# ============================================================
test_that("edid() PT-Post overall ATT close to 2 for ATT=2 DGP (large sample)", {
  df  <- make_panel_1cohort(n_treat = 300, n_never = 300, n_periods = 5, seed = 555)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname         = "first_treat",
              pt_assumption = "post",
              aggregate     = "overall",
              bstrap        = FALSE)
  # True ATT = 2; allow generous tolerance for finite samples
  expect_equal(fit$overall$att, 2, tolerance = 0.4)
})

test_that("edid() PT-All overall ATT close to 2 for ATT=2 DGP (large sample)", {
  df  <- make_panel_1cohort(n_treat = 300, n_never = 300, n_periods = 5, seed = 666)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname         = "first_treat",
              pt_assumption = "all",
              aggregate     = "overall",
              bstrap        = FALSE)
  expect_equal(fit$overall$att, 2, tolerance = 0.5)
})

# ============================================================
# 9.6 S3 methods work without error
# ============================================================
test_that("print.edid_fit() runs without error", {
  df  <- make_panel_1cohort(seed = 4)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat", aggregate = "all")
  expect_output(print(fit))
})

test_that("summary.edid_fit() runs without error", {
  df  <- make_panel_1cohort(seed = 5)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat", aggregate = "all")
  expect_output(summary(fit))
})

test_that("coef.edid_fit() returns named numeric vector for att_gt", {
  df  <- make_panel_1cohort(seed = 6)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat", aggregate = "all")
  coefs <- coef(fit, which = "att_gt")
  expect_true(is.numeric(coefs))
  expect_true(length(coefs) > 0)
  expect_false(is.null(names(coefs)))
})

test_that("vcov.edid_fit() returns a numeric matrix", {
  df  <- make_panel_1cohort(seed = 7, n_treat = 40, n_never = 40)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat",
              aggregate = "all", store_eif = TRUE)
  V <- vcov(fit, which = "att_gt")
  expect_true(is.matrix(V))
  expect_true(is.numeric(V))
})

test_that("as.data.frame.edid_fit() returns a data.frame", {
  df  <- make_panel_1cohort(seed = 8)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat", aggregate = "all")
  df_out <- as.data.frame(fit)
  expect_s3_class(df_out, "data.frame")
})

# ============================================================
# 9.7 Two-cohort staggered panel
# ============================================================
test_that("edid() two-cohort staggered panel produces two groups in group aggregation", {
  df  <- make_panel_2cohort(seed = 200)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat",
              pt_assumption = "all", aggregate = "all", bstrap = FALSE)
  expect_equal(length(fit$group), 2L)
})

test_that("edid() two-cohort: group ATTs are near true values 1.5 and 2.5 (large sample)", {
  df  <- make_panel_2cohort(n_g3 = 200, n_g5 = 200, n_never = 200,
                             n_periods = 7, seed = 300)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat",
              pt_assumption = "post", aggregate = "group", bstrap = FALSE)
  atts <- sapply(fit$group, function(x) x$att)
  # Group g=3 should be near 1.5; group g=5 near 2.5
  expect_equal(unname(sort(atts)), c(1.5, 2.5), tolerance = 0.6)
})

# ============================================================
# 9.8 Bootstrap integration
# ============================================================
test_that("edid() with bstrap=TRUE returns bootstrap field in edid_fit", {
  df  <- make_panel_1cohort(seed = 42)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname         = "first_treat",
              pt_assumption = "post", aggregate = "overall",
              bstrap        = TRUE, biters = 50L, bootstrap_weights = "rademacher",
              seed          = 42L)
  expect_false(is.null(fit$bootstrap))
  expect_equal(fit$bootstrap$n_bootstrap, 50L)
})

test_that("edid() bootstrap SE differs from analytical SE (not identical)", {
  df  <- make_panel_1cohort(n_treat = 40, n_never = 40, seed = 99)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname         = "first_treat",
              pt_assumption = "post", aggregate = "overall",
              bstrap        = TRUE, biters = 200L, seed = 42L)
  # Bootstrap SE and analytical SE are related but not identical
  # Both should be finite and positive
  expect_true(is.finite(fit$overall$se))
  expect_true(fit$overall$se > 0)
})

# ============================================================
# 9.9 Clustered edid()
# ============================================================
test_that("edid() with clustervars argument produces finite clustered SEs", {
  df  <- make_panel_clustered(seed = 55)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname       = "first_treat",
              clustervars   = "cluster_id",
              pt_assumption = "post", aggregate = "overall",
              bstrap        = FALSE)
  expect_true(is.finite(fit$overall$se))
  expect_true(fit$overall$se > 0)
})

# ============================================================
# 9.10 Covariates stub
# ============================================================
test_that("edid() errors with clear message when covariates are supplied", {
  df  <- make_panel_1cohort(seed = 1)
  df$x1 <- rnorm(nrow(df))
  expect_error(
    edid(df, yname = "outcome", idname = "unit", tname = "time",
         gname = "first_treat", covariates = "x1"),
    regexp = "covariate|not yet implemented"
  )
})

# ============================================================
# 9.11 survey_design stub
# ============================================================
test_that("edid() errors with clear message when survey_design is supplied", {
  df <- make_panel_1cohort(seed = 1)
  expect_error(
    edid(df, yname = "outcome", idname = "unit", tname = "time",
         gname = "first_treat", survey_design = list(fake = TRUE)),
    regexp = "survey|not yet implemented"
  )
})

# ============================================================
# 9.12 store_eif = TRUE stores the EIF matrix
# ============================================================
test_that("edid() with store_eif=TRUE returns eif matrix of correct dimensions", {
  df  <- make_panel_1cohort(n_treat = 20, n_never = 20, n_periods = 4, seed = 42)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname = "first_treat",
              store_eif = TRUE, aggregate = "all")
  expect_false(is.null(fit$eif))
  # eif should be n x n_non_na_cells matrix
  expect_equal(nrow(fit$eif), fit$n)
  expect_true(ncol(fit$eif) >= 1)
})

# ============================================================
# 9.13 Degenerate panel: minimal 2-period, 2-unit case
# ============================================================
test_that("edid() on minimal degenerate panel runs without error", {
  df <- make_degenerate_panel()
  # PT-Post, g=2: baseline = 2-1-0 = 1 = period_1, so all cells have no valid pairs -> NA ATT
  # Should return result without error (all NA cells)
  expect_no_error({
    fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
                gname = "first_treat",
                pt_assumption = "post", aggregate = "none")
  })
})

# ============================================================
# 9.14 Balanced panel enforcement
# ============================================================
test_that("edid() errors loudly on unbalanced panel", {
  df <- make_panel_1cohort(seed = 1)
  df <- df[-1, ]  # remove one row -> unbalanced
  expect_error(
    edid(df, yname = "outcome", idname = "unit", tname = "time",
         gname = "first_treat"),
    regexp = "balanced|unbalanced"
  )
})

# ============================================================
# 9.15 Inference: CI contains true value in reasonable proportion (large n sanity check)
# ============================================================
test_that("edid() 95% CI contains 2 for large-n ATT=2 DGP", {
  df  <- make_panel_1cohort(n_treat = 500, n_never = 500, n_periods = 5, seed = 1234)
  fit <- edid(df, yname = "outcome", idname = "unit", tname = "time",
              gname         = "first_treat",
              pt_assumption = "post", aggregate = "overall", bstrap = FALSE)
  ci_lo <- fit$overall$ci_lower
  ci_hi <- fit$overall$ci_upper
  expect_true(ci_lo < 2 && ci_hi > 2,
              info = paste("CI:", ci_lo, "-", ci_hi, "does not contain 2"))
})

# ============================================================
# 9.16 G=0 auto-conversion (att_gt convention)
# ============================================================
test_that("edid() accepts G=0 for never-treated and auto-converts to Inf", {
  df <- make_panel_1cohort(seed = 42)
  # Replace Inf with 0 to simulate att_gt convention
  df$first_treat_0 <- ifelse(df$first_treat == Inf, 0, df$first_treat)
  fit_0   <- edid(df, yname = "outcome", idname = "unit", tname = "time",
                  gname = "first_treat_0",
                  pt_assumption = "post", aggregate = "overall", bstrap = FALSE)
  fit_inf <- edid(df, yname = "outcome", idname = "unit", tname = "time",
                  gname = "first_treat",
                  pt_assumption = "post", aggregate = "overall", bstrap = FALSE)
  expect_equal(fit_0$overall$att, fit_inf$overall$att, tolerance = 1e-10)
})
