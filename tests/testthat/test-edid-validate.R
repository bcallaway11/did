library(testthat)
# helper-edid.R is auto-loaded

# ============================================================
# 3.1 Valid inputs pass without error
# ============================================================
test_that("validate_edid_inputs() passes on valid one-cohort panel", {
  df <- make_panel_1cohort()
  expect_silent(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    )
  )
})

test_that("validate_edid_inputs() passes on two-cohort panel", {
  df <- make_panel_2cohort()
  expect_silent(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "post",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 100L, anticipation = 0L, survey_design = NULL
    )
  )
})

# ============================================================
# 3.2 Missing column names
# ============================================================
test_that("validate_edid_inputs() errors on missing yname column", {
  df <- make_panel_1cohort()
  expect_error(
    validate_edid_inputs(
      data = df, yname = "y_outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "y_outcome"
  )
})

test_that("validate_edid_inputs() errors on missing tname column", {
  df <- make_panel_1cohort()
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "t_var", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "t_var"
  )
})

# ============================================================
# 3.3 Non-numeric outcome
# ============================================================
test_that("validate_edid_inputs() errors on character outcome column", {
  df <- make_panel_1cohort()
  df$outcome <- as.character(df$outcome)
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    )
  )
})

# ============================================================
# 3.4 Non-finite outcomes
# ============================================================
test_that("validate_edid_inputs() errors on Inf outcome", {
  df <- make_panel_1cohort()
  df$outcome[1] <- Inf
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "finite|non-finite|Inf"
  )
})

test_that("validate_edid_inputs() errors on NA outcome", {
  df <- make_panel_1cohort()
  df$outcome[5] <- NA_real_
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    )
  )
})

# ============================================================
# 3.5 Unbalanced panel
# ============================================================
test_that("validate_edid_inputs() errors on unbalanced panel", {
  df <- make_panel_1cohort()
  df_unbal <- df[-1, ]  # drop one row -> unit 1 missing period 1
  expect_error(
    validate_edid_inputs(
      data = df_unbal, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "balanced|unbalanced"
  )
})

# ============================================================
# 3.6 Duplicate (unit, time) rows
# ============================================================
test_that("validate_edid_inputs() errors on duplicate (idname, tname) rows", {
  df <- make_panel_1cohort()
  df_dup <- rbind(df, df[1, ])  # duplicate first row
  expect_error(
    validate_edid_inputs(
      data = df_dup, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "[Dd]uplicate"
  )
})

# ============================================================
# 3.7 Non-absorbing treatment
# ============================================================
test_that("validate_edid_inputs() errors on non-absorbing treatment", {
  df <- make_panel_1cohort()
  # Make first_treat time-varying within unit 1
  df$first_treat[df$unit == 1 & df$time == 2] <- 4L
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "absorbing|time-varying|constant"
  )
})

# ============================================================
# 3.8 No never-treated units
# ============================================================
test_that("validate_edid_inputs() errors when no never-treated units and control_group='nevertreated'", {
  df <- make_panel_1cohort()
  # relabel all never-treated as cohort 4
  df$first_treat[df$first_treat == Inf] <- 4L
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "never.treated|nevertreated"
  )
})

# ============================================================
# 3.9 Covariate stub
# ============================================================
test_that("validate_edid_inputs() errors when covariates supplied (stub)", {
  df <- make_panel_1cohort()
  df$x1 <- rnorm(nrow(df))
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = "x1", pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "covariate|not yet implemented"
  )
})

# ============================================================
# 3.10 Survey stub
# ============================================================
test_that("validate_edid_inputs() errors when survey_design supplied (stub)", {
  df <- make_panel_1cohort()
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L,
      survey_design = list(strata = "fake")
    ),
    regexp = "survey|not yet implemented"
  )
})

# ============================================================
# 3.11 Invalid alp
# ============================================================
test_that("validate_edid_inputs() errors on alp outside (0,1)", {
  df <- make_panel_1cohort()
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 1.5, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    )
  )
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0, clustervars = NULL, control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    )
  )
})

# ============================================================
# 3.12 Cluster column time-varying check
# ============================================================
test_that("validate_edid_inputs() errors on time-varying cluster variable", {
  df <- make_panel_clustered()
  # Make cluster_id time-varying for unit 1
  df$cluster_id[df$unit == 1 & df$time == 2] <- 999L
  expect_error(
    validate_edid_inputs(
      data = df, yname = "outcome", idname = "unit",
      tname = "time", gname = "first_treat",
      covariates = NULL, pt_assumption = "all",
      alp = 0.05, clustervars = "cluster_id", control_group = "nevertreated",
      biters = 0L, anticipation = 0L, survey_design = NULL
    ),
    regexp = "cluster|time.invariant|time-invariant"
  )
})
