# =============================================================================
# Systematic faster_mode=TRUE vs FALSE consistency tests
# =============================================================================

# Shared setup
set.seed(20260401)
sp <- reset.sim()
data_fm <- build_sim_dataset(sp)

# Unbalanced version
data_ub <- data_fm[-c(1, 5, 10), ]

# Helper to compare two att_gt results
compare_modes <- function(res_slow, res_fast, label) {
  expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10, label = paste(label, "ATT"))
  expect_equal(res_slow$group, res_fast$group, label = paste(label, "group"))
  expect_equal(res_slow$t, res_fast$t, label = paste(label, "t"))
}

# =============================================================================
# Core grid: est_method x panel_type x control_group x base_period
# =============================================================================

# --- Balanced panel ---

for (em in c("dr", "ipw", "reg")) {
  for (cg in c("nevertreated", "notyettreated")) {
    for (bp in c("varying", "universal")) {
      label <- paste(em, "balanced", cg, bp)
      test_that(paste("consistency:", label), {
        res_slow <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 faster_mode = FALSE, bstrap = FALSE)
        ))
        res_fast <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 faster_mode = TRUE, bstrap = FALSE)
        ))
        compare_modes(res_slow, res_fast, label)
      })
    }
  }
}

# --- Unbalanced panel ---

for (em in c("dr", "ipw", "reg")) {
  for (cg in c("nevertreated", "notyettreated")) {
    for (bp in c("varying", "universal")) {
      label <- paste(em, "unbalanced", cg, bp)
      test_that(paste("consistency:", label), {
        res_slow <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_ub, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 allow_unbalanced_panel = TRUE,
                 faster_mode = FALSE, bstrap = FALSE)
        ))
        res_fast <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_ub, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 allow_unbalanced_panel = TRUE,
                 faster_mode = TRUE, bstrap = FALSE)
        ))
        compare_modes(res_slow, res_fast, label)
      })
    }
  }
}

# --- Repeated cross sections ---

for (em in c("dr", "ipw", "reg")) {
  for (cg in c("nevertreated", "notyettreated")) {
    for (bp in c("varying", "universal")) {
      label <- paste(em, "RC", cg, bp)
      test_that(paste("consistency:", label), {
        res_slow <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 panel = FALSE,
                 faster_mode = FALSE, bstrap = FALSE)
        ))
        res_fast <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 panel = FALSE,
                 faster_mode = TRUE, bstrap = FALSE)
        ))
        compare_modes(res_slow, res_fast, label)
      })
    }
  }
}

# =============================================================================
# Additional consistency tests beyond the core grid
# =============================================================================

test_that("consistency with anticipation > 0", {
  res_slow <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           anticipation = 1, faster_mode = FALSE, bstrap = FALSE)
  ))
  res_fast <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           anticipation = 1, faster_mode = TRUE, bstrap = FALSE)
  ))
  compare_modes(res_slow, res_fast, "anticipation=1")
})

test_that("consistency without covariates", {
  res_slow <- att_gt(yname = "Y", data = data_fm, tname = "period",
                     idname = "id", gname = "G", est_method = "reg",
                     faster_mode = FALSE, bstrap = FALSE)
  res_fast <- att_gt(yname = "Y", data = data_fm, tname = "period",
                     idname = "id", gname = "G", est_method = "reg",
                     faster_mode = TRUE, bstrap = FALSE)
  compare_modes(res_slow, res_fast, "no covariates")
})

# =============================================================================
# aggte consistency across modes
# =============================================================================

test_that("aggte consistency across faster_mode for all types", {
  mp_slow <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           faster_mode = FALSE, bstrap = FALSE)
  ))
  mp_fast <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           faster_mode = TRUE, bstrap = FALSE)
  ))

  for (tp in c("simple", "dynamic", "group", "calendar")) {
    agg_slow <- suppressWarnings(aggte(mp_slow, type = tp))
    agg_fast <- suppressWarnings(aggte(mp_fast, type = tp))
    expect_equal(agg_slow$overall.att, agg_fast$overall.att,
                 tolerance = 1e-8, label = paste("aggte", tp))
  }
})
