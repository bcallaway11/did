# =============================================================================
# Tests for ggdid.MP and ggdid.AGGTEobj plotting functions
# =============================================================================

# Shared setup
set.seed(20260401)
sp <- reset.sim()
data_gg <- build_sim_dataset(sp)

mp_gg <- suppressWarnings(suppressMessages(
  att_gt(yname = "Y", xformla = ~X, data = data_gg, tname = "period",
         idname = "id", gname = "G", est_method = "dr",
         bstrap = FALSE)
))

agg_dyn <- suppressWarnings(aggte(mp_gg, type = "dynamic"))
agg_grp <- suppressWarnings(aggte(mp_gg, type = "group"))
agg_cal <- suppressWarnings(aggte(mp_gg, type = "calendar"))
agg_sim <- suppressWarnings(aggte(mp_gg, type = "simple"))

# =============================================================================
# ggdid.MP tests
# =============================================================================

test_that("ggdid.MP returns a ggplot object", {
  p <- ggdid(mp_gg)
  expect_s3_class(p, "gg")
})

test_that("ggdid.MP accepts group parameter", {
  groups <- unique(mp_gg$group)
  p <- ggdid(mp_gg, group = groups[1])
  expect_s3_class(p, "gg")
})

test_that("ggdid.MP warns for non-existent group values", {
  expect_warning(
    ggdid(mp_gg, group = c(9999)),
    "do not exist"
  )
})

test_that("ggdid.MP accepts custom labels", {
  p <- ggdid(mp_gg, xlab = "Time", ylab = "Effect", title = "Test Plot")
  expect_s3_class(p, "gg")
})

# =============================================================================
# ggdid.AGGTEobj tests
# =============================================================================

test_that("ggdid works for type = dynamic", {
  p <- ggdid(agg_dyn)
  expect_s3_class(p, "gg")
})

test_that("ggdid works for type = group", {
  p <- suppressWarnings(ggdid(agg_grp))
  expect_s3_class(p, "gg")
})

test_that("ggdid works for type = calendar", {
  p <- ggdid(agg_cal)
  expect_s3_class(p, "gg")
})

test_that("ggdid errors for type = simple", {
  expect_error(
    ggdid(agg_sim),
    "not available"
  )
})

test_that("ggdid.AGGTEobj accepts custom labels and theme settings", {
  p <- ggdid(agg_dyn, xlab = "Event Time", ylab = "ATT",
             title = "Event Study", theming = TRUE, legend = TRUE)
  expect_s3_class(p, "gg")
})

test_that("ggdid.AGGTEobj works with theming=FALSE", {
  p <- ggdid(agg_dyn, theming = FALSE)
  expect_s3_class(p, "gg")
})

test_that("ggdid.AGGTEobj works with ref_line=NULL", {
  p <- ggdid(agg_dyn, ref_line = NULL)
  expect_s3_class(p, "gg")
})

test_that("splot works for group-type and renders without deprecation warnings", {
  # splot() is the path that uses the ggplot2 version-gated errorbar layer.
  # It is called via ggdid.AGGTEobj() for type="group". Verify the output is
  # a valid ggplot object and no deprecation warnings leak through.
  expect_warning(p <- ggdid(agg_grp), regexp = NA)
  expect_s3_class(p, "gg")
  # Verify the plot actually builds (catches any layer construction errors)
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})
