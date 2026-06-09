# =============================================================================
# Slow-path (faster_mode = FALSE) panel precompute: the per-period-block fast path
# must be BIT-IDENTICAL to the original get_wide_data() reshape. The escape hatch
# options(did.disable_precompute = TRUE) forces the original path, so we can compare
# the two implementations directly. Also guards the precompute precondition (the
# panel is re-sorted by (id, period) in preprocessing, so input row order is
# irrelevant).
# =============================================================================

test_that("panel precompute is bit-identical to the original get_wide_data path", {
  old_opt <- getOption("did.disable_precompute")
  withr::defer(options(did.disable_precompute = old_opt))
  set.seed(101)
  sp <- did::reset.sim(time.periods = 5)
  d <- did::build_sim_dataset(sp)
  for (cg in c("nevertreated", "notyettreated"))
  for (em in c("dr", "reg", "ipw"))
  for (bp in c("varying", "universal"))
  for (ant in c(0, 1)) {
    options(did.disable_precompute = TRUE)
    ref <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
      tname = "period", idname = "id", gname = "G", faster_mode = FALSE, bstrap = FALSE,
      control_group = cg, est_method = em, base_period = bp, anticipation = ant)))
    options(did.disable_precompute = FALSE)
    new <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
      tname = "period", idname = "id", gname = "G", faster_mode = FALSE, bstrap = FALSE,
      control_group = cg, est_method = em, base_period = bp, anticipation = ant)))
    lab <- paste(cg, em, bp, "ant", ant)
    expect_equal(ref$att, new$att, tolerance = 0, label = paste(lab, "ATT"))
    expect_equal(as.matrix(ref$inffunc), as.matrix(new$inffunc), tolerance = 0,
                 label = paste(lab, "IF"))
    expect_equal(is.na(ref$att), is.na(new$att), label = paste(lab, "NA pattern"))
  }
})

test_that("panel precompute also bit-identical for factor / transform covariates", {
  old_opt <- getOption("did.disable_precompute")
  withr::defer(options(did.disable_precompute = old_opt))
  set.seed(202)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  d$fac <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  for (xf in c(~X, ~fac, ~I(X^2), ~X + fac)) {
    options(did.disable_precompute = TRUE)
    ref <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = xf, data = d,
      tname = "period", idname = "id", gname = "G", faster_mode = FALSE, bstrap = FALSE)))
    options(did.disable_precompute = FALSE)
    new <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = xf, data = d,
      tname = "period", idname = "id", gname = "G", faster_mode = FALSE, bstrap = FALSE)))
    expect_equal(ref$att, new$att, tolerance = 0, label = deparse(xf))
    expect_equal(as.matrix(ref$inffunc), as.matrix(new$inffunc), tolerance = 0,
                 label = deparse(xf))
  }
})

test_that("input row order does not affect slow-path results (precompute precondition)", {
  set.seed(7)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  ref <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
    tname = "period", idname = "id", gname = "G", faster_mode = FALSE, bstrap = FALSE)))
  shf <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
    data = d[sample(nrow(d)), ], tname = "period", idname = "id", gname = "G",
    faster_mode = FALSE, bstrap = FALSE)))
  expect_equal(ref$att, shf$att, tolerance = 1e-12)
  expect_equal(ref$se, shf$se, tolerance = 1e-12)
})
