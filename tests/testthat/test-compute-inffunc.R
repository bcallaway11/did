# =============================================================================
# Tests for compute_inffunc = FALSE (point-estimates-only mode): skips influence
# functions, standard errors, bands, the pre-test, and the bootstrap. Point
# estimates must be identical to a full run.
# =============================================================================

test_that("compute_inffunc = FALSE gives point estimates identical to a full run", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  for (fm in c(TRUE, FALSE)) for (est in c("dr", "reg", "ipw")) for (pn in c(TRUE, FALSE)) {
    full <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", est_method = est, panel = pn, bstrap = FALSE, faster_mode = fm)))
    pe <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", est_method = est, panel = pn, faster_mode = fm,
             compute_inffunc = FALSE)))
    lab <- paste(est, "panel", pn, "fm", fm)
    expect_equal(full$att, pe$att, tolerance = 1e-12, label = paste(lab, "ATT"))
    expect_true(all(is.na(pe$se)), label = paste(lab, "se all NA"))
    expect_null(pe$inffunc, label = paste(lab, "inffunc NULL"))
    expect_null(pe$W, label = paste(lab, "Wald NULL"))
  }
})

test_that("compute_inffunc = FALSE forces bstrap/cband off and yields a valid MP object", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  pe <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
           gname = "G", bstrap = TRUE, cband = TRUE, compute_inffunc = FALSE)))
  expect_s3_class(pe, "MP")
  expect_false(pe$DIDparams$bstrap)
  expect_false(pe$DIDparams$cband)
  # print / summary / tidy / glance must work with NA SEs and no influence functions
  expect_output(print(pe))
  expect_output(summary(pe))
  td <- tidy(pe)
  expect_equal(nrow(td), length(pe$att))
  expect_true(all(is.na(td$std.error)))
  expect_s3_class(glance(pe), "data.frame")
})

test_that("aggte() errors clearly on a compute_inffunc = FALSE object", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  pe <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
           gname = "G", compute_inffunc = FALSE)))
  expect_error(suppressWarnings(suppressMessages(aggte(pe, type = "group"))),
               "compute_inffunc = FALSE")
})

test_that("compute_inffunc must be a single logical", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  expect_error(
    att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
           gname = "G", compute_inffunc = "yes"),
    "compute_inffunc must be a single logical")
})
