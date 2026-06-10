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

# -----------------------------------------------------------------------------
# Non-default option matrix: the do_inf skips on the varying-weights, unbalanced,
# universal-base, clustered, and custom-est_method branches must leave point
# estimates identical to a full run in both modes.
# -----------------------------------------------------------------------------

test_that("compute_inffunc = FALSE matches a full run with fix_weights = 'varying'", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  data$tvw <- data$period + runif(nrow(data), -0.1, 0.1)
  for (fm in c(TRUE, FALSE)) {
    full <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", weightsname = "tvw", fix_weights = "varying",
             bstrap = FALSE, faster_mode = fm)))
    pe <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", weightsname = "tvw", fix_weights = "varying",
             faster_mode = fm, compute_inffunc = FALSE)))
    lab <- paste("fix_weights varying fm", fm)
    expect_equal(full$att, pe$att, tolerance = 1e-12, label = paste(lab, "ATT"))
    expect_true(identical(is.na(full$att), is.na(pe$att)), label = paste(lab, "NA pattern"))
    expect_null(pe$inffunc, label = paste(lab, "inffunc NULL"))
  }
})

test_that("compute_inffunc = FALSE matches a full run with anticipation, universal base period, not-yet-treated controls, and clustering", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  for (fm in c(TRUE, FALSE)) {
    full <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", anticipation = 1, base_period = "universal",
             control_group = "notyettreated", clustervars = "cluster",
             bstrap = FALSE, faster_mode = fm)))
    pe <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", anticipation = 1, base_period = "universal",
             control_group = "notyettreated", clustervars = "cluster",
             faster_mode = fm, compute_inffunc = FALSE)))
    lab <- paste("anticipation/universal/nyt/cluster fm", fm)
    expect_equal(full$att, pe$att, tolerance = 1e-12, label = paste(lab, "ATT"))
    expect_true(identical(is.na(full$att), is.na(pe$att)), label = paste(lab, "NA pattern"))
    expect_null(pe$inffunc, label = paste(lab, "inffunc NULL"))
  }
})

test_that("compute_inffunc = FALSE matches a full run on an unbalanced panel", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  data_ub <- data[-sample(nrow(data), 150), ]
  for (fm in c(TRUE, FALSE)) {
    full <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data_ub, tname = "period", idname = "id",
             gname = "G", allow_unbalanced_panel = TRUE, bstrap = FALSE,
             faster_mode = fm)))
    pe <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data_ub, tname = "period", idname = "id",
             gname = "G", allow_unbalanced_panel = TRUE, faster_mode = fm,
             compute_inffunc = FALSE)))
    lab <- paste("unbalanced panel fm", fm)
    expect_equal(full$att, pe$att, tolerance = 1e-12, label = paste(lab, "ATT"))
    expect_true(identical(is.na(full$att), is.na(pe$att)), label = paste(lab, "NA pattern"))
    expect_null(pe$inffunc, label = paste(lab, "inffunc NULL"))
  }
})

test_that("custom est_method receives inffunc = FALSE and may return a NULL influence function", {
  set.seed(20260609)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  seen <- new.env(parent = emptyenv())
  my_panel_est <- function(y1, y0, D, covariates, i.weights, inffunc, ...) {
    seen$inffunc <- c(seen$inffunc, isTRUE(inffunc))
    delta <- y1 - y0
    att <- weighted.mean(delta[D == 1], i.weights[D == 1])
    inf <- if (isTRUE(inffunc)) (D / mean(D)) * (delta - att) else NULL
    list(ATT = att, att.inf.func = inf)
  }
  for (fm in c(TRUE, FALSE)) {
    seen$inffunc <- logical(0)
    full <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", est_method = my_panel_est, bstrap = FALSE,
             faster_mode = fm)))
    expect_true(length(seen$inffunc) > 0 && all(seen$inffunc),
                label = paste("full run passes inffunc = TRUE, fm", fm))
    seen$inffunc <- logical(0)
    pe <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~X, data = data, tname = "period", idname = "id",
             gname = "G", est_method = my_panel_est, faster_mode = fm,
             compute_inffunc = FALSE)))
    expect_true(length(seen$inffunc) > 0 && !any(seen$inffunc),
                label = paste("point-estimate run passes inffunc = FALSE, fm", fm))
    lab <- paste("custom est_method fm", fm)
    expect_s3_class(pe, "MP")
    expect_equal(full$att, pe$att, tolerance = 1e-12, label = paste(lab, "ATT"))
    expect_null(pe$inffunc, label = paste(lab, "inffunc NULL"))
  }
})
