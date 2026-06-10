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

test_that("panel precompute also bit-identical with time-varying sampling weights", {
  # the precompute retains per-period weights (period_w[[min(t + tfac, pret)]]);
  # time-varying weights are the only sim quantity that can expose a wrong-period index
  old_opt <- getOption("did.disable_precompute")
  withr::defer(options(did.disable_precompute = old_opt))
  set.seed(303)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  d$tvw <- d$period + runif(nrow(d), -0.1, 0.1)
  for (bp in c("varying", "universal")) {
    options(did.disable_precompute = TRUE)
    ref <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
      tname = "period", idname = "id", gname = "G", weightsname = "tvw",
      faster_mode = FALSE, bstrap = FALSE, est_method = "dr", base_period = bp)))
    options(did.disable_precompute = FALSE)
    new <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
      tname = "period", idname = "id", gname = "G", weightsname = "tvw",
      faster_mode = FALSE, bstrap = FALSE, est_method = "dr", base_period = bp)))
    expect_equal(ref$att, new$att, tolerance = 0, label = paste("tvw", bp, "ATT"))
    expect_equal(as.matrix(ref$inffunc), as.matrix(new$inffunc), tolerance = 0,
                 label = paste("tvw", bp, "IF"))
  }
})

test_that("RC/unbalanced precompute is bit-identical to the legacy subset path", {
  # the positional assembly (per-period row indices + plain vectors) must equal
  # the legacy rightids/%in%/data.table-subset construction exactly; the escape
  # hatch options(did.disable_precompute = TRUE) forces the legacy path
  old_opt <- getOption("did.disable_precompute")
  withr::defer(options(did.disable_precompute = old_opt))
  set.seed(404)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  set.seed(99)
  d_ub <- d[-sample(nrow(d), floor(nrow(d) * 0.10)), ]
  configs <- list(
    rc_nev_var_dr = list(data = d, args = list(panel = FALSE,
      control_group = "nevertreated", base_period = "varying", est_method = "dr")),
    rc_nyt_uni_dr = list(data = d, args = list(panel = FALSE,
      control_group = "notyettreated", base_period = "universal", est_method = "dr")),
    rc_ipw = list(data = d, args = list(panel = FALSE,
      control_group = "nevertreated", base_period = "varying", est_method = "ipw")),
    rc_nyt_ant = list(data = d, args = list(panel = FALSE,
      control_group = "notyettreated", base_period = "varying", est_method = "dr",
      anticipation = 1)),
    unb_nyt_var = list(data = d_ub, args = list(allow_unbalanced_panel = TRUE,
      control_group = "notyettreated", base_period = "varying", est_method = "dr")),
    unb_nev_uni = list(data = d_ub, args = list(allow_unbalanced_panel = TRUE,
      control_group = "nevertreated", base_period = "universal", est_method = "dr"))
  )
  for (nm in names(configs)) {
    cfg <- configs[[nm]]
    common <- c(list(yname = "Y", xformla = ~X, data = cfg$data, tname = "period",
                     idname = "id", gname = "G", faster_mode = FALSE, bstrap = FALSE),
                cfg$args)
    options(did.disable_precompute = TRUE)
    ref <- suppressWarnings(suppressMessages(do.call(att_gt, common)))
    options(did.disable_precompute = FALSE)
    new <- suppressWarnings(suppressMessages(do.call(att_gt, common)))
    expect_equal(ref$att, new$att, tolerance = 0, label = paste(nm, "ATT"))
    expect_equal(ref$se,  new$se,  tolerance = 0, label = paste(nm, "se"))
    expect_equal(as.matrix(ref$inffunc), as.matrix(new$inffunc), tolerance = 0,
                 label = paste(nm, "IF"))
  }
})

test_that("RC/unbalanced precompute bit-identical with fix_weights base/first, incl. dropped units", {
  # the per-period .rowid/.w lookup vectors must reproduce the legacy per-cell
  # data[t_col == target_period, ] subset + match() exactly, including the
  # dropped-units branch (units missing from the target period) and its warnings
  old_opt <- getOption("did.disable_precompute")
  withr::defer(options(did.disable_precompute = old_opt))
  set.seed(505)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  d$tvw <- d$period + runif(nrow(d), -0.1, 0.1)
  set.seed(99)
  d_ub <- d[-sample(nrow(d), floor(nrow(d) * 0.10)), ]
  collect <- function(fw) {
    ws <- character(0)
    res <- withCallingHandlers(
      suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d_ub, tname = "period",
        idname = "id", gname = "G", allow_unbalanced_panel = TRUE, weightsname = "tvw",
        fix_weights = fw, est_method = "dr", faster_mode = FALSE, bstrap = FALSE)),
      warning = function(w) { ws[[length(ws) + 1]] <<- conditionMessage(w); invokeRestart("muffleWarning") })
    list(att = res$att, se = res$se, inffunc = as.matrix(res$inffunc), warns = ws)
  }
  for (fw in c("base_period", "first_period")) {
    options(did.disable_precompute = TRUE)
    ref <- collect(fw)
    options(did.disable_precompute = FALSE)
    new <- collect(fw)
    expect_equal(ref$att, new$att, tolerance = 0, label = paste(fw, "ATT"))
    expect_equal(ref$se,  new$se,  tolerance = 0, label = paste(fw, "se"))
    expect_equal(ref$inffunc, new$inffunc, tolerance = 0, label = paste(fw, "IF"))
    expect_identical(ref$warns, new$warns, label = paste(fw, "warnings"))
    # the dropped-units branch must actually be exercised
    expect_true(any(grepl("^Dropped", new$warns)), label = paste(fw, "drop branch hit"))
  }
})

test_that("RC/unbalanced precompute keeps slow-fast parity", {
  set.seed(606)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  set.seed(99)
  d_ub <- d[-sample(nrow(d), floor(nrow(d) * 0.10)), ]
  for (cfg in list(list(data = d, extra = list(panel = FALSE)),
                   list(data = d_ub, extra = list(allow_unbalanced_panel = TRUE)))) {
    common <- c(list(yname = "Y", xformla = ~X, data = cfg$data, tname = "period",
                     idname = "id", gname = "G", est_method = "dr", bstrap = FALSE),
                cfg$extra)
    res_slow <- suppressWarnings(suppressMessages(do.call(att_gt, c(common, list(faster_mode = FALSE)))))
    res_fast <- suppressWarnings(suppressMessages(do.call(att_gt, c(common, list(faster_mode = TRUE)))))
    lab <- paste(names(cfg$extra)[1])
    expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10, label = paste(lab, "ATT"))
    expect_equal(res_slow$se,  res_fast$se,  tolerance = 1e-10, label = paste(lab, "se"))
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
