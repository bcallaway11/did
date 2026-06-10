# =============================================================================
# Tests for covariate handling in compute.att_gt.R (the slow path). The design
# matrix is built once over the full data and row-subset per (g,t) cell, giving:
#   * numeric formulae      -> bit-identical to a per-cell model.matrix()
#   * factor covariates     -> identical to manually-expanded dummy columns
#   * transform formulae    -> work (previously errored); slow == fast
# =============================================================================

# helper: add treatment-contrast dummy columns for a factor, named f_<level>
add_dummies <- function(d, fac = "fac") {
  lv <- levels(d[[fac]])
  for (l in lv[-1L]) d[[paste0("f_", l)]] <- as.numeric(d[[fac]] == l)
  d
}
dummy_formula <- function(d) {
  stats::as.formula(paste("~", paste(grep("^f_", names(d), value = TRUE), collapse = " + ")))
}

test_that("global design slice is bit-identical to per-subset model.matrix for numeric formulae", {
  set.seed(20260608)
  n <- 80
  df <- data.frame(X1 = rnorm(n), X2 = rnorm(n), Xp = runif(n, 1, 5))
  rows <- sample(n, 40)
  for (f in list(~X1, ~X1 + X2, ~X1 * X2, ~X1:X2, ~I(X1^2) + X2, ~log(Xp))) {
    a <- model.matrix(f, df[rows, ])
    b <- model.matrix(f, df)[rows, , drop = FALSE]
    expect_identical(as.vector(a), as.vector(b), label = deparse(f))
    expect_identical(colnames(a), colnames(b), label = deparse(f))
  }
})

test_that("slow path matches faster_mode for numeric / transform formulae", {
  set.seed(20260401)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  set.seed(13)
  data$X2 <- rnorm(nrow(data))
  data$Xpos <- runif(nrow(data), 1, 5)

  # bare numeric, interactions, and transform formulae (the latter previously errored)
  for (f in list(~X + X2, ~X * X2, ~X:X2, ~I(X^2) + X2, ~log(Xpos), ~poly(X, 2))) {
    res_slow <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = f, data = data, tname = "period",
             idname = "id", gname = "G", est_method = "dr",
             bstrap = FALSE, faster_mode = FALSE)))
    res_fast <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = f, data = data, tname = "period",
             idname = "id", gname = "G", est_method = "dr",
             bstrap = FALSE, faster_mode = TRUE)))
    expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10, label = paste(deparse(f), "ATT"))
    expect_equal(res_slow$se,  res_fast$se,  tolerance = 1e-10, label = paste(deparse(f), "se"))
  }
})

test_that("RC-path (positional global slice) matches faster_mode", {
  set.seed(20260401)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  set.seed(13)
  data$X2 <- rnorm(nrow(data))
  for (f in list(~X, ~X + X2)) {
    res_slow <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = f, data = data, tname = "period", idname = "id",
             gname = "G", est_method = "dr", bstrap = FALSE, panel = FALSE, faster_mode = FALSE)))
    res_fast <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = f, data = data, tname = "period", idname = "id",
             gname = "G", est_method = "dr", bstrap = FALSE, panel = FALSE, faster_mode = TRUE)))
    expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10, label = paste("RC", deparse(f), "ATT"))
    expect_equal(res_slow$se,  res_fast$se,  tolerance = 1e-10, label = paste("RC", deparse(f), "se"))
  }
})

test_that("a factor covariate equals manually-expanded dummies EXACTLY (dense levels)", {
  set.seed(11)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  data$fac <- factor(sample(c("a", "b", "c"), nrow(data), TRUE))
  data <- add_dummies(data)
  dform <- dummy_formula(data)
  for (est in c("dr", "reg", "ipw")) for (panel in c(TRUE, FALSE)) for (fm in c(TRUE, FALSE)) {
    rf <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = ~fac, data = data, tname = "period", idname = "id",
             gname = "G", est_method = est, panel = panel, bstrap = FALSE, faster_mode = fm)))
    rd <- suppressWarnings(suppressMessages(
      att_gt(yname = "Y", xformla = dform, data = data, tname = "period", idname = "id",
             gname = "G", est_method = est, panel = panel, bstrap = FALSE, faster_mode = fm)))
    lab <- paste(est, "panel", panel, "fm", fm)
    expect_equal(rf$att, rd$att, tolerance = 1e-12, label = paste(lab, "ATT"))
    expect_equal(rf$se,  rd$se,  tolerance = 1e-12, label = paste(lab, "se"))
    expect_equal(as.matrix(rf$inffunc), as.matrix(rd$inffunc), tolerance = 1e-12,
                 label = paste(lab, "inffunc"))
  }
})

test_that("transform formulae that evaluate to NaN drop those rows instead of crashing", {
  # Regression: the non-finite-row check must keep NA/NaN rows visible (model.frame
  # with na.pass) so they are dropped; model.matrix(na.pass) silently drops NaN rows,
  # which used to leave them in the data and crash att_gt with 'subscript out of bounds'.
  set.seed(404)
  sp <- did::reset.sim(time.periods = 4)
  data <- did::build_sim_dataset(sp)
  ids <- unique(data$id)
  xv <- stats::setNames(runif(length(ids), -0.5, 3), ids)   # some <= 0  -> log() is NaN
  data$Xneg <- xv[as.character(data$id)]
  logXneg <- suppressWarnings(log(data$Xneg))
  stopifnot(any(!is.finite(logXneg)))                        # the pathological rows exist

  rs <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~log(Xneg), data = data, tname = "period",
           idname = "id", gname = "G", est_method = "dr", bstrap = FALSE, faster_mode = FALSE)))
  rf <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~log(Xneg), data = data, tname = "period",
           idname = "id", gname = "G", est_method = "dr", bstrap = FALSE, faster_mode = TRUE)))
  expect_true(any(is.finite(rs$att)))
  expect_equal(rs$att, rf$att, tolerance = 1e-10)
  # equals manually pre-cleaning the NaN units
  d2 <- data[is.finite(logXneg), ]
  rref <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~log(Xneg), data = d2, tname = "period",
           idname = "id", gname = "G", est_method = "dr", bstrap = FALSE, faster_mode = FALSE)))
  expect_equal(rs$att, rref$att, tolerance = 1e-12)
})

test_that("a globally-empty factor level is dropped instead of NA-failing every cell", {
  # Regression: factor(levels = c('a','b','c')) where 'c' never occurs in the data
  # (common after users subset their data, since R keeps empty levels) used to emit
  # an all-zero dummy column in EVERY (g,t) cell, NA-failing the whole estimation
  # with a misleading 'Not enough control units' warning. Globally-empty levels are
  # now dropped once up front on both paths, so results equal the droplevels()-ed
  # data exactly. Levels absent only from particular cells keep their documented
  # NA-that-cell-only behavior (see the sparse-factor test below).
  set.seed(22)
  nunit <- 80
  G <- rep(c(0, 2), each = 40)
  funit <- factor(sample(c("a", "b"), nunit, TRUE), levels = c("a", "b", "c"))
  dat <- expand.grid(id = 1:nunit, period = 1:3)
  dat$G <- G[dat$id]
  dat$flev <- funit[dat$id]
  eta <- rnorm(nunit)
  dat$Y <- eta[dat$id] + dat$period + 1 * (dat$G == 2 & dat$period >= 2) +
    rnorm(nrow(dat), sd = 0.1)
  dat2 <- dat
  dat2$flev <- droplevels(dat2$flev)

  run <- function(d, fm, est) suppressWarnings(suppressMessages(
    att_gt(yname = "Y", tname = "period", idname = "id", gname = "G", data = d,
           xformla = ~flev, est_method = est, bstrap = FALSE, faster_mode = fm)))

  for (fm in c(FALSE, TRUE)) for (est in c("dr", "ipw", "reg")) {
    lab <- paste("fm", fm, est)
    res <- run(dat, fm, est)    # empty level 'c' declared
    ref <- run(dat2, fm, est)   # empty level pre-dropped by the user
    expect_false(anyNA(res$att), label = paste(lab, "no NA"))
    expect_equal(res$att, ref$att, tolerance = 1e-12, label = paste(lab, "ATT"))
    expect_equal(res$se,  ref$se,  tolerance = 1e-12, label = paste(lab, "se"))
    expect_equal(as.matrix(res$inffunc), as.matrix(ref$inffunc), tolerance = 1e-12,
                 label = paste(lab, "inffunc"))
  }

  # known values (master's slow path estimated these with the empty level present)
  res_slow <- run(dat, FALSE, "dr")
  res_fast <- run(dat, TRUE, "dr")
  expect_equal(res_slow$att, c(1.0275, 0.9826), tolerance = 1e-3)
  expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10)
  expect_equal(res_slow$se,  res_fast$se,  tolerance = 1e-10)
})

test_that("a sparse factor (level absent from some cells) matches manual dummies, incl. warnings", {
  # 'b' appears only in group 2, so the never-treated control never has it; every cell
  # is rank-deficient and must return NA -- exactly as manual dummies do, with the same
  # warnings, and WITHOUT crashing (the old per-cell droplevels() raised a contrast error).
  set.seed(11)
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
  fac <- rep("a", nrow(data)); fac[data$G == 2] <- sample(c("a", "b"), sum(data$G == 2), TRUE)
  data$fac <- factor(fac)
  data$f_b <- as.numeric(data$fac == "b")

  collect <- function(f, fm) {
    ws <- character(0)
    res <- withCallingHandlers(
      suppressMessages(att_gt(yname = "Y", xformla = f, data = data, tname = "period",
                              idname = "id", gname = "G", est_method = "reg",
                              panel = TRUE, bstrap = FALSE, faster_mode = fm)),
      warning = function(w) { ws[[length(ws) + 1]] <<- conditionMessage(w); invokeRestart("muffleWarning") })
    list(att = res$att, se = res$se, warns = sort(unique(ws)))
  }
  for (fm in c(TRUE, FALSE)) {
    rf <- collect(~fac, fm)
    rd <- collect(~f_b, fm)
    expect_equal(rf$att, rd$att, label = paste("fm", fm, "ATT"))     # both all-NA, same cells
    expect_equal(rf$se,  rd$se,  label = paste("fm", fm, "se"))
    expect_identical(rf$warns, rd$warns, label = paste("fm", fm, "warnings"))
  }
})
