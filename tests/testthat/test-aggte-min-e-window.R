# =============================================================================
# min_e now restricts the event-time window of the "simple" and "group"
# aggregations, symmetrically with the long-standing max_e. The defaults
# (min_e = -Inf, max_e = Inf) are unchanged, negative min_e is a no-op for these
# two types (only post-treatment cells enter), and a window that empties the
# selection gives a clear error. balance_e, which only applies to "dynamic",
# now warns for "simple" and "group".
# =============================================================================

# unit-level group shares (the pg weights used by the aggregations; there are no
# sampling weights in mpdta, so weights.ind = 1)
mpdta_group_share <- function() {
  data(mpdta, package = "did")
  units <- unique(mpdta[, c("countyreal", "first.treat")])
  tab <- table(units$first.treat) / nrow(units)
  tab
}

mpdta_mp <- function(fm) {
  data(mpdta, package = "did")
  suppressWarnings(suppressMessages(att_gt(yname = "lemp", tname = "year",
    idname = "countyreal", gname = "first.treat", data = mpdta,
    faster_mode = fm, bstrap = FALSE, cband = FALSE)))
}

test_that("aggte(type='simple', min_e=) equals the hand-computed restricted weighted mean", {
  share <- mpdta_group_share()
  for (fm in c(TRUE, FALSE)) {
    mp <- mpdta_mp(fm)
    for (me in c(1, 2)) {
      keep <- which(mp$group <= mp$t & (mp$t - mp$group) >= me)
      pgv <- as.numeric(share[as.character(mp$group[keep])])
      manual <- sum(mp$att[keep] * pgv) / sum(pgv)
      agg <- suppressWarnings(suppressMessages(aggte(mp, type = "simple", min_e = me)))
      expect_equal(agg$overall.att, manual, tolerance = 1e-10)
    }
    # the restriction actually bites
    agg_all <- suppressWarnings(suppressMessages(aggte(mp, type = "simple")))
    agg_min <- suppressWarnings(suppressMessages(aggte(mp, type = "simple", min_e = 1)))
    expect_false(isTRUE(all.equal(agg_all$overall.att, agg_min$overall.att)))
  }
})

test_that("the default min_e = -Inf leaves simple and group aggregations unchanged", {
  for (fm in c(TRUE, FALSE)) {
    mp <- mpdta_mp(fm)
    for (ty in c("simple", "group")) {
      a_default <- suppressWarnings(suppressMessages(aggte(mp, type = ty)))
      a_minf <- suppressWarnings(suppressMessages(aggte(mp, type = ty, min_e = -Inf)))
      # a negative min_e is also a no-op: only post-treatment cells (e >= 0) enter
      a_neg <- suppressWarnings(suppressMessages(aggte(mp, type = ty, min_e = -5)))
      expect_equal(a_minf$overall.att, a_default$overall.att, tolerance = 1e-10)
      expect_equal(a_neg$overall.att, a_default$overall.att, tolerance = 1e-10)
      expect_equal(a_minf$overall.se, a_default$overall.se, tolerance = 1e-10)
      expect_equal(a_neg$overall.se, a_default$overall.se, tolerance = 1e-10)
      if (ty == "group") {
        expect_equal(a_minf$att.egt, a_default$att.egt, tolerance = 1e-10)
        expect_equal(a_neg$att.egt, a_default$att.egt, tolerance = 1e-10)
      }
    }
  }
})

test_that("aggte(type='group', min_e=) restricts each group's average to the window", {
  for (fm in c(TRUE, FALSE)) {
    mp <- mpdta_mp(fm)
    # mpdta's last cohort (2007) has only its e = 0 cell, so min_e = 1 leaves it
    # with an empty window; na.rm = TRUE drops it, as it does for all-NA groups
    agg <- suppressWarnings(suppressMessages(aggte(mp, type = "group", min_e = 1,
                                                   na.rm = TRUE)))
    for (g in agg$egt) {
      whichg <- which(mp$group == g & mp$group <= mp$t & (mp$t - mp$group) >= 1)
      expect_equal(agg$att.egt[agg$egt == g], mean(mp$att[whichg]), tolerance = 1e-10)
    }
    # hand-check the earliest group explicitly (it has the most post-treatment cells)
    g1 <- min(agg$egt)
    manual_g1 <- mean(mp$att[mp$group == g1 & mp$t >= g1 + 1])
    expect_equal(agg$att.egt[agg$egt == g1], manual_g1, tolerance = 1e-10)
    # and it differs from the unrestricted group average
    agg_all <- suppressWarnings(suppressMessages(aggte(mp, type = "group")))
    expect_false(isTRUE(all.equal(agg_all$att.egt[agg_all$egt == g1],
                                  agg$att.egt[agg$egt == g1])))
  }
})

test_that("a min_e window that excludes every post-treatment cell errors clearly", {
  for (fm in c(TRUE, FALSE)) {
    mp <- mpdta_mp(fm)
    expect_error(suppressWarnings(suppressMessages(aggte(mp, type = "simple", min_e = 100))),
                 "No group-time average treatment effects fall within the requested window")
    # the group aggregation names the groups whose window came out empty
    expect_error(suppressWarnings(suppressMessages(aggte(mp, type = "group", min_e = 100))),
                 "No group-time average treatment effects fall within the requested window")
    expect_error(suppressWarnings(suppressMessages(aggte(mp, type = "group", min_e = 100,
                                                         na.rm = TRUE))),
                 "min_e")
  }
})

test_that("a window that empties only some groups errors clearly, naming them", {
  for (fm in c(TRUE, FALSE)) {
    mp <- mpdta_mp(fm)
    # only the 2007 cohort loses every cell at min_e = 1
    expect_error(suppressWarnings(suppressMessages(aggte(mp, type = "group", min_e = 1))),
                 "for group\\(s\\) 2007")
  }
})

test_that("min_e/max_e for simple and group are in original time units on irregularly spaced panels", {
  # periods 1, 2, 4: for the g = 2 cohort the cell (g = 2, t = 4) sits at event
  # time e = 2 in original units but only one recoded index step past treatment.
  # The window must follow original units, matching type = "dynamic" and the
  # documented meaning of e = t - g.
  set.seed(20260814)
  n <- 400
  d <- expand.grid(id = seq_len(n), period = c(1, 2, 4))
  d$G <- ifelse(d$id <= n / 2, 0, 2)
  d$X <- stats::rnorm(nrow(d))
  d$y <- 0.1 * d$X + 0.05 * d$period +
    (d$G == 2 & d$period >= 2) * (d$period - 1) +   # effect grows with exposure
    stats::rnorm(nrow(d), 0, 0.1)
  for (fm in c(TRUE, FALSE)) {
    mp <- suppressWarnings(suppressMessages(att_gt(yname = "y", tname = "period",
      idname = "id", gname = "G", xformla = ~X, data = d, faster_mode = fm,
      bstrap = FALSE, cband = FALSE)))
    k0 <- which(mp$group == 2 & mp$t == 2)   # e = 0
    k2 <- which(mp$group == 2 & mp$t == 4)   # e = 2 (one index step)
    expect_false(isTRUE(all.equal(mp$att[k0], mp$att[k2])))
    # max_e = 1 must keep only the e = 0 cell (index-step counting would also
    # keep (2, 4), which is one recoded step past treatment)
    a_max <- suppressWarnings(suppressMessages(aggte(mp, type = "simple", max_e = 1)))
    expect_equal(a_max$overall.att, mp$att[k0], tolerance = 1e-10)
    g_max <- suppressWarnings(suppressMessages(aggte(mp, type = "group", max_e = 1)))
    expect_equal(g_max$att.egt[g_max$egt == 2], mp$att[k0], tolerance = 1e-10)
    # min_e = 2 must keep the e = 2 cell (index-step counting would find no
    # cell two steps past treatment and error)
    a_min <- suppressWarnings(suppressMessages(aggte(mp, type = "simple", min_e = 2)))
    expect_equal(a_min$overall.att, mp$att[k2], tolerance = 1e-10)
    # consistent with the dynamic aggregation's original-unit event times
    a_dyn <- suppressWarnings(suppressMessages(aggte(mp, type = "dynamic", max_e = 1)))
    expect_true(all(a_dyn$egt <= 1))
    expect_equal(a_dyn$att.egt[a_dyn$egt == 0], mp$att[k0], tolerance = 1e-10)
  }
})

test_that("the min_e/max_e window boundary is inclusive: an existing e = max_e cell is kept", {
  # periods 1, 2, 3, 5: the g = 2 cohort has cells at e = 0, 1, and 3, so
  # max_e = 1 must keep BOTH e = 0 and e = 1 (the window is e <= max_e,
  # inclusive) and drop e = 3
  set.seed(20260814)
  n <- 400
  d <- expand.grid(id = seq_len(n), period = c(1, 2, 3, 5))
  d$G <- ifelse(d$id <= n / 2, 0, 2)
  d$X <- stats::rnorm(nrow(d))
  d$y <- 0.1 * d$X + 0.05 * d$period +
    (d$G == 2 & d$period >= 2) * (d$period - 1) +
    stats::rnorm(nrow(d), 0, 0.1)
  for (fm in c(TRUE, FALSE)) {
    mp <- suppressWarnings(suppressMessages(att_gt(yname = "y", tname = "period",
      idname = "id", gname = "G", xformla = ~X, data = d, faster_mode = fm,
      bstrap = FALSE, cband = FALSE)))
    k0 <- which(mp$group == 2 & mp$t == 2)   # e = 0
    k1 <- which(mp$group == 2 & mp$t == 3)   # e = 1
    a <- suppressWarnings(suppressMessages(aggte(mp, type = "simple", max_e = 1)))
    # single cohort -> equal pg weights -> the simple average of the two cells
    expect_equal(a$overall.att, mean(mp$att[c(k0, k1)]), tolerance = 1e-10)
    g <- suppressWarnings(suppressMessages(aggte(mp, type = "group", max_e = 1)))
    expect_equal(g$att.egt[g$egt == 2], mean(mp$att[c(k0, k1)]), tolerance = 1e-10)
  }
})

test_that("balance_e warns that it is ignored for the simple and group aggregations", {
  mp <- mpdta_mp(TRUE)
  for (ty in c("simple", "group")) {
    ws <- testthat::capture_warnings(suppressMessages(aggte(mp, type = ty, balance_e = 1)))
    expect_true(any(grepl("`balance_e` is ignored for type", ws)))
    expect_true(any(grepl(paste0("type = \"", ty, "\""), ws)))
    # no spurious warning when balance_e is not set
    ws0 <- testthat::capture_warnings(suppressMessages(aggte(mp, type = ty)))
    expect_false(any(grepl("`balance_e` is ignored for type", ws0)))
  }
  # the dynamic aggregation still honors balance_e without that warning
  ws_dyn <- testthat::capture_warnings(suppressMessages(aggte(mp, type = "dynamic", balance_e = 1)))
  expect_false(any(grepl("`balance_e` is ignored for type", ws_dyn)))
  a_bal <- suppressWarnings(suppressMessages(aggte(mp, type = "dynamic", balance_e = 1)))
  a_unb <- suppressWarnings(suppressMessages(aggte(mp, type = "dynamic")))
  expect_false(isTRUE(all.equal(a_bal$overall.att, a_unb$overall.att)))
})
