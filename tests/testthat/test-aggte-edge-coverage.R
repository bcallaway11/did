# =============================================================================
# aggte() edge-case coverage: a single treated group exercises the single-column
# matrix-collapse paths in compute.aggte (where group/event-time matrices reduce to
# one column). All four aggregation types must still return a valid AGGTEobj.
# =============================================================================

test_that("aggte handles a single treated group for every aggregation type", {
  set.seed(5)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  g1 <- sort(unique(d$G[d$G > 0]))[1]
  d1 <- d[d$G == 0 | d$G == g1, ]   # one treated group + never-treated controls
  mp <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d1,
    tname = "period", idname = "id", gname = "G", bstrap = FALSE)))
  for (ty in c("simple", "group", "dynamic", "calendar")) {
    agg <- suppressWarnings(aggte(mp, type = ty))
    expect_s3_class(agg, "AGGTEobj")
    expect_true(is.finite(agg$overall.att))
    expect_true(is.finite(agg$overall.se))
  }
})

test_that("aggte dynamic with min_e / max_e windows returns a consistent event-time set", {
  set.seed(6)
  sp <- did::reset.sim(time.periods = 5)
  d <- did::build_sim_dataset(sp)
  mp <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
    tname = "period", idname = "id", gname = "G", bstrap = FALSE)))
  agg <- suppressWarnings(aggte(mp, type = "dynamic", min_e = -1, max_e = 1))
  expect_true(all(agg$egt >= -1 & agg$egt <= 1))
  expect_equal(length(agg$att.egt), length(agg$egt))
})
