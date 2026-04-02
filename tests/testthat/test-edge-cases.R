# =============================================================================
# Tests for edge cases and boundary conditions
# =============================================================================

test_that("single treated group produces valid att_gt", {
  set.seed(20260401)
  data_sg <- data.frame(
    id = rep(1:200, each = 5),
    period = rep(1:5, 200),
    G = rep(c(rep(4, 50), rep(0, 150)), each = 5),
    Y = rnorm(1000)
  )
  data_sg$Y[data_sg$G == 4 & data_sg$period >= 4] <-
    data_sg$Y[data_sg$G == 4 & data_sg$period >= 4] + 1

  result <- att_gt(yname = "Y", data = data_sg, tname = "period", idname = "id",
                   gname = "G", bstrap = FALSE)
  expect_s3_class(result, "MP")
  expect_true(any(!is.na(result$att)))
})

test_that("single treated group works with all 4 aggte types", {
  set.seed(20260401)
  data_sg <- data.frame(
    id = rep(1:200, each = 5),
    period = rep(1:5, 200),
    G = rep(c(rep(4, 50), rep(0, 150)), each = 5),
    Y = rnorm(1000)
  )
  data_sg$Y[data_sg$G == 4 & data_sg$period >= 4] <-
    data_sg$Y[data_sg$G == 4 & data_sg$period >= 4] + 1

  mp <- att_gt(yname = "Y", data = data_sg, tname = "period", idname = "id",
               gname = "G", bstrap = FALSE)

  for (tp in c("simple", "dynamic", "group", "calendar")) {
    agg <- suppressWarnings(aggte(mp, type = tp))
    expect_s3_class(agg, "AGGTEobj")
    expect_false(is.na(agg$overall.att), label = tp)
  }
})

test_that("two-period data with universal base period", {
  set.seed(20260401)
  data_2p <- data.frame(
    id = rep(1:200, each = 2),
    period = rep(1:2, 200),
    G = rep(c(rep(2, 50), rep(0, 150)), each = 2),
    Y = rnorm(400)
  )
  data_2p$Y[data_2p$G == 2 & data_2p$period == 2] <-
    data_2p$Y[data_2p$G == 2 & data_2p$period == 2] + 1

  result <- suppressWarnings(
    att_gt(yname = "Y", data = data_2p, tname = "period", idname = "id",
           gname = "G", base_period = "universal", bstrap = FALSE)
  )
  expect_s3_class(result, "MP")
  expect_true(any(!is.na(result$att)))
})

test_that("data with no never-treated group works with notyettreated", {
  set.seed(20260401)
  # All units are eventually treated (groups 3 and 5)
  data_nnt <- data.frame(
    id = rep(1:200, each = 6),
    period = rep(1:6, 200),
    G = rep(c(rep(3, 100), rep(5, 100)), each = 6),
    Y = rnorm(1200)
  )
  result <- suppressWarnings(
    att_gt(yname = "Y", data = data_nnt, tname = "period", idname = "id",
           gname = "G", control_group = "notyettreated", bstrap = FALSE)
  )
  expect_s3_class(result, "MP")
})

test_that("data with no never-treated group warns with nevertreated", {
  set.seed(20260401)
  data_nnt <- data.frame(
    id = rep(1:200, each = 6),
    period = rep(1:6, 200),
    G = rep(c(rep(3, 100), rep(5, 100)), each = 6),
    Y = rnorm(1200)
  )
  expect_warning(
    att_gt(yname = "Y", data = data_nnt, tname = "period", idname = "id",
           gname = "G", control_group = "nevertreated", bstrap = FALSE),
    "No never-treated group|never-treated"
  )
})

test_that("groups treated in first period are dropped", {
  set.seed(20260401)
  sp <- did::reset.sim()
  data_fp <- did::build_sim_dataset(sp)
  # Add units treated in the very first period
  first_per <- min(data_fp$period)
  extra <- data_fp[data_fp$G == sort(unique(data_fp$G[data_fp$G > 0]))[1], ]
  extra$G <- first_per
  extra$id <- extra$id + max(data_fp$id)
  data_fp <- rbind(data_fp, extra)

  expect_warning(
    att_gt(yname = "Y", data = data_fp, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE),
    "already treated|Dropped"
  )
})

test_that("non-consecutive time periods work", {
  set.seed(20260401)
  data_nc <- data.frame(
    id = rep(1:200, each = 4),
    period = rep(c(2000, 2003, 2007, 2010), 200),
    G = rep(c(rep(2007, 50), rep(0, 150)), each = 4),
    Y = rnorm(800)
  )
  data_nc$Y[data_nc$G == 2007 & data_nc$period >= 2007] <-
    data_nc$Y[data_nc$G == 2007 & data_nc$period >= 2007] + 1

  result <- suppressWarnings(
    att_gt(yname = "Y", data = data_nc, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE)
  )
  expect_s3_class(result, "MP")
  expect_true(any(!is.na(result$att)))

  agg <- suppressWarnings(aggte(result, type = "dynamic"))
  expect_s3_class(agg, "AGGTEobj")
})

test_that("non-consecutive group values work", {
  set.seed(20260401)
  data_ng <- data.frame(
    id = rep(1:200, each = 5),
    period = rep(1:5, 200),
    G = rep(c(rep(3, 50), rep(5, 50), rep(0, 100)), each = 5),
    Y = rnorm(1000)
  )
  result <- att_gt(yname = "Y", data = data_ng, tname = "period", idname = "id",
                   gname = "G", bstrap = FALSE)
  expect_s3_class(result, "MP")
  expect_true(length(unique(result$group)) >= 2)
})

test_that("allow_unbalanced_panel=TRUE with balanced data proceeds normally", {
  set.seed(20260401)
  sp <- did::reset.sim()
  data_bal <- did::build_sim_dataset(sp)

  result <- suppressMessages(
    att_gt(yname = "Y", data = data_bal, tname = "period", idname = "id",
           gname = "G", allow_unbalanced_panel = TRUE, bstrap = FALSE)
  )
  expect_s3_class(result, "MP")
  expect_true(any(!is.na(result$att)))
})

test_that("allow_unbalanced_panel=TRUE with truly unbalanced data", {
  set.seed(20260401)
  sp <- did::reset.sim()
  data_ub <- did::build_sim_dataset(sp)
  # Drop a few rows to make it unbalanced
  data_ub <- data_ub[-c(1, 5, 10), ]

  result <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", data = data_ub, tname = "period", idname = "id",
           gname = "G", allow_unbalanced_panel = TRUE, bstrap = FALSE)
  ))
  expect_s3_class(result, "MP")
  expect_true(any(!is.na(result$att)))
})

test_that("single post-treatment period produces valid results", {
  set.seed(20260401)
  data_spt <- data.frame(
    id = rep(1:200, each = 3),
    period = rep(1:3, 200),
    G = rep(c(rep(3, 50), rep(0, 150)), each = 3),
    Y = rnorm(600)
  )
  data_spt$Y[data_spt$G == 3 & data_spt$period == 3] <-
    data_spt$Y[data_spt$G == 3 & data_spt$period == 3] + 1

  result <- att_gt(yname = "Y", data = data_spt, tname = "period", idname = "id",
                   gname = "G", bstrap = FALSE)
  expect_s3_class(result, "MP")
  # Should have at least one post-treatment ATT
  post_atts <- result$att[result$group <= result$t]
  expect_true(any(!is.na(post_atts)))
})
