# =============================================================================
# panel = FALSE declares genuine repeated cross sections: every observation is
# its own sampling unit, so a supplied idname must be unique across ALL rows.
# Both code paths must reject a repeating idname with the same message, while
# panel = TRUE (including allow_unbalanced_panel = TRUE) is unaffected.
# =============================================================================

# long panel data: id repeats across periods
make_rc_panel <- function(seed = 20260814, n = 240, time.periods = 4) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(n), period = seq_len(time.periods))
  d$G <- ifelse(d$id <= n / 3, 0, ifelse(d$id <= 2 * n / 3, 3, 4))
  d$X <- stats::rnorm(nrow(d))
  d$Y <- 0.5 * d$X + 0.2 * d$period + (d$G > 0 & d$period >= d$G) +
    stats::rnorm(nrow(d), 0, 0.2)
  d
}

rc_msg <- "must be unique when panel = FALSE"

test_that("panel = FALSE rejects an idname that repeats across periods, in both modes", {
  d <- make_rc_panel()
  for (fm in c(TRUE, FALSE)) {
    expect_error(suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
      data = d, tname = "period", idname = "id", gname = "G", panel = FALSE,
      faster_mode = fm, bstrap = FALSE))), rc_msg)
  }
})

test_that("panel = FALSE rejects an idname that repeats within a period, in both modes", {
  d <- make_rc_panel()
  # one row per observation -> a valid repeated cross section, then duplicate a
  # single row so its obsid appears twice in the SAME period
  d$obsid <- seq_len(nrow(d))
  d_dup <- rbind(d, d[1, , drop = FALSE])
  for (fm in c(TRUE, FALSE)) {
    expect_error(suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
      data = d_dup, tname = "period", idname = "obsid", gname = "G", panel = FALSE,
      faster_mode = fm, bstrap = FALSE))), rc_msg)
  }
})

test_that("panel = FALSE runs when the supplied idname is unique, in both modes", {
  d <- make_rc_panel()
  d$obsid <- seq_len(nrow(d))
  for (fm in c(TRUE, FALSE)) {
    res <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
      data = d, tname = "period", idname = "obsid", gname = "G", panel = FALSE,
      faster_mode = fm, bstrap = FALSE)))
    expect_s3_class(res, "MP")
    expect_true(any(is.finite(res$att)))
  }
})

test_that("panel = FALSE runs when idname is omitted even though an id column is present", {
  d <- make_rc_panel()
  for (fm in c(TRUE, FALSE)) {
    res <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
      data = d, tname = "period", gname = "G", panel = FALSE,
      faster_mode = fm, bstrap = FALSE)))
    expect_s3_class(res, "MP")
    expect_true(any(is.finite(res$att)))
  }
})

test_that("panel = TRUE with allow_unbalanced_panel = TRUE is unaffected by the new check", {
  d <- make_rc_panel()
  set.seed(20260815)
  d_ub <- d[-sample(seq_len(nrow(d)), 40), ]   # genuinely unbalanced
  for (fm in c(TRUE, FALSE)) {
    res <- suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X,
      data = d_ub, tname = "period", idname = "id", gname = "G", panel = TRUE,
      allow_unbalanced_panel = TRUE, faster_mode = fm, bstrap = FALSE)))
    expect_s3_class(res, "MP")
    expect_true(any(is.finite(res$att)))
  }
})
