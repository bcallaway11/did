.balance_support_data <- function(periods = 1:5, groups = c(0, 2, 3)) {
  set.seed(711)
  n <- 180L
  d <- expand.grid(id = seq_len(n), t = periods)
  d$g <- rep(rep(groups, each = n / length(groups)), length(periods))
  d$y <- 0.1 * d$t + ifelse(d$g > 0 & d$t >= d$g, d$g, 0) +
    stats::rnorm(nrow(d), sd = 0.1)
  d
}

.balance_support_fit <- function(d, faster_mode) {
  suppressWarnings(suppressMessages(att_gt(
    "y", "t", "id", "g", data = d, faster_mode = faster_mode,
    est_method = "reg", allow_unbalanced_panel = TRUE,
    bstrap = FALSE, cband = FALSE
  )))
}

test_that("balance_e retains the same cohorts after missing cells are removed", {
  d <- .balance_support_data()
  d <- d[!(d$g == 2 & d$t == 3), ]
  for (fm in c(TRUE, FALSE)) {
    mp <- .balance_support_fit(d, fm)
    expect_true(is.na(mp$att[mp$group == 2 & mp$t == 3]))
    a <- aggte(mp, "dynamic", balance_e = 1, min_e = 0, na.rm = TRUE)
    k <- vapply(0:1, function(e) which(mp$group == 3 & mp$t == 3 + e), integer(1))
    expect_equal(a$egt, 0:1)
    expect_equal(unname(a$att.egt), mp$att[k], tolerance = 1e-10)
    expect_equal(unname(a$se.egt), mp$se[k], tolerance = 1e-10)
    expect_equal(unname(a$inf.function$dynamic.inf.func.e),
                 unname(as.matrix(mp$inffunc[, k])), tolerance = 1e-10)

    # Restricting the displayed window does not redefine which cohorts have
    # complete support through balance_e: g=2 is still excluded at e=0.
    a0 <- aggte(mp, "dynamic", balance_e = 1, min_e = 0, max_e = 0,
                na.rm = TRUE)
    expect_equal(a0$att.egt, mp$att[k[1]], tolerance = 1e-10)
  }
})

test_that("balance_e requires common event-time support on irregular grids", {
  d <- .balance_support_data(c(1, 2, 4, 5), c(0, 2, 4))
  for (fm in c(TRUE, FALSE)) {
    mp <- .balance_support_fit(d, fm)
    a <- aggte(mp, "dynamic", balance_e = 1, min_e = 0)
    # g=2 reaches the horizon but has no e=1 cell. Only g=4 contributes.
    k <- vapply(0:1, function(e) which(mp$group == 4 & mp$t == 4 + e), integer(1))
    expect_equal(a$egt, 0:1)
    expect_equal(unname(a$att.egt), mp$att[k], tolerance = 1e-10)
  }
})

test_that("balance_e uses observed support without inventing absent event times", {
  d <- .balance_support_data(c(1, 3, 5), c(0, 3, 5))
  for (fm in c(TRUE, FALSE)) {
    mp <- .balance_support_fit(d, fm)
    a <- aggte(mp, "dynamic", balance_e = 2, min_e = 0)
    k <- vapply(c(0, 2), function(e) which(mp$group == 3 & mp$t == 3 + e), integer(1))
    expect_equal(a$egt, c(0, 2))
    expect_equal(unname(a$att.egt), mp$att[k], tolerance = 1e-10)
  }
})

test_that("balance_e errors when no cohort has complete observed support", {
  d <- .balance_support_data()
  d <- d[!(d$g == 2 & d$t == 3) & !(d$g == 3 & d$t == 3), ]
  for (fm in c(TRUE, FALSE)) {
    mp <- .balance_support_fit(d, fm)
    expect_error(aggte(mp, "dynamic", balance_e = 1, min_e = 0, na.rm = TRUE),
                 "No groups have non-missing ATT")
    expect_error(aggte(mp, "dynamic", balance_e = 10, na.rm = TRUE),
                 "No groups have non-missing ATT")
  }
})

test_that("an event time with all missing estimates still has to be balanced", {
  d <- .balance_support_data()
  d <- d[!(d$g == 2 & d$t == 3) & !(d$g == 3 & d$t == 4), ]
  for (fm in c(TRUE, FALSE)) {
    mp <- .balance_support_fit(d, fm)
    k <- which(mp$t - mp$group == 1)
    expect_equal(length(k), 2)
    expect_true(all(is.na(mp$att[k])))
    expect_error(aggte(mp, "dynamic", balance_e = 1, min_e = 0, na.rm = TRUE),
                 "No groups have non-missing ATT")
  }
})

test_that("balance_e preserves the complete regularly spaced aggregation", {
  d <- .balance_support_data()
  for (fm in c(TRUE, FALSE)) {
    mp <- .balance_support_fit(d, fm)
    a <- aggte(mp, "dynamic", balance_e = 1, min_e = 0)
    for (e in 0:1) {
      k <- which(mp$t - mp$group == e)
      # The two treated cohorts have equal population shares.
      expect_equal(a$att.egt[a$egt == e], mean(mp$att[k]), tolerance = 1e-10)
    }
  }
})
