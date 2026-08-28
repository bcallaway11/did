# Balancing the panel drops units whole. att_gt() checked for a never-treated group
# before balancing, so when balancing removed every never-treated unit (e.g. a
# covariate missing in one period for all of them) every ATT(g,t) came back NA --
# silently under faster_mode = TRUE, as "overlap condition violated" under
# faster_mode = FALSE. A treated cohort removed the same way crashed the fast path.

make_panel <- function(cohorts = c(0, 3, 4), n_per = 30L, periods = 1:5, seed = 1) {
  set.seed(seed)
  ids <- seq_len(n_per * length(cohorts))
  d <- data.frame(id = rep(ids, each = length(periods)),
                  time = rep(periods, times = length(ids)),
                  g = rep(rep(cohorts, each = n_per), each = length(periods)))
  d$x1 <- rnorm(nrow(d))
  d$y <- 1 + d$x1 + 0.3 * d$time + (d$g > 0 & d$time >= d$g) + rnorm(nrow(d))
  d
}
# x1 missing in one period for every unit of cohort(s) gsel
kill_x1 <- function(d, gsel, tsel = 2) { d$x1[d$g %in% gsel & d$time %in% tsel] <- NA; d }

fit_gt <- function(d, ...) {
  att_gt(yname = "y", tname = "time", idname = "id", gname = "g", xformla = ~x1,
         data = d, bstrap = FALSE, cband = FALSE, ...)
}
quiet_fit <- function(d, ...) suppressWarnings(suppressMessages(fit_gt(d, ...)))
expect_same_fit <- function(fit, ref) {
  expect_identical(fit$group, ref$group); expect_identical(fit$t, ref$t); expect_identical(fit$n, ref$n)
  expect_equal(fit$att, ref$att, tolerance = 1e-10); expect_equal(fit$se, ref$se, tolerance = 1e-10)
}
RECHECK <- "No never-treated group is available after converting to balanced panel"

test_that("never-treated group removed by balancing is re-checked", {
  d <- kill_x1(make_panel(), gsel = 0)
  d_bal <- d[d$g != 0, ]   # what balancing leaves behind
  for (cg in c("nevertreated", "notyettreated")) {
    msgs <- list()
    for (fm in c(TRUE, FALSE)) {
      w <- capture_warnings(suppressMessages(fit <- fit_gt(d, control_group = cg, faster_mode = fm)))
      msgs[[as.character(fm)]] <- grep(RECHECK, w, fixed = TRUE, value = TRUE)
      expect_length(msgs[[as.character(fm)]], 1L)
      expect_false(any(grepl("overlap", w)))
      expect_false(anyNA(fit$att))
      expect_false(4 %in% fit$group)   # last cohort is the comparison group only
      expect_same_fit(fit, quiet_fit(d_bal, control_group = cg, faster_mode = fm))
    }
    expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
    slow <- quiet_fit(d, control_group = cg, faster_mode = FALSE)
    fast <- quiet_fit(d, control_group = cg, faster_mode = TRUE)
    expect_same_fit(slow, fast)
    expect_length(fast$att, 2L)
  }
})

test_that("re-check honors anticipation", {
  d <- kill_x1(make_panel(), gsel = 0)
  for (cg in c("nevertreated", "notyettreated")) for (fm in c(TRUE, FALSE)) {
    expect_warning(fit <- suppressMessages(fit_gt(d, control_group = cg, anticipation = 1, faster_mode = fm)), RECHECK, fixed = TRUE)
    expect_identical(fit$t, 2)   # cutoff 4 - 1 = 3 leaves periods 1-2
    expect_same_fit(fit, quiet_fit(d[d$g != 0, ], control_group = cg, anticipation = 1, faster_mode = fm))
  }
})

test_that("re-check also covers the comparison cohort coerced from the raw data", {
  # no never-treated group: cohort 5 becomes the comparison group, then balancing removes it
  d <- kill_x1(make_panel(cohorts = c(3, 4, 5)), gsel = 5)
  for (cg in c("nevertreated", "notyettreated")) for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d, control_group = cg, faster_mode = fm)))
    expect_length(grep(RECHECK, w, fixed = TRUE), 1L)
    expect_false(anyNA(fit$att))
    expect_same_fit(fit, quiet_fit(d[d$g != 5, ], control_group = cg, faster_mode = fm))
  }
})

test_that("a treated cohort removed by balancing is dropped with a warning, not an internal error", {
  cases <- list(list(d = make_panel(cohorts = c(0, 3, 4, 5)), kill = 4),
                list(d = make_panel(cohorts = c(0, 3, 4, 5)), kill = 5),
                list(d = make_panel(cohorts = c(3, 4, 5, 6), periods = 1:7), kill = 4))
  for (cs in cases) for (cg in c("nevertreated", "notyettreated")) for (fm in c(TRUE, FALSE)) {
    d <- kill_x1(cs$d, gsel = cs$kill)
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d, control_group = cg, faster_mode = fm)))
    expect_length(grep(paste0("Dropped cohort(s) ", cs$kill, " while converting to balanced panel"), w, fixed = TRUE), 1L)
    expect_false(any(grepl("Internal error", w)))
    expect_false(cs$kill %in% fit$group)
    expect_false(anyNA(fit$att))
    expect_same_fit(fit, quiet_fit(cs$d[cs$d$g != cs$kill, ], control_group = cg, faster_mode = fm))
  }
})

test_that("no re-check when never-treated units survive balancing or the panel is left unbalanced", {
  d <- make_panel()
  clean <- quiet_fit(d)
  d_half <- d; d_half$x1[d_half$g == 0 & d_half$time == 2 & d_half$id <= 15] <- NA
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d_half, faster_mode = fm)))
    expect_false(any(grepl("never-treated|Dropped cohort", w)))
    expect_identical(fit$group, clean$group); expect_identical(fit$t, clean$t); expect_identical(fit$n, 75L)
    w <- capture_warnings(suppressMessages(fit <- fit_gt(kill_x1(d, 0), allow_unbalanced_panel = TRUE, faster_mode = fm)))
    expect_false(any(grepl("never-treated|Dropped cohort", w)))
    late4 <- fit$att[fit$group == 4 & fit$t %in% c(4, 5)]
    expect_length(late4, 2L); expect_false(anyNA(late4))
  }
})

test_that("'notyettreated' without a never-treated group in the raw data stays silent", {
  d <- make_panel(cohorts = c(3, 4))
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d, control_group = "notyettreated", faster_mode = fm)))
    expect_false(any(grepl("never-treated", w)))
    expect_false(4 %in% fit$group)
  }
})

test_that("nothing left to estimate gives an error that names the cause, identically in both modes", {
  cases <- list(
    list(d = kill_x1(make_panel(cohorts = c(0, 3)), gsel = 0), extra = list(),
         msg = "the last treated cohort (first treated at 3) is used only as a comparison group and no other treated cohort remains"),
    list(d = make_panel(cohorts = 1), extra = list(), msg = "No periods before the last treated cohort's treatment date"),
    list(d = make_panel(cohorts = 1, periods = 2010:2014), extra = list(), msg = "same scale as 'tname'"),
    list(d = kill_x1(make_panel(cohorts = c(0, 3, 5), periods = c(1, 4, 7)), gsel = 0, tsel = 4), extra = list(anticipation = 1),
         msg = "Only one time period remains"),
    list(d = kill_x1(make_panel(), gsel = c(3, 4)), extra = list(), msg = "No valid groups")
  )
  for (cs in cases) for (cg in c("nevertreated", "notyettreated")) {
    e <- lapply(c(TRUE, FALSE), function(fm)
      tryCatch(do.call(quiet_fit, c(list(cs$d, control_group = cg, faster_mode = fm), cs$extra)), error = function(e) conditionMessage(e)))
    expect_type(e[[1]], "character")
    expect_match(e[[1]], cs$msg, fixed = TRUE)
    expect_identical(e[[1]], e[[2]])
  }
})

test_that("a sample emptied by the missing-data drop stops with the missing-data message", {
  d <- make_panel(); d$y <- NA_real_
  for (fm in c(TRUE, FALSE)) {
    w <- character()
    e <- tryCatch(withCallingHandlers(suppressMessages(fit_gt(d, faster_mode = fm)),
                                      warning = function(w0) { w <<- c(w, conditionMessage(w0)); invokeRestart("muffleWarning") }),
                  error = function(e) conditionMessage(e))
    expect_match(e, "All observations were dropped due to missing data", fixed = TRUE)
    expect_false(any(grepl("-Inf|no non-missing arguments", c(e, w))))
  }
})

test_that("a NaN 2x2 estimate is announced (both modes)", {
  d <- make_panel(); d$w0 <- ifelse(d$g == 4, 0, 1)   # cohort 4 has no effective units
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d, weightsname = "w0", faster_mode = fm)))
    na_cells <- which(is.na(fit$att))
    expect_identical(fit$group[na_cells], rep(4, 4))
    for (i in na_cells) expect_true(any(grepl(paste0("(g,t) = (4,", fit$t[i], ") is NaN"), w, fixed = TRUE)))
  }
})

test_that("faster_mode = TRUE names a cell with no control units instead of returning NA silently", {
  set.seed(3)
  d <- expand.grid(rep = 1:40, time = 1:5, g = c(0, 3))
  d$y <- rnorm(nrow(d)) + 0.2 * d$time + (d$g > 0 & d$time >= d$g)
  d <- d[!(d$g == 0 & d$time %in% 2:4), ]   # never-treated observed only in periods 1 and 5
  d$id <- seq_len(nrow(d))
  w <- capture_warnings(suppressMessages(
    fit <- att_gt(yname = "y", tname = "time", idname = "id", gname = "g", data = d, panel = FALSE, bstrap = FALSE, cband = FALSE)
  ))
  expect_true(all(is.na(fit$att)))
  for (tt in c(3, 4)) expect_true(any(grepl(paste0("No treated or control units available for group 3 in time period ", tt), w, fixed = TRUE)))
})

test_that("aggte() explains an object with no post-treatment cells", {
  d <- make_panel(cohorts = c(4, 6), periods = 1:6)   # anticipation 2 drops periods >= 4
  for (fm in c(TRUE, FALSE)) {
    fit <- quiet_fit(d, anticipation = 2, faster_mode = fm)
    expect_true(all(fit$t < fit$group))
    for (ty in c("simple", "dynamic", "group", "calendar"))
      expect_error(suppressWarnings(aggte(fit, type = ty)), "No post-treatment ATT(g,t) estimates are available", fixed = TRUE)
  }
})
