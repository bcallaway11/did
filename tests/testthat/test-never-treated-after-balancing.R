# Regression tests: a never-treated group (or a whole cohort) that the balanced-panel
# coercion removes must be re-checked. att_gt() decided whether a never-treated group
# exists BEFORE coercing the data to a balanced panel; balancing drops units whole, so
# when it removed every never-treated unit (one covariate missing in one period for
# all of them is enough) nothing noticed, no control units remained, and every
# ATT(g,t) came back NA -- silently under faster_mode = TRUE and as a misleading
# "overlap condition violated" under faster_mode = FALSE. The fix re-runs the existing
# "no never-treated group" fallback on the balanced sample, in both code paths, and
# reconciles the list of treated cohorts with the balanced sample (a treated cohort
# removed whole used to crash the fast path with an "Internal error").

make_panel <- function(cohorts = c(0, 3, 4), n_per = 30L, periods = 1:5, seed = 1) {
  set.seed(seed)
  ids <- seq_len(n_per * length(cohorts))
  gvec <- rep(cohorts, each = n_per)
  d <- data.frame(id = rep(ids, each = length(periods)),
                  time = rep(periods, times = length(ids)),
                  g = rep(gvec, each = length(periods)))
  d$x1 <- rnorm(nrow(d))
  d$y <- 1 + d$x1 + 0.3 * d$time + (d$g > 0 & d$time >= d$g) + rnorm(nrow(d))
  d
}

# one covariate missing in ONE period for every unit of cohort(s) `gsel`
kill_x1 <- function(d, gsel, tsel = 2) {
  d$x1[d$g %in% gsel & d$time %in% tsel] <- NA
  d
}

fit_gt <- function(d, ...) {
  att_gt(yname = "y", tname = "time", idname = "id", gname = "g", xformla = ~x1,
         data = d, bstrap = FALSE, cband = FALSE, ...)
}

quiet_fit <- function(d, ...) suppressWarnings(suppressMessages(fit_gt(d, ...)))

expect_same_fit <- function(fit, ref, tol = 1e-10, info = NULL) {
  expect_identical(fit$group, ref$group, info = info)
  expect_identical(fit$t, ref$t, info = info)
  expect_identical(fit$n, ref$n, info = info)
  expect_equal(fit$att, ref$att, tolerance = tol, info = info)
  expect_equal(fit$se, ref$se, tolerance = tol, info = info)
}

RECHECK_MSG <- "No never-treated group is available after converting to a balanced panel"

test_that("never-treated group removed by balancing is re-checked (nevertreated, both modes)", {
  d <- kill_x1(make_panel(), gsel = 0)
  d_bal <- d[d$g != 0, ]   # exactly what balancing leaves behind
  msgs <- list()
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d, control_group = "nevertreated", faster_mode = fm)
    ))
    m <- grep(RECHECK_MSG, w, fixed = TRUE, value = TRUE)
    expect_length(m, 1L)
    # the warning says what was lost and how to keep it
    expect_match(m, "all never-treated units (30 remained after dropping rows with missing data) were removed", fixed = TRUE)
    expect_match(m, "allow_unbalanced_panel = TRUE", fixed = TRUE)
    expect_match(m, "first treated at 4", fixed = TRUE)
    msgs[[as.character(fm)]] <- m
    # the slow path used to report every cell as an overlap violation
    expect_false(any(grepl("overlap", w)), info = paste("faster_mode =", fm))
    expect_false(anyNA(fit$att), info = paste("faster_mode =", fm))
    expect_identical(fit$n, 60L)

    # the raw data must now reach the rule the package already applies when handed
    # the balanced sample directly
    ref <- quiet_fit(d_bal, control_group = "nevertreated", faster_mode = fm)
    expect_same_fit(fit, ref, info = paste("faster_mode =", fm))
  }
  # byte-identical across modes
  expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
})

test_that("never-treated group removed by balancing is re-checked (notyettreated, both modes)", {
  d <- kill_x1(make_panel(), gsel = 0)
  d_bal <- d[d$g != 0, ]
  msgs <- list()
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d, control_group = "notyettreated", faster_mode = fm)
    ))
    m <- grep(RECHECK_MSG, w, fixed = TRUE, value = TRUE)
    expect_length(m, 1L)
    expect_match(m, "used only as a not-yet-treated comparison group", fixed = TRUE)
    msgs[[as.character(fm)]] <- m
    expect_false(any(grepl("overlap", w)), info = paste("faster_mode =", fm))
    expect_false(anyNA(fit$att), info = paste("faster_mode =", fm))

    ref <- quiet_fit(d_bal, control_group = "notyettreated", faster_mode = fm)
    expect_same_fit(fit, ref, info = paste("faster_mode =", fm))
    # same rule as without a never-treated group in the raw data: the last cohort
    # is a comparison group only and gets no ATT(g,t) of its own
    expect_false(4 %in% fit$group)
  }
  expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
})

test_that("fast and slow modes agree on the re-checked sample", {
  d <- kill_x1(make_panel(), gsel = 0)
  for (cg in c("nevertreated", "notyettreated")) {
    slow <- quiet_fit(d, control_group = cg, faster_mode = FALSE)
    fast <- quiet_fit(d, control_group = cg, faster_mode = TRUE)
    expect_same_fit(slow, fast, info = cg)
    # not vacuous: on the unfixed code both paths agreed on an all-NA grid
    expect_false(anyNA(fast$att), info = cg)
    expect_length(fast$att, 2L)
  }
})

test_that("re-check honors anticipation (both control groups)", {
  d <- kill_x1(make_panel(), gsel = 0)
  d_bal <- d[d$g != 0, ]
  for (cg in c("nevertreated", "notyettreated")) for (fm in c(TRUE, FALSE)) {
    expect_warning(
      fit <- suppressMessages(fit_gt(d, control_group = cg, anticipation = 1, faster_mode = fm)),
      RECHECK_MSG, fixed = TRUE
    )
    ref <- quiet_fit(d_bal, control_group = cg, anticipation = 1, faster_mode = fm)
    expect_false(anyNA(fit$att))
    # cutoff = 4 - 1 = 3: only periods 1-2 remain, cohort 3's single pre-treatment cell
    expect_identical(fit$group, 3)
    expect_identical(fit$t, 2)
    expect_same_fit(fit, ref, info = paste(cg, "faster_mode =", fm))
  }
})

test_that("re-check also covers the comparison cohort created by the raw-data fallback", {
  # No never-treated group in the raw data: the fallback turns the latest cohort (5)
  # into the comparison group. Balancing then removes that whole cohort, so the
  # rule must be applied again -- to cohort 4 -- exactly as when the package is
  # handed cohorts 3 and 4 directly.
  d <- kill_x1(make_panel(cohorts = c(3, 4, 5)), gsel = 5)
  d_bal <- d[d$g != 5, ]
  msgs <- list()
  for (cg in c("nevertreated", "notyettreated")) for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d, control_group = cg, faster_mode = fm)
    ))
    # the raw-data fallback warns under "nevertreated" only (routine design under
    # "notyettreated", kept silent on purpose)
    expect_identical(sum(grepl("No never-treated group is available. ", w, fixed = TRUE)),
                     if (cg == "nevertreated") 1L else 0L,
                     info = paste(cg, "faster_mode =", fm))
    m <- grep(RECHECK_MSG, w, fixed = TRUE, value = TRUE)
    expect_length(m, 1L)
    expect_match(m, "every unit of the comparison cohort (first treated at 5) was removed", fixed = TRUE)
    msgs[[paste(cg, fm)]] <- m
    expect_false(anyNA(fit$att))
    ref <- quiet_fit(d_bal, control_group = cg, faster_mode = fm)
    expect_same_fit(fit, ref, info = paste(cg, "faster_mode =", fm))
  }
  for (cg in c("nevertreated", "notyettreated")) expect_identical(msgs[[paste(cg, TRUE)]], msgs[[paste(cg, FALSE)]], info = cg)
})

test_that("a treated cohort removed whole by balancing is dropped with a warning, not an internal error", {
  # Same mechanism, applied to a treated cohort: glist was fixed before balancing, so
  # faster_mode = TRUE (the default) stopped with "Internal error: treated group g not
  # found in cohort_vec" and faster_mode = FALSE returned NA for every cell of the
  # cohort with no explanation.
  cases <- list(
    list(d = make_panel(cohorts = c(0, 3, 4, 5)), kill = 4),                 # middle cohort, never-treated present
    list(d = make_panel(cohorts = c(0, 3, 4, 5)), kill = 5),                 # latest cohort, never-treated present
    list(d = make_panel(cohorts = c(3, 4, 5, 6), periods = 1:7), kill = 4)   # middle cohort, no never-treated group
  )
  msgs <- list()
  for (ci in seq_along(cases)) for (cg in c("nevertreated", "notyettreated")) for (fm in c(TRUE, FALSE)) {
    cs <- cases[[ci]]
    info <- paste("kill", cs$kill, cg, "faster_mode =", fm)
    d <- kill_x1(cs$d, gsel = cs$kill)
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d, control_group = cg, faster_mode = fm)
    ))
    m <- grep("removed every unit of the cohort(s) first treated at", w, fixed = TRUE, value = TRUE)
    expect_length(m, 1L)
    expect_match(m, paste0("first treated at ", cs$kill), fixed = TRUE)
    msgs[[paste(ci, cg, fm)]] <- m
    expect_false(any(grepl("Internal error", w)), info = info)
    expect_false(any(grepl(RECHECK_MSG, w, fixed = TRUE)), info = info)
    expect_false(cs$kill %in% fit$group, info = info)
    expect_false(anyNA(fit$att), info = info)
    ref <- quiet_fit(cs$d[cs$d$g != cs$kill, ], control_group = cg, faster_mode = fm)
    expect_same_fit(fit, ref, info = info)
  }
  for (ci in seq_along(cases)) for (cg in c("nevertreated", "notyettreated"))
    expect_identical(msgs[[paste(ci, cg, TRUE)]], msgs[[paste(ci, cg, FALSE)]], info = paste(ci, cg))
})

test_that("the re-check's count refers to the never-treated units left after the row drop", {
  # 10 never-treated units have x1 missing in EVERY period (gone at the row drop),
  # the other 20 only in period 2 (gone at balancing)
  d <- make_panel()
  d$x1[d$g == 0 & d$id <= 10] <- NA
  d$x1[d$g == 0 & d$id > 10 & d$time == 2] <- NA
  msgs <- list()
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d, faster_mode = fm)))
    m <- grep(RECHECK_MSG, w, fixed = TRUE, value = TRUE)
    expect_length(m, 1L)
    expect_match(m, "(20 remained after dropping rows with missing data)", fixed = TRUE)
    expect_false(anyNA(fit$att))
    msgs[[as.character(fm)]] <- m
  }
  expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
})

test_that("balancing that removes every treated cohort errors informatively", {
  cases <- list(
    list(d = kill_x1(make_panel(), gsel = c(3, 4)), msg = "removed every unit of every treated cohort (first treated at 3, 4)"),
    # no never-treated group: the raw-data fallback uses cohort 4 as the comparison, then
    # balancing removes the only estimable cohort
    list(d = kill_x1(make_panel(cohorts = c(3, 4)), gsel = 3), msg = "removed every unit of every treated cohort (first treated at 3)",
         extra = "(There is no never-treated group, so the last treated cohort, first treated at 4, serves only as the comparison group and is not counted.)"),
    list(d = kill_x1(make_panel(cohorts = c(3, 4, 5)), gsel = c(3, 4)), msg = "removed every unit of every treated cohort (first treated at 3, 4)",
         extra = "first treated at 5, serves only as the comparison group")
  )
  for (cs in cases) for (cg in c("nevertreated", "notyettreated")) {
    msgs <- list()
    for (fm in c(TRUE, FALSE)) {
      e <- tryCatch(quiet_fit(cs$d, control_group = cg, faster_mode = fm), error = function(e) conditionMessage(e))
      expect_type(e, "character")
      expect_match(e, cs$msg, fixed = TRUE)
      expect_match(e, "allow_unbalanced_panel = TRUE", fixed = TRUE)
      if (!is.null(cs$extra)) expect_match(e, cs$extra, fixed = TRUE)
      msgs[[as.character(fm)]] <- e
    }
    expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
  }
  # every unit removed (each cohort misses a different period): the balanced-panel
  # remedy, not "panel = FALSE" / "revisiting idname"
  d <- make_panel()
  d$x1[d$g == 0 & d$time == 2] <- NA
  d$x1[d$g == 3 & d$time == 3] <- NA
  d$x1[d$g == 4 & d$time == 4] <- NA
  msgs <- list()
  for (fm in c(TRUE, FALSE)) {
    e <- tryCatch(quiet_fit(d, faster_mode = fm), error = function(e) conditionMessage(e))
    expect_type(e, "character")
    expect_match(e, "Converting to a balanced panel removed every unit", fixed = TRUE)
    expect_match(e, "allow_unbalanced_panel = TRUE", fixed = TRUE)
    msgs[[as.character(fm)]] <- e
  }
  expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
})

test_that("the no-never-treated error also names a cohort that balancing removed", {
  # balancing removes the never-treated group AND cohort 3; the re-check then uses
  # cohort 4 as the comparison, leaving nothing to estimate: the error must mention
  # cohort 3 and the way to keep it, not just "supply two cohorts"
  cases <- list(
    list(d = kill_x1(make_panel(), gsel = c(0, 3)),
         need = c("first treated at 4", "removed every unit of the cohort(s) first treated at 3", "all never-treated units (30 remained")),
    # no never-treated group: balancing removes cohort 4 and the comparison cohort 5
    list(d = kill_x1(make_panel(cohorts = c(3, 4, 5)), gsel = c(4, 5)),
         need = c("first treated at 3", "every unit of the cohort(s) first treated at 4", "every unit of the comparison cohort (first treated at 5)"))
  )
  for (cs in cases) for (cg in c("nevertreated", "notyettreated")) {
    msgs <- list()
    for (fm in c(TRUE, FALSE)) {
      e <- tryCatch(quiet_fit(cs$d, control_group = cg, faster_mode = fm), error = function(e) conditionMessage(e))
      expect_type(e, "character")
      expect_match(e, "no never-treated group in the estimation sample", fixed = TRUE)
      for (k in cs$need) expect_match(e, k, fixed = TRUE)
      expect_match(e, "allow_unbalanced_panel = TRUE", fixed = TRUE)
      msgs[[as.character(fm)]] <- e
    }
    expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
  }
})

test_that("a single treated cohort left without any comparison group errors informatively", {
  d <- kill_x1(make_panel(cohorts = c(0, 3)), gsel = 0)
  for (cg in c("nevertreated", "notyettreated")) {
    msgs <- list()
    for (fm in c(TRUE, FALSE)) {
      w <- capture_warnings(
        res <- tryCatch(suppressMessages(fit_gt(d, control_group = cg, faster_mode = fm)),
                        error = function(e) e)
      )
      expect_s3_class(res, "error")
      e <- conditionMessage(res)
      expect_match(e, "no never-treated group in the estimation sample", fixed = TRUE)
      expect_match(e, "first treated at 3", fixed = TRUE)
      # the error itself (printed before the deferred warnings) says what balancing
      # removed and how to keep it
      expect_match(e, "all never-treated units (30 remained after dropping rows with missing data)", fixed = TRUE)
      expect_match(e, "allow_unbalanced_panel = TRUE", fixed = TRUE)
      expect_true(any(grepl(RECHECK_MSG, w, fixed = TRUE)), info = paste(cg, "faster_mode =", fm))
      msgs[[as.character(fm)]] <- e
    }
    expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
  }
})

test_that("a fallback that leaves no period, or a single period, errors identically on both paths", {
  # every period at or after the cutoff: a single cohort treated in the first period
  d1 <- make_panel(cohorts = 1)
  # gname on a different scale than tname (every unit coded 1 in a gname meant as a
  # 0/1 treatment indicator, against calendar years)
  d2 <- make_panel(cohorts = 1, periods = 2010:2014)
  # exactly one period left after the re-check (anticipation shifts the cutoff)
  d3 <- kill_x1(make_panel(cohorts = c(0, 3, 5), periods = c(1, 4, 7)), gsel = 0, tsel = 4)
  for (cg in c("nevertreated", "notyettreated")) {
    msgs <- list()
    for (fm in c(TRUE, FALSE)) {
      e1 <- tryCatch(quiet_fit(d1, control_group = cg, faster_mode = fm), error = function(e) conditionMessage(e))
      e2 <- tryCatch(quiet_fit(d2, control_group = cg, faster_mode = fm), error = function(e) conditionMessage(e))
      e3 <- tryCatch(quiet_fit(d3, control_group = cg, anticipation = 1, faster_mode = fm), error = function(e) conditionMessage(e))
      expect_type(e1, "character"); expect_type(e2, "character"); expect_type(e3, "character")
      expect_match(e1, "no period remains once that cohort is used as the comparison group", fixed = TRUE)
      expect_match(e2, "coded as the period of first treatment on the same scale as 'tname'", fixed = TRUE)
      expect_match(e3, "Only one time period remains", fixed = TRUE)
      msgs[[as.character(fm)]] <- c(e1, e2, e3)
    }
    expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
  }
})

test_that("re-check is a no-op when never-treated units survive balancing or no balancing happens", {
  d <- make_panel()
  clean <- quiet_fit(d)
  # only half of the never-treated units lose a period
  d_half <- d
  d_half$x1[d_half$g == 0 & d_half$time == 2 & d_half$id <= 15] <- NA
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d_half, faster_mode = fm)))
    expect_false(any(grepl("No never-treated group|removed every unit", w)), info = paste("faster_mode =", fm))
    expect_false(anyNA(fit$att))
    # same cell set as the clean run, 15 fewer units, no period dropped
    expect_identical(fit$group, clean$group)
    expect_identical(fit$t, clean$t)
    expect_identical(fit$n, 75L)
    # unbalanced panel allowed: the never-treated units are kept (minus their
    # period-2 rows), so there is nothing to re-check; cells that do not touch
    # period 2 are estimated as usual
    d_all <- kill_x1(d, gsel = 0)
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d_all, allow_unbalanced_panel = TRUE, faster_mode = fm)
    ))
    expect_false(any(grepl("No never-treated group|removed every unit", w)), info = paste("faster_mode =", fm))
    late4 <- fit$att[fit$group == 4 & fit$t %in% c(4, 5)]
    expect_length(late4, 2L)
    expect_false(anyNA(late4))
  }
})

test_that("the re-check applies to the panel balanced over all periods (documented rule)", {
  # A cohort-3 unit missing only in period 5 -- a period the fallback then drops -- is
  # still removed by balancing, because the panel is balanced over every period of the
  # data before the fallback runs. The reference is the same balanced sample.
  d <- kill_x1(make_panel(), gsel = 0)
  d$x1[d$id == 31 & d$time == 5] <- NA
  d_bal <- d[d$g != 0 & d$id != 31, ]
  for (fm in c(TRUE, FALSE)) {
    fit <- quiet_fit(d, faster_mode = fm)
    expect_identical(fit$n, 59L)
    expect_same_fit(fit, quiet_fit(d_bal, faster_mode = fm), info = paste("faster_mode =", fm))
  }
})

test_that("the re-check reproduces the balanced-sample fit under every estimator and base period", {
  d <- kill_x1(make_panel(), gsel = 0)
  d_bal <- d[d$g != 0, ]
  for (m in c("dr", "reg", "ipw")) for (bp in c("varying", "universal")) for (fm in c(TRUE, FALSE)) {
    info <- paste(m, bp, "faster_mode =", fm)
    fit <- quiet_fit(d, est_method = m, base_period = bp, faster_mode = fm)
    ref <- quiet_fit(d_bal, est_method = m, base_period = bp, faster_mode = fm)
    expect_same_fit(fit, ref, info = info)
    expect_false(anyNA(fit$att), info = info)
  }
})

test_that("'notyettreated' without a never-treated group in the raw data stays silent", {
  # Using the last cohort as a comparison group only is the routine design under
  # "notyettreated" (common in practice), so it must not warn on every call; the
  # fallback is announced only when balancing is what removed the never-treated
  # group (tests above).
  d <- make_panel(cohorts = c(3, 4))
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d, control_group = "notyettreated", faster_mode = fm)
    ))
    expect_false(any(grepl("No never-treated group", w)), info = paste("faster_mode =", fm))
    expect_false(4 %in% fit$group)
    expect_false(anyNA(fit$att))
  }
})

test_that("the 'nevertreated' fallback warning is precise and identical across modes", {
  d <- make_panel(cohorts = c(3, 4))
  msgs <- list()
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d, control_group = "nevertreated", faster_mode = fm)))
    m <- grep("No never-treated group is available. ", w, fixed = TRUE, value = TRUE)
    expect_length(m, 1L)
    expect_match(m, "first treated at 4", fixed = TRUE)
    expect_match(m, "at or after its treatment date (net of anticipation) is being filtered out", fixed = TRUE)
    msgs[[as.character(fm)]] <- m
  }
  expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
  # integer gname with large codes: the slow path used to label cohorts from the
  # integer ("500000") and the fast path from a double ("5e+05")
  di <- make_panel(cohorts = c(300000, 400000, 500000), periods = 100000L * 1:6)
  di$g <- as.integer(di$g); di$time <- as.integer(di$time)
  msgs <- list(); fits <- list()
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fits[[as.character(fm)]] <- fit_gt(di, control_group = "nevertreated", faster_mode = fm)))
    msgs[[as.character(fm)]] <- grep("No never-treated group is available. ", w, fixed = TRUE, value = TRUE)
  }
  expect_length(msgs[["TRUE"]], 1L)
  expect_identical(msgs[["TRUE"]], msgs[["FALSE"]])
  expect_identical(fits[["TRUE"]]$group, fits[["FALSE"]]$group)
})

test_that("an estimation sample emptied by the missing-data row drop stops with the missing-data message", {
  d <- make_panel()
  d$y <- NA_real_   # e.g. the wrong outcome column
  d_rc <- d; d_rc$id <- seq_len(nrow(d_rc))
  for (fm in c(TRUE, FALSE)) {
    for (dd in list(list(d = d, panel = TRUE), list(d = d_rc, panel = FALSE))) {
      w <- character()
      e <- tryCatch(withCallingHandlers(
        suppressMessages(fit_gt(dd$d, panel = dd$panel, faster_mode = fm)),
        warning = function(w0) { w <<- c(w, conditionMessage(w0)); invokeRestart("muffleWarning") }),
        error = function(e) conditionMessage(e))
      expect_type(e, "character")
      expect_match(e, "All observations were dropped due to missing data", fixed = TRUE)
      expect_false(grepl("-Inf", e, fixed = TRUE))
      expect_false(any(grepl("-Inf|no non-missing arguments", w)), info = paste("panel =", dd$panel, "faster_mode =", fm))
    }
    # an all-NA weights column empties the sample too, and the message says so
    dw <- make_panel(); dw$w <- NA_real_
    e <- tryCatch(quiet_fit(dw, weightsname = "w", faster_mode = fm), error = function(e) conditionMessage(e))
    expect_match(e, "All observations were dropped due to missing data", fixed = TRUE)
    expect_match(e, "weights", fixed = TRUE)
    # a zero-row input is not "dropped"
    e0 <- tryCatch(quiet_fit(make_panel()[0, ], faster_mode = fm), error = function(e) conditionMessage(e))
    expect_type(e0, "character")
    expect_false(grepl("-Inf", e0, fixed = TRUE))
  }
})

test_that("print_details labels each cell by its own period on the fast path", {
  d <- make_panel()
  for (bp in c("varying", "universal")) {
    out <- capture.output(fit <- quiet_fit(d, base_period = bp, faster_mode = TRUE, print_details = TRUE))
    labels <- regmatches(out, regexpr("Evaluating \\(g,t\\) = \\([0-9]+,[0-9]+\\)", out))
    labels <- sub("Evaluating \\(g,t\\) = ", "", labels)
    expect_identical(sort(labels), sort(paste0("(", fit$group, ",", fit$t, ")")), info = bp)
  }
})

test_that("aggte() explains an object with no post-treatment cells", {
  # no never-treated group + anticipation: the fallback drops periods >= 6 - 2, so
  # only cohort 4's anticipation-window cells (4,2), (4,3) remain -- none NA
  d <- make_panel(cohorts = c(4, 6), periods = 1:6)
  for (fm in c(TRUE, FALSE)) {
    fit <- quiet_fit(d, anticipation = 2, faster_mode = fm)
    expect_true(all(fit$t < fit$group))
    expect_false(anyNA(fit$att))
    for (ty in c("simple", "dynamic", "group", "calendar")) {
      expect_error(suppressWarnings(aggte(fit, type = ty)),
                   "No post-treatment ATT(g,t) cells are available to aggregate", fixed = TRUE)
    }
  }
})

test_that("cohort and period codes in messages are printed as written, not in scientific notation", {
  d <- make_panel(cohorts = c(300000, 400000, 500000), periods = 100000 * 1:6)
  d$w0 <- ifelse(d$g == 400000, 0, 1)
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(fit <- fit_gt(d, control_group = "nevertreated", weightsname = "w0", faster_mode = fm)))
    expect_true(any(grepl("first treated at 500000)", w, fixed = TRUE)), info = paste("faster_mode =", fm))
    expect_true(any(grepl("(g, t) = (400000, 200000) is NaN", w, fixed = TRUE)), info = paste("faster_mode =", fm))
    expect_false(any(grepl("e+05", w, fixed = TRUE)), info = paste("faster_mode =", fm))
  }
  # a lost cohort and the informative error use the same formatting
  d2 <- kill_x1(make_panel(cohorts = c(0, 300000, 400000), periods = 100000 * 1:5), gsel = c(0, 300000), tsel = 200000)
  e <- tryCatch(quiet_fit(d2), error = function(e) conditionMessage(e))
  expect_match(e, "first treated at 400000", fixed = TRUE)
  expect_match(e, "cohort(s) first treated at 300000", fixed = TRUE)
  expect_false(grepl("e+05", e, fixed = TRUE))
})

test_that("a NaN 2x2 estimate is announced, not silently set to NA (both modes)", {
  # every unit of cohort 4 has weight zero: its cells have no effective treated
  # observations and the 2x2 estimator returns NaN
  d <- make_panel()
  d$w0 <- ifelse(d$g == 4, 0, 1)
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d, weightsname = "w0", faster_mode = fm)
    ))
    na_cells <- which(is.na(fit$att))
    expect_identical(fit$group[na_cells], rep(4, 4))
    nan_w <- grep("is NaN", w, fixed = TRUE, value = TRUE)
    expect_length(nan_w, length(na_cells))
    for (i in na_cells) {
      expect_true(any(grepl(paste0("(g, t) = (4, ", fit$t[i], ")"), nan_w, fixed = TRUE)),
                  info = paste("cell (4,", fit$t[i], ") faster_mode =", fm))
    }
  }
  slow <- quiet_fit(d, weightsname = "w0", faster_mode = FALSE)
  fast <- quiet_fit(d, weightsname = "w0", faster_mode = TRUE)
  expect_equal(slow$att, fast$att, tolerance = 1e-10)
  # repeated cross sections reach the slow path's other NaN site
  d_rc <- d; d_rc$id <- seq_len(nrow(d_rc))
  for (m in c("reg", "ipw")) for (fm in c(TRUE, FALSE)) {
    info <- paste("RC", m, "faster_mode =", fm)
    w <- capture_warnings(suppressMessages(
      fit <- fit_gt(d_rc, panel = FALSE, est_method = m, weightsname = "w0", faster_mode = fm)
    ))
    na_cells <- which(is.na(fit$att))
    expect_identical(fit$group[na_cells], rep(4, 4), info = info)
    nan_w <- grep("is NaN", w, fixed = TRUE, value = TRUE)
    expect_length(nan_w, length(na_cells))
    for (i in na_cells) {
      expect_true(any(grepl(paste0("(g, t) = (4, ", fit$t[i], ")"), nan_w, fixed = TRUE)), info = info)
    }
  }
  slow <- quiet_fit(d_rc, panel = FALSE, est_method = "reg", weightsname = "w0", faster_mode = FALSE)
  fast <- quiet_fit(d_rc, panel = FALSE, est_method = "reg", weightsname = "w0", faster_mode = TRUE)
  expect_equal(slow$att, fast$att, tolerance = 1e-10)
})

test_that("faster_mode = TRUE names a cell with no control units instead of returning NA silently", {
  # Repeated cross sections with never-treated units observed only in periods 1 and
  # 5: cells (3,3) and (3,4) have no control observations in either of their two
  # periods. The slow path has always warned per period; the fast path returned NA
  # for these cells without any warning.
  set.seed(3)
  d <- expand.grid(rep = 1:40, time = 1:5, g = c(0, 3))
  d$y <- rnorm(nrow(d)) + 0.2 * d$time + (d$g > 0 & d$time >= d$g)
  d <- d[!(d$g == 0 & d$time %in% 2:4), ]
  d$id <- seq_len(nrow(d))
  for (fm in c(TRUE, FALSE)) {
    w <- capture_warnings(suppressMessages(
      fit <- att_gt(yname = "y", tname = "time", idname = "id", gname = "g", data = d,
                    panel = FALSE, bstrap = FALSE, cband = FALSE, faster_mode = fm)
    ))
    expect_true(all(is.na(fit$att)))
    # at least one control-unit warning per NA cell, on both paths
    expect_gte(sum(grepl("No available control units for group 3", w, fixed = TRUE)), sum(is.na(fit$att)))
    if (fm) {
      for (tt in c(3, 4)) {
        expect_true(any(grepl(paste0("No available control units for group 3 in time period ", tt,
                                     "; the ATT for this cell is set to NA"), w, fixed = TRUE)),
                    info = paste("cell (3,", tt, ")"))
      }
    }
  }
  # the other branch: the treated group itself is absent from both periods of a cell
  set.seed(4)
  d <- expand.grid(rep = 1:40, time = 1:5, g = c(0, 3))
  d$y <- rnorm(nrow(d)) + 0.2 * d$time + (d$g > 0 & d$time >= d$g)
  d <- d[!(d$g == 3 & d$time %in% 2:3), ]
  d$id <- seq_len(nrow(d))
  w <- capture_warnings(suppressMessages(
    fit <- att_gt(yname = "y", tname = "time", idname = "id", gname = "g", data = d,
                  panel = FALSE, bstrap = FALSE, cband = FALSE, faster_mode = TRUE)
  ))
  expect_true(is.na(fit$att[fit$group == 3 & fit$t == 3]))
  expect_true(any(grepl("No units in group 3 in time period 3; the ATT for this cell is set to NA", w, fixed = TRUE)))
})
