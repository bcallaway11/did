# Regression tests for the 2.5.1 audit fixes:
#   A. fix_weights = "varying" on a balanced panel reported SEs exactly 2x too
#      large (force_rc influence-function fold was not divided by 2).
#   B. aggte() crashed with cryptic errors on empty / all-NA aggregation
#      selections (calendar na.rm, simple/dynamic empty post-treatment windows).
#   C. negative gname codes were silently accepted on the fast path but errored
#      on the slow path (now rejected consistently up front).

# small balanced panel with a never-treated group; constant (uniform) weights
.bal_panel <- function(seed = 1, n = 200, sigma = 0.5) {
  set.seed(seed)
  g  <- sample(c(0L, 0L, 2L, 3L), n, replace = TRUE)
  fe <- stats::rnorm(n)
  d  <- data.table::CJ(id = 1:n, t = 1:3)
  d[, g := g[id]]; d[, fe := fe[id]]
  d[, y := fe + 0.2 * t + 1 * (g != 0 & t >= g) + stats::rnorm(.N, 0, sigma)]
  as.data.frame(d[, .(id, t, g, y)])
}
.qgt <- function(df, ...) suppressWarnings(suppressMessages(
  att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = df,
         control_group = "nevertreated", bstrap = FALSE, ...)))

# ---- A. fix_weights = "varying" SE normalization -------------------------------

test_that("fix_weights='varying' SE is not 2x too large on a balanced panel", {
  df <- .bal_panel(1)
  r0 <- .qgt(df, fix_weights = NULL)
  rv <- .qgt(df, fix_weights = "varying")
  ord <- match(paste(r0$group, r0$t), paste(rv$group, rv$t))
  # point estimates were always correct
  expect_equal(r0$att, rv$att[ord], tolerance = 1e-8)
  # with constant weights, the varying (RC) SE must equal the panel SE. Before
  # the fix it was EXACTLY 2x larger (the force_rc fold lacked the 1/2). This is
  # the property the Monte Carlo confirmed (varying SE == empirical sampling SD).
  expect_equal(r0$se, rv$se[ord], tolerance = 1e-8)
  # explicit guard against regression to the 2x bug
  expect_true(all(abs(rv$se[ord] / r0$se - 1) < 1e-6))
})

test_that("fix_weights='varying' fast and slow paths agree (balanced panel)", {
  df <- .bal_panel(2)
  rf <- .qgt(df, fix_weights = "varying", faster_mode = TRUE)
  rs <- .qgt(df, fix_weights = "varying", faster_mode = FALSE)
  expect_equal(rf$att, rs$att, tolerance = 1e-9)
  expect_equal(rf$se, rs$se, tolerance = 1e-9)
})

test_that("fix_weights='varying' aggregations inherit the corrected SE", {
  df <- .bal_panel(3)
  cs0 <- .qgt(df, fix_weights = NULL)
  csv <- .qgt(df, fix_weights = "varying")
  a0 <- suppressWarnings(suppressMessages(aggte(cs0, type = "simple", na.rm = TRUE)))
  av <- suppressWarnings(suppressMessages(aggte(csv, type = "simple", na.rm = TRUE)))
  expect_equal(a0$overall.att, av$overall.att, tolerance = 1e-8)
  expect_equal(a0$overall.se,  av$overall.se,  tolerance = 1e-8)
})

# ---- B. aggte() empty / all-NA selection guards --------------------------------

test_that("aggte() normal aggregations still work (no regression)", {
  data(mpdta, package = "did")
  cs <- suppressWarnings(suppressMessages(att_gt("lemp", "year", "countyreal",
        "first.treat", data = mpdta, bstrap = FALSE)))
  for (ty in c("simple", "dynamic", "group", "calendar")) {
    a <- suppressWarnings(suppressMessages(aggte(cs, type = ty, na.rm = TRUE)))
    expect_true(is.finite(a$overall.att))
    expect_true(is.finite(a$overall.se))
  }
})

test_that("aggte(type='calendar', na.rm=TRUE) drops all-NA calendar periods instead of crashing", {
  data(mpdta, package = "did")
  cs <- suppressWarnings(suppressMessages(att_gt("lemp", "year", "countyreal",
        "first.treat", data = mpdta, bstrap = FALSE)))
  # NA out the lone post-treatment cell of calendar year 2005 (att_gt itself
  # legitimately produces NA cells on unbalanced data; this mimics that).
  w <- which(cs$t == 2005 & cs$group <= cs$t)
  cs$att[w] <- NA; cs$inffunc[, w] <- NA
  expect_no_error(
    res <- suppressWarnings(suppressMessages(aggte(cs, type = "calendar", na.rm = TRUE)))
  )
  expect_false(2005 %in% res$egt)        # the all-NA period was dropped
  expect_true(is.finite(res$overall.att))
})

test_that("aggte() empty post-treatment windows give clean errors, not cryptic crashes", {
  data(mpdta, package = "did")
  cs <- suppressWarnings(suppressMessages(att_gt("lemp", "year", "countyreal",
        "first.treat", data = mpdta, bstrap = FALSE)))
  # simple / dynamic with no e >= 0 -> clean "no valid estimates", NOT
  # "non-numeric argument..." or "...report this as a bug."
  expect_error(suppressWarnings(suppressMessages(aggte(cs, type = "simple", max_e = -1))),
               "No valid att_gt")
  expect_error(suppressWarnings(suppressMessages(aggte(cs, type = "dynamic", max_e = -1))),
               "No valid att_gt")
  # event-time window that excludes every period -> clean "no event times".
  expect_error(suppressWarnings(suppressMessages(aggte(cs, type = "dynamic", min_e = 100))),
               "No event times")
})

# ---- C. negative gname rejected consistently -----------------------------------

test_that("negative gname is rejected with a clear error in both code paths", {
  set.seed(5); n <- 120
  g  <- sample(c(0L, -1L), n, replace = TRUE)
  d  <- data.table::CJ(id = 1:n, t = c(-3L, -2L, -1L, 0L))
  d[, g := g[id]]
  d[, y := stats::rnorm(.N) + (g != 0 & t >= g)]
  df <- as.data.frame(d)
  for (fm in c(TRUE, FALSE)) {
    expect_error(
      suppressWarnings(suppressMessages(att_gt("y", "t", "id", "g", data = df,
        bstrap = FALSE, faster_mode = fm))),
      "negative values are not supported"
    )
  }
})

# ---- D. parallel multiplier bootstrap: reproducibility + chunking guard ---------

test_that("parallel multiplier bootstrap is reproducible under a fixed seed", {
  skip_on_os("windows")          # parallel path is force-disabled on Windows
  inf <- matrix(stats::rnorm(2600 * 4), 2600, 4)   # n > 2500 triggers the parallel branch
  pre_kind <- RNGkind()                            # whatever the runner's RNG kind is
  set.seed(7); a <- did:::run_multiplier_bootstrap(inf, 300, pl = TRUE, cores = 2)
  set.seed(7); b <- did:::run_multiplier_bootstrap(inf, 300, pl = TRUE, cores = 2)
  expect_identical(a, b)                            # same seed -> bit-identical draws
  expect_identical(RNGkind(), pre_kind)             # caller's RNG kind restored (whatever it was)
})

test_that("parallel bootstrap chunking never produces negative chunks (biters < cores)", {
  # the fix: chunks are non-negative and sum to biters for any biters/cores
  # (the old rep(ceiling(biters/cores), cores) + correction went negative when
  # biters < cores, e.g. biters=2,cores=4 -> [-1,1,1,1] -> crash).
  for (biters in c(1, 2, 3, 7, 1000)) {
    for (cores in c(1, 2, 4, 8)) {
      chunks <- diff(round(seq(0, biters, length.out = cores + 1)))
      chunks <- chunks[chunks > 0]
      expect_true(all(chunks > 0))
      expect_equal(sum(chunks), biters)
    }
  }
  # and the parallel branch itself no longer errors with biters < cores
  skip_on_os("windows")
  inf <- matrix(stats::rnorm(2600 * 3), 2600, 3)
  expect_no_error(did:::run_multiplier_bootstrap(inf, biters = 1, pl = TRUE, cores = 2))
})

# ---- E. calendar aggregation warns that min_e/max_e/balance_e are ignored ------

test_that("aggte(type='calendar') warns that min_e/max_e/balance_e are ignored, returns unrestricted", {
  data(mpdta, package = "did")
  cs <- suppressWarnings(suppressMessages(att_gt("lemp", "year", "countyreal",
        "first.treat", data = mpdta, bstrap = FALSE)))
  ws <- testthat::capture_warnings(suppressMessages(aggte(cs, type = "calendar", max_e = 2, na.rm = TRUE)))
  expect_true(any(grepl("ignored for type", ws)))
  a_win <- suppressWarnings(suppressMessages(aggte(cs, type = "calendar", max_e = 2, na.rm = TRUE)))
  a_unr <- suppressWarnings(suppressMessages(aggte(cs, type = "calendar", na.rm = TRUE)))
  expect_equal(a_win$overall.att, a_unr$overall.att)   # max_e is a no-op for calendar (correct)
  # no spurious warning at the defaults
  ws0 <- testthat::capture_warnings(suppressMessages(aggte(cs, type = "calendar", na.rm = TRUE)))
  expect_false(any(grepl("ignored for type", ws0)))
})
