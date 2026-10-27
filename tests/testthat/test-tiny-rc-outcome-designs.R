# Issue #64: the default RC DR estimator needs treated outcome regressions.
# Its old regression test lost the default-DR call after a native fastglm crash
# (a811d46). Use deterministic data and verify rejection BEFORE that backend.
tiny_rc_data <- function(n_treated = 3L) {
  cube <- as.matrix(expand.grid(x1 = c(-1, 1), x2 = c(-1, 1),
                               x3 = c(-1, 1)))
  regular <- cube[rep(1:8, 2), ]
  tiny <- if (n_treated == 3L) {
    rbind(c(-1, -1, -1), c(1, 1, -1), c(0, 0, 1))
  } else if (n_treated == 4L) {
    rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  } else {
    regular
  }
  cells <- lapply(1:3, function(t) {
    do.call(rbind, lapply(c(0, 2, 3), function(g) {
      x <- if (g == 2) tiny else regular
      data.frame(period = t, G = g, X = x[, 1], X2 = x[, 2], X3 = x[, 3])
    }))
  })
  d <- do.call(rbind, cells)
  d$id <- seq_len(nrow(d))
  d$w <- 1
  d$Y <- 0.2 * d$period + 0.3 * d$X - 0.1 * d$X2 + 0.2 * d$X3 +
    (d$G > 0 & d$period >= d$G) + 0.1 * sin(d$id)
  d
}

fit_tiny_rc <- function(d, mode, ...) {
  suppressMessages(att_gt(yname = "Y", tname = "period", idname = "id",
    gname = "G", xformla = ~X + X2 + X3, data = d, panel = FALSE,
    faster_mode = mode, bstrap = FALSE, cband = FALSE, ...))
}

test_that("small RC treated groups are rejected before the DR backend", {
  d <- tiny_rc_data()
  original <- DRDID::drdid_rc
  calls <- integer()
  backend <- function(y, post, D, ...) {
    calls <<- c(calls, sum(D))
    # A failed preflight must never reach an unsafe native fit, even on an old
    # dependency. Fail safely here so this regression test cannot segfault R.
    if (sum(D) == 6L) stop("undersized treated regression reached the backend")
    original(y = y, post = post, D = D, ...)
  }
  local_mocked_bindings(drdid_rc = backend, .package = "DRDID")
  local_mocked_bindings(drdid_rc = backend, .package = "did")

  results <- lapply(c(FALSE, TRUE), function(mode) {
    w <- capture_warnings(r <- fit_tiny_rc(d, mode, est_method = "dr"))
    failures <- grep("fewer positive-weight observations than covariates", w,
                     value = TRUE)
    expect_length(failures, 2L)
    expect_true(any(grepl("group 2 in time period 2", failures)))
    expect_true(any(grepl("group 2 in time period 3", failures)))
    expect_false(any(grepl("reached the backend", w)))
    expect_identical(is.na(r$att), r$group == 2)
    expect_identical(is.na(r$se), r$group == 2)
    expect_true(all(is.na(as.matrix(r$inffunc)[, r$group == 2])))
    expect_true(all(is.finite(r$att[r$group == 3])))
    expect_true(all(is.finite(r$se[r$group == 3])))
    expect_true(all(is.finite(as.matrix(r$inffunc)[, r$group == 3])))
    r
  })
  expect_equal(calls, rep(32L, 4))
  expect_equal(results[[1]]$att, results[[2]]$att, tolerance = 1e-10)
  expect_equal(results[[1]]$se, results[[2]]$se, tolerance = 1e-10)
})

test_that("RC methods without treated outcome regressions retain small groups", {
  d <- tiny_rc_data()
  for (mode in c(FALSE, TRUE)) {
    for (method in list("reg", "ipw", DRDID::drdid_rc1)) {
      w <- capture_warnings(r <- fit_tiny_rc(d, mode, est_method = method))
      expect_false(any(grepl("fewer positive-weight observations", w)))
      expect_true(all(is.finite(r$att)))
      expect_true(all(is.finite(r$se)))
      expect_true(all(is.finite(r$inffunc)))
    }
  }
})

test_that("zero-weight rows cannot identify RC treated outcome regressions", {
  d <- tiny_rc_data(16L)
  for (t in 1:3) {
    rows <- which(d$G == 2 & d$period == t)
    d$w[rows[-(1:3)]] <- 0
  }
  for (mode in c(FALSE, TRUE)) {
    w <- capture_warnings(r <- fit_tiny_rc(d, mode, est_method = "dr",
                                          weightsname = "w"))
    expect_equal(sum(grepl("fewer positive-weight observations", w)), 2L)
    expect_identical(is.na(r$att), r$group == 2)
    expect_identical(is.na(r$se), r$group == 2)
    expect_true(all(is.finite(r$att[r$group == 3])))
  }
})

test_that("square full-rank RC treated regressions are allowed", {
  d <- tiny_rc_data(4L)
  for (mode in c(FALSE, TRUE)) {
    w <- capture_warnings(r <- fit_tiny_rc(d, mode, est_method = "dr"))
    expect_false(any(grepl("fewer positive-weight observations", w)))
    expect_true(all(is.finite(r$att)))
    expect_true(all(is.finite(r$se)))
    expect_true(all(is.finite(r$inffunc)))

    # A direct 2x2 reference fixes both the estimator and IF normalization for
    # an estimable cell, rather than only comparing fast and slow code paths.
    cell <- d$G %in% c(0, 2) & d$period %in% 1:2
    ref <- DRDID::drdid_rc(d$Y[cell], as.numeric(d$period[cell] == 2),
      as.numeric(d$G[cell] == 2), model.matrix(~X + X2 + X3, d[cell, ]),
      inffunc = TRUE)
    col <- which(r$group == 2 & r$t == 2)
    expect_equal(r$att[col], ref$ATT, tolerance = 1e-10)
    expected_if <- numeric(nrow(d))
    expected_if[cell] <- nrow(d) / sum(cell) * ref$att.inf.func
    # The modes use different row orders; compare by sampling-unit id.
    aligned_if <- expected_if[match(rownames(r$inffunc), as.character(d$id))]
    expect_equal(as.numeric(r$inffunc[, col]), aligned_if, tolerance = 1e-10)
    expect_equal(r$se[col], sqrt(mean(expected_if^2) / nrow(d)),
                 tolerance = 1e-10)
  }
})

test_that("undersized control-period regressions fail only affected RC cells", {
  d <- tiny_rc_data(16L)
  rows <- which(d$G == 0 & d$period == 1)
  d$w[rows[-(1:3)]] <- 0
  for (mode in c(FALSE, TRUE)) {
    for (method in c("dr", "reg")) {
      w <- capture_warnings(r <- fit_tiny_rc(d, mode, est_method = method,
                                            weightsname = "w"))
      # Group 2 uses period 1 as its base in both post-treatment cells;
      # group 3 uses it only in its period-2 placebo cell.
      failed <- r$group == 2 | r$t == 2
      expect_equal(sum(grepl("fewer positive-weight observations", w)), 3L)
      expect_identical(is.na(r$att), failed)
      expect_identical(is.na(r$se), failed)
      expect_true(all(is.finite(r$att[!failed])))
      expect_true(all(is.finite(r$se[!failed])))
    }
  }
})

test_that("unbalanced and varying-weight panel RC routes check treated regressions", {
  d <- tiny_rc_data()
  d$id <- (d$G + 1) * 1000 + ave(d$id, d$G, d$period, FUN = seq_along)
  d$w <- 1 + d$period / 10
  for (mode in c(FALSE, TRUE)) {
    for (route in c("unbalanced", "varying")) {
      dat <- if (route == "unbalanced") {
        d[-which(d$G == 3 & d$period == 1)[1], ]
      } else d
      for (do_inf in c(TRUE, FALSE)) {
        w <- capture_warnings(r <- suppressMessages(att_gt(yname = "Y",
          tname = "period", idname = "id", gname = "G", xformla = ~X + X2 + X3,
          data = dat, panel = TRUE, faster_mode = mode, bstrap = FALSE,
          cband = FALSE, weightsname = "w", est_method = "dr",
          allow_unbalanced_panel = route == "unbalanced",
          fix_weights = if (route == "varying") "varying" else NULL,
          compute_inffunc = do_inf)))
        expect_equal(sum(grepl("fewer positive-weight observations", w)), 2L)
        expect_identical(is.na(r$att), r$group == 2)
        expect_true(all(is.finite(r$att[r$group == 3])))
        if (do_inf) {
          expect_identical(is.na(r$se), r$group == 2)
          expect_true(all(is.finite(r$se[r$group == 3])))
          expect_true(all(is.na(as.matrix(r$inffunc)[, r$group == 2])))
        } else expect_null(r$inffunc)
      }
    }

    # Ordinary balanced-panel DR only fits a control outcome regression.
    w <- capture_warnings(r <- suppressMessages(att_gt(yname = "Y",
      tname = "period", idname = "id", gname = "G", xformla = ~X + X2 + X3,
      data = d, panel = TRUE, faster_mode = mode, bstrap = FALSE,
      cband = FALSE, weightsname = "w", est_method = "dr")))
    expect_false(any(grepl("fewer positive-weight observations", w)))
    expect_true(all(is.finite(r$att)))
    expect_true(all(is.finite(r$se)))
  }
})
