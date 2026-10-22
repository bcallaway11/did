test_that("RC support counts use the weights after DRDID normalization", {
  D <- c(0, 0, 1, 1)
  post <- c(0, 1, 0, 1)
  for (rescale in c(1, 2^-900)) {
    w <- c(2^1000, 2^1000, 2^-100, 2^-100) * rescale
    expect_true(all(is.finite(w) & w > 0))
    expect_identical(w / mean(w), c(2, 2, 0, 0))
    expect_true(rc_outcome_has_too_few_obs(D, post, w, 1L))
    # Regression-only needs the two controls, which still have positive mass.
    expect_false(rc_outcome_has_too_few_obs(D, post, w, 1L, FALSE))
    expect_true(rc_outcome_has_too_few_obs(D, post, rev(w), 1L, FALSE))
  }
  # Do not impose a positive-weight tolerance: the smallest subnormal still
  # counts, while an exactly rounded zero does not. No native fit is invoked.
  retained <- c(2^1000, 2^1000, 2^-75, 2^-75)
  lost <- c(2^1000, 2^1000, 2^-76, 2^-76)
  expect_identical((retained / mean(retained))[3], 2^-1074)
  expect_false(rc_outcome_has_too_few_obs(D, post, retained, 1L))
  expect_identical((lost / mean(lost))[3], 0)
  expect_true(rc_outcome_has_too_few_obs(D, post, lost, 1L))
  for (level in c(2^-1000, 1, 2^1000)) {
    expect_false(rc_outcome_has_too_few_obs(D, post, rep(level, 4), 1L))
  }
  for (w in list(numeric(), rep(0, 4), rep(Inf, 4), rep(NA_real_, 4))) {
    expect_true(rc_outcome_has_too_few_obs(D[seq_along(w)], post[seq_along(w)],
                                         w, 1L))
  }
})

test_that("normalization loss is rejected before DR/reg backend entry", {
  # Enough physical rows for the existing global cohort checks and every OR.
  d <- expand.grid(j = 1:16, G = c(0, 2, 3), period = 1:3)
  d$id <- seq_len(nrow(d))
  d$Y <- 0.2 * d$period + (d$G > 0 & d$period >= d$G) + 0.1 * sin(d$id)
  originals <- list(dr = DRDID::drdid_rc, reg = DRDID::reg_did_rc)
  entries <- list()
  make_backend <- function(method) {
    force(method)
    function(y, post, D, covariates, i.weights, ...) {
      normalized <- i.weights / mean(i.weights)
      required <- if (method == "dr") rep(TRUE, length(D)) else D == 0
      lost <- any(!is.finite(normalized[required]) | normalized[required] <= 0)
      entries[[length(entries) + 1L]] <<- lost
      # Safely fail on a missing preflight. Never run a vulnerable native fit.
      if (lost) stop("normalized support loss reached the backend")
      originals[[method]](y = y, post = post, D = D,
        covariates = covariates, i.weights = i.weights, ...)
    }
  }
  dr_backend <- make_backend("dr")
  reg_backend <- make_backend("reg")
  local_mocked_bindings(drdid_rc = dr_backend, reg_did_rc = reg_backend,
                        .package = "DRDID")
  local_mocked_bindings(drdid_rc = dr_backend, reg_did_rc = reg_backend,
                        .package = "did")

  for (small_group in c(2, 0)) {
    for (method in if (small_group == 2) "dr" else c("dr", "reg")) {
      for (rescale in c(1, 2^-900)) {
        d$w <- ifelse(d$G == small_group, 2^-100, 2^1000) * rescale
        expect_true(all(is.finite(d$w) & d$w > 0))
        for (mode in c(FALSE, TRUE)) {
          for (do_inf in c(FALSE, TRUE)) {
            warnings <- capture_warnings(r <- suppressMessages(att_gt(
              yname = "Y", tname = "period", idname = "id", gname = "G",
              xformla = ~1, data = d, panel = FALSE, weightsname = "w",
              est_method = method, faster_mode = mode, compute_inffunc = do_inf,
              bstrap = FALSE, cband = FALSE)))
            failed <- if (small_group == 2) r$group == 2 else rep(TRUE, length(r$att))
            expect_equal(sum(grepl("fewer positive-weight observations", warnings)),
                         sum(failed))
            expect_false(any(grepl("reached the backend", warnings)))
            expect_identical(is.na(r$att), failed)
            expect_true(all(is.finite(r$att[!failed])))
            if (do_inf) {
              expect_identical(is.na(r$se), failed)
              expect_true(all(is.finite(r$se[!failed])))
              expect_true(all(is.na(as.matrix(r$inffunc)[, failed])))
            } else expect_null(r$inffunc)
          }
        }
      }
    }
  }
  expect_length(entries, 16L)
  expect_false(any(unlist(entries)))
})

test_that("a second local normalization is checked on every RC route", {
  d <- expand.grid(j = 1:16, period = 1:4, D = 0:1)
  d$G <- 2 * d$D
  d$Y <- d$period
  d$w <- ifelse(d$D == 0 & d$period %in% 1:2, 6, 2^-1074)
  # Global mean1.5 leaves the smallest subnormal positive. In the period1/2
  # cell the normalized weights have mean2, so local division rounds it to0.
  expect_identical(mean(d$w), 1.5)
  expect_true(all(d$w / mean(d$w) > 0))
  seen_times <- numeric()
  backend <- function(y, post, ...) {
    seen_times <<- c(seen_times, unique(y[post == 1]))
    # All cells stay behind a safe probe; this fixture estimates no ATT.
    stop("safe non-target backend probe")
  }
  local_mocked_bindings(drdid_rc = backend, .package = "DRDID")
  local_mocked_bindings(drdid_rc = backend, .package = "did")
  for (route in c("rc", "unbalanced", "varying")) {
    dat <- d
    dat$id <- if (route == "rc") seq_len(nrow(dat)) else 1000 * (dat$D + 1) + dat$j
    if (route == "unbalanced") {
      extra <- dat[dat$D == 1 & dat$period == 4 & dat$j == 16, ]
      extra$id <- 3000
      dat <- rbind(dat, extra)
    }
    for (mode in c(FALSE, TRUE)) {
      for (do_inf in c(FALSE, TRUE)) {
        seen_times <- numeric()
        warnings <- capture_warnings(r <- suppressMessages(att_gt(
          yname = "Y", tname = "period", idname = "id", gname = "G",
          xformla = ~1, data = dat, panel = route != "rc", weightsname = "w",
          est_method = "dr", faster_mode = mode, compute_inffunc = do_inf,
          allow_unbalanced_panel = route == "unbalanced",
          fix_weights = if (route == "varying") "varying" else NULL,
          bstrap = FALSE, cband = FALSE)))
        expect_equal(seen_times, c(3, 4))
        expect_true(any(grepl(paste0("fewer positive-weight observations.*",
          "group 2 in time period 2"), warnings)))
        col <- which(r$t == 2)
        expect_true(is.na(r$att[col]))
        if (do_inf) {
          expect_true(is.na(r$se[col]))
          expect_true(all(is.na(r$inffunc[, col])))
          expect_equal(nrow(r$inffunc), if (route == "rc") nrow(dat) else
                         length(unique(dat$id)))
        } else expect_null(r$inffunc)
      }
    }
  }
})
