# Overlap/rcond guard short-circuit and per-(g, covariate-period) cache:
# the intercept-only closed forms and the cached guard booleans must reproduce
# the per-cell preliminary-fit decisions, results, and warnings bit-identically.

test_that("intercept-only overlap short-circuit matches the real fit's decision, including the 0.999 knife edge", {
  decisions_match <- function(D) {
    X <- as.matrix(rep(1, length(D)))
    full <- suppressWarnings(max(did:::overlap_logit_fit(X, D)$fitted.values) >= 0.999)
    shortcut <- suppressWarnings(did:::overlap_check_fail(X, D, intercept_only = TRUE))
    identical(full, shortcut)
  }
  # interior, just-below-cutoff, the exact 0.999 knife edge (where the IRLS
  # iterate lands just BELOW the cutoff and the cell must keep passing),
  # just-above-cutoff, and the degenerate no-control / no-treated designs
  expect_true(decisions_match(rep(c(1, 0), c(500, 500))))   # mean = 0.5
  expect_true(decisions_match(rep(c(1, 0), c(998, 2))))     # mean = 0.998
  expect_true(decisions_match(rep(c(1, 0), c(999, 1))))     # mean = 0.999 exactly
  expect_false(did:::overlap_check_fail(as.matrix(rep(1, 1000)),
                                        rep(c(1, 0), c(999, 1)),
                                        intercept_only = TRUE))
  expect_true(decisions_match(rep(c(1, 0), c(9991, 9))))    # mean = 0.9991
  expect_true(decisions_match(rep(c(1, 0), c(1500, 1))))    # mean ~ 0.99933
  expect_true(decisions_match(rep(1, 100)))                 # mean = 1 (no controls)
  expect_true(decisions_match(rep(0, 100)))                 # mean = 0 (no treated)
})

test_that("intercept-only rcond short-circuit matches the Gram-matrix test", {
  X <- as.matrix(rep(1, 50))
  ref <- function(D) {
    cc <- X[D == 0, , drop = FALSE]
    rcond(crossprod(cc)) < .Machine$double.eps
  }
  D_some_controls <- rep(c(1, 0), c(40, 10))
  D_no_controls <- rep(1, 50)
  expect_identical(did:::rcond_check_fail(X, D_some_controls, intercept_only = TRUE),
                   ref(D_some_controls))
  expect_identical(did:::rcond_check_fail(X, D_no_controls, intercept_only = TRUE),
                   ref(D_no_controls))
  expect_false(did:::rcond_check_fail(X, D_some_controls, intercept_only = TRUE))
  expect_true(did:::rcond_check_fail(X, D_no_controls, intercept_only = TRUE))
})

test_that("guard cache returns bit-identical results to per-cell checks (fast and slow paths)", {
  set.seed(20260610)
  sp <- did::reset.sim(time.periods = 5, n = 500)
  d <- did::build_sim_dataset(sp)
  run_one <- function(faster, disable_cache) {
    old <- options(did.disable_check_cache = disable_cache)
    on.exit(options(old), add = TRUE)
    suppressWarnings(suppressMessages(att_gt(yname = "Y", xformla = ~X, data = d,
      tname = "period", idname = "id", gname = "G", est_method = "dr",
      control_group = "nevertreated", faster_mode = faster, bstrap = FALSE)))
  }
  for (faster in c(TRUE, FALSE)) {
    cached <- run_one(faster, disable_cache = FALSE)
    uncached <- run_one(faster, disable_cache = TRUE)
    expect_identical(cached$att, uncached$att)
    expect_identical(cached$se, uncached$se)
    expect_identical(as.matrix(cached$inffunc), as.matrix(uncached$inffunc))
  }
})

test_that("cached guard failures still warn once per cell with unchanged text", {
  set.seed(20260610)
  d <- expand.grid(id = 1:200, period = 1:4)
  gvals <- c(0, 2, 3)
  d$G <- gvals[(d$id - 1) %% 3 + 1]
  d$Xsep <- 1 * (d$G > 0)  # separates treated from never-treated: overlap fails in all 6 cells
  d$Y <- 0.1 * d$period + (d$G > 0 & d$period >= d$G) + rnorm(nrow(d), 0, 0.5)
  run_one <- function(faster, disable_cache) {
    old <- options(did.disable_check_cache = disable_cache)
    on.exit(options(old), add = TRUE)
    capture_warnings(suppressMessages(att_gt(yname = "Y", xformla = ~Xsep,
      data = d, tname = "period", idname = "id", gname = "G", est_method = "dr",
      faster_mode = faster, bstrap = FALSE)))
  }
  # cache on/off must produce identical warning vectors (count, text, order)
  w_fast <- run_one(TRUE, disable_cache = FALSE)
  expect_identical(w_fast, run_one(TRUE, disable_cache = TRUE))
  w_slow <- run_one(FALSE, disable_cache = FALSE)
  expect_identical(w_slow, run_one(FALSE, disable_cache = TRUE))
  expect_identical(w_fast, w_slow)
  # one overlap warning per failed cell, even when the boolean came from the cache
  expect_identical(sum(grepl("overlap condition violated for group", w_fast)), 6L)
})
