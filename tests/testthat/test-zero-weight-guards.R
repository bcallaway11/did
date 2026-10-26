# A zero-weight observation does not constrain a weighted 2x2 estimator's
# overlap or regression design, even when its covariates are extreme.
zero_weight_guard_data <- function() {
  set.seed(50)
  units <- data.frame(id = 1:100, G = rep(c(0, 3), each = 50),
                      X = rnorm(100), w = 1, a = rnorm(100))
  d <- merge(units, data.frame(period = 1:4))
  d$Y <- d$a + d$period * (0.2 + 0.3 * d$X) + rnorm(nrow(d)) +
    2 * (d$G == 3 & d$period >= 3)
  d
}

test_that("zero-weight covariates cannot cause false overlap or rank failures", {
  d <- zero_weight_guard_data()
  for (method in c("reg", "dr", "ipw")) {
    # A zero-weight control spoils the unweighted control Gram matrix; a
    # zero-weight treated observation spoils the preliminary overlap logit.
    extra <- data.frame(id = 101, G = if (method == "reg") 0 else 3,
                        X = if (method == "reg") 1e12 else 1000,
                        w = 0, a = 0, period = 1:4, Y = 0)
    for (structure in c("panel", "varying", "rc", "unbalanced")) {
      original <- d
      if (structure == "unbalanced") original <- original[-seq(17, nrow(d), 17), ]
      augmented <- rbind(original, extra)
      if (structure == "rc") {
        original$id <- seq_len(nrow(original))
        augmented$id <- seq_len(nrow(augmented))
      }
      for (faster in c(TRUE, FALSE)) {
        run <- function(dat) suppressMessages(att_gt(
          yname = "Y", tname = "period", idname = "id", gname = "G",
          xformla = ~X, data = dat, weightsname = "w", est_method = method,
          panel = structure != "rc", allow_unbalanced_panel = structure == "unbalanced",
          fix_weights = if (structure == "varying") "varying" else NULL,
          faster_mode = faster, bstrap = FALSE, cband = FALSE))
        expected <- run(original)
        expect_warning(actual <- run(augmented), NA)
        expect_true(all(is.finite(actual$att)))
        expect_equal(actual$att, expected$att, tolerance = 1e-10,
                     info = paste(method, structure, faster))
      }
    }
  }
})

test_that("zero-weight guard fix also applies without positional precompute or guard cache", {
  d <- zero_weight_guard_data()
  extra <- data.frame(id = 101, G = 3, X = 1000, w = 0,
                      a = 0, period = 1:4, Y = 0)
  withr::local_options(list(did.disable_precompute = TRUE, did.disable_check_cache = TRUE))
  for (faster in c(TRUE, FALSE)) {
    run <- function(dat) att_gt(yname = "Y", tname = "period", idname = "id",
      gname = "G", xformla = ~X, data = dat, weightsname = "w",
      faster_mode = faster, bstrap = FALSE, cband = FALSE)
    expected <- run(d)
    expect_warning(actual <- run(rbind(d, extra)), NA)
    expect_equal(actual$att, expected$att, tolerance = 1e-10)
    # Independently assemble the post-treatment (g=3,t=3) 2x2 estimator.
    z <- rbind(d, extra)
    pre <- z[z$period == 2, ]
    post <- z[z$period == 3, ]
    pre <- pre[order(pre$id), ]
    post <- post[order(post$id), ]
    reference <- DRDID::drdid_panel(y1 = post$Y, y0 = pre$Y,
      D = as.numeric(pre$G == 3), covariates = model.matrix(~X, pre),
      i.weights = pre$w, inffunc = TRUE)
    expect_equal(actual$att[actual$t == 3], reference$ATT, tolerance = 1e-10)
  }
})

test_that("guard support follows normalization and still rejects unidentified designs", {
  X <- cbind(1, c(-1, 0, 1, 1000))
  D <- c(0, 1, 0, 1)
  w <- c(1e200, 1e200, 1e200, 1e-200)
  # The last positive raw weight becomes zero in the estimator's normalization.
  expect_identical(did:::overlap_check_fail(X, D, i.weights = w),
                   did:::overlap_check_fail(X[1:3, ], D[1:3]))
  expect_true(did:::overlap_check_fail(matrix(1, 3, 1), c(0, 1, 1),
                                      intercept_only = TRUE, i.weights = c(0, 1, 1)))
  expect_true(did:::rcond_check_fail(matrix(1, 3, 1), c(0, 1, 1),
                                    intercept_only = TRUE, i.weights = c(0, 1, 1)))
  expect_true(did:::overlap_check_fail(matrix(1, 3, 1), c(0, 1, 1),
                                      intercept_only = TRUE, i.weights = rep(0, 3)))
  expect_true(did:::rcond_check_fail(matrix(1, 3, 1), c(0, 1, 1),
                                    intercept_only = TRUE, i.weights = rep(0, 3)))
})

test_that("regression feasibility uses the actual weighted control design", {
  d <- zero_weight_guard_data()
  extra <- data.frame(id = 101, G = 0, X = 1e12, w = 1e-24,
                      a = 0, period = 1:4, Y = 0)
  z <- rbind(d, extra)
  pre <- z[z$period == 2, ]
  post <- z[z$period == 3, ]
  pre <- pre[order(pre$id), ]
  post <- post[order(post$id), ]
  reference <- DRDID::reg_did_panel(y1 = post$Y, y0 = pre$Y,
    D = as.numeric(pre$G == 3), covariates = model.matrix(~X, pre),
    i.weights = pre$w, inffunc = TRUE)
  for (faster in c(TRUE, FALSE)) {
    expect_warning(actual <- att_gt(yname = "Y", tname = "period",
      idname = "id", gname = "G", xformla = ~X, data = z,
      weightsname = "w", est_method = "reg", faster_mode = faster,
      bstrap = FALSE, cband = FALSE), NA)
    expect_equal(actual$att[actual$t == 3], reference$ATT, tolerance = 1e-10)
  }

  # Conversely, the unweighted control design is full-rank but its second
  # covariate direction has negligible weight in the actual regression.
  X <- cbind(1, c(-1, 1, -1, 1, 0, 0))
  D <- c(0, 0, 0, 0, 1, 1)
  w <- c(1, 1e-20, 1, 1e-20, 1, 1)
  expect_false(did:::rcond_check_fail(X, D))
  expect_true(did:::rcond_check_fail(X, D, i.weights = w))
  for (scale in c(1e-100, 1, 1e100)) {
    expect_identical(did:::rcond_check_fail(X, D, i.weights = w * scale), TRUE)
    expect_identical(did:::rcond_check_fail(X, D, i.weights = rep(scale, 6)), FALSE)
  }
  # The weighted regression itself has the same numerical feasibility limit.
  expect_error(DRDID::reg_did_panel(y1 = 1:6, y0 = rep(0, 6), D = D,
    covariates = X, i.weights = w, inffunc = TRUE),
    "singular|coefficients have NA")
})
