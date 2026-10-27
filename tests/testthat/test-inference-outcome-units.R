.inference_units_data <- function() {
  set.seed(730)
  n <- 120L
  d <- expand.grid(id = seq_len(n), period = seq_len(4L))
  d$G <- rep(rep(c(0, 3, 4), each = n / 3), 4L)
  d$cluster <- (d$id - 1L) %% 20L + 1L
  d$X <- stats::rnorm(n)[d$id]
  d$weight <- exp(0.2 * d$X)
  shock <- matrix(stats::rnorm(20L * 4L), 20L, 4L)
  d$Y <- stats::rnorm(n)[d$id] + 0.3 * d$period + 0.2 * d$X * d$period +
    shock[cbind(d$cluster, d$period)] + stats::rnorm(nrow(d)) +
    (d$G > 0 & d$period >= d$G)
  d
}

.inference_units_fit <- function(d, fm, bootstrap, clustered = FALSE, ...) {
  set.seed(731)
  suppressWarnings(suppressMessages(att_gt(
    "Y", "period", "id", "G", data = d, xformla = ~X,
    weightsname = "weight", faster_mode = fm, bstrap = bootstrap,
    cband = bootstrap, biters = 199L,
    clustervars = if (clustered) "cluster" else NULL, ...
  )))
}

.inference_units_aggregate <- function(mp, type) {
  set.seed(732)
  suppressWarnings(aggte(mp, type = type))
}

test_that("outcome units preserve SEs, Wald tests, and simultaneous bands", {
  d <- .inference_units_data()
  for (fm in c(TRUE, FALSE)) for (bootstrap in c(FALSE, TRUE)) {
    for (clustered in c(FALSE, TRUE)) {
      ref <- .inference_units_fit(d, fm, bootstrap, clustered)
      expect_true(all(is.finite(ref$se) & ref$se > 0))
      if (!bootstrap || !clustered) expect_false(is.null(ref$W))
      if (bootstrap) {
        set.seed(731)
        direct <- mboot(ref$inffunc, ref$DIDparams, return_V = FALSE)
        expect_equal(ref$c, direct$crit.val, tolerance = 1e-12)
      }
      types <- c("simple", "dynamic", "group", "calendar")
      ref_agg <- lapply(types, function(type) .inference_units_aggregate(ref, type))

      for (scale in c(1e-12, 1e-6, 1e6)) {
        scaled_data <- d
        scaled_data$Y <- d$Y * scale
        out <- .inference_units_fit(scaled_data, fm, bootstrap, clustered)
        expect_equal(out$att / scale, ref$att, tolerance = 1e-9)
        expect_equal(out$se / scale, ref$se, tolerance = 1e-9)
        expect_equal(out$c, ref$c, tolerance = 1e-9)
        expect_equal(out$W, ref$W, tolerance = 1e-9)
        expect_equal(out$Wpval, ref$Wpval)
        expect_equal(out$DIDparams$cband, ref$DIDparams$cband)

        for (j in seq_along(types)) {
          a <- .inference_units_aggregate(out, types[j])
          expect_true(is.finite(a$overall.se) && a$overall.se > 0)
          expect_equal(a$overall.att / scale, ref_agg[[j]]$overall.att,
                       tolerance = 1e-9)
          expect_equal(a$overall.se / scale, ref_agg[[j]]$overall.se,
                       tolerance = 1e-9)
          if (types[j] != "simple") {
            expect_equal(a$att.egt / scale, ref_agg[[j]]$att.egt,
                         tolerance = 1e-9)
            expect_equal(a$se.egt / scale, ref_agg[[j]]$se.egt,
                         tolerance = 1e-9)
            expect_equal(a$crit.val.egt, ref_agg[[j]]$crit.val.egt,
                         tolerance = 1e-9)
          }
        }
      }
    }
  }
})

test_that("one low-variance cohort is not filtered relative to other cohorts", {
  d <- .inference_units_data()
  set.seed(733)
  noise_scale <- ifelse(d$G == 4, 1, 1e-11)
  d$Y <- stats::rnorm(120L)[d$id] + noise_scale * stats::rnorm(nrow(d)) +
    (d$G > 0 & d$period >= d$G)
  for (fm in c(TRUE, FALSE)) {
    out <- .inference_units_fit(d, fm, FALSE)
    small <- out$se[out$group == 3]
    large <- out$se[out$group == 4]
    expect_true(all(is.finite(small) & small > 0))
    expect_lt(max(small), 1e-10)
    expect_gt(min(large), 1e-2)
  }
})

test_that("exactly zero influence functions remain unavailable in panels and RC", {
  d <- .inference_units_data()
  d$Y <- 0
  zero_estimator <- function(D, ...) {
    list(ATT = 0, att = 0, att.inf.func = rep(0, length(D)))
  }
  for (fm in c(TRUE, FALSE)) for (panel in c(TRUE, FALSE)) {
    dd <- d
    if (!panel) dd$id <- seq_len(nrow(dd))
    for (method in list("dr", "reg", "ipw", zero_estimator)) {
      out <- .inference_units_fit(dd, fm, TRUE, panel = panel, est_method = method)
      expect_true(all(out$att == 0))
      expect_true(all(as.matrix(out$inffunc) == 0))
      expect_true(all(is.na(out$se)))
      expect_null(out$W)
      expect_false(out$DIDparams$cband)
      expect_equal(out$c, stats::qnorm(0.975))
      for (type in c("simple", "dynamic", "group", "calendar")) {
        a <- suppressWarnings(aggte(out, type = type, cband = TRUE))
        expect_equal(a$overall.att, 0)
        expect_true(is.na(a$overall.se))
        if (type != "simple") {
          expect_true(all(is.na(a$se.egt)))
          expect_false(a$DIDparams$cband)
          expect_equal(a$crit.val.egt, stats::qnorm(0.975))
        }
      }
    }
  }
})

test_that("non-finite custom influence functions do not produce finite inference", {
  d <- .inference_units_data()
  invalid_if <- function(D, ...) {
    list(ATT = 0, att = 0, att.inf.func = rep(Inf, length(D)))
  }
  for (fm in c(TRUE, FALSE)) for (bootstrap in c(FALSE, TRUE)) {
    out <- .inference_units_fit(d, fm, bootstrap, est_method = invalid_if)
    expect_true(all(is.na(out$se)))
    expect_null(out$W)
    expect_true(is.finite(out$c))
    if (bootstrap) expect_false(out$DIDparams$cband)
    a <- suppressWarnings(aggte(out, type = "simple"))
    expect_true(is.na(a$overall.se))
  }
})

test_that("zero-IQR bootstrap draws are not reported as zero standard errors", {
  dp <- list(idname = "id", clustervars = NULL, biters = 1L,
             tname = "period", alp = 0.05, panel = TRUE,
             allow_unbalanced_panel = FALSE, faster_mode = FALSE)
  set.seed(734)
  out <- mboot(matrix(stats::rnorm(20L), ncol = 1L), dp)
  expect_true(is.na(out$se))
  expect_true(is.na(out$crit.val))
})

test_that("a single bootstrap draw retains columns and cannot estimate their scales", {
  dp <- list(idname = "id", clustervars = NULL, biters = 1L,
             tname = "period", alp = 0.05, panel = TRUE,
             allow_unbalanced_panel = FALSE, faster_mode = FALSE)
  x <- cbind(seq_len(20), 2 * seq_len(20), 3 * seq_len(20), 0, Inf, NA_real_)
  set.seed(834)
  out <- mboot(x, dp)
  expect_identical(dim(out$bres), c(1L, 3L))
  expect_length(out$se, ncol(x))
  expect_true(all(is.na(out$se)))
  expect_true(is.na(out$crit.val))

  d <- .inference_units_data()
  for (fm in c(TRUE, FALSE)) {
    set.seed(835)
    out <- suppressWarnings(att_gt(
      "Y", "period", "id", "G", data = d, xformla = ~X,
      faster_mode = fm, bstrap = TRUE, cband = TRUE, biters = 1L
    ))
    expect_true(all(is.na(out$se)))
    expect_false(out$DIDparams$cband)
    expect_equal(out$c, stats::qnorm(0.975))
    # The analytical covariance still supports the separate Wald pre-test.
    expect_true(is.finite(as.numeric(out$W)))
  }
})
