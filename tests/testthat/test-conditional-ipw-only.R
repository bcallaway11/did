test_that("conditional pre-test rejects unsupported methods before processing data", {
  methods <- list("reg", "dr", function(...) stop("must not be called"),
                  "unknown", c("ipw", "dr"), NA_character_, NULL)
  for (method in methods) {
    expect_error(conditional_did_pretest(
      yname = "Y", tname = "period", idname = "id", gname = "G",
      data = stop("must not process data for an unsupported method"),
      est_method = method
    ), 'only supports est_method = "ipw"', fixed = TRUE)
  }
})

test_that("correct original outcome regression cannot validate the legacy regression test", {
  # Conditional trends are exactly X in both groups. A correct linear model
  # for that change becomes misspecified after multiplication by 1{X <= 0}.
  x <- c(-1, 0, 1)
  controls <- c(20, 80, 20)
  treated <- c(80, 20, 80)
  X <- cbind(1, x)
  beta <- stats::lm.wfit(X, x, controls)$coefficients
  h <- as.numeric(x <= 0)
  original_residual_moment <- stats::weighted.mean(h * (x - X %*% beta), treated)
  transformed_beta <- stats::lm.wfit(X, h * x, controls)$coefficients
  transformed_moment <- stats::weighted.mean(h * x - X %*% transformed_beta, treated)
  expect_equal(as.numeric(beta), c(0, 1), tolerance = 1e-12)
  expect_equal(original_residual_moment, 0, tolerance = 1e-12)
  expect_equal(transformed_moment, -5 / 18, tolerance = 1e-12)

  units <- data.frame(X = rep(x, each = 100),
    G = c(rep(4, 80), rep(0, 20), rep(4, 20), rep(0, 80),
          rep(4, 80), rep(0, 20)), id = seq_len(300))
  d <- merge(units, data.frame(period = 1:4))
  d$Y <- d$period * d$X + (d$G == 4 & d$period >= 4)
  for (method in c("reg", "dr")) {
    # Those methods remain supported by att_gt(), where the ORIGINAL linear
    # outcome model correctly yields zero pre-effects and a unit post-effect.
    fit <- suppressWarnings(att_gt(yname = "Y", tname = "period",
      idname = "id", gname = "G", data = d, xformla = ~X,
      est_method = method, bstrap = FALSE, cband = FALSE))
    expect_equal(fit$att, c(0, 0, 1), tolerance = 1e-12)
    expect_error(conditional_did_pretest(yname = "Y", tname = "period",
      idname = "id", gname = "G", data = d, xformla = ~X,
      est_method = method), 'only supports est_method = "ipw"', fixed = TRUE)
  }
})

test_that("default and explicit IPW conditional pre-tests remain identical", {
  set.seed(490)
  d <- expand.grid(id = seq_len(45), period = seq_len(4))
  d$G <- rep(c(0, 3, 4), each = 15)[d$id]
  d$X <- rnorm(45)[d$id]
  d$Y <- d$period + d$X + rnorm(nrow(d)) +
    (d$G > 0 & d$period >= d$G)
  run <- function(panel, method = NULL) {
    set.seed(491)
    args <- list(yname = "Y", tname = "period", gname = "G", data = d,
      panel = panel, idname = if (panel) "id" else NULL,
      xformla = ~X, biters = 31)
    if (!is.null(method)) args$est_method <- method
    invisible(capture.output(answer <- suppressMessages(
      do.call(conditional_did_pretest, args))))
    answer
  }
  for (panel in c(TRUE, FALSE)) {
    default <- run(panel)
    explicit <- run(panel, "ipw")
    expect_identical(default$CvM, explicit$CvM)
    expect_identical(default$CvMb, explicit$CvMb)
    expect_identical(default$CvMcval, explicit$CvMcval)
    expect_identical(default$CvMpval, explicit$CvMpval)
  }
})

test_that("conditional pre-test uses the resolved panel structure", {
  set.seed(490)
  d <- expand.grid(id = seq_len(45), period = seq_len(4))
  d$G <- rep(c(0, 3, 4), each = 15)[d$id]
  d$X <- rnorm(45)[d$id]
  d$Y <- d$period + d$X + rnorm(nrow(d)) +
    (d$G > 0 & d$period >= d$G)
  run <- function(panel, allow) {
    set.seed(491)
    invisible(capture.output(answer <- suppressMessages(
      conditional_did_pretest(yname = "Y", tname = "period", gname = "G",
        data = d, panel = panel, idname = if (panel) "id" else NULL,
        allow_unbalanced_panel = allow, xformla = ~X, biters = 31))))
    answer
  }
  for (panel in c(TRUE, FALSE)) {
    # A balanced panel resolves TRUE to FALSE; genuine cross sections ignore
    # it. Neither case should be mistaken for an unsupported incomplete panel.
    expected <- run(panel, FALSE)
    actual <- run(panel, TRUE)
    for (field in c("CvM", "CvMb", "CvMcval", "CvMpval")) {
      expect_identical(actual[[field]], expected[[field]])
    }
  }
  d <- d[-1L, ]
  expect_error(suppressMessages(conditional_did_pretest(
    yname = "Y", tname = "period", idname = "id", gname = "G",
    data = d, allow_unbalanced_panel = TRUE, xformla = ~X, biters = 31)),
    "conditional pre-test is not supported for unbalanced panel data")
})
