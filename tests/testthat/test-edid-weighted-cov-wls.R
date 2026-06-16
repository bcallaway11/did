library(testthat)

# ============================================================
# Phase 1 gate for the weighted-covariate rollout: the sieve /
# Riesz nuisance fits gained an obs-weight (`weights=`) slot.
#
# Two invariants, asserted at the FITTER level (the end-to-end
# weighted-covariate path is still guarded off in edid() until a
# later phase, so these call the nuisance estimators directly):
#
#   (A) WLS(w == 1) is byte-identical to OLS (weights = NULL),
#       including the full M-estimator aux (score_mat, H_inv).
#       This is the dominant regression guard: the new code path
#       must collapse to the verbatim unweighted computation.
#   (B) WLS(w == c) for a CONSTANT c != 1 leaves the fitted
#       nuisance (predictions) invariant -- the loss scales by c,
#       so the minimizer is unchanged. (The aux scales with c and
#       is not asserted equal here; (A) pins the w==1 case.)
# ============================================================

make_wls_panel <- function(n = 240, seed = 7) {
  set.seed(seed)
  x1   <- runif(n, -2, 2)
  X    <- matrix(x1, ncol = 1)
  eta  <- 1.1 * x1 + 0.7 * x1^2 - 0.5
  P    <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  G    <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  y    <- rnorm(n, 0.5 * x1, 1)
  list(X = X, G = G, y = y, n = n)
}

# Compact comparison helpers ------------------------------------------------
expect_byte_identical <- function(a, b, info = NULL) {
  expect_equal(a, b, tolerance = 1e-12, info = info)
}

test_that("propensity ratio (direct LS): WLS(w==1)==OLS incl. aux; WLS(w==c) pred-invariant", {
  d <- make_wls_panel()
  args0 <- list(X_train = d$X, G_train = d$G, X_test = d$X, g = 2, gp = Inf,
                bs_df = 4L, return_aux = TRUE)
  f0 <- do.call(estimate_propensity_ratio_edid, args0)
  f1 <- do.call(estimate_propensity_ratio_edid, c(args0, list(weights = rep(1, d$n))))
  fc <- do.call(estimate_propensity_ratio_edid, c(args0, list(weights = rep(3, d$n))))

  expect_byte_identical(f1$pred,      f0$pred,      "ratio pred w==1")
  expect_byte_identical(f1$score_mat, f0$score_mat, "ratio score w==1")
  expect_byte_identical(f1$H_inv,     f0$H_inv,     "ratio H_inv w==1")
  expect_equal(fc$pred, f0$pred, tolerance = 1e-8)            # scale-invariance of the fit
})

test_that("inverse propensity (direct LS): WLS(w==1)==OLS incl. aux; WLS(w==c) pred-invariant", {
  d <- make_wls_panel()
  args0 <- list(X_train = d$X, G_train = d$G, X_test = d$X, gp = 2,
                bs_df = 4L, return_aux = TRUE)
  f0 <- do.call(estimate_inverse_propensity_edid, args0)
  f1 <- do.call(estimate_inverse_propensity_edid, c(args0, list(weights = rep(1, d$n))))
  fc <- do.call(estimate_inverse_propensity_edid, c(args0, list(weights = rep(3, d$n))))

  expect_byte_identical(f1$s_hat,     f0$s_hat,     "invp s_hat w==1")
  expect_byte_identical(f1$score_mat, f0$score_mat, "invp score w==1")
  expect_byte_identical(f1$H_inv,     f0$H_inv,     "invp H_inv w==1")
  expect_equal(fc$s_hat, f0$s_hat, tolerance = 1e-8)
})

test_that("conditional mean (OLS): WLS(w==1)==OLS incl. aux; WLS(w==c) pred-invariant", {
  d <- make_wls_panel()
  args0 <- list(X_train = d$X, Y_delta_train = d$y, G_train = d$G, X_test = d$X,
                gp = 2, bs_df = 4L, return_aux = TRUE)
  f0 <- do.call(estimate_conditional_mean_edid, args0)
  f1 <- do.call(estimate_conditional_mean_edid, c(args0, list(weights = rep(1, d$n))))
  fc <- do.call(estimate_conditional_mean_edid, c(args0, list(weights = rep(3, d$n))))

  expect_byte_identical(f1$pred,      f0$pred,      "cmean pred w==1")
  expect_byte_identical(f1$score_mat, f0$score_mat, "cmean score w==1")
  expect_byte_identical(f1$H_inv,     f0$H_inv,     "cmean H_inv w==1")
  expect_equal(fc$pred, f0$pred, tolerance = 1e-8)
})

test_that("propensity ratio (exp Riesz): WLS(w==1)==OLS incl. aux; WLS(w==c) pred-invariant", {
  d <- make_wls_panel()
  # finite comparison cohort -> exercises the exponential-link Riesz solver
  args0 <- list(X_train = d$X, G_train = d$G, X_test = d$X, g = 2, gp = 3,
                bs_df = 4L, return_aux = TRUE)
  f0 <- do.call(estimate_propensity_ratio_exp_edid, args0)
  f1 <- do.call(estimate_propensity_ratio_exp_edid, c(args0, list(weights = rep(1, d$n))))
  fc <- do.call(estimate_propensity_ratio_exp_edid, c(args0, list(weights = rep(3, d$n))))

  expect_byte_identical(f1$pred,      f0$pred,      "exp ratio pred w==1")
  expect_byte_identical(f1$beta,      f0$beta,      "exp ratio beta w==1")
  expect_byte_identical(f1$score_mat, f0$score_mat, "exp ratio score w==1")
  expect_byte_identical(f1$H_inv,     f0$H_inv,     "exp ratio H_inv w==1")
  expect_equal(fc$pred, f0$pred, tolerance = 1e-6)            # tailored loss scales by c -> same minimizer
})

test_that("inverse propensity (exp Riesz): WLS(w==1)==OLS incl. aux; WLS(w==c) pred-invariant", {
  d <- make_wls_panel()
  args0 <- list(X_train = d$X, G_train = d$G, X_test = d$X, gp = 3,
                bs_df = 4L, return_aux = TRUE)
  f0 <- do.call(estimate_inverse_propensity_exp_edid, args0)
  f1 <- do.call(estimate_inverse_propensity_exp_edid, c(args0, list(weights = rep(1, d$n))))
  fc <- do.call(estimate_inverse_propensity_exp_edid, c(args0, list(weights = rep(3, d$n))))

  expect_byte_identical(f1$s_hat,     f0$s_hat,     "exp invp s_hat w==1")
  expect_byte_identical(f1$beta,      f0$beta,      "exp invp beta w==1")
  expect_byte_identical(f1$score_mat, f0$score_mat, "exp invp score w==1")
  expect_byte_identical(f1$H_inv,     f0$H_inv,     "exp invp H_inv w==1")
  expect_equal(fc$s_hat, f0$s_hat, tolerance = 1e-6)
})

test_that("a genuinely dispersed weight moves the fit (sanity: weights are not inert)", {
  d  <- make_wls_panel()
  set.seed(99)
  wd <- exp(0.8 * rnorm(d$n)); wd <- wd / mean(wd)            # dispersed, mean-1
  f0 <- estimate_conditional_mean_edid(d$X, d$y, d$G, d$X, gp = 2, return_aux = TRUE)
  fw <- estimate_conditional_mean_edid(d$X, d$y, d$G, d$X, gp = 2, return_aux = TRUE, weights = wd)
  expect_false(isTRUE(all.equal(fw$pred, f0$pred, tolerance = 1e-6)))
})
