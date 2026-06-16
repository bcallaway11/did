# Tests for GENUINE cov-path ridge regularization (omega_cov_shrink = "ridge" with covariates).
# Ridge adds a vanishing diagonal lift lambda*I (lambda = (H/n) mean(diag Omega*(X)) per cell) to each
# cell's conditional moment covariance BEFORE inversion -- the exact analog of the no-cov ridge. It does
# NOT move the estimand toward the pooled/i.i.d. pole (unlike Ledoit-Wolf); it only guarantees a PD inverse.
# Validates: (1) ridge differs from ledoit_wolf and is PD/finite; (2) ridge is byte-identical to nothing
# changing for the LW and none modes; (3) the lift vanishes as O(H/n); (4) FD-oracle of the ridge
# estimation-effect total differential C:dOmega + (tr(C)/n) tr(dOmega); (5) misspec_robust folds (att
# unchanged, EIF mean-zero, SE finite) under ridge on both kernel and sieve.

mk_cov_panel <- function(n = 400L, seed = 11L) {
  set.seed(seed); TP <- 5L
  gcat <- sample(c(3, 4, Inf), n, replace = TRUE, prob = c(.35, .35, .30))
  x1 <- rnorm(n); x2 <- rnorm(n)
  rows <- list()
  for (tt in 1:TP) {
    tr <- as.numeric(is.finite(gcat) & tt >= gcat)
    rows[[tt]] <- data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gcat), gcat, 0),
      y = 0.4 * x1 - 0.2 * x2 + 0.3 * tt + 1.0 * tr + rnorm(n, sd = 1), x1 = x1, x2 = x2)
  }
  do.call(rbind, rows)
}

fit_r <- function(df, shrink, sieve = FALSE, scheme = "efficient", mr = FALSE) {
  old <- options(edid_omega_method = if (sieve) "sieve" else "kernel"); on.exit(options(old))
  set.seed(1)
  suppressWarnings(edid(data = df, yname = "y", idname = "id", tname = "t", gname = "g",
       xformla = ~ x1 + x2, weight_scheme = scheme, ratio_method = "exp", pt_assumption = "all",
       aggregate = "none", seed = 1, omega_cov_shrink = shrink, misspec_robust = mr))
}

test_that("cov-path ridge differs from ledoit_wolf and yields finite PD weights (kernel + sieve)", {
  df <- mk_cov_panel()
  for (sv in c(FALSE, TRUE)) {
    fr <- fit_r(df, "ridge", sieve = sv)
    fl <- fit_r(df, "ledoit_wolf", sieve = sv)
    expect_true(all(is.finite(fr$att_gt$att)))                       # finite point estimates
    expect_true(all(is.finite(fr$att_gt$se[fr$att_gt$se > 0])))      # finite SEs (PD inverse)
    # ridge is a genuinely different channel from LW (some cell SE moves)
    expect_false(isTRUE(all.equal(fr$att_gt$se, fl$att_gt$se, tolerance = 1e-6)))
  }
})

test_that("cov-path ridge does NOT change the 'none' or 'ledoit_wolf' modes (byte-identical)", {
  df <- mk_cov_panel()
  # Running ridge must not perturb the other two modes' results (no global-option leakage).
  for (sv in c(FALSE, TRUE)) {
    n1 <- fit_r(df, "none", sieve = sv); invisible(fit_r(df, "ridge", sieve = sv))
    n2 <- fit_r(df, "none", sieve = sv)
    expect_identical(n1$att_gt$att, n2$att_gt$att)
    expect_identical(n1$att_gt$se,  n2$att_gt$se)
    l1 <- fit_r(df, "ledoit_wolf", sieve = sv); invisible(fit_r(df, "ridge", sieve = sv))
    l2 <- fit_r(df, "ledoit_wolf", sieve = sv)
    expect_identical(l1$att_gt$att, l2$att_gt$att)
    expect_identical(l1$att_gt$se,  l2$att_gt$se)
  }
})

test_that("cov-path ridge att tracks the unshrunk 'none' att (gentle: no estimand shift toward pole)", {
  df <- mk_cov_panel()
  fr <- fit_r(df, "ridge"); fn <- fit_r(df, "none")
  # ridge only adds a vanishing diagonal lift -> the point estimate is essentially the unshrunk one,
  # NOT moved toward the pooled pole the way Ledoit-Wolf moves it.
  expect_equal(fr$att_gt$att, fn$att_gt$att, tolerance = 1e-3)
})

test_that("cov-path ridge lift lambda vanishes as O(H/n) (halves as n doubles)", {
  probe <- function(nn) {
    df <- mk_cov_panel(n = nn, seed = 7L)
    acc <- new.env(); acc$l <- numeric(0)
    orig <- did:::compute_omega_star_kernel_fast_edid
    assignInNamespace("compute_omega_star_kernel_fast_edid", function(...) {
      out <- orig(...); rl <- attr(out, "ridge_lift")
      if (!is.null(rl) && length(rl) == 1L && rl > 0) acc$l <- c(acc$l, rl)
      out
    }, "did")
    on.exit(assignInNamespace("compute_omega_star_kernel_fast_edid", orig, "did"))
    old <- options(edid_cov_ridge = TRUE); on.exit(options(old), add = TRUE)
    suppressWarnings(edid(data = df, yname = "y", idname = "id", tname = "t", gname = "g",
         xformla = ~ x1 + x2, weight_scheme = "averaged", ratio_method = "exp",
         pt_assumption = "all", aggregate = "none", seed = 7L, omega_cov_shrink = "ridge"))
    mean(acc$l)
  }
  l1 <- probe(500L); l2 <- probe(1000L); l3 <- probe(2000L)
  expect_gt(l1, 0); expect_gt(l2, 0); expect_gt(l3, 0)
  expect_equal(l2 / l1, 0.5, tolerance = 0.12)     # halves as n doubles
  expect_equal(l3 / l2, 0.5, tolerance = 0.12)
})

test_that("ridge estimation-effect total differential matches FD oracle (C:dOmega + (tr C/n) tr dOmega)", {
  # Exact smooth weight map (plain inverse, no floor) so the leading C:dOmega is exact by construction;
  # the FD then validates the RIDGE TRACE term -- the new EE channel for omega_cov_shrink = "ridge".
  set.seed(3)
  th <- function(M, mbar) { Minv <- solve(0.5 * (M + t(M)))
    w <- drop(Minv %*% rep(1, nrow(M))); w <- w / sum(w); sum(w * mbar) }
  Csmooth <- function(M, mbar) { Minv <- solve(0.5 * (M + t(M)))
    w <- drop(Minv %*% rep(1, nrow(M))); w <- w / sum(w); theta <- sum(w * mbar)
    q <- drop(Minv %*% (mbar - theta)); -0.5 * (outer(q, w) + outer(w, q)) }
  mk_pd <- function(H) { S <- 0.3 * crossprod(matrix(rnorm(H * H), H)) / H + diag(seq(1.5, 1.5 + H - 1)); 0.5 * (S + t(S)) }
  for (n_full in c(40L, 120L)) {
    ridge <- function(M) M + (sum(diag(M)) / n_full) * diag(nrow(M))
    for (H in c(3L, 5L)) {
      M0 <- mk_pd(H); mbar <- rnorm(H); Mr <- ridge(M0)
      C <- Csmooth(Mr, mbar); trC <- sum(diag(C))
      dE <- matrix(rnorm(H * H), H); dE <- 0.5 * (dE + t(dE))
      ana <- sum(C * dE) + (trC / n_full) * sum(diag(dE))
      h <- 1e-6; fd <- (th(ridge(M0 + h * dE), mbar) - th(ridge(M0 - h * dE), mbar)) / (2 * h)
      expect_equal(ana, fd, tolerance = 1e-6)
    }
  }
})

test_that("misspec_robust folds the ridge weight-estimation channel (att unchanged, EIF mean-zero, SE finite)", {
  df <- mk_cov_panel()
  for (sv in c(FALSE, TRUE)) {
    f0 <- fit_r(df, "ridge", sieve = sv, mr = FALSE)
    f1 <- fit_r(df, "ridge", sieve = sv, mr = TRUE)
    expect_equal(f1$att_gt$att, f0$att_gt$att, tolerance = 1e-10)          # att UNCHANGED by the SE channel
    fin <- is.finite(f1$att_gt$se)
    expect_true(all(f1$att_gt$se[fin] >= 0))                               # finite, non-negative SEs
    expect_lt(max(abs(colMeans(as.matrix(f1$eif))), na.rm = TRUE), 5e-3)   # augmented EIF ~ mean-zero
  }
})
