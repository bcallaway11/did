library(testthat)

# ============================================================
# Phase 3 gate for the weighted-covariate rollout: end-to-end
# edid() on the COVARIATE path with observation weights.
#
# Validates the obs-weighted plug-in moment / EIF (Hajek) and the
# obs-weighted ACH estimation-effect correction, now that the
# weightsname x xformla guard is relaxed.
#
#   (A) A CONSTANT weight column normalizes to all-ones, so the entire
#       weighted-covariate pipeline (weighted nuisance WLS + weighted
#       Omega*(X) + Hajek plug-in + obs-weighted ACH) must collapse to
#       the unweighted covariate fit -- end-to-end, including SEs.
#   (B) Under a GENUINELY dispersed weight, the analytic ACH reproduces
#       the finite-difference oracle (the obs-weighted moment's exact
#       nuisance sensitivity) to 1e-5 -- the EE channel is correct.
#   (C) The ACH is a real, non-negligible correction here.
#   (D) Dispersed weights MOVE the estimate/SE (weights are not inert).
# ============================================================

make_wcov_panel <- function(n = 320, seed = 21) {
  set.seed(seed)
  Tn  <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P   <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  wu   <- exp(0.8 * rnorm(n)); wu <- wu / mean(wu)         # dispersed, mean-1 (Kish n_eff << n)
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

fit_se <- function(df, wt = NULL, ee = TRUE, ach = "analytic", misspec = FALSE) {
  op <- options(edid_ach = ach); on.exit(options(op), add = TRUE)
  suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient",
       aggregate = "none", bstrap = FALSE, seed = 1L, misspec_robust = misspec,
       estimation_effect = ee, weightsname = wt))
}

test_that("constant weight column => weighted-cov fit == unweighted-cov fit (plug-in + ACH)", {
  df <- make_wcov_panel(); df$w <- 5                        # constant column -> unit_weights == 1
  for (ee in c(FALSE, TRUE)) {
    f_unw <- fit_se(df, wt = NULL, ee = ee)
    f_w1  <- fit_se(df, wt = "w",  ee = ee)
    ok <- is.finite(f_unw$att_gt$att) & is.finite(f_w1$att_gt$att)
    expect_equal(f_w1$att_gt$att[ok], f_unw$att_gt$att[ok], tolerance = 1e-8,
                 info = sprintf("att, estimation_effect=%s", ee))
    expect_equal(f_w1$att_gt$se[ok],  f_unw$att_gt$se[ok],  tolerance = 1e-8,
                 info = sprintf("se, estimation_effect=%s", ee))
  }
})

test_that("constant weight column => weighted-cov == unweighted-cov under misspec_robust", {
  df <- make_wcov_panel(); df$w <- 5
  f_unw <- fit_se(df, wt = NULL, ee = TRUE, misspec = TRUE)
  f_w1  <- fit_se(df, wt = "w",  ee = TRUE, misspec = TRUE)
  ok <- is.finite(f_unw$att_gt$se) & is.finite(f_w1$att_gt$se)
  expect_equal(f_w1$att_gt$se[ok], f_unw$att_gt$se[ok], tolerance = 1e-8)
})

test_that("weighted-covariate ACH reproduces the finite-difference oracle (analytic == FD, 1e-5)", {
  df <- make_wcov_panel(seed = 7)
  se_noee <- fit_se(df, wt = "w", ee = FALSE, ach = "analytic")$att_gt$se   # weighted plug-in SE
  se_an   <- fit_se(df, wt = "w", ee = TRUE,  ach = "analytic")$att_gt$se   # weighted ACH, analytic Gamma
  se_fd   <- fit_se(df, wt = "w", ee = TRUE,  ach = "fd")$att_gt$se         # weighted ACH, FD oracle
  ok <- is.finite(se_noee) & is.finite(se_an) & is.finite(se_fd)
  skip_if(sum(ok) < 2L, "too few non-degenerate weighted cells")
  expect_gt(max(abs(se_an[ok] - se_noee[ok])), 1e-3)   # (C) ACH is a real correction under weights
  expect_equal(se_an[ok], se_fd[ok], tolerance = 1e-5) # (B) analytic reproduces the weighted FD oracle
})

test_that("dispersed observation weights move the covariate-path estimate (not inert)", {
  df <- make_wcov_panel(seed = 3)
  f_unw <- fit_se(df, wt = NULL, ee = TRUE)
  f_w   <- fit_se(df, wt = "w",  ee = TRUE)
  ok <- is.finite(f_unw$att_gt$att) & is.finite(f_w$att_gt$att)
  expect_false(isTRUE(all.equal(f_w$att_gt$att[ok], f_unw$att_gt$att[ok], tolerance = 1e-4)))
})
