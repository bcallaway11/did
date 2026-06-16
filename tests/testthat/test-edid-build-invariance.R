# Regression guards for the two biggest landed kernel optimizations, which the rest of the suite does not
# pin tightly (its end-to-end tolerances are statistical, ~0.4-0.6, and would survive a 1e-3..1e-2 drift):
#
#   (1) BUILD-INVARIANCE. The default fast BLAS Nadaraya-Watson Omega build (edid_omega_method = "kernel")
#       is a hand-derived BLAS-3 reformulation of the exact original per-pair build ("kernel_orig"). It is
#       meant to be bit-identical; this asserts att AND se agree, so a future edit to the fast build that
#       silently shifts the kernel cannot pass CI.
#   (2) GOLDEN VALUES. A pinned att/se snapshot on a fixed shipped fixture (mpdta + ~lpop) guards the
#       kp_cache / m_eff micro-optimizations, which BOTH kernel variants share -- so (1) cannot catch a
#       cache regression, only a pinned snapshot can. The values were recorded on the clean optimization
#       HEAD; att and se are RNG-independent (the seed only affects the sup-t band crit, not se).

data(mpdta, package = "did")

fit_mpdta_bi <- function(omega = "kernel", misspec_robust = TRUE) {
  old <- options(edid_omega_method = omega)
  on.exit(options(old))
  set.seed(123)
  edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
       gname = "first.treat", xformla = ~ lpop,
       misspec_robust = misspec_robust, aggregate = "none")
}

test_that("fast BLAS kernel build is invariant to the exact original per-pair build (att + se)", {
  # default path (misspec_robust = TRUE: Hessian + ACH + psi channel all active)
  fk  <- fit_mpdta_bi("kernel")
  fko <- fit_mpdta_bi("kernel_orig")
  expect_equal(fk$att_gt$att, fko$att_gt$att, tolerance = 1e-9)
  expect_equal(fk$att_gt$se,  fko$att_gt$se,  tolerance = 1e-9)
  # plug-in path (misspec_robust = FALSE)
  gk  <- fit_mpdta_bi("kernel",      misspec_robust = FALSE)
  gko <- fit_mpdta_bi("kernel_orig", misspec_robust = FALSE)
  expect_equal(gk$att_gt$att, gko$att_gt$att, tolerance = 1e-9)
  expect_equal(gk$att_gt$se,  gko$att_gt$se,  tolerance = 1e-9)
})

test_that("edid() reproduces the pinned golden att/se on mpdta + ~lpop (guards kp_cache / m_eff)", {
  fk <- fit_mpdta_bi("kernel")
  # att + se RE-PINNED for the GENUINE cov-path ridge round (2026-06-15). The covariate-path DEFAULT
  # omega_cov_shrink = "ridge" previously FELL BACK to "ledoit_wolf" on the covariate path; it now does a
  # genuine ridge (a vanishing diagonal lift lambda*I, lambda = (H/n) mean(diag Omega*(X)) per cell, added
  # to each cell's conditional moment covariance before inversion -- the exact analog of the no-cov ridge).
  # Ridge does NOT move the estimand toward the pooled/i.i.d. pole (unlike Ledoit-Wolf), so the with-X
  # default att/se move from the LW-fallback numbers to the (gentler) ridge numbers; the estimation-effect
  # channel carries the FD-oracled ridge EE (C:dOmega + (tr C/n) tr dOmega). The earlier exp-default re-pin
  # note still holds (exp-link Riesz ratios with full estimation-effect aux). Structural notes remain in
  # force: pooled/pointwise eigen floors on the pooled-diagonal scale; degenerate self-pair moments excluded
  # (pinv-style, weight 0); cell-common overlap-trim estimand; eigen-floor-aware Daleckii-Krein psi. (Prior
  # LW-fallback golden att[1]/se[6] were -0.0228661807787743 / 0.0262026852334394.)
  golden_att <- c(-0.0243245745421897, -0.0847227022585466, -0.147109182022892, -0.112359638295074,
                  -0.0125883974483623, -0.0116946248568905, -0.00648970279305543, -0.0485840675412391,
                   0.00805951068627602, 0.0159797568874081, -0.0269735268359635, -0.0460664656465983)
  golden_se  <- c(0.0222542941868945, 0.0285030682244642, 0.0346143085890438, 0.0325485794739207,
                  0.021198086086665, 0.0240138541201071, 0.023112371360369, 0.0192390338814384,
                  0.00995375824804451, 0.0104266223437504, 0.0175624581131783, 0.0135426516524408)
  expect_equal(fk$att_gt$att, golden_att, tolerance = 1e-7)
  expect_equal(fk$att_gt$se,  golden_se,  tolerance = 1e-7)
})
