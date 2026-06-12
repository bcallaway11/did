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
  # att + se RE-PINNED for the coherent-nuisance + pooled-scale-floor round (2026-06, the with-X
  # cross-cohort repair): (i) ratio_method = "coherent" replaces the per-pair LS cross-cohort ratio and
  # finite-cohort inverse-propensity sieves with the multinomial-logit system (mpdta has cross pairs in
  # every cohort-2006/2007 cell, so their moments and Omega prefactors move); (ii) the pooled/pointwise
  # eigen floors act on the pooled-diagonal scale (the d-dependent raw-scale cap erased the per-moment
  # variance ordering); (iii) STRUCTURALLY DEGENERATE moments (the tpre == t self pair of PRE-treatment
  # cells, an exact-zero moment) now get weight 0 (pinv-style exclusion) instead of diluting pre-cell
  # placebos toward 0 -- the pre cells (5, 6, 9, 10, 11) move the most for that reason.
  # (Earlier re-pin notes remain in force: cell-common overlap-trim estimand; the eigen-floor-aware
  # Daleckii-Krein psi coupling.)
  #
  # SE (only) RE-PINNED again for the coherent FULL estimation-effect round (2026-06): the
  # cross-cohort ratio and finite-cohort inverse-propensity entries now carry complete
  # joint-stacked-system M-estimator aux (no fallback-marking), so the default misspec_robust
  # bundle's ACH / higher-order / inv-p channels cover them -- SEs move by up to ~3e-4
  # (FD-oracled in test-edid-coherent-ee.R). The ATT estimates are BYTE-IDENTICAL (the
  # corrections enter the influence function only), so golden_att is unchanged -- asserted here.
  golden_att <- c(-0.022135707672074, -0.083097479920636, -0.144928755380812, -0.113505421220554,
                   0.002630097380151, -0.008934558404687, -0.004837910447134, -0.047712053995789,
                   0.016118304720655,  0.010700359324765, -0.026537068859659, -0.041435874184193)
  golden_se  <- c(0.02172263305400, 0.02838394621347, 0.03450642588284, 0.03290100269279,
                  0.01187899032301, 0.01975408280814, 0.01835050209244, 0.02013058035024,
                  0.00991035133013, 0.01168187818878, 0.01790787579469, 0.01417812336181)
  expect_equal(fk$att_gt$att, golden_att, tolerance = 1e-7)
  expect_equal(fk$att_gt$se,  golden_se,  tolerance = 1e-7)
})
