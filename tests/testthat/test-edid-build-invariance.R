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
  # att + se RE-PINNED for the exp-default + coherent-removal round (2026-06-12). The covariate-path
  # default flipped from ratio_method = "coherent" (now removed; see ratio_method_thinshares_mc.md --
  # coherent hard-failed ~6.7% of thin-share draws and undercovered) to ratio_method = "exp": every
  # cross-cohort ratio r_{g,g'} AND the never-treated ratio r_{g,Inf} are now per-target exp-link Riesz
  # fits (mpdta has cross pairs in every cohort-2006/2007 cell, and r_{g,Inf} enters every cell, so
  # both att and se move). The exp engine carries FULL estimation-effect aux (FD-oracled in
  # test-edid-exp-ratio.R), so the default misspec_robust bundle's ACH / higher-order / inv-p channels
  # cover the cross-cohort channels. Earlier structural re-pin notes remain in force: pooled/pointwise
  # eigen floors on the pooled-diagonal scale; structurally degenerate self-pair moments excluded
  # (pinv-style, weight 0); cell-common overlap-trim estimand; eigen-floor-aware Daleckii-Krein psi.
  golden_att <- c(-0.0228661807787743, -0.084640492698104, -0.147602286243642, -0.113826874446173,
                  -0.0124328665332939, -0.011218170244851, -0.00644463130881531, -0.0407414526270808,
                   0.00822268853901105, 0.0160788895431695, -0.0270681259670889, -0.0461158513996921)
  golden_se  <- c(0.0224055257979494, 0.0283967569471087, 0.0345608877673499, 0.0325131158273884,
                  0.0214105300899993, 0.0262026852334394, 0.0232546409208079, 0.0208387910133422,
                  0.00994222279586637, 0.0104250894785987, 0.0175349826093005, 0.0135277283813801)
  expect_equal(fk$att_gt$att, golden_att, tolerance = 1e-7)
  expect_equal(fk$att_gt$se,  golden_se,  tolerance = 1e-7)
})
