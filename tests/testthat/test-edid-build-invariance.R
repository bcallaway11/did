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
  golden_att <- c(-0.021325657732767, -0.081566921128498, -0.144734309317491, -0.111693456073102,
                  -0.000489304301903, -0.006978740079083, -0.005584399663392, -0.048704071009899,
                   0.011343970337863,  0.006903644769504, -0.014167439733459, -0.047166800055300)
  golden_se  <- c(0.02243280187454, 0.02894131655886, 0.03480847736478, 0.03337624223599,
                  0.00849702796138, 0.01250901383451, 0.02021629886772, 0.02171703477515,
                  0.00732151564137, 0.00812716405792, 0.01127193490819, 0.01573256470791)
  expect_equal(fk$att_gt$att, golden_att, tolerance = 1e-7)
  expect_equal(fk$att_gt$se,  golden_se,  tolerance = 1e-7)
})
