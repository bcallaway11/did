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
  # att + se RE-PINNED for the cell-common overlap-trim estimand (2026-06): the default trim_level = 200
  # BINDS on mpdta + ~lpop (max |att(200) - att(Inf)| = 0.0069), and in 8 of the 12 cells the comparison
  # pairs carry DIFFERENT overlap masks, so the common-mask renormalization (every moment masked on the
  # intersection of the surviving pairs' masks -- ONE estimand per cell) moves those cells by 3e-4..7e-3.
  # No pair is dropped here (n_pairs_dropped == 0 in every cell); the 4 unmoved cells are those whose
  # pair masks coincide. trim_level = Inf remains byte-identical to the pre-change build.
  golden_att <- c(-0.020876019051455, -0.081264904800779, -0.144247924115520, -0.110234369322343,
                  -0.000489304301903, -0.006978740079083, -0.005584399663392, -0.048704071009899,
                   0.007680747334852,  0.012341556720145, -0.016242058275798, -0.053993457185514)
  # (Earlier se re-pin note, still in force: the kernel misspec_robust weight channel uses the
  # eigen-floor-aware Daleckii-Krein coupling + (1-lambda) shrinkage factor + Eq.(3.12) Term-1
  # contribution -- the same construction the sieve channel already used.)
  golden_se  <- c(0.02255763635885, 0.02891300644631, 0.03471181464718, 0.03302420914630,
                  0.00888452687748, 0.01250881585930, 0.02013082743995, 0.02140868178139,
                  0.00702085767154, 0.00770758312701, 0.01108837878524, 0.01586968653181)
  expect_equal(fk$att_gt$att, golden_att, tolerance = 1e-7)
  expect_equal(fk$att_gt$se,  golden_se,  tolerance = 1e-7)
})
