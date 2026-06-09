# Tests for the sieve Omega smoother (options(edid_omega_method = "sieve")) and, in particular, the fix for
# the silent mixed-estimator: the weight-estimation channel is NOT implemented for the sieve (its psi pass
# would rebuild Omega with the kernel), so under the default misspec_robust = TRUE it must be DISABLED with a
# warning and the plug-in efficient SE reported for it -- never silently combined into a sieve-weights +
# kernel-Sigma_Omega mix. estimation_effect (ACH) and higher_order do not rebuild Omega, so they stay active.

data(mpdta, package = "did")

test_that("the sieve smoother runs and yields a mean-zero EIF with finite, positive SEs", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve"); on.exit(options(old))
  f <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
            gname = "first.treat", xformla = ~ lpop, misspec_robust = FALSE, aggregate = "none")
  expect_true(is.matrix(f$eif))                          # the influence-function slot is $eif (not $inffunc)
  expect_lt(max(abs(colMeans(f$eif))), 1e-8)             # mean-zero per cell
  expect_true(all(is.finite(f$att_gt$se)) && all(f$att_gt$se > 0))
})

test_that("sieve + misspec_robust = TRUE warns and disables ONLY the weight channel (ACH / higher_order stay on)", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve"); on.exit(options(old))
  # The weight channel is not implemented for the sieve -> must WARN (not silently mix in a kernel Sigma_Omega).
  w <- testthat::capture_warnings(
    f_mr <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, misspec_robust = TRUE, aggregate = "none"))
  expect_true(any(grepl("weight-estimation channel is not implemented", w)))
  expect_true(all(is.finite(f_mr$att_gt$se)))
  # Disabling ONLY the weight channel means sieve + misspec_robust = TRUE must be IDENTICAL to sieve with the
  # bundle's other two effects on and the weight channel off: estimation_effect + higher_order, mr = FALSE.
  f_eq <- suppressWarnings(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, misspec_robust = FALSE,
                 estimation_effect = TRUE, higher_order = TRUE, aggregate = "none"))
  expect_equal(f_mr$att_gt$att, f_eq$att_gt$att, tolerance = 1e-10)
  expect_equal(f_mr$att_gt$se,  f_eq$att_gt$se,  tolerance = 1e-10)
})
