# Tests for the sieve Omega smoother (options(edid_omega_method = "sieve")) and its weight-estimation channel.
# EFFICIENT + sieve + misspec_robust now runs the proper sieve Sigma_Omega (OLS-projection IF with the
# eigen-floor-aware Daleckii-Krein coupling; validated to nominal coverage + jackknife slope ~1). AVERAGED +
# sieve degrades (the pooled-Omega-bar psi is numerically unstable). estimation_effect / higher_order are
# smoother-agnostic (frozen weights) and stay active under either scheme.

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

test_that("sieve EFFICIENT + misspec_robust runs the weight channel (no warning, mean-zero EIF, finite SEs)", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve"); on.exit(options(old))
  w <- testthat::capture_warnings(
    f_mr <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "efficient",
                 misspec_robust = TRUE, aggregate = "none"))
  # the efficient sieve weight channel IS implemented -> NO "not implemented" warning
  expect_false(any(grepl("not implemented", w)))
  expect_true(is.matrix(f_mr$eif) && max(abs(colMeans(f_mr$eif))) < 1e-8)   # EIF still mean-zero with the channel folded
  expect_true(all(is.finite(f_mr$att_gt$se)) && all(f_mr$att_gt$se > 0))
  # folding the weight channel changes the SE vs the plug-in (it is a genuine, non-degenerate contribution)
  f_pl <- suppressWarnings(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "efficient",
                 misspec_robust = FALSE, aggregate = "none"))
  expect_true(max(abs(f_mr$att_gt$se - f_pl$att_gt$se)) > 1e-8)
  expect_equal(f_mr$att_gt$att, f_pl$att_gt$att, tolerance = 1e-10)          # point estimates unchanged
})

test_that("sieve AVERAGED + misspec_robust degrades with a warning (pooled-Omega-bar channel unstable)", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve"); on.exit(options(old))
  w <- testthat::capture_warnings(
    f_mr <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "averaged",
                 misspec_robust = TRUE, aggregate = "none"))
  expect_true(any(grepl("'averaged' is not implemented", w)))
  expect_true(all(is.finite(f_mr$att_gt$se)))
  # degraded == averaged with the bundle's other two effects on and the weight channel off
  f_eq <- suppressWarnings(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "averaged",
                 misspec_robust = FALSE, estimation_effect = TRUE, higher_order = TRUE, aggregate = "none"))
  expect_equal(f_mr$att_gt$se, f_eq$att_gt$se, tolerance = 1e-10)
})
