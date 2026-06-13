#-----------------------------------------------------------------------------
#
# These are tests for the conditional pre-tests in the conditional_did_pretest
# function.
#
#-----------------------------------------------------------------------------


test_that("conditional did pre-test", {
  skip_on_cran()  # bootstrap-based CvM test is time-consuming
  sp <- did::reset.sim(time.periods=3, n=1000)
  data <- did::build_sim_dataset(sp)

  cdp <- conditional_did_pretest(yname="Y",
                                 tname="period",
                                 idname="id",
                                 gname="G",
                                 xformla=~X,
                                 data=data)


  # check that test statistic and critical value take reasonable values
  expect_true((cdp$CvM > 0) & (cdp$CvM < 20))
  expect_true((cdp$CvMcval > 0) & (cdp$CvMcval < 100))
})

test_that("pretest setup-bundle/y-override path is bit-identical to the legacy data-copy loop", {
  skip_on_cran()
  # The pretest's per-unit loop reuses a y-invariant setup bundle and passes
  # only a modified outcome vector (dp$setup_precomp / dp$y_override), and
  # compute.att_gt skips post-treatment cells (dp$pretreatment_cells_only).
  # options(did.disable_precompute = TRUE) forces the legacy data-copy loop
  # (the bundle is only consumed on the panel precompute path), so the two
  # runs must agree EXACTLY: same kept cells, same estimator inputs, and no
  # RNG use in compute.att_gt, hence identical(), not just a tolerance.
  old_opt <- getOption("did.disable_precompute")
  withr::defer(options(did.disable_precompute = old_opt))
  set.seed(20260401)
  sp <- did::reset.sim(time.periods = 3, n = 150)
  data <- did::build_sim_dataset(sp)
  run_pt <- function() {
    set.seed(42)  # same multiplier-bootstrap draws in both runs
    suppressWarnings(suppressMessages(
      conditional_did_pretest(yname = "Y", tname = "period", idname = "id",
                              gname = "G", xformla = ~X, data = data,
                              est_method = "ipw", biters = 100)))
  }
  options(did.disable_precompute = NULL)
  pt_bundle <- run_pt()
  options(did.disable_precompute = TRUE)
  pt_legacy <- run_pt()
  expect_identical(pt_bundle$CvM, pt_legacy$CvM)
  expect_identical(pt_bundle$CvMpval, pt_legacy$CvMpval)
  expect_identical(as.numeric(pt_bundle$CvMcval), as.numeric(pt_legacy$CvMcval))
  expect_identical(pt_bundle$CvMb, pt_legacy$CvMb)
})

test_that("conditional pre-test CvM is on the bootstrap scale (R>=4.0 orientation regression)", {
  skip_on_cran()
  # Regression guard for the class() == "matrix" length-2 bug (R >= 4.0): with MORE
  # THAN ONE pre-treatment (g,t) cell the observed CvM statistic was left in
  # (n_gt x nX) orientation and therefore scaled by n / n_gt relative to the bootstrap
  # null distribution (which is correctly (nX x n_gt) via test.mboot's t(apply(...))),
  # driving the p-value spuriously to ~0 so the pre-test almost always rejected.
  # mpdta (conditioning on lpop) has several pre-treatment cells and plausibly
  # satisfies conditional parallel trends, so the *fixed* test should NOT reject.
  data(mpdta, package = "did")
  set.seed(42)
  keep <- sample(unique(mpdta$countyreal), 150)
  d <- mpdta[mpdta$countyreal %in% keep, ]
  pt <- suppressWarnings(suppressMessages(
    conditional_did_pretest(yname = "lemp", tname = "year", idname = "countyreal",
                            gname = "first.treat", xformla = ~lpop, data = d,
                            est_method = "ipw", biters = 199)))
  expect_true(is.finite(pt$CvM))
  # the bug inflated CvM ~(n / n_gt)-fold above the critical value -> p ~ 0
  expect_lt(pt$CvM, pt$CvMcval)
  expect_gt(pt$CvMpval, 0.05)
})

