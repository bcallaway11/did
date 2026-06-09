#-----------------------------------------------------------------------------
#
# These are tests for the conditional pre-tests in the conditional_did_pretest
# function.
#
#-----------------------------------------------------------------------------

library(DRDID)
library(BMisc)
#library(ggplot2)
#library(ggpubr)


test_that("conditional did pre-test", {
  # this code is time-consuming, and we are just going to test that the code runs.
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

