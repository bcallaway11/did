#-----------------------------------------------------------------------------
#
# These are tests for the conditional pre-tests in the conditional_did_pretest
# function.
#
#-----------------------------------------------------------------------------

library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)


test_that("conditional did pre-test", {
  # this code is time-consuming, and we are just going to test that the code runs.
  sp <- reset.sim(time.periods=3, n=1000)
  data <- build_sim_dataset(sp)

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

