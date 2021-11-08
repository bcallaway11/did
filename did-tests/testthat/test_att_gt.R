library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)


## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test each estimation method with panel data
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, ipw model is incorectly specified here
#-----------------------------------------------------------------------------
test_that("att_gt works w/o dynamics, time effects, or group effects", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)

  # dr
  res_dr <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="dr")
  # reg
  res_reg <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="reg")


  res_ipw <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
              gname="G", est_method="ipw")
  res_ipw

  expect_equal(res_dr$att[1], 1, tol=.5)
  expect_equal(res_reg$att[1], 1, tol=.5)
  expect_equal(res_ipw$att[1], 1, tol=.5)
})
