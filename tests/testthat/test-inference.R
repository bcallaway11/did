#-----------------------------------------------------------------------------
#
# These are tests for the inference procedure for computing ATT(g,t)'s and
# different aggregation methods.  These tests are meant to run fast ---
# they only confirm that standard errors are not changing across different
# iterations of the code.  More extensive/direct tests are available at
# did/tests/att_gt_inference_tests.Rmd.
#
# IMPORTANT: For these tests to run, version 2.0.0 of the did package should
# be available.  The path of to the old version should be set in the
# following variable.
old_did_path <- "~/R/old_packages/did"
#
# The following code will install version 2.0.0 of the did package into
# the directory above:
# install.packages("did", version="2.0.0", lib=old_did_path)
# 
#-----------------------------------------------------------------------------

library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)


test_that("inference with balanced panel data and aggregations", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)

  tryCatch(detach("package:did"), error=function(e) "")

  library(did, lib.loc="~/R/old_packages/")

  # first: make sure that we are using right version of package for these estimates
  expect_true(packageVersion("did") == "2.0.0", "wrong version of package")
  
  set.seed(1234)
  # dr
  dr_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr")
  # reg
  reg_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg")
  # reg
  ipw_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw")

  # aggregations
  dyn_2.0 <- aggte(ipw_2.0, type="dynamic")
  group_2.0 <- aggte(reg_2.0, type="group")
  cal_2.0 <- aggte(dr_2.0, type="calendar")
  detach("package:did")

  
  devtools::load_all("~/Dropbox/did")
  expect_true(packageVersion("did") != "2.0.0", "wrong version of package")
  set.seed(1234)
  #dr 
  dr_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr")
  # reg
  reg_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg")
  # reg
  ipw_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw")

  # aggregations
  dyn_new <- aggte(ipw_new, type="dynamic")
  group_new <- aggte(reg_new, type="group")
  cal_new <- aggte(dr_new, type="calendar")
  detach("package:did")

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
  
})  


test_that("inference with clustering", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)

  tryCatch(detach("package:did"), error=function(e) "")

  library(did, lib.loc="~/R/old_packages/")

  # first: make sure that we are using right version of package for these estimates
  expect_true(packageVersion("did") == "2.0.0", "wrong version of package")
  
  set.seed(1234)
  # dr
  dr_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", clustervars="cluster")
  # reg
  reg_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", clustervars="cluster")
  # reg
  ipw_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", clustervars="cluster")

  # aggregations
  dyn_2.0 <- aggte(ipw_2.0, type="dynamic")
  group_2.0 <- aggte(reg_2.0, type="group")
  cal_2.0 <- aggte(dr_2.0, type="calendar")
  detach("package:did")

  
  devtools::load_all("~/Dropbox/did")
  expect_true(packageVersion("did") != "2.0.0", "wrong version of package")
  set.seed(1234)
  #dr 
  dr_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", clustervars="cluster")
  # reg
  reg_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", clustervars="cluster")
  # reg
  ipw_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", clustervars="cluster")

  # aggregations
  dyn_new <- aggte(ipw_new, type="dynamic")
  group_new <- aggte(reg_new, type="group")
  cal_new <- aggte(dr_new, type="calendar")
  detach("package:did")

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
  
})  

test_that("same inference with unbalanced panel and panel data", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  
  res_factor <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", 
                        gname="G", est_method="dr", clustervars="cluster")

  #-----------------------------------------------------------------------------
  # clustered standard errors with unbalanced panel
  res_ub <- att_gt(yname="Y",
              tname="period",
              idname="id",
              gname="G",
              xformla=~X,
              data=data,
              panel=TRUE,
              allow_unbalanced_panel=TRUE,
              clustervars="cluster")

  # check that influence function is the same if we use unbalanced panel approach
  # vs. balanced panel approach
  # Note: differences are showing up in never-treated units (that are towards end of sample)
  expect_true(all(res_factor$inffunc == res_ub$inffunc))
})


test_that("inference with repeated cross sections", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp, panel=FALSE)

  tryCatch(detach("package:did"), error=function(e) "")

  library(did, lib.loc="~/R/old_packages/")

  # first: make sure that we are using right version of package for these estimates
  expect_true(packageVersion("did") == "2.0.0", "wrong version of package")
  
  set.seed(1234)
  # dr
  dr_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=FALSE)
  # reg
  reg_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", panel=FALSE)
  # reg
  ipw_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", panel=FALSE)

  # aggregations
  dyn_2.0 <- aggte(ipw_2.0, type="dynamic")
  group_2.0 <- aggte(reg_2.0, type="group")
  cal_2.0 <- aggte(dr_2.0, type="calendar")
  detach("package:did")

  
  devtools::load_all("~/Dropbox/did")
  expect_true(packageVersion("did") != "2.0.0", "wrong version of package")
  set.seed(1234)
  #dr 
  dr_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=FALSE)
  # reg
  reg_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", panel=FALSE)
  # reg
  ipw_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", panel=FALSE)

  # aggregations
  dyn_new <- aggte(ipw_new, type="dynamic")
  group_new <- aggte(reg_new, type="group")
  cal_new <- aggte(dr_new, type="calendar")
  detach("package:did")

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
  
})


test_that("inference with repeated cross sections and clustering", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp, panel=FALSE)

  tryCatch(detach("package:did"), error=function(e) "")

  library(did, lib.loc="~/R/old_packages/")

  # first: make sure that we are using right version of package for these estimates
  expect_true(packageVersion("did") == "2.0.0", "wrong version of package")
  
  set.seed(1234)
  # dr
  dr_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", clustervars="cluster", panel=FALSE)
  # reg
  reg_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", clustervars="cluster", panel=FALSE)
  # reg
  ipw_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", clustervars="cluster", panel=FALSE)

  # aggregations
  dyn_2.0 <- aggte(ipw_2.0, type="dynamic")
  group_2.0 <- aggte(reg_2.0, type="group")
  cal_2.0 <- aggte(dr_2.0, type="calendar")
  detach("package:did")

  
  devtools::load_all("~/Dropbox/did")
  expect_true(packageVersion("did") != "2.0.0", "wrong version of package")
  set.seed(1234)
  #dr 
  dr_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", clustervars="cluster", panel=FALSE)
  # reg
  reg_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", clustervars="cluster", panel=FALSE)
  # reg
  ipw_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", clustervars="cluster", panel=FALSE)

  # aggregations
  dyn_new <- aggte(ipw_new, type="dynamic")
  group_new <- aggte(reg_new, type="group")
  cal_new <- aggte(dr_new, type="calendar")
  detach("package:did")

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
  
})  



test_that("inference with unbalanced panel", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data <- data[-3,]

  tryCatch(detach("package:did"), error=function(e) "")

  library(did, lib.loc="~/R/old_packages/")

  # first: make sure that we are using right version of package for these estimates
  expect_true(packageVersion("did") == "2.0.0", "wrong version of package")
  
  set.seed(1234)
  # dr
  dr_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=TRUE, allow_unbalanced_panel=TRUE)
  # reg
  reg_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", panel=TRUE, allow_unbalanced_panel=TRUE)
  # reg
  ipw_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", panel=TRUE, allow_unbalanced_panel=TRUE)

  # aggregations
  dyn_2.0 <- aggte(ipw_2.0, type="dynamic")
  group_2.0 <- aggte(reg_2.0, type="group")
  cal_2.0 <- aggte(dr_2.0, type="calendar")
  detach("package:did")

  
  devtools::load_all("~/Dropbox/did")
  expect_true(packageVersion("did") != "2.0.0", "wrong version of package")
  set.seed(1234)
  #dr 
  dr_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", panel=TRUE, allow_unbalanced_panel=TRUE)
  # reg
  reg_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", panel=TRUE, allow_unbalanced_panel=TRUE)
  # reg
  ipw_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", panel=TRUE, allow_unbalanced_panel=TRUE)

  # aggregations
  dyn_new <- aggte(ipw_new, type="dynamic")
  group_new <- aggte(reg_new, type="group")
  cal_new <- aggte(dr_new, type="calendar")
  detach("package:did")

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
  
})

test_that("inference with unbalanced panel and clustering", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data <- data[-3,]

  tryCatch(detach("package:did"), error=function(e) "")

  library(did, lib.loc="~/R/old_packages/")

  # first: make sure that we are using right version of package for these estimates
  expect_true(packageVersion("did") == "2.0.0", "wrong version of package")
  
  set.seed(1234)
  # dr
  dr_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", clustervars="cluster", allow_unbalanced_panel=TRUE)
  # reg
  reg_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", clustervars="cluster", allow_unbalanced_panel=TRUE)
  # reg
  ipw_2.0 <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", clustervars="cluster", allow_unbalanced_panel=TRUE)

  # aggregations
  dyn_2.0 <- aggte(ipw_2.0, type="dynamic")
  group_2.0 <- aggte(reg_2.0, type="group")
  cal_2.0 <- aggte(dr_2.0, type="calendar")
  detach("package:did")

  
  devtools::load_all("~/Dropbox/did")
  expect_true(packageVersion("did") != "2.0.0", "wrong version of package")
  set.seed(1234)
  #dr 
  dr_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                   gname="G", est_method="dr", clustervars="cluster", allow_unbalanced_panel=TRUE)
  # reg
  reg_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="reg", clustervars="cluster", allow_unbalanced_panel=TRUE)
  # reg
  ipw_new <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
                    gname="G", est_method="ipw", clustervars="cluster", allow_unbalanced_panel=TRUE)

  # aggregations
  dyn_new <- aggte(ipw_new, type="dynamic")
  group_new <- aggte(reg_new, type="group")
  cal_new <- aggte(dr_new, type="calendar")
  detach("package:did")

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
  
})  
