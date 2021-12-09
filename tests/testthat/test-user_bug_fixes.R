library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)


#-----------------------------------------------------------------------------
#
# These are tests (primarily coming from github) related to bugs
# encountered by users.
#
#-----------------------------------------------------------------------------

test_that("having column named t1 causes code to crash", {
  data(mpdta, package="did")
  out <- att_gt(yname = "lemp",
                gname = "first.treat",
                idname = "countyreal",
                tname = "year",
                xformla = ~1,
                data = mpdta,
                est_method = "reg",
                control_group="notyettreated"
                )
  mpdta$t1 <- 1

  out <- att_gt(yname = "lemp",
                gname = "first.treat",
                idname = "countyreal",
                tname = "year",
                xformla = ~1,
                data = mpdta,
                est_method = "reg",
                control_group="notyettreated"
                )
  expect_false(is.null(out), "code crashed due to stange variable names")
})

test_that("missing covariates", {
  # should warn about missing data but otherwise run
  data(mpdta, package="did")

  mpdta[1, "lpop"] <- NA
  expect_warning(out <- att_gt(yname = "lemp",
                gname = "first.treat",
                idname = "countyreal",
                tname = "year",
                xformla = ~lpop,
                data = mpdta,
                est_method = "reg",
                control_group="notyettreated"
                ))

})
