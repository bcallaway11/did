# =============================================================================
# With only two time periods, uniform confidence bands coincide with pointwise
# confidence intervals. Pre-processing overrides cband to FALSE; the returned MP
# object must be self-consistent with that override (pointwise critical value,
# DIDparams$cband = FALSE, "Pointwise" summary label) and the override must be
# announced. Objects with three or more periods are unaffected.
# =============================================================================

pointwise_cval <- stats::qnorm(1 - 0.05 / 2)

test_that("two periods + cband = TRUE announces the override and returns a pointwise critical value", {
  set.seed(20260814)
  sp <- did::reset.sim(time.periods = 2)
  sp$n <- 2000
  d2 <- did::build_sim_dataset(sp)

  for (fm in c(TRUE, FALSE)) {
    msgs <- testthat::capture_messages(
      res <- suppressWarnings(att_gt(yname = "Y", xformla = ~X, data = d2,
        tname = "period", idname = "id", gname = "G", faster_mode = fm,
        bstrap = TRUE, biters = 200, cband = TRUE))
    )
    expect_true(any(grepl("Only two time periods", msgs)))
    expect_equal(unname(res$c), pointwise_cval)
    expect_false(res$DIDparams$cband)
  }
})

test_that("two periods without cband is silent about the override", {
  set.seed(20260814)
  sp <- did::reset.sim(time.periods = 2)
  sp$n <- 2000
  d2 <- did::build_sim_dataset(sp)

  for (fm in c(TRUE, FALSE)) {
    msgs <- testthat::capture_messages(
      res <- suppressWarnings(att_gt(yname = "Y", xformla = ~X, data = d2,
        tname = "period", idname = "id", gname = "G", faster_mode = fm,
        bstrap = TRUE, biters = 200, cband = FALSE))
    )
    expect_false(any(grepl("Only two time periods", msgs)))
    expect_false(res$DIDparams$cband)
  }
})

test_that("three or more periods still get the bootstrap uniform band, with no message", {
  set.seed(20260814)
  sp <- did::reset.sim(time.periods = 4)
  sp$n <- 2000
  d4 <- did::build_sim_dataset(sp)

  for (fm in c(TRUE, FALSE)) {
    msgs <- testthat::capture_messages(
      res <- suppressWarnings(att_gt(yname = "Y", xformla = ~X, data = d4,
        tname = "period", idname = "id", gname = "G", faster_mode = fm,
        bstrap = TRUE, biters = 200, cband = TRUE))
    )
    expect_false(any(grepl("Only two time periods", msgs)))
    expect_true(res$DIDparams$cband)
    expect_gt(unname(res$c), pointwise_cval)
  }
})
