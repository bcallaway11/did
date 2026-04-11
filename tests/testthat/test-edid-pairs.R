library(testthat)

# Helper: consistent args for enumerate_valid_pairs_edid
default_treatment_groups <- c(3L, 5L)
default_time_periods     <- 1:7
default_period_1         <- 1L

# ============================================================
# 4.1 PT-Post: exactly one pair
# ============================================================
test_that("enumerate_valid_pairs_edid() returns 1 pair under PT-Post for post-treatment period", {
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = default_treatment_groups,
    time_periods     = default_time_periods,
    period_1         = default_period_1,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(pairs), 1L)
  expect_equal(pairs$gp[1], Inf)
  expect_equal(pairs$tpre[1], 2L)  # g - 1 - anticipation = 3 - 1 - 0 = 2
})

test_that("enumerate_valid_pairs_edid() returns 1 pair under PT-Post with anticipation=1", {
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 1L,
    never_treated_val = Inf
  )
  # baseline = g - 1 - anticipation = 3 - 1 - 1 = 1 = period_1
  # Since period_1 is excluded, should return 0 pairs
  expect_equal(nrow(pairs), 0L)
})

test_that("enumerate_valid_pairs_edid() PT-Post baseline = period_1 returns 0 pairs", {
  # When g - 1 - anticipation equals period_1, no valid pairs
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 2L,
    treatment_groups = c(2L),
    time_periods     = 1:4,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  # g - 1 - 0 = 1 = period_1 -> 0 pairs
  expect_equal(nrow(pairs), 0L)
})

# ============================================================
# 4.2 PT-All: multiple pairs including same-cohort
# ============================================================
test_that("enumerate_valid_pairs_edid() PT-All includes same-cohort comparisons", {
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L, 5L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_true(nrow(pairs) > 1L)
  # Same-cohort comparison gp=3 must be present
  expect_true(any(pairs$gp == 3L))
  # Never-treated comparison gp=Inf must be present
  expect_true(any(is.infinite(pairs$gp)))
  # period_1 must NOT appear as tpre
  expect_false(any(pairs$tpre == 1L))
})

test_that("enumerate_valid_pairs_edid() PT-All excludes period_1 as tpre", {
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_false(1L %in% pairs$tpre)
})

test_that("enumerate_valid_pairs_edid() PT-All with anticipation=1 adjusts gp's effective treatment", {
  # gp=3 with anticipation=1: effective treatment = 3-1=2; valid tpre < 2, excl period_1=1 -> none
  # gp=Inf: all periods except period_1 -> tpre = 2,3,4
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 1L,
    never_treated_val = Inf
  )
  # gp=3 has no valid tpre (effective treatment = 2, only tpre < 2 excl 1 = none)
  expect_false(any(pairs$gp == 3L))
  # never-treated pairs should still exist
  expect_true(any(is.infinite(pairs$gp)))
})

# ============================================================
# 4.3 Never-treated pairs retained
# ============================================================
test_that("enumerate_valid_pairs_edid() PT-All includes all never-treated periods except period_1", {
  time_periods <- 1:6
  period_1     <- 1L
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L),
    time_periods     = time_periods,
    period_1         = period_1,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  inf_pairs <- pairs[is.infinite(pairs$gp), ]
  # All periods 2..6 should appear as tpre for gp=Inf (5 periods)
  expect_equal(sort(inf_pairs$tpre), 2:6)
})

# ============================================================
# 4.4 Empty pairs cases
# ============================================================
test_that("enumerate_valid_pairs_edid() returns 0 rows when target_g is first-ever cohort with period_1 baseline", {
  # g=2, periods=1:4, period_1=1: PT-Post baseline = g-1-0=1 = period_1, so 0 pairs
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 2L,
    treatment_groups = c(2L),
    time_periods     = 1:4,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(pairs), 0L)
  expect_true(is.data.frame(pairs))
  expect_true("gp" %in% names(pairs))
  expect_true("tpre" %in% names(pairs))
})

# ============================================================
# 4.5 Return type
# ============================================================
test_that("enumerate_valid_pairs_edid() always returns a data.frame with gp and tpre columns", {
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_s3_class(pairs, "data.frame")
  expect_named(pairs, c("gp", "tpre"))
})
