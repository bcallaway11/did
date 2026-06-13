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
  # When tpre == period_1, the standard 2x2 DiD is still valid
  expect_equal(nrow(pairs), 1L)
  expect_equal(pairs$tpre, 1L)
})

test_that("enumerate_valid_pairs_edid() PT-Post baseline = period_1 returns 1 pair", {
  # When g - 1 - anticipation equals period_1, this is the standard 2x2 DiD
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 2L,
    treatment_groups = c(2L),
    time_periods     = 1:4,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  # g - 1 - 0 = 1 = period_1 -> valid pair (Inf, 1)
  expect_equal(nrow(pairs), 1L)
  expect_equal(pairs$tpre, 1L)
})

# ============================================================
# 4.2 PT-All: multiple pairs including same-cohort
# Updated 2026-04-13: gp=Inf is no longer included in PT-All;
# period_1 IS valid as tpre for self-pairs (gp == target_g).
# ============================================================
test_that("enumerate_valid_pairs_edid() PT-All includes same-cohort comparisons but no gp=Inf", {
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
  # Never-treated comparison gp=Inf must NOT be present (PT-All uses only treated cohorts)
  expect_false(any(is.infinite(pairs$gp)))
  # period_1 IS valid as tpre for self-pair (gp=3, tpre=1)
  expect_true(any(pairs$gp == 3L & pairs$tpre == 1L))
})

test_that("enumerate_valid_pairs_edid() PT-All includes period_1 as tpre for self-pair", {
  # Single cohort: gp=target_g is the only comparison cohort (self-pair).
  # Self-pair includes period_1 as a valid tpre (degenerate CS DiD moment).
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  # period_1=1 MUST appear as tpre for the self-pair
  expect_true(1L %in% pairs$tpre)
  # All pairs have finite gp (no gp=Inf)
  expect_true(all(is.finite(pairs$gp)))
})

test_that("enumerate_valid_pairs_edid() PT-All with anticipation=1 adjusts effective treatment", {
  # target_g=3, anticipation=1: eff_start(3) = 2
  # Self-pair (gp=3): valid tpre < 2, includes period_1=1 -> {1} = 1 row
  # No gp=Inf in PT-All
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 1L,
    never_treated_val = Inf
  )
  # gp=3 self-pair with eff_start=2 has tpre=1 (period_1 is valid)
  expect_true(any(pairs$gp == 3L))
  expect_equal(nrow(pairs[pairs$gp == 3L, ]), 1L)
  expect_equal(pairs$tpre[pairs$gp == 3L], 1L)
  # No never-treated pairs in PT-All
  expect_false(any(is.infinite(pairs$gp)))
})

# ============================================================
# 4.3 Self-pair structure in PT-All
# Updated 2026-04-13: gp=Inf no longer exists in PT-All;
# the correct invariant is that only treated cohorts appear as gp.
# ============================================================
test_that("enumerate_valid_pairs_edid() PT-All has only treated-cohort gp values", {
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
  # Only gp=3 (self-pair) should appear; no gp=Inf
  expect_true(all(pairs$gp == 3L))
  expect_true(all(is.finite(pairs$gp)))
  # Self-pair: tpre < 3, includes period_1=1 -> {1, 2}
  expect_equal(sort(pairs$tpre), c(1L, 2L))
})

# ============================================================
# 4.4 Empty pairs cases
# ============================================================
test_that("enumerate_valid_pairs_edid() returns 1 pair when target_g is first-ever cohort with period_1 baseline", {
  # g=2, periods=1:4, period_1=1: PT-Post baseline = g-1-0=1 = period_1
  # This is the standard 2x2 DiD with the first period as base — valid
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 2L,
    treatment_groups = c(2L),
    time_periods     = 1:4,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(pairs), 1L)
  expect_true(is.data.frame(pairs))
  expect_true("gp" %in% names(pairs))
  expect_true("tpre" %in% names(pairs))
  expect_equal(pairs$tpre, 1L)
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
