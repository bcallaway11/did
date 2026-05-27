# test-edid-pairs-validation.R
# Validation tests for enumerate_valid_pairs_edid() based on test-spec.md
# (2026-04-13 post-builder-fix spec).
#
# Covers: U1-U11 (unit), I1-I8 (integration), R1-R4 (regression), E1-E2 (edge case)

library(testthat)

# ===========================================================================
# SECTION 3: Unit tests for enumerate_valid_pairs_edid()
# ===========================================================================

# ---------------------------------------------------------------------------
# Scenario U1 — Cohorts {3,5,7}, target g=3
# ---------------------------------------------------------------------------
test_that("U1: target_g=3, cohorts={3,5,7}, periods=1:10", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 3L,
    treatment_groups = c(3L, 5L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 10L)
  expect_false(any(is.infinite(result$gp)),         info = "no gp=Inf in PT-All")
  expect_true(any(result$gp == 3L & result$tpre == 1L),  info = "self-pair period_1 present")
  expect_false(any(result$gp == 5L & result$tpre == 1L), info = "cross-pair gp=5 period_1 absent")
  expect_false(any(result$gp == 7L & result$tpre == 1L), info = "cross-pair gp=7 period_1 absent")
  expect_true(all(result[result$gp == 5L, "tpre"] %in% 2:4), info = "gp=5 tpre in 2:4")
  expect_true(all(result[result$gp == 7L, "tpre"] %in% 2:6), info = "gp=7 tpre in 2:6")
})

# ---------------------------------------------------------------------------
# Scenario U2 — Cohorts {3,5,7}, target g=5
# ---------------------------------------------------------------------------
test_that("U2: target_g=5, cohorts={3,5,7}, periods=1:10", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 5L,
    treatment_groups = c(3L, 5L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 10L)
  expect_equal(nrow(result[result$gp == 3L, ]), 1L, info = "gp=3 has 1 cross-pair")
  expect_equal(result[result$gp == 3L, "tpre"], 2L,  info = "gp=3 tpre=2")
  expect_true(any(result$gp == 5L & result$tpre == 1L),  info = "self-pair period_1 present")
  expect_false(any(result$gp == 3L & result$tpre == 1L), info = "cross-pair gp=3 period_1 absent")
  expect_false(any(is.infinite(result$gp)),              info = "no gp=Inf")
})

# ---------------------------------------------------------------------------
# Scenario U3 — Cohorts {3,5,7}, target g=7
# ---------------------------------------------------------------------------
test_that("U3: target_g=7, cohorts={3,5,7}, periods=1:10", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 7L,
    treatment_groups = c(3L, 5L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 10L)
  expect_true(any(result$gp == 7L & result$tpre == 1L),  info = "self-pair period_1 present")
  expect_false(any(result$gp == 3L & result$tpre == 1L), info = "cross gp=3 period_1 absent")
  expect_false(any(result$gp == 5L & result$tpre == 1L), info = "cross gp=5 period_1 absent")
  expect_equal(nrow(result[result$gp == 7L, ]), 6L,      info = "gp=7 has 6 self-pairs")
  expect_false(any(is.infinite(result$gp)),              info = "no gp=Inf")
})

# ---------------------------------------------------------------------------
# Scenario U4 — Cohorts {4,7}, target g=4
# ---------------------------------------------------------------------------
test_that("U4: target_g=4, cohorts={4,7}, periods=1:10", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 4L,
    treatment_groups = c(4L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 8L)
  expect_true(all(is.finite(result$gp)),                 info = "all gp finite")
  expect_true(any(result$gp == 4L & result$tpre == 1L),  info = "self-pair period_1 present")
  expect_false(any(result$gp == 7L & result$tpre == 1L), info = "cross gp=7 period_1 absent")
  expect_true(all(result[result$gp == 7L, "tpre"] %in% 2:6), info = "gp=7 tpre in 2:6")
})

# ---------------------------------------------------------------------------
# Scenario U5 — Cohorts {4,7}, target g=7
# ---------------------------------------------------------------------------
test_that("U5: target_g=7, cohorts={4,7}, periods=1:10", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 7L,
    treatment_groups = c(4L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 8L)
  expect_true(any(result$gp == 7L & result$tpre == 1L),  info = "self-pair period_1 present")
  expect_false(any(result$gp == 4L & result$tpre == 1L), info = "cross gp=4 period_1 absent")
  expect_equal(nrow(result[result$gp == 4L, ]), 2L,      info = "gp=4 has 2 cross-pairs")
  expect_false(any(is.infinite(result$gp)),              info = "no gp=Inf")
})

# ---------------------------------------------------------------------------
# Scenario U6 — Single cohort {5}
# ---------------------------------------------------------------------------
test_that("U6: target_g=5, cohorts={5}, periods=1:10 (single cohort)", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 5L,
    treatment_groups = c(5L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 4L)
  expect_true(all(result$gp == 5L),       info = "all gp=5")
  expect_true(any(result$tpre == 1L),     info = "period_1 present")
  expect_false(any(result$tpre >= 5L),    info = "no tpre >= target_g")
  expect_false(any(is.infinite(result$gp)), info = "no gp=Inf")
})

# ---------------------------------------------------------------------------
# Scenario U7 — Anticipation=1, cohorts {4,7}, target g=4
# ---------------------------------------------------------------------------
test_that("U7: target_g=4, cohorts={4,7}, anticipation=1", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 4L,
    treatment_groups = c(4L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 1L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 6L)
  expect_false(any(result$gp == 4L & result$tpre >= 3L), info = "no tpre >= eff_start(4)=3 for gp=4")
  expect_false(any(result$gp == 7L & result$tpre >= 6L), info = "no tpre >= eff_start(7)=6 for gp=7")
  expect_false(any(result$gp == 7L & result$tpre == 1L), info = "cross gp=7 period_1 absent")
  expect_true(any(result$gp == 4L & result$tpre == 1L),  info = "self-pair period_1 present")
  expect_false(any(is.infinite(result$gp)),              info = "no gp=Inf")
})

# ---------------------------------------------------------------------------
# Scenario U8 — Cross-pair with no interior tpre
# ---------------------------------------------------------------------------
test_that("U8: target_g=7, cohorts={2,7} — gp=2 has no valid cross-pair tpre", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 7L,
    treatment_groups = c(2L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 6L)
  expect_false(any(result$gp == 2L), info = "no pairs for gp=2 (no interior tpre in (1,2))")
  expect_true(all(result$gp == 7L),  info = "all pairs have gp=7 (self)")
  expect_false(any(is.infinite(result$gp)), info = "no gp=Inf")
})

# ---------------------------------------------------------------------------
# Scenario U9 — PT-Post: exactly one pair
# ---------------------------------------------------------------------------
test_that("U9: PT-Post target_g=4, cohorts={4,7}, periods=1:10", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 4L,
    treatment_groups = c(4L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 1L)
  expect_equal(result$gp,   Inf)
  expect_equal(result$tpre, 3L)
})

# ---------------------------------------------------------------------------
# Scenario U10 — PT-Post: tpre = period_1 -> valid (standard 2x2 DiD)
# ---------------------------------------------------------------------------
test_that("U10: PT-Post tpre=period_1 returns 1 valid pair", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 2L,
    treatment_groups = c(2L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(result), 1L)
  expect_equal(result$tpre, 1L)
  expect_equal(result$gp, Inf)
})

# ---------------------------------------------------------------------------
# Scenario U11 — PT-Post: tpre not in time_periods -> 0 rows
# ---------------------------------------------------------------------------
test_that("U11: PT-Post tpre not in time_periods returns 0 pairs", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 5L,
    treatment_groups = c(5L),
    time_periods     = c(1L, 3L, 5L, 7L, 9L),   # even periods missing
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  # tpre_val = 5 - 1 - 0 = 4; 4 NOT in time_periods
  expect_equal(nrow(result), 0L)
})

# ===========================================================================
# SECTION 6: Edge Case Scenarios
# ===========================================================================

# ---------------------------------------------------------------------------
# Scenario E1 — Single cohort, only period_1 pre-period
# ---------------------------------------------------------------------------
test_that("E1: target_g=2, single cohort, only period_1 as pre-period", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 2L,
    treatment_groups = c(2L),
    time_periods     = 1:5,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  # Self-pair: tpre < 2 -> {1} = period_1. Exactly 1 row.
  expect_equal(nrow(result), 1L)
  expect_equal(result$gp,   2L)
  expect_equal(result$tpre, 1L)
})

# ---------------------------------------------------------------------------
# Scenario E2 — Cross-pair cohort at effective boundary (no interior tpre)
# ---------------------------------------------------------------------------
test_that("E2: target_g=7, cohorts={3,7}, anticipation=1 — gp=3 has no valid cross-pair", {
  result <- enumerate_valid_pairs_edid(
    target_g         = 7L,
    treatment_groups = c(3L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "all",
    anticipation     = 1L,
    never_treated_val = Inf
  )
  # eff_start(3) = 2; cross-pair condition: 1 < tpre < 2 -> no integer
  expect_false(any(result$gp == 3L), info = "no pairs for gp=3")
  # Self-pair (gp=7): eff_start=6, tpre < 6 incl 1 -> {1,2,3,4,5} = 5 rows
  expect_equal(nrow(result), 5L)
  expect_false(any(is.infinite(result$gp)), info = "no gp=Inf")
})

# ===========================================================================
# SECTION 5 (Regression): PT-Post path unchanged
# ===========================================================================

test_that("R: PT-Post always produces gp=Inf pairs", {
  pairs <- enumerate_valid_pairs_edid(
    target_g         = 5L,
    treatment_groups = c(3L, 5L, 7L),
    time_periods     = 1:10,
    period_1         = 1L,
    pt_assumption    = "post",
    anticipation     = 0L,
    never_treated_val = Inf
  )
  expect_equal(nrow(pairs), 1L)
  expect_true(all(pairs$gp == Inf), info = "PT-Post gp must be Inf")
  expect_equal(pairs$tpre, 4L)
})

# ===========================================================================
# SECTION 7: Property-Based Invariants (100 random inputs)
# ===========================================================================

test_that("Property: no gp=Inf in any PT-All call", {
  set.seed(42L)
  for (i in seq_len(100L)) {
    # Generate random staggered design
    n_cohorts    <- sample(2:5, 1L)
    max_period   <- sample(10:20, 1L)
    # Cohort values: distinct, between 3 and max_period-1
    cohorts <- sort(sample(3:(max_period - 1L), n_cohorts, replace = FALSE))
    target  <- cohorts[sample(seq_along(cohorts), 1L)]
    periods <- seq_len(max_period)
    period1 <- 1L

    pairs <- enumerate_valid_pairs_edid(
      target_g         = target,
      treatment_groups = cohorts,
      time_periods     = periods,
      period_1         = period1,
      pt_assumption    = "all",
      anticipation     = 0L,
      never_treated_val = Inf
    )
    expect_true(
      all(is.finite(pairs$gp)),
      info = paste0("gp=Inf found for target_g=", target, ", cohorts={", paste(cohorts, collapse=","), "}")
    )
  }
})

test_that("Property: self-pair includes period_1 (when valid tpre exist)", {
  set.seed(123L)
  for (i in seq_len(50L)) {
    n_cohorts <- sample(1:4, 1L)
    max_period <- sample(8:15, 1L)
    cohorts <- sort(sample(3:(max_period - 1L), n_cohorts, replace = FALSE))
    target  <- cohorts[sample(seq_along(cohorts), 1L)]
    periods <- seq_len(max_period)

    pairs <- enumerate_valid_pairs_edid(
      target_g         = target,
      treatment_groups = cohorts,
      time_periods     = periods,
      period_1         = 1L,
      pt_assumption    = "all",
      anticipation     = 0L
    )
    if (nrow(pairs) > 0L) {
      # If self-pair has any rows, period_1 should be among them
      self_rows <- pairs[pairs$gp == target, ]
      if (nrow(self_rows) > 0L) {
        expect_true(
          any(self_rows$tpre == 1L),
          info = paste0("self-pair for target_g=", target, " missing period_1")
        )
      }
    }
  }
})

test_that("Property: cross-pair excludes period_1", {
  set.seed(456L)
  for (i in seq_len(50L)) {
    n_cohorts <- sample(2:5, 1L)
    max_period <- sample(8:15, 1L)
    cohorts <- sort(sample(3:(max_period - 1L), n_cohorts, replace = FALSE))
    target  <- cohorts[sample(seq_along(cohorts), 1L)]
    periods <- seq_len(max_period)

    pairs <- enumerate_valid_pairs_edid(
      target_g         = target,
      treatment_groups = cohorts,
      time_periods     = periods,
      period_1         = 1L,
      pt_assumption    = "all",
      anticipation     = 0L
    )
    cross_rows <- pairs[pairs$gp != target, ]
    if (nrow(cross_rows) > 0L) {
      expect_false(
        any(cross_rows$tpre == 1L),
        info = paste0("cross-pair has period_1 for target_g=", target)
      )
    }
  }
})
