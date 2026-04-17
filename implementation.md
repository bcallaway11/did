# Implementation Record — PT-All Pair Enumeration Fix

**Date**: 2026-04-13
**Branch**: feature-efficient-DiD-estimator
**Commit**: 10af087

---

## Files Modified

### `R/edid-pairs.R`

- **Roxygen `@details`**: Replaced the inaccurate description ("ranges over all
  cohorts including never-treated") with the correct description of the two-branch
  tpre rule: self-pairs include `period_1`, cross-pairs exclude it.

- **PT-All loop (lines 56–87)**: Replaced the old loop that iterated over
  `c(treatment_groups, never_treated_val)` with a new loop that iterates over
  `treatment_groups` only.
  - Self-pair (`gp == target_g`): `valid_tpre = {t : t < eff_start}` — includes
    `period_1` (degenerate CS DiD moment; comparison-cohort EIF is identically zero
    because `delta_gp = Y_g(period_1) - Y_g(period_1) = 0`).
  - Cross-pair (`gp != target_g`): `valid_tpre = {t : period_1 < t < eff_start}` —
    excludes `period_1` (non-degenerate moments only).
  - PT-Post branch: unchanged byte-for-byte.

### `R/edid-nocov.R`

- **`compute_omega_star_nocov_edid` — `delta_gp_cache` loop**: Removed
  `if (is.finite(gp_val)) { ... } else { ... }` wrapper; kept only the finite
  branch. `gp_val` is always finite for PT-All after the fix.

- **`compute_omega_star_nocov_edid` — `n_gp_j` lookup**: Replaced
  `if (is.finite(gp_j)) { sum(...) } else { n_inf }` with direct
  `sum(panel_obj$cohort_masks[[as.character(gp_j)]])`.

- **`compute_omega_star_nocov_edid` — `n_gp_k` lookup**: Same simplification
  as `n_gp_j`.

- **`compute_generated_outcomes_nocov_edid` — `gp_j` branch**: Removed
  `if (is.finite(gp_j)) { ... } else { ... }` wrapper; kept only the finite
  branch (direct lookup of comparison cohort mask and mean).

- **`compute_eif_nocov_edid` — comparison cohort branch**: Removed the
  `else` branch that handled `gp_j == Inf` by adding an extra
  `E[Y_inf(tpre) - Y_inf(1)]` EIF term. This was dead code for PT-All since
  `Inf` is no longer a comparison cohort.

All PT-Post branches in all four functions are unchanged.

---

## Root Cause Fixed

`candidate_gps <- c(treatment_groups, never_treated_val)` in the old
PT-All loop included `Inf` as a comparison cohort. For a panel with T time
periods, this produced T-1 `(gp=Inf, tpre=t)` pairs that all collapse to
the same CS DiD moment `E[Y_g(t)-Y_g(1)] - E[Y_inf(t)-Y_inf(tpre)]`, creating
near-singular Omega and a misspecified analytical covariance matrix.

The author's reference (`benchmark/edid_sim_original.R`, line 35) limits
`g_prime` to `g_treated` (finite cohorts only), which matches the new
implementation.

---

## Verification Results

Pair enumeration for `g=3, treatment_groups=c(3,5,7), time_periods=0:9, period_1=0`:
- Pair count: 13 (correct)
- Any Inf gp: FALSE (correct)
- Degenerate `(gp=3, tpre=0)` pair present: TRUE (correct)

Pairs: gp=3 with tpre in {0,1,2}; gp=5 with tpre in {1,2,3,4}; gp=7 with
tpre in {1,2,3,4,5,6}.

ATT comparison against `efficient_did_unc_stagg()` on benchmark data:
- Matched cells: 15
- Max ATT diff: 2.2e-15 (machine precision)
- All ATT within 1e-10: TRUE

SE comparison:
- Max SE diff: 0.000165
- This small discrepancy is a pre-existing finite-sample design difference:
  the reference uses R's `cov()` (n-1 denominator) while our `cov_nn_edid`
  uses n. Not a regression introduced by this fix.

---

## Known Limitations / Deferred Items

- Covariate path is still a stub; no changes needed here.
- SE denominator choice (n vs n-1) is pre-existing; tester should confirm
  this is within acceptable tolerance.

---

## Design Choices

- The degenerate `(gp=target_g, tpre=period_1)` pair is included via the
  self-pair branch and requires no special downstream handling: the zero
  `delta_gp` vector produces zero Omega Terms C and D contributions and
  zero EIF contribution automatically via the existing formula.
- The `is.finite()` guards on Term C_j and C_k checks in the Omega inner
  loop (`if (is.finite(gp_j) && gp_j == target_g)`) were left in place —
  they are now vacuously satisfied (always finite) but harmless and do not
  affect correctness.
