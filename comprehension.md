# Comprehension Record — PT-All Pair Enumeration Fix

**Date**: 2026-04-13
**HOLD rounds used**: 0
**Verdict**: FULLY UNDERSTOOD

---

## Input Materials Read

1. `R/edid-pairs.R` — the function under repair (`enumerate_valid_pairs_edid`)
2. `R/edid-nocov.R` — Omega, EIF, generated outcomes; four functions affected
3. `R/edid-data.R` — panel preparation; no changes needed
4. `R/edid-utils.R` — constants; no changes needed
5. `R/edid-linalg.R` — pseudoinverse/condition helpers; no changes needed
6. `R/edid-fit.R` — outer loop; no changes needed
7. `benchmark/edid_sim_original.R` — author's reference; primary source for pair rule

---

## Restated Core Requirement

Under `pt_assumption = "all"`, `enumerate_valid_pairs_edid()` currently includes the never-treated cohort (`gp = Inf`) as a comparison cohort, generating 9 redundant moments (in 10-period data) that all collapse to the same CS DiD estimate. The correct design from the author's reference excludes the never-treated as a comparison cohort entirely, and instead uses only treated cohorts (`g' ∈ G_treated`) as comparison cohorts, with a tpre rule that depends on whether `g' == target_g` or `g' != target_g`. This produces a non-redundant, well-conditioned moment system. The `pt_assumption = "post"` path is unaffected.

---

## All Formulas Restated and Verified

### Identifying Moment for pair `(g', t')`

For target cell `(g, t)` with `t >= g` (post-treatment):

```
psi_{g,t}^{(g', t')} = E[Y_g(t) - Y_g(1)]              [treated group change from period_1]
                      - E[Y_0(t) - Y_0(t')]              [never-treated time control]
                      - E[Y_{g'}(t') - Y_{g'}(1)]         [comparison cohort pre-trend baseline]
```

Where:
- `Y_g(s)` = potential outcome of cohort g at time s
- `Y_0(s)` = potential outcome of never-treated at time s
- `period_1 = min(time_periods)` = universal first period
- `g' ∈ G_treated` (never the never-treated cohort)

### Degenerate CS DiD case (`g' = g`, `t' = period_1`)

When `t' = period_1`:
```
E[Y_{g'}(t') - Y_{g'}(1)] = E[Y_g(period_1) - Y_g(period_1)] = 0
```
So the moment reduces to:
```
psi = E[Y_g(t) - Y_g(1)] - E[Y_0(t) - Y_0(t')]
```
This is the CS DiD estimand (with `t'` chosen as the pre-period for the never-treated).

**Cross-check against author's code (benchmark lines 47–52)**:
```r
filter((g == g_prime) & (t_prime < g_prime))         # includes t_prime = min(tlist)
filter((g != g_prime) & (min(t_prime) < t_prime) & (t_prime < g_prime))  # excludes min
```
This matches exactly. The author's `g_prime` ranges over `g_treated` only (line 35: `g_prime = g_treated`).

### Valid tpre Rules (with anticipation)

Effective treatment start of comparison cohort: `eff_start(g') = g' - anticipation`

| Condition | Valid tpre |
|-----------|-----------|
| `g' == target_g` | `{t : t < eff_start(g')}` — includes `period_1` |
| `g' != target_g` | `{t : period_1 < t < eff_start(g')}` — excludes `period_1` |

Source: author's code `t_prime < g_prime` (with `anticipation=0`) for `g'=g`, and `min(t_prime) < t_prime < g_prime` for `g'≠g`.

### Omega Entry for Pair `(j, k)` — PT-All, Finite gp only

```
Omega[j,k] = Cov(delta_g_t_1, delta_g_t_1)/n_g          [Term A: always present]
           + Cov(delta_inf_j, delta_inf_k)/n_inf           [Term B: never-treated cross-cov]
           - 1(gp_j == target_g) * Cov(delta_g_t_1, delta_gp_j)/n_g   [Term C_j]
           - 1(gp_k == target_g) * Cov(delta_g_t_1, delta_gp_k)/n_g   [Term C_k]
           + 1(gp_j == gp_k) * Cov(delta_gp_j, delta_gp_k)/n_gp       [Term D]
```

Where:
- `delta_g_t_1 = Y_g(t) - Y_g(period_1)` (vector over g-units)
- `delta_inf_j = Y_inf(t) - Y_inf(tpre_j)` (vector over never-treated units)
- `delta_gp_j = Y_{gp_j}(tpre_j) - Y_{gp_j}(period_1)` (vector over gp_j units)

With the new pair rule, `gp` is always finite — the `else` branch for `gp = Inf` in the Omega loop (line 91-93 in edid-nocov.R) is dead code for PT-All and should be removed.

### EIF for pair `(g', t')` — PT-All, Finite gp only

Three components per unit i:

1. **Treated group** (if unit i is in cohort g):
   `(1/pi_g) * [Y_i(t) - Y_i(period_1) - E[Y_g(t) - Y_g(period_1)]]`

2. **Never-treated** (if unit i is never-treated):
   `-(1/pi_inf) * [Y_i(t) - Y_i(t') - E[Y_0(t) - Y_0(t')]]`

3. **Comparison cohort** (if unit i is in cohort g'):
   `-(1/pi_{g'}) * [Y_i(t') - Y_i(period_1) - E[Y_{g'}(t') - Y_{g'}(period_1)]]`

With the new rule, component 3 always uses a finite `g'`. The `else` branch in the EIF loop (lines 341-351 in edid-nocov.R) handling `gp_j == Inf` is dead code for PT-All.

### Generated Outcome for pair `(g', t')` — PT-All

```
y_hat_j = E[Y_g(t) - Y_g(1)] - E[Y_inf(t) - Y_inf(t')] - E[Y_{g'}(t') - Y_{g'}(1)]
```

With the new rule, the `else` branch (lines 245-247 in edid-nocov.R) for `is.infinite(gp_j)` is dead code for PT-All.

---

## Comprehension Self-Test Answers

1. **Core requirement restated**: See above. ✓
2. **Formulas written from memory**: See above. Verified against source. ✓
3. **Undefined symbols**: None. All symbols defined. ✓
4. **Judgment calls not in source**: None. The tpre rule is explicit in author's filter logic. ✓
5. **Mathematical intuition**: The `(g'=g, t'=period_1)` pair is the degenerate case where the comparison cohort provides no pre-trend information (its baseline change is 0). Including it exactly once via the `g'==g` branch gives one clean CS DiD moment. The `g'≠g` pairs provide additional identifying variation from other treated cohorts' pre-treatment differences from period_1, which is only meaningful for `t' > period_1`. Never-treated as comparison cohort is redundant because it is already the time control in every moment. ✓
6. **Implicit assumptions**: The fix assumes `never_treated_val = Inf` throughout (already the default). The `last_cohort` control group relabels the last cohort as Inf in `panel_obj$unit_cohorts` and removes it from `treatment_groups`, so `enumerate_valid_pairs_edid` naturally handles it: the relabeled cohort is absent from `treatment_groups`, and `control_group = "last_cohort"` works without any additional changes. ✓

---

## No HOLD Raised — Proceeding to Spec Production
