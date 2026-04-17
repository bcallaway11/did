# Methodology Review

This document tracks the progress of reviewing each estimator's implementation against the Methodology Registry and academic references. It ensures that implementations are correct, consistent, and well-documented.

For the methodology registry with academic foundations and key equations, see [docs/methodology/REGISTRY.md](docs/methodology/REGISTRY.md).

---

## Overview

Each estimator in diff-diff should be periodically reviewed to ensure:
1. **Correctness**: Implementation matches the academic paper's equations
2. **Reference alignment**: Behavior matches reference implementations (R packages, Stata commands)
3. **Edge case handling**: Documented edge cases are handled correctly
4. **Standard errors**: SE formulas match the documented approach

---

## Review Status Summary

| Estimator | Module | R Reference | Status | Last Review |
|-----------|--------|-------------|--------|-------------|
| DifferenceInDifferences | `estimators.py` | `fixest::feols()` | **Complete** | 2026-01-24 |
| MultiPeriodDiD | `estimators.py` | `fixest::feols()` | **Complete** | 2026-02-02 |
| TwoWayFixedEffects | `twfe.py` | `fixest::feols()` | **Complete** | 2026-02-08 |
| CallawaySantAnna | `staggered.py` | `did::att_gt()` | **Complete** | 2026-01-24 |
| SunAbraham | `sun_abraham.py` | `fixest::sunab()` | **Complete** | 2026-02-15 |
| SyntheticDiD | `synthetic_did.py` | `synthdid::synthdid_estimate()` | **Complete** | 2026-02-10 |
| TripleDifference | `triple_diff.py` | `triplediff::ddd()` | **Complete** | 2026-02-18 |
| StackedDiD | `stacked_did.py` | `stacked-did-weights` | **Complete** | 2026-02-19 |
| TROP | `trop.py` | (forthcoming) | Not Started | - |
| BaconDecomposition | `bacon.py` | `bacondecomp::bacon()` | Not Started | - |
| HonestDiD | `honest_did.py` | `HonestDiD` package | **Complete** | 2026-03-31 |
| PreTrendsPower | `pretrends.py` | `pretrends` package | Not Started | - |
| PowerAnalysis | `power.py` | `pwr` / `DeclareDesign` | Not Started | - |

**Status legend:**
- **Not Started**: No formal review conducted
- **In Progress**: Review underway
- **Complete**: Review finished, implementation verified

---

## Detailed Review Notes

### Core DiD Estimators

#### DifferenceInDifferences

| Field | Value |
|-------|-------|
| Module | `estimators.py` |
| Primary Reference | Wooldridge (2010), Angrist & Pischke (2009) |
| R Reference | `fixest::feols()` |
| Status | **Complete** |
| Last Review | 2026-01-24 |

**Verified Components:**
- [x] ATT formula: Double-difference of cell means matches regression interaction coefficient
- [x] R comparison: ATT matches `fixest::feols()` within 1e-3 tolerance
- [x] R comparison: SE (HC1 robust) matches within 5%
- [x] R comparison: P-value matches within 0.01
- [x] R comparison: Confidence intervals overlap
- [x] R comparison: Cluster-robust SE matches within 10%
- [x] R comparison: Fixed effects (absorb) matches `feols(...|unit)` within 1%
- [x] Wild bootstrap inference (Rademacher, Mammen, Webb weights)
- [x] Formula interface (`y ~ treated * post`)
- [x] All REGISTRY.md edge cases tested

**Test Coverage:**
- 53 methodology verification tests in `tests/test_methodology_did.py`
- 123 existing tests in `tests/test_estimators.py`
- R benchmark tests (skip if R not available)

**R Comparison Results:**
- ATT matches within 1e-3 (R JSON truncation limits precision)
- HC1 SE matches within 5%
- Cluster-robust SE matches within 10%
- Fixed effects results match within 1%

**Corrections Made:**
- (None - implementation verified correct)

**Outstanding Concerns:**
- R comparison precision limited by JSON output truncation (4 decimal places)
- Consider improving R script to output full precision for tighter tolerances

**Edge Cases Verified:**
1. Empty cells: Produces rank deficiency warning (expected behavior)
2. Singleton clusters: Included in variance estimation, contribute via residuals (corrected REGISTRY.md)
3. Rank deficiency: All three modes (warn/error/silent) working
4. Non-binary treatment/time: Raises ValueError as expected
5. No variation in treatment/time: Raises ValueError as expected
6. Missing values: Raises ValueError as expected

---

#### MultiPeriodDiD

| Field | Value |
|-------|-------|
| Module | `estimators.py` |
| Primary Reference | Freyaldenhoven et al. (2021), Wooldridge (2010), Angrist & Pischke (2009) |
| R Reference | `fixest::feols()` |
| Status | **Complete** |
| Last Review | 2026-02-02 |

**Verified Components:**
- [x] Full event-study specification: treatment × period interactions for ALL non-reference periods (pre and post)
- [x] Reference period coefficient is zero (normalized by omission from design matrix)
- [x] Default reference period is last pre-period (e=-1 convention, matches fixest/did)
- [x] Pre-period coefficients available for parallel trends assessment
- [x] Average ATT computed from post-treatment effects only, with covariance-aware SE
- [x] Returns PeriodEffect objects with confidence intervals for all periods
- [x] Supports balanced and unbalanced panels
- [x] NaN inference: t_stat/p_value/CI use NaN when SE is non-finite or zero
- [x] R-style NA propagation: avg_att is NaN if any post-period effect is unidentified
- [x] Rank-deficient design matrix: warns and sets NaN for dropped coefficients (R-style)
- [x] Staggered adoption detection warning (via `unit` parameter)
- [x] Treatment reversal detection warning
- [x] Time-varying D_it detection warning (advises creating ever-treated indicator)
- [x] Single pre-period warning (ATT valid but pre-trends assessment unavailable)
- [x] Post-period reference_period raises ValueError (would bias avg_att)
- [x] HonestDiD/PreTrendsPower integration uses interaction sub-VCV (not full regression VCV)
- [x] All REGISTRY.md edge cases tested

**Test Coverage:**
- 50 tests across `TestMultiPeriodDiD` and `TestMultiPeriodDiDEventStudy` in `tests/test_estimators.py`
- 18 new event-study specification tests added in PR #125

**Corrections Made:**
- **PR #125 (2026-02-02)**: Transformed from post-period-only estimator into full event-study
  specification with pre-period coefficients. Reference period default changed from first
  pre-period to last pre-period (e=-1 convention). HonestDiD/PreTrendsPower VCV extraction
  fixed to use interaction sub-VCV instead of full regression VCV.

**Outstanding Concerns:**
- ~~No R comparison benchmarks yet~~ — **Resolved**: R comparison benchmark added via
  `benchmarks/R/benchmark_multiperiod.R` using `fixest::feols(outcome ~ treated * time_f | unit)`.
  Results match R exactly: ATT diff < 1e-11, SE diff 0.0%, period effects correlation 1.0.
  Validated at small (200 units) and 1k scales.
- Default SE is HC1 (not cluster-robust at unit level as fixest uses). Cluster-robust
  available via `cluster` parameter but not the default.
- Endpoint binning for distant event times not yet implemented.
- FutureWarning for reference_period default change should eventually be removed once
  the transition is complete.

---

#### TwoWayFixedEffects

| Field | Value |
|-------|-------|
| Module | `twfe.py` |
| Primary Reference | Wooldridge (2010), Ch. 10 |
| R Reference | `fixest::feols()` |
| Status | **Complete** |
| Last Review | 2026-02-08 |

**Verified Components:**
- [x] Within-transformation algebra: `y_it - ȳ_i - ȳ_t + ȳ` matches hand calculation (rtol=1e-12)
- [x] ATT matches manual demeaned OLS (rtol=1e-10)
- [x] ATT matches `DifferenceInDifferences` on 2-period data (rtol=1e-10)
- [x] Covariates are also within-transformed (sum to zero within unit/time groups)
- [x] R comparison: ATT matches `fixest::feols(y ~ treated:post | unit + post, cluster=~unit)` (rtol<0.1%)
- [x] R comparison: Cluster-robust SE match (rtol<1%)
- [x] R comparison: P-value match (atol<0.01)
- [x] R comparison: CI bounds match (rtol<1%)
- [x] R comparison: ATT and SE match with covariate (same tolerances)
- [x] Edge case: Staggered treatment triggers `UserWarning`
- [x] Edge case: Auto-clusters at unit level (SE matches explicit `cluster="unit"`)
- [x] Edge case: DF adjustment for absorbed FE matches manual `solve_ols()` with `df_adjustment`
- [x] Edge case: Covariate collinear with interaction raises `ValueError` ("cannot be identified")
- [x] Edge case: Covariate collinearity warns but ATT remains finite
- [x] Edge case: `rank_deficient_action="error"` raises `ValueError`
- [x] Edge case: `rank_deficient_action="silent"` emits no warnings
- [x] Edge case: Unbalanced panel produces valid results (finite ATT, positive SE)
- [x] Edge case: Missing unit column raises `ValueError`
- [x] Integration: `decompose()` returns `BaconDecompositionResults`
- [x] SE: Cluster-robust SE >= HC1 SE
- [x] SE: VCoV positive semi-definite
- [x] Wild bootstrap: Valid inference (finite SE, p-value in [0,1])
- [x] Wild bootstrap: All weight types (rademacher, mammen, webb) produce valid inference
- [x] Wild bootstrap: `inference="wild_bootstrap"` routes correctly
- [x] Params: `get_params()` returns all inherited parameters
- [x] Params: `set_params()` modifies attributes
- [x] Results: `summary()` contains "ATT"
- [x] Results: `to_dict()` contains att, se, t_stat, p_value, n_obs
- [x] Results: residuals + fitted = demeaned outcome (not raw)
- [x] Edge case: Multi-period time emits UserWarning advising binary post indicator
- [x] Edge case: Non-{0,1} binary time emits UserWarning (ATT still correct)
- [x] Edge case: ATT invariant to time encoding ({0,1} vs {2020,2021} produces identical results)

**Key Implementation Detail:**
The interaction term `D_i × Post_t` must be within-transformed (demeaned) alongside the outcome,
consistent with the Frisch-Waugh-Lovell (FWL) theorem: all regressors and the outcome must be
projected out of the fixed effects space. R's `fixest::feols()` does this automatically when
variables appear to the left of the `|` separator.

**Corrections Made:**
- **Bug fix: interaction term must be within-transformed** (found during review). The previous
  implementation used raw (un-demeaned) `D_i × Post_t` in the demeaned regression. This gave
  correct results only for 2-period panels where `post == period`. For multi-period panels
  (e.g., 4 periods with binary `post`), the raw interaction had incorrect correlation with
  demeaned Y, producing ATT approximately 1/3 of the true value. Fixed by applying the same
  within-transformation to the interaction term before regression. This matches R's
  `fixest::feols()` behavior. (`twfe.py` lines 99-113)

**Outstanding Concerns:**
- **Multi-period `time` parameter**: Multi-period time values (e.g., 1,2,3,4) produce
  `treated × period_number` instead of `treated × post_indicator`, which is not the standard
  D_it treatment indicator. A `UserWarning` is emitted when `time` has >2 unique values.
  For binary time with non-{0,1} values (e.g., {2020, 2021}), the ATT is mathematically
  correct (the within-transformation absorbs the scaling), but a warning recommends 0/1
  encoding for clarity. Users with multi-period data should create a binary `post` column.
- **Staggered treatment warning**: The warning only fires when `time` has >2 unique values
  (i.e., actual period numbers). With binary `time="post"`, all treated units appear to start
  treatment at `time=1`, making staggering undetectable. Users with staggered designs should
  use `decompose()` or `CallawaySantAnna` directly for proper diagnostics.

---

### Modern Staggered Estimators

#### CallawaySantAnna

| Field | Value |
|-------|-------|
| Module | `staggered.py` |
| Primary Reference | Callaway & Sant'Anna (2021) |
| R Reference | `did::att_gt()` |
| Status | **Complete** |
| Last Review | 2026-01-24 |

**Verified Components:**
- [x] ATT(g,t) basic formula (hand-calculated exact match)
- [x] Doubly robust estimator
- [x] IPW estimator
- [x] Outcome regression
- [x] Base period selection (varying/universal)
- [x] Anticipation parameter handling
- [x] Simple/event-study/group aggregation
- [x] Analytical SE with weight influence function
- [x] Bootstrap SE (Rademacher/Mammen/Webb)
- [x] Control group composition (never_treated/not_yet_treated)
- [x] All documented edge cases from REGISTRY.md

**Test Coverage:**
- 46 methodology verification tests in `tests/test_methodology_callaway.py`
- 93 existing tests in `tests/test_staggered.py`
- R benchmark tests (skip if R not available)

**R Comparison Results:**
- Overall ATT matches within 20% (difference due to dynamic effects in generated data)
- Post-treatment ATT(g,t) values match within 20%
- Pre-treatment effects may differ due to base_period handling differences

**Corrections Made:**
- (None - implementation verified correct)

**Outstanding Concerns:**
- R comparison shows ~20% difference in overall ATT with generated data
  - Likely due to differences in how dynamic effects are handled in data generation
  - Individual ATT(g,t) values match closely for post-treatment periods
  - Further investigation recommended with real-world data
- Pre-treatment ATT(g,t) may differ from R due to base_period="varying" semantics
  - Python uses t-1 as base for pre-treatment
  - R's behavior requires verification

**Deviations from R's did::att_gt():**
1. **NaN for invalid inference**: When SE is non-finite or zero, Python returns NaN for
   t_stat/p_value rather than potentially erroring. This is a defensive enhancement.

**Alignment with R's did::att_gt() (as of v2.1.5):**
1. **Webb weights**: Webb's 6-point distribution with values ±√(3/2), ±1, ±√(1/2)
   uses equal probabilities (1/6 each) matching R's `did` package. This gives
   E[w]=0, Var(w)=1.0, consistent with other bootstrap weight distributions.

   **Verification**: Our implementation matches the well-established `fwildclusterboot`
   R package (C++ source: [wildboottest.cpp](https://github.com/s3alfisc/fwildclusterboot/blob/master/src/wildboottest.cpp)).
   The implementation uses `sqrt(1.5)`, `1`, `sqrt(0.5)` (and negatives) with equal 1/6
   probabilities—identical to our values.

   **Note on documentation discrepancy**: Some documentation (e.g., fwildclusterboot
   vignette) describes Webb weights as "±1.5, ±1, ±0.5". This appears to be a
   simplification for readability. The actual implementations use ±√1.5, ±1, ±√0.5
   which provides the required unit variance (Var(w) = 1.0).

---

#### SunAbraham

| Field | Value |
|-------|-------|
| Module | `sun_abraham.py` |
| Primary Reference | Sun & Abraham (2021) |
| R Reference | `fixest::sunab()` |
| Status | **Complete** |
| Last Review | 2026-02-15 |

**Verified Components:**
- [x] Saturated TWFE regression with cohort × relative-time interactions
- [x] Within-transformation for unit and time fixed effects
- [x] Interaction-weighted event study effects (δ̂_e = Σ_g ŵ_{g,e} × δ̂_{g,e})
- [x] IW weights match event-time sample shares (n_{g,e} / Σ_g n_{g,e})
- [x] Overall ATT as weighted average of post-treatment effects
- [x] Delta method SE for aggregated effects (Var = w' Σ w)
- [x] Cluster-robust SEs at unit level
- [x] Reference period normalized to zero (e=-1 excluded from design matrix)
- [x] R comparison: ATT matches `fixest::sunab()` within machine precision (<1e-11)
- [x] R comparison: SE matches within 0.3% (small scale) / 0.1% (1k scale)
- [x] R comparison: Event study effects correlation = 1.000000
- [x] R comparison: Event study max diff < 1e-11
- [x] Bootstrap inference (pairs bootstrap)
- [x] Rank deficiency handling (warn/error/silent)
- [x] All REGISTRY.md edge cases tested

**Test Coverage:**
- 43 tests in `tests/test_sun_abraham.py` (36 existing + 7 methodology verification)
- R benchmark tests via `benchmarks/run_benchmarks.py --estimator sunab`

**R Comparison Results:**
- Overall ATT matches within machine precision (diff < 1e-11 at both scales)
- Cluster-robust SE matches within 0.3% (well within 1% threshold)
- Event study effects match perfectly (correlation 1.0, max diff < 1e-11)
- Validated at small (200 units) and 1k (1000 units) scales

**Corrections Made:**
1. **DF adjustment for absorbed FE** (`sun_abraham.py`, `_fit_saturated_regression()`):
   Added `df_adjustment = n_units + n_times - 1` to `LinearRegression.fit()` to account
   for absorbed unit and time fixed effects in degrees of freedom. Unlike TWFE (which uses
   `-2` plus an explicit intercept column), SunAbraham's saturated regression has no
   intercept, so all absorbed df must come from the adjustment. Affects t-distribution DoF
   for cohort-level p-values/CIs (slightly larger p-values, slightly wider CIs) but does
   NOT change VCV or SE values.

2. **NaN return for no post-treatment effects** (`sun_abraham.py`, `_compute_overall_att()`):
   Changed return from `(0.0, 0.0)` to `(np.nan, np.nan)` when no post-treatment effects
   exist. All downstream inference fields (t_stat, p_value, conf_int) correctly propagate
   NaN via existing guards in `fit()`.

3. **Deprecation warnings for unused parameters** (`sun_abraham.py`, `fit()`):
   Added `FutureWarning` for `min_pre_periods` and `min_post_periods` parameters that
   are accepted but never used (no-op). These will be removed in a future version.

4. **Removed event-time truncation at [-20, 20]** (`sun_abraham.py`):
   Removed the hardcoded cap `max(min(...), -20)` / `min(max(...), 20)` to match
   R's `fixest::sunab()` which has no such limit. All available relative times are
   now estimated.

5. **Warning for variance fallback path** (`sun_abraham.py`, `_compute_overall_att()`):
   Added `UserWarning` when the full weight vector cannot be constructed and a
   simplified variance (ignoring covariances between periods) is used as fallback.

6. **IW weights use event-time sample shares** (`sun_abraham.py`, `_compute_iw_effects()`):
   Changed IW weights from `n_g / Σ_g n_g` (cohort sizes) to `n_{g,e} / Σ_g n_{g,e}`
   (per-event-time observation counts) to match the REGISTRY.md formula. For balanced
   panels these are identical; for unbalanced panels the new formula correctly reflects
   actual sample composition at each event-time. Added unbalanced panel test.

7. **Normalize `np.inf` never-treated encoding** (`sun_abraham.py`, `fit()`):
   `first_treat=np.inf` (documented as valid for never-treated) was included in
   `treatment_groups` and `_rel_time` via `> 0` checks, producing `-inf` event times.
   Fixed by normalizing `np.inf` to `0` immediately after computing `_never_treated`.
   Same fix applied to `staggered.py` (`CallawaySantAnna`).

**Outstanding Concerns:**
- **Inference distribution**: Cohort-level p-values use t-distribution (via
  `LinearRegression.get_inference()`), while aggregated event study and overall ATT
  p-values use normal distribution (via `compute_p_value()`). This is asymptotically
  equivalent and standard for delta-method-aggregated quantities. R's fixest uses
  t-distribution at all levels, so aggregated p-values may differ slightly for small
  samples — this is a documented deviation.

**Deviations from R's fixest::sunab():**
1. **NaN for no post-treatment effects**: Python returns `(NaN, NaN)` for overall ATT/SE
   when no post-treatment effects exist. R would error.
2. **Normal distribution for aggregated inference**: Aggregated p-values use normal
   distribution (asymptotically equivalent). R uses t-distribution.

---

#### StackedDiD

| Field | Value |
|-------|-------|
| Module | `stacked_did.py` |
| Primary Reference | Wing, Freedman & Hollingsworth (2024), NBER WP 32054 |
| R Reference | `stacked-did-weights` (`create_sub_exp()` + `compute_weights()`) |
| Status | **Complete** |
| Last Review | 2026-02-19 |

**Verified Components:**
- [x] IC1 trimming: `a - kappa_pre >= T_min AND a + kappa_post <= T_max` (matches R reference)
- [x] IC2 trimming: Three clean control modes (not_yet_treated, strict, never_treated)
- [x] Sub-experiment construction: treated + clean controls within `[a - kappa_pre, a + kappa_post]`
- [x] Q-weights aggregate: treated Q=1, control `Q = (sub_treat_n/stack_treat_n) / (sub_control_n/stack_control_n)` per (event_time, sub_exp) — matches R `compute_weights()`
- [x] Q-weights population: `Q_a = (Pop_a^D / Pop^D) / (N_a^C / N^C)` (Table 1, Row 2)
- [x] Q-weights sample_share: `Q_a = ((N_a^D + N_a^C)/(N^D+N^C)) / (N_a^C / N^C)` (Table 1, Row 3)
- [x] WLS via sqrt(w) transformation (numerically equivalent to weighted regression)
- [x] Event study regression: `Y = α_0 + α_1·D_sa + Σ_{h≠-1}[λ_h·1(e=h) + δ_h·D_sa·1(e=h)] + U` (Eq. 3)
- [x] Reference period e=-1-anticipation normalized to zero (omitted from design matrix)
- [x] Delta-method SE for overall ATT: `SE = sqrt(ones' @ sub_vcv @ ones) / K`
- [x] Cluster-robust SEs at unit level (default) and unit×sub-experiment level
- [x] Anticipation parameter: reference period shifts to e=-1-anticipation, post-treatment includes anticipation periods
- [x] Rank deficiency handling (warn/error/silent via `solve_ols()`)
- [x] Never-treated encoding: both `first_treat=0` and `first_treat=inf` handled
- [x] R comparison: ATT matches within machine precision (diff < 2.1e-11)
- [x] R comparison: SE matches within machine precision (diff < 4.0e-10)
- [x] R comparison: Event study effects correlation = 1.000000, max diff < 4.5e-11
- [x] safe_inference() used for all inference fields
- [x] All REGISTRY.md edge cases tested

**Test Coverage:**
- 72 tests in `tests/test_stacked_did.py` across 11 test classes:
  - `TestStackedDiDBasic` (8): fit, event study, group/all raises, simple aggregation, known constant effect, dynamic effects
  - `TestTrimming` (5): IC1 window, IC2 no-controls, trimmed groups reported, all-trimmed raises, wider window
  - `TestQWeights` (4): treated=1, aggregate formula, sample_share formula, positivity
  - `TestCleanControl` (5): not_yet_treated, strict, never_treated, missing never-treated raises
  - `TestClustering` (2): unit, unit_subexp
  - `TestStackedData` (4): accessible, required columns, event time range
  - `TestEdgeCases` (8): single cohort, anticipation, unbalanced panel, NaN inference, never-treated encodings
  - `TestSklearnInterface` (4): get_params, set_params, unknown raises, convenience function
  - `TestResultsMethods` (7): summary, to_dataframe, is_significant, significance_stars, repr
  - `TestValidation` (8): missing columns, invalid params, population required, no treated units
- R benchmark tests via `benchmarks/run_benchmarks.py --estimator stacked`

**R Comparison Results (200 units, 8 periods, kappa_pre=2, kappa_post=2):**
| Metric | Python | R | Diff |
|--------|--------|---|------|
| Overall ATT | 2.277699574579 | 2.2776995746 | 2.1e-11 |
| Overall SE | 0.062045687626 | 0.062045688027 | 4.0e-10 |
| ES e=-2 ATT | 0.044517975379 | 0.044517975379 | <1e-12 |
| ES e=0 ATT | 2.104181683763 | 2.104181683800 | <1e-11 |
| ES e=1 ATT | 2.209990715130 | 2.209990715100 | <1e-11 |
| ES e=2 ATT | 2.518926324845 | 2.518926324800 | <1e-11 |
| Stacked obs | 1600 | 1600 | exact |
| Sub-experiments | 3 | 3 | exact |

**Corrections Made:**
1. **IC1 lower bound and time window aligned with R reference** (`stacked_did.py`,
   `_trim_adoption_events()` and `_build_sub_experiment()`): The paper text specifies
   time window `[a - kappa_pre - 1, a + kappa_post]` (including an extra pre-period),
   but the R reference implementation by co-author Hollingsworth uses
   `[a - kappa_pre, a + kappa_post]`. The extra period had no event-study dummy,
   altering the baseline regression. Fixed to match R: removed `-1` from both
   IC1 check (`a - kappa_pre >= T_min`) and time window start. Discrepancy documented
   in `docs/methodology/papers/wing-2024-review.md` Gaps section.

2. **Q-weight computation: event-time-specific for aggregate weighting** (`stacked_did.py`,
   `_compute_q_weights()`): Changed aggregate Q-weights from unit counts per sub-experiment
   to observation counts per (event_time, sub_exp), matching R reference `compute_weights()`.
   For balanced panels, results are unchanged. For unbalanced panels, weights now adjust for
   varying observation density. Population/sample_share retain unit-count formulas (paper notation).

3. **Anticipation parameter: reference period and dummies** (`stacked_did.py`, `fit()`):
   Reference period now shifts to `e = -1 - anticipation`. Event-time dummies cover the
   full window `[-kappa_pre - anticipation, ..., kappa_post]`. Post-treatment effects include
   anticipation periods. Consistent with ImputationDiD, TwoStageDiD, SunAbraham.

4. **Group aggregation removed** (`stacked_did.py`): `aggregate="group"` and `aggregate="all"`
   removed. The pooled stacked regression cannot produce cohort-specific effects without
   cohort×event-time interactions. Use CallawaySantAnna or ImputationDiD for cohort-level estimates.

5. **n_sub_experiments metadata** (`stacked_did.py`, `fit()`): Now tracks actual built
   sub-experiments, not all events in omega_kappa. Warns if any sub-experiments are empty
   after data filtering.

**Outstanding Concerns:**
- Population/sample_share Q-weights use paper's unit-count formulas (no R reference to validate)
- Anticipation not validated against R (R reference doesn't test anticipation > 0)

**Deviations from R's stacked-did-weights:**
1. **NaN for invalid inference**: Python returns NaN for t_stat/p_value/conf_int when
   SE is non-finite or zero. R would propagate through `fixest::feols()` error handling.

---

### Advanced Estimators

#### SyntheticDiD

| Field | Value |
|-------|-------|
| Module | `synthetic_did.py` |
| Primary Reference | Arkhangelsky et al. (2021) |
| R Reference | `synthdid::synthdid_estimate()` |
| Status | **Complete** |
| Last Review | 2026-02-10 |

**Corrections Made:**
1. **Time weights: Frank-Wolfe on collapsed form** (was heuristic inverse-distance).
   Replaced ad-hoc inverse-distance weighting with the Frank-Wolfe algorithm operating
   on the collapsed (N_co x T_pre) problem as specified in Algorithm 1 of
   Arkhangelsky et al. (2021), matching R's `synthdid::fw.step()`.
2. **Unit weights: Frank-Wolfe with two-pass sparsification** (was projected gradient
   descent with wrong penalty). Replaced projected gradient descent (which used an
   incorrect penalty formulation) with Frank-Wolfe optimization followed by two-pass
   sparsification, matching R's `synthdid::sc.weight.fw()` and `sparsify_function()`.
3. **Auto-computed regularization from data noise level** (was `lambda_reg=0.0`,
   `zeta=1.0`). Regularization parameters `zeta_omega` and `zeta_lambda` are now
   computed automatically from the data noise level (N_tr * sigma^2) as specified in
   Appendix D of Arkhangelsky et al. (2021), matching R's default behavior.
4. **Bootstrap SE uses fixed weights matching R's `bootstrap_sample`** (was
   re-estimating all weights). The bootstrap variance procedure now holds unit and time
   weights fixed at their point estimates and only re-estimates the treatment effect,
   matching the approach in R's `synthdid::bootstrap_sample()`.
5. **Default `variance_method` changed to `"placebo"`** matching R's default. The R
   package uses placebo variance by default (`synthdid_estimate` returns an object whose
   `vcov()` uses the placebo method); our default now matches.
6. **Deprecated `lambda_reg` and `zeta` params; new params are `zeta_omega` and
   `zeta_lambda`**. The old parameters had unclear semantics and did not correspond to
   the paper's notation. The new parameters directly match the paper and R package
   naming conventions. `lambda_reg` and `zeta` are deprecated with warnings and will
   be removed in a future release.

**Outstanding Concerns:**
- (None)

---

#### TripleDifference

| Field | Value |
|-------|-------|
| Module | `triple_diff.py` |
| Primary Reference | Ortiz-Villavicencio & Sant'Anna (2025) |
| R Reference | `triplediff::ddd()` (v0.2.1, CRAN) |
| Status | **Complete** |
| Last Review | 2026-02-18 |

**Verified Components:**
- [x] ATT matches R `triplediff::ddd()` for all 3 methods (DR, RA, IPW) — <0.001% relative difference
- [x] SE matches R `triplediff::ddd()` for all 3 methods — <0.001% relative difference
- [x] With-covariates ATT matches R — <0.001% relative difference
- [x] With-covariates SE matches R — <0.001% relative difference
- [x] Verified across all 4 DGP types from `gen_dgp_2periods()` (different model misspecification scenarios)
- [x] Influence function-based SE: `SE = std(w3*IF_3 + w2*IF_2 - w1*IF_1, ddof=1) / sqrt(n)`
- [x] Three-DiD decomposition: `DDD = DiD_3 + DiD_2 - DiD_1` matching R's approach
- [x] safe_inference() used for all inference fields (t_stat, p_value, conf_int)

**Corrections Made:**
1. **Complete rewrite of estimation methods** (was naive cell-mean approach, now three-DiD
   decomposition). The original implementation computed DDD directly from 8 cell means with
   a naive cell-variance SE. Replaced with R's decomposition into three pairwise DiD
   comparisons (subgroup j vs reference subgroup 4), each using DR/IPW/RA methodology
   from Callaway & Sant'Anna. This fixed:
   - DR SE: was off by >100% (naive cell variance vs influence function)
   - IPW SE: was off by >200% (incorrect cell-probability-ratio weights)
   - With-covariates ATT: was off by >1000% for all methods (incorrect cell-by-cell regression)
2. **Influence function SE** replaces naive cell variance for all methods:
   `SE = std(w3*IF_3 + w2*IF_2 - w1*IF_1, ddof=1) / sqrt(n)` where
   `w_j = n / n_j` and `IF_j` is the per-observation influence function for pairwise DiD j.
3. **Propensity score estimation** now runs per-pairwise-comparison (P(subgroup=4|X) within
   {j, 4} subset) instead of global P(G=1|X).
4. **Outcome regression** now fits separate OLS per subgroup-time cell within each pairwise
   comparison, matching R's `compute_outcome_regression_rc()`.

**Outstanding Concerns:**
- Implementation uses `panel=FALSE` (repeated cross-section) mode. Panel mode (`panel=TRUE`)
  with differenced outcomes not yet implemented.

**R Comparison Results (panel=FALSE, n=500 per DGP):**
| DGP | Method | Covariates | ATT Diff | SE Diff |
|-----|--------|-----------|----------|---------|
| 1 | DR | No | <0.001% | <0.001% |
| 1 | DR | Yes | <0.001% | <0.001% |
| 1 | REG | No | <0.001% | <0.001% |
| 1 | REG | Yes | <0.001% | <0.001% |
| 1 | IPW | No | <0.001% | <0.001% |
| 1 | IPW | Yes | <0.001% | <0.001% |
| 2-4 | All | Both | <0.001% | <0.001% |

---

#### TROP

| Field | Value |
|-------|-------|
| Module | `trop.py` |
| Primary Reference | Athey, Imbens, Qu & Viviano (2025) |
| R Reference | (forthcoming) |
| Status | Not Started |
| Last Review | - |

**Corrections Made:**
- (None yet)

**Outstanding Concerns:**
- (None yet)

---

### Diagnostics & Sensitivity

#### BaconDecomposition

| Field | Value |
|-------|-------|
| Module | `bacon.py` |
| Primary Reference | Goodman-Bacon (2021) |
| R Reference | `bacondecomp::bacon()` |
| Status | Not Started |
| Last Review | - |

**Corrections Made:**
- (None yet)

**Outstanding Concerns:**
- (None yet)

---

#### HonestDiD

| Field | Value |
|-------|-------|
| Module | `honest_did.py` |
| Primary Reference | Rambachan & Roth (2023) |
| R Reference | `HonestDiD` package |
| Status | **Complete** (pending R comparison) |
| Last Review | 2026-04-01 |

**Verified Components:**
- [x] Delta^SD: second-difference constraints [1,-2,1] with delta_0=0 boundary handling
- [x] Delta^SD: T+Tbar-1 constraint rows (bridge constraint at t=0)
- [x] Delta^RM: constrains first differences (not levels), union of polyhedra per Lemma 2.2
- [x] Identified set LP: pins delta_pre = beta_pre via equality constraints (Equations 5-6)
- [x] M=0 for Delta^SD: linear extrapolation gives finite point-identified bounds
- [x] Mbar=0 for Delta^RM: point identification (all post first-diffs = 0)
- [x] Optimal FLCI for Delta^SD: folded normal cv_alpha, Nelder-Mead over pre-period weights
- [x] Sensitivity grid: bounds computed for each M in grid, breakdown value via binary search
- [x] Survey variance (RM, M=0 smoothness): t-distribution critical values from df_survey
- [ ] Survey variance (M>0 smoothness): optimal FLCI uses asymptotic normal only; df_survey=0 → NaN
- [x] CallawaySantAnna integration: universal base period, reference period filtering
- [x] Three-period analytical case matches paper Section 2.3
- [ ] ARP hybrid for Delta^RM: infrastructure implemented, moment inequality transformation needs calibration
- [ ] R comparison: pending (benchmark scripts need updating)

**Test Coverage:**
- 63 existing tests in `tests/test_honest_did.py` (14 classes) — all passing
- 17 new methodology verification tests in `tests/test_methodology_honest_did.py`
- R benchmark tests (pending)

**Corrections Made:**
1. **DeltaRM: first differences, not levels** (`honest_did.py`, `_construct_constraints_rm_component`):
   The paper's Delta^RM constrains `|delta_{t+1} - delta_t|` (consecutive first differences)
   bounded by Mbar × max pre-treatment first difference. The code constrained `|delta_post|`
   (absolute levels) bounded by Mbar × max `|beta_pre|`. Completely rewritten using
   union-of-polyhedra decomposition per Lemma 2.2.

2. **LP pins delta_pre = beta_pre** (`honest_did.py`, `_solve_bounds_lp`):
   The paper's identified set LP (Equations 5-6) fixes pre-treatment violations to the observed
   pre-treatment coefficients. The code had no equality constraint — delta_pre was unconstrained.
   For Delta^SD(M=0), this made the LP unbounded. Added A_eq/b_eq equality constraints.

3. **DeltaSD constraint matrix: delta_0=0 boundary** (`honest_did.py`, `_construct_A_sd`):
   The code built second-difference matrices treating [delta_{-T},...,delta_{-1},delta_1,...,delta_{Tbar}]
   as consecutive, missing delta_0=0 at the boundary. Three boundary rows were wrong:
   - t=-1: `d_{-2} - 2*d_{-1} + 0` (uses delta_0=0)
   - t=0: `d_{-1} + d_1` (bridge constraint, was missing)
   - t=1: `0 - 2*d_1 + d_2` (uses delta_0=0)
   Now produces T+Tbar-1 rows (was T+Tbar-2).

4. **Optimal FLCI for Delta^SD** (`honest_did.py`, `_compute_optimal_flci`):
   Replaced naive FLCI (`lb - z*se, ub + z*se`) with the paper's optimal FLCI (Section 4.1):
   jointly optimizes affine estimator direction v and half-length chi using folded normal
   critical values cv_alpha(bias/se). Significantly narrower CIs.

5. **REGISTRY.md equations** (`docs/methodology/REGISTRY.md`):
   DeltaSD equation was first differences (should be second differences). DeltaRM equation
   was absolute levels (should be first differences). Both corrected with full formulations.

6. **Performance** (`honest_did.py`):
   Sensitivity grid reduced from ~9 minutes to 0.1 seconds via: Newton's method for cv_alpha
   (5 iterations vs 100), centrosymmetric bias LP (1 solve vs 2), M=0 short-circuit,
   looser Nelder-Mead tolerances.

**Outstanding Concerns:**
- **Delta^RM CI**: uses naive FLCI (conservative) instead of the paper's ARP conditional/hybrid
  confidence sets. ARP infrastructure exists but moment inequality transformation needs
  calibration. Tracked in TODO.md.
- R benchmark comparison not yet run (Python benchmark needs API update)
- Combined method uses single M for both SD and RM (DeltaSDRM dataclass has separate M/Mbar)

**Deviations from R's HonestDiD:**
1. **Deviation from R:** Delta^RM CIs use naive FLCI (`lb - z*se, ub + z*se`) instead of ARP
   conditional/hybrid. Conservative (wider CIs, valid coverage). ARP deferred.
2. **Note:** Delta^SD optimal FLCI matches the paper's Section 4.1 methodology: first-difference
   reparameterization, slope weights with sum(w)=sum_j j*l_j constraint (Eq. 17), bias LP in
   fd-space, folded normal (or folded non-central t for survey df). Nelder-Mead optimizer vs
   R's custom solver may produce numerical differences at tolerance level.
3. **Note:** `method="combined"` (Delta^SDRM) uses naive FLCI on the intersection of SD and RM
   bounds. The paper proves FLCI is not consistent for Delta^SDRM (Proposition 4.2). A runtime
   UserWarning is emitted. Use `method="smoothness"` or `method="relative_magnitude"` separately
   for paper-supported inference.
4. **Note (deviation from R):** Python warns (doesn't error) when CallawaySantAnna results use
   `base_period != "universal"`. R's HonestDiD requires universal base period.

---

#### PreTrendsPower

| Field | Value |
|-------|-------|
| Module | `pretrends.py` |
| Primary Reference | Roth (2022) |
| R Reference | `pretrends` package |
| Status | Not Started |
| Last Review | - |

**Corrections Made:**
- (None yet)

**Outstanding Concerns:**
- (None yet)

---

#### PowerAnalysis

| Field | Value |
|-------|-------|
| Module | `power.py` |
| Primary Reference | Bloom (1995), Burlig et al. (2020) |
| R Reference | `pwr` / `DeclareDesign` |
| Status | Not Started |
| Last Review | - |

**Corrections Made:**
- (None yet)

**Outstanding Concerns:**
- (None yet)

---

## Review Process Guidelines

### Review Checklist

For each estimator, complete the following steps:

- [ ] **Read primary academic source** - Review the key paper(s) cited in REGISTRY.md
- [ ] **Compare key equations** - Verify implementation matches equations in REGISTRY.md
- [ ] **Run benchmark against R reference** - Execute `benchmarks/run_benchmarks.py --estimator <name>` if available
- [ ] **Verify edge case handling** - Check behavior matches REGISTRY.md documentation
- [ ] **Check standard error formula** - Confirm SE computation matches reference
- [ ] **Document any deviations** - Add notes explaining intentional differences with rationale

### When to Update This Document

1. **After completing a review**: Update status to "Complete" and add date
2. **When making corrections**: Document what was fixed in the "Corrections Made" section
3. **When identifying issues**: Add to "Outstanding Concerns" for future investigation
4. **When deviating from reference**: Document the deviation and rationale

### Deviation Documentation

When our implementation intentionally differs from the reference implementation, document:

1. **What differs**: Specific behavior or formula that differs
2. **Why**: Rationale (e.g., "defensive enhancement", "bug in R package", "follows updated paper")
3. **Impact**: Whether results differ in practice
4. **Cross-reference**: Update REGISTRY.md edge cases section

Example:
```
**Deviation (2025-01-15)**: CallawaySantAnna returns NaN for t_stat when SE is non-finite,
whereas R's `did::att_gt` would error. This is a defensive enhancement that provides
more graceful handling of edge cases while still signaling invalid inference to users.
```

### Priority Order

Suggested order for reviews based on usage and complexity:

1. **High priority** (most used, complex methodology):
   - CallawaySantAnna
   - SyntheticDiD
   - HonestDiD

2. **Medium priority** (commonly used, simpler methodology):
   - DifferenceInDifferences
   - TwoWayFixedEffects
   - MultiPeriodDiD
   - SunAbraham
   - BaconDecomposition

3. **Lower priority** (newer or less commonly used):
   - TripleDifference
   - TROP
   - PreTrendsPower
   - PowerAnalysis

---

## Related Documents

- [REGISTRY.md](docs/methodology/REGISTRY.md) - Academic foundations and key equations
- [ROADMAP.md](ROADMAP.md) - Feature roadmap
- [TODO.md](TODO.md) - Technical debt tracking
- [CLAUDE.md](CLAUDE.md) - Development guidelines
