# IMPLEMENTATION_PLAN

## 1. High-Level Overview

### What the Python library does

The relevant Python implementation is the `EfficientDiD` estimator in `diff_diff`, built around the papers by Chen, Sant'Anna, and Xie (2025). It estimates staggered-adoption group-time treatment effects `ATT(g, t)` and then optionally aggregates them into:

- an overall ATT
- event-study effects
- treatment-cohort averages
- a Hausman-style PT-All vs PT-Post pretest

It supports two estimation regimes:

- no covariates: closed-form efficient DiD using unconditional covariance matrices
- covariates: doubly robust efficient DiD using outcome regression, sieve propensity ratios, inverse propensities, and conditional covariance estimation

It also supports:

- analytical EIF-based standard errors
- cluster-robust analytical SEs
- multiplier bootstrap
- survey-design-based inference

### Core features of the DiD method implemented

The estimator is built around one core idea: for each target cell `(g, t)`, construct all valid DiD moments and optimally combine them using inverse-covariance weights.

Core method features:

- supports `PT-All` and `PT-Post`
- enumerates valid comparison moments `(g', t_pre)` per target `(g, t)`
- computes efficient weights `w = 1' Omega^{-1} / (1' Omega^{-1} 1)`
- uses efficient influence functions (EIFs) as the main inference primitive
- aggregates cell-level effects with weight-influence-function (WIF) corrections
- uses a plug-in doubly robust score in the covariate path

### Key design patterns used in the Python code

- thin orchestration layer in `efficient_did.py`
- pure mathematical helpers in dedicated modules
- cache-heavy nuisance estimation keyed by tuples
- matrix-first internal representation after validation
- results stored in immutable dataclass-like containers
- inference and bootstrap both driven from stored EIF vectors

For R, the key architectural takeaway is: keep a functional core, separate estimator math from orchestration, and represent fitted results as an S3 object with generics.

## 2. Architecture Mapping (Python -> R)

### Recommended R object model

Do not mirror the mutable Python class directly. The idiomatic R design should be:

- a primary fitting function `edid()`
- an S3 fitted object class `edid_fit`
- optional S3 helper class `edid_bootstrap`
- optional S3 helper class `edid_hausman`

Why S3 instead of R6/S4:

- the estimator is naturally "fit-once, return-result"
- methods like `print()`, `summary()`, `coef()`, `vcov()`, `as.data.frame()`, and `plot()` are easy with S3
- the Python class mostly stores configuration plus a single fit result, which is more naturally expressed in R as a function plus a returned object

Use R6 only if the package wants sklearn-like mutable estimators. That is not necessary here.

### Python module/class/function -> R mapping

| Python | Role | Recommended R equivalent |
|---|---|---|
| `EfficientDiD` class | main estimator | `edid()` fitting function + `edid_fit` S3 object |
| `EfficientDiDResults` | fitted result container | `edid_fit` list with class + S3 methods |
| `EDiDBootstrapResults` | bootstrap result container | `edid_bootstrap` list or nested list inside `edid_fit` |
| `efficient_did.py` | orchestration | `R/edid-fit.R`, `R/edid-aggregate.R`, `R/edid-hausman.R` |
| `efficient_did_weights.py` | no-covariate math | `R/edid-nocov.R` |
| `efficient_did_covariates.py` | DR covariate math | `R/edid-covariates.R` |
| `efficient_did_bootstrap.py` | multiplier bootstrap | `R/edid-bootstrap.R` |
| `efficient_did_results.py` | output formatting | `R/edid-methods.R` |
| `utils.safe_inference()` | scalar inference helper | `safe_inference_edid()` in `R/edid-inference.R` |
| `linalg.solve_ols()` | weighted OLS helper | `solve_ols_edid()` in `R/edid-linalg.R` |
| survey helpers | design-based variance | wrapper functions in `R/edid-survey.R` using `survey` package where possible |

### Recommended R package file structure

```text
R/
  edid.R
  edid-fit.R
  edid-validate.R
  edid-data.R
  edid-pairs.R
  edid-nocov.R
  edid-covariates.R
  edid-linalg.R
  edid-inference.R
  edid-aggregate.R
  edid-bootstrap.R
  edid-survey.R
  edid-hausman.R
  edid-methods.R
  utils.R

tests/testthat/
  test-edid-validate.R
  test-edid-pairs.R
  test-edid-nocov.R
  test-edid-covariates.R
  test-edid-aggregate.R
  test-edid-bootstrap.R
  test-edid-survey.R
  test-edid-hausman.R
  test-edid-parity.R

inst/extdata/
  parity-fixtures/

vignettes/
  efficient-did.Rmd
```

## 3. Core Components Breakdown

### A. Data preprocessing

Purpose:

- validate panel structure
- normalize treatment cohorts
- reshape long data into unit-level matrices
- build cohort masks and period mappings

Key Python functions/classes:

- `EfficientDiD.fit()`
- `_validate_and_build_cluster_mapping()`
- internal pivoting and cohort-fraction logic in `efficient_did.py`

Planned equivalent in R:

- `validate_edid_inputs()`
- `prepare_edid_panel()`
- `build_cluster_index()`

Notes on translation challenges:

- Python uses sorted unique values and `pivot`; R should use `data.table::dcast()` or base reshape plus explicit sorting
- Python uses `0` or `Inf` internally for never-treated; R should standardize on `Inf` at the API boundary, but may internally map to `0` or keep `Inf` consistently
- unit-level ordering must be deterministic and reused everywhere

### B. Valid-pair enumeration

Purpose:

- construct the identification set `H_gt = {(g', t_pre)}` for each target `(g, t)`

Key Python functions/classes:

- `enumerate_valid_triples()`

Planned equivalent in R:

- `enumerate_valid_pairs_edid(target_g, treatment_groups, time_periods, period_1, pt_assumption, anticipation = 0, never_treated_val = Inf)`

Notes on translation challenges:

- same-group comparisons `g' = g` are valid under `PT-All`
- all never-treated pairs are retained under `PT-All`, even when redundant
- `period_1` is excluded as a valid `t_pre`

### C. No-covariate estimator logic

Purpose:

- compute closed-form efficient DiD when no covariates are supplied

Key Python functions/classes:

- `compute_omega_star_nocov()`
- `compute_efficient_weights()`
- `compute_generated_outcomes_nocov()`
- `compute_eif_nocov()`

Planned equivalent in R:

- `compute_omega_star_nocov_edid()`
- `compute_efficient_weights_edid()`
- `compute_generated_outcomes_nocov_edid()`
- `compute_eif_nocov_edid()`

Notes on translation challenges:

- weighted covariance formulas must match Python exactly when survey weights are used
- pseudoinverse fallback should be explicit and deterministic
- `MASS::ginv()` is acceptable, but condition-number checks must be standardized

### D. Covariate / doubly robust estimator logic

Purpose:

- estimate DR efficient DiD with nuisance models and conditional efficient weights

Key Python functions/classes:

- `estimate_outcome_regression()`
- `_polynomial_sieve_basis()`
- `estimate_propensity_ratio_sieve()`
- `estimate_inverse_propensity_sieve()`
- `compute_generated_outcomes_cov()`
- `compute_omega_star_conditional()`
- `compute_per_unit_weights()`
- `compute_eif_cov()`

Planned equivalent in R:

- `estimate_outcome_regression_edid()`
- `build_polynomial_sieve_basis_edid()`
- `estimate_propensity_ratio_sieve_edid()`
- `estimate_inverse_propensity_sieve_edid()`
- `compute_generated_outcomes_cov_edid()`
- `compute_omega_star_conditional_edid()`
- `compute_per_unit_weights_edid()`
- `compute_eif_cov_edid()`

Notes on translation challenges:

- Python caches nuisance objects by tuple keys; in R use named lists keyed by strings like `"group|a|b"`
- kernel-smoothed conditional covariance is potentially expensive; block-wise computation may be needed
- preserve DR structure even if the underlying nuisance estimators are later swapped out

### E. Model fitting orchestration

Purpose:

- iterate over cohorts and periods
- dispatch to no-covariate or covariate path
- compute per-cell effects and inference
- store EIFs for later aggregation/bootstrap

Key Python functions/classes:

- `EfficientDiD.fit()`

Planned equivalent in R:

- `edid()`
- `fit_edid_cells()`

Notes on translation challenges:

- R should avoid deeply nested mutable state; return explicit lists from helpers
- the cell loop is still appropriate; vectorize inside each helper, not necessarily across all cells

### F. Inference / standard errors

Purpose:

- compute analytical, clustered, and survey-design-based SEs
- gate inference when SEs are invalid

Key Python functions/classes:

- `_compute_se_from_eif()`
- `_cluster_aggregate()`
- `safe_inference()`
- survey variance helpers in `survey.py`

Planned equivalent in R:

- `compute_eif_se_edid()`
- `cluster_aggregate_edid()`
- `safe_inference_edid()`
- `compute_survey_se_edid()`

Notes on translation challenges:

- R's `survey` package should replace most custom survey variance code
- cluster and survey handling should be mutually exclusive unless explicit design rules are defined

### G. Aggregation

Purpose:

- aggregate `ATT(g, t)` into overall, event-study, and group-level summaries
- apply WIF correction where aggregation weights are estimated

Key Python functions/classes:

- `_aggregate_overall()`
- `_aggregate_event_study()`
- `_aggregate_by_group()`
- `_compute_wif_contribution()`

Planned equivalent in R:

- `aggregate_overall_edid()`
- `aggregate_event_study_edid()`
- `aggregate_group_edid()`
- `compute_wif_contribution_edid()`

Notes on translation challenges:

- WIF correction is essential for matching inference parity
- overall and event-study aggregation use cohort-share weights; group aggregation uses equal weights over time

### H. Bootstrap

Purpose:

- generate multiplier perturbations from stored EIFs
- reaggregate perturbed cell effects
- produce bootstrap SEs, CIs, and p-values

Key Python functions/classes:

- `EfficientDiDBootstrapMixin`
- `_run_multiplier_bootstrap()`
- `generate_bootstrap_weights_batch()`
- `compute_effect_bootstrap_stats()`

Planned equivalent in R:

- `run_multiplier_bootstrap_edid()`
- `generate_multiplier_weights_edid()`
- `compute_bootstrap_stats_edid()`

Notes on translation challenges:

- cluster-level multipliers should be generated once and expanded to units
- survey bootstrap should be deferred unless the package explicitly targets design-based parity in the first release

### I. Result containers and interfaces

Purpose:

- provide stable access to fitted values and diagnostics
- support summary and tidy output

Key Python functions/classes:

- `EfficientDiDResults`
- `HausmanPretestResult`
- `EDiDBootstrapResults`

Planned equivalent in R:

- `edid_fit` S3 object
- `edid_hausman` S3 object
- nested `bootstrap` list with class `edid_bootstrap`

Notes on translation challenges:

- keep the returned object flat enough to inspect easily
- expose diagnostics like efficient weights, condition numbers, and stored EIFs only when requested

## 4. Data Structures & Interfaces

### Input data format in Python

The Python implementation expects a balanced long panel with columns for:

- outcome
- unit id
- time
- first treatment period
- optional time-invariant covariates
- optional cluster id
- optional survey-design columns

Internal numeric structures:

- long `DataFrame` for validation
- wide numeric outcome matrix `(n_units x n_periods)`
- logical cohort masks
- unit-level covariate matrix
- named dictionaries keyed by tuples

### Recommended R equivalents

External API:

- accept `data.frame`, `data.table`, or tibble

Internal representation:

- coerce to `data.table` for validation and grouping
- use base numeric matrices for:
  - wide outcomes
  - covariates
  - kernel weight matrices
  - covariance matrices
- use logical vectors for cohort masks
- use named lists for caches

Recommendation:

- use `data.table` internally for reshaping and grouping
- do not use tibble internally for performance-critical steps

### Recommended R function signatures

Primary fitting function:

```r
edid <- function(
  data,
  outcome,
  unit,
  time,
  first_treat,
  covariates = NULL,
  pt_assumption = c("all", "post"),
  alpha = 0.05,
  cluster = NULL,
  control_group = c("never_treated", "last_cohort"),
  n_bootstrap = 0L,
  bootstrap_weights = c("rademacher", "mammen", "webb"),
  seed = NULL,
  anticipation = 0L,
  sieve_k_max = NULL,
  sieve_criterion = c("bic", "aic"),
  ratio_clip = 20,
  kernel_bandwidth = NULL,
  aggregate = c("none", "simple", "event_study", "group", "all"),
  balance_e = NULL,
  survey_design = NULL,
  store_eif = FALSE
) { }
```

Core helpers:

```r
prepare_edid_panel <- function(data, outcome, unit, time, first_treat, covariates = NULL) { }
enumerate_valid_pairs_edid <- function(target_g, treatment_groups, time_periods, period_1, pt_assumption, anticipation = 0L, never_treated_val = Inf) { }
compute_omega_star_nocov_edid <- function(...) { }
compute_generated_outcomes_nocov_edid <- function(...) { }
estimate_propensity_ratio_sieve_edid <- function(covariate_matrix, mask_g, mask_gp, k_max = NULL, criterion = "bic", ratio_clip = 20, unit_weights = NULL) { }
compute_omega_star_conditional_edid <- function(...) { }
run_multiplier_bootstrap_edid <- function(group_time_effects, eif_by_gt, n_units, aggregate, balance_e, treatment_groups, cohort_fractions, cluster_indices = NULL, survey_design = NULL, unit_level_weights = NULL) { }
hausman_pretest_edid <- function(data, outcome, unit, time, first_treat, covariates = NULL, cluster = NULL, anticipation = 0L, control_group = "never_treated", alpha = 0.05, ...) { }
```

S3 methods:

- `print.edid_fit`
- `summary.edid_fit`
- `coef.edid_fit`
- `vcov.edid_fit`
- `as.data.frame.edid_fit`
- `plot.edid_fit`

## 5. Step-by-Step Implementation Plan

### Step 1: Set up package structure

Concrete tasks:

- create package skeleton with `DESCRIPTION`, `NAMESPACE`, `R/`, `tests/testthat/`
- choose imports: `data.table`, `MASS`, `survey`, `stats`
- add `testthat` edition and CI config
- define exported API surface: `edid()`, `hausman_pretest_edid()`

Expected outputs:

- installable R package skeleton
- passing package load tests

Dependencies on previous steps:

- none

### Step 2: Implement core data handling

Concrete tasks:

- implement input validation
- implement balanced-panel check
- implement absorbing-treatment check
- implement duplicate `(unit, time)` check
- implement `last_cohort` control handling
- reshape long panel to wide outcome matrix
- extract unit-level covariates and weights
- build cohort masks, period maps, and cluster indices

Expected outputs:

- `prepare_edid_panel()` returning a deterministic internal structure
- tests covering valid and invalid input cases

Dependencies on previous steps:

- Step 1

### Step 3: Implement no-covariate estimator logic

Concrete tasks:

- implement valid-pair enumeration
- implement weighted/unweighted sample covariance helper
- implement `Omega*` for no-covariate path
- implement efficient weights with condition-number and pseudoinverse fallback
- implement generated outcomes
- implement no-covariate EIF
- build cell-level fit loop using only no-covariate path

Expected outputs:

- correct cell-level `ATT(g, t)` estimates without covariates
- stored cell-level EIFs
- tests for PT-All, PT-Post, and single-valid-pair cases

Dependencies on previous steps:

- Step 2

### Step 4: Implement inference

Concrete tasks:

- implement scalar inference helper
- implement standard analytical EIF-based SEs
- implement cluster-robust EIF SEs
- implement aggregation helpers
- implement WIF correction
- return overall, event-study, and group summaries

Expected outputs:

- usable `edid_fit` object with summary statistics
- tests for analytical and clustered inference

Dependencies on previous steps:

- Step 3

### Step 5: Implement multiplier bootstrap

Concrete tasks:

- implement multiplier weight generator
- implement per-cell perturbation from stored EIFs
- implement bootstrap reaggregation
- implement bootstrap SE, percentile CI, and bootstrap p-value logic
- merge bootstrap inference back into fitted result

Expected outputs:

- bootstrap-enhanced `edid_fit`
- tests for unit-level and cluster-level bootstrap

Dependencies on previous steps:

- Step 4

### Step 6: Implement covariate / doubly robust path

Concrete tasks:

- implement outcome regression helper
- implement polynomial sieve basis construction
- implement sieve ratio estimation and degree selection
- implement inverse propensity estimation
- implement DR generated outcomes
- implement conditional `Omega*(X)`
- implement per-unit efficient weights
- implement DR EIF
- integrate covariate path into main fit loop

Expected outputs:

- cell-level and aggregate results with covariates
- diagnostics for clipping, singularity, and bandwidth
- parity tests against Python fixtures

Dependencies on previous steps:

- Step 4
- Step 5 is not strictly required, but helps verify inference parity once DR fit works

### Step 7: Implement survey support

Concrete tasks:

- design an R survey interface, ideally accepting a `survey::svydesign` or metadata used to build one
- route point estimation through unit-level survey weights
- route variance estimation through `survey` package design-based variance machinery
- explicitly reject unsupported combinations until implemented

Expected outputs:

- survey-weighted point estimates
- design-based SEs
- clear warnings/errors for unsupported bootstrap-survey combinations

Dependencies on previous steps:

- Step 4
- Step 6 if covariate+survey parity is required

### Step 8: Implement Hausman pretest

Concrete tasks:

- fit under `PT-All` and `PT-Post`
- aggregate shared post-treatment horizons to event-study vectors
- build covariance difference from EIF matrices
- compute pseudoinverse-based test statistic and p-value

Expected outputs:

- `edid_hausman` result object
- tests for common-horizon and degenerate-covariance cases

Dependencies on previous steps:

- Step 4

### Step 9: Documentation and examples

Concrete tasks:

- write `roxygen2` docs for all exported functions
- add a vignette for no-covariate usage
- add a vignette for covariate usage
- add an advanced vignette for cluster/bootstrap and survey options

Expected outputs:

- documented exported API
- reproducible examples

Dependencies on previous steps:

- all prior implementation steps

## 6. Key Algorithm Translation

### Step-by-step description

For each treatment cohort `g` and target period `t`:

1. Determine the baseline convention.
   - `PT-All`: universal earliest period
   - `PT-Post`: cohort-specific baseline `g - 1 - anticipation`
2. Enumerate valid comparison pairs `(g', t_pre)`.
3. If no covariates:
   - compute scalar generated outcomes for each pair
   - compute unconditional covariance matrix `Omega*`
   - compute efficient weights
   - combine moments to get `ATT(g, t)`
   - compute the no-covariate EIF
4. If covariates:
   - estimate outcome regression nuisances
   - estimate propensity ratio nuisances
   - estimate inverse propensity nuisances
   - compute per-unit DR generated outcomes
   - compute per-unit conditional covariance matrices
   - compute per-unit efficient weights
   - average the per-unit efficient score to get `ATT(g, t)`
   - compute the DR EIF
5. Compute analytical SE from the EIF.
6. Store cell effect, SE, inference fields, and EIF.
7. After all cells are fit, aggregate them and apply WIF correction to aggregate EIFs.
8. If bootstrap is requested, perturb cell effects using multiplier weights and recompute aggregation statistics.

### Pseudocode

```text
prepare panel objects

for g in treatment_groups:
    base_ref = if PT-Post then g - 1 - anticipation else first_period

    for t in time_periods excluding universal first period:
        pairs = enumerate_valid_pairs(g, time_periods, treatment_groups, pt_assumption)

        if pairs is empty:
            store NA cell
            continue

        if covariates is NULL:
            y_hat = generated_outcomes_nocov(g, t, pairs, panel_objects)
            omega = omega_star_nocov(g, t, pairs, panel_objects)
            weights = efficient_weights(omega)
            att_gt = sum(weights * y_hat)
            eif_gt = eif_nocov(g, t, pairs, weights, panel_objects, att_gt)
        else:
            fill nuisance caches needed for this cell
            gen_out = generated_outcomes_cov(g, t, pairs, nuisance_cache, panel_objects)
            omega_x = omega_star_conditional(g, t, pairs, nuisance_cache, panel_objects)
            weights_x = per_unit_weights(omega_x)
            score_i = row_sum(weights_x * gen_out)
            att_gt = weighted_or_unweighted_mean(score_i)
            eif_gt = score_i - att_gt

        se_gt = se_from_eif(eif_gt, cluster_or_survey_info)
        store cell result
        store eif_gt

overall = aggregate_overall(cell_results, eifs, cohort_fractions, wif = TRUE)
event_study = aggregate_event_study(cell_results, eifs, cohort_fractions, wif = TRUE)
group = aggregate_group(cell_results, eifs, equal_time_weights)

if n_bootstrap > 0:
    bootstrap_results = run_multiplier_bootstrap(stored_eifs, aggregation_info)
    replace analytical inference with bootstrap inference

return fitted object
```

### How it should be implemented in R

- keep the outer `(g, t)` loop explicit
- keep heavy operations matrix-based inside helper functions
- cache nuisance estimates in named lists
- represent panel objects as a single internal list passed across helpers
- keep inference as post-estimation logic based on stored EIF vectors

Recommended internal object:

```r
panel_obj <- list(
  outcome_wide = ...,
  covariate_matrix = ...,
  all_units = ...,
  time_periods = ...,
  period_to_col = ...,
  period_1 = ...,
  unit_cohorts = ...,
  cohort_masks = ...,
  never_treated_mask = ...,
  cohort_fractions = ...,
  cluster_indices = ...,
  survey_design = ...,
  unit_level_weights = ...
)
```

## 7. Numerical & Statistical Considerations

### Precision issues

- matrix inversion and pseudoinversion are the main numerical risk points
- the condition-number threshold should be centralized in one helper
- ratio and inverse-propensity clipping thresholds must be identical across code paths

### Differences between Python and R numerical behavior

- Python `np.linalg.solve` and R `solve()` can behave slightly differently near singularity
- `MASS::ginv()` may not match NumPy pseudoinverse exactly; define acceptable tolerances in tests
- weighted means and covariances in base R may not match Python's explicit formulas unless implemented manually

### Handling of missing data

- the Python implementation rejects non-finite outcomes and covariates early
- the R port should do the same rather than silently dropping rows
- missing cluster or survey identifiers should error, not recycle or coerce

### Performance considerations

- use vectorized matrix operations inside helpers
- keep loops at the `(g, t)` level rather than trying to vectorize the whole estimator
- kernel covariance estimation is `O(n^2 * H^2)` in the current design and is the main performance bottleneck
- if performance becomes limiting, optimize `compute_omega_star_conditional_edid()` first
- delay Rcpp optimization until parity is established in pure R

## 8. Dependencies & Libraries

### Python libraries used

- `numpy` for arrays and linear algebra
- `pandas` for panel wrangling
- `scipy` for distances and inference distributions
- custom `diff_diff.linalg`
- custom survey utilities

### Recommended R equivalents

- `data.table` for panel preprocessing and reshaping
- `stats` for linear algebra basics, `lm.wfit`, `pnorm`, `qt`, `qnorm`
- `MASS` for pseudoinverse via `ginv`
- `survey` for survey-design objects and variance estimation
- base matrices for kernel and covariance work

### Justification

- `data.table` is the most natural equivalent to `pandas` for deterministic, fast grouped operations
- `survey` is far more idiomatic and robust in R than reimplementing design-based variance logic from scratch
- S3 + base matrices keeps dependencies light while still matching the Python design

Avoid depending on `fixest` for the core estimator logic. It is excellent for regression-based DiD, but this estimator is not a fixed-effects regression workflow.

## 9. Testing Strategy

### How to ensure parity between Python and R implementations

Use a three-layer testing strategy:

1. unit tests for helper math
2. integration tests for complete `edid()` fits
3. parity tests against Python-generated fixtures

### Suggested unit tests

- valid-pair enumeration under `PT-All` and `PT-Post`
- weighted/unweighted covariance helper
- efficient weights sum to 1
- single valid pair returns weight 1
- singular `Omega*` triggers pseudoinverse or uniform fallback
- DR generated outcomes shape and finiteness
- cluster aggregation formula
- bootstrap weight generators

### Golden test cases

Create Python fixture files for:

- simple one-cohort balanced panel, no covariates
- staggered panel, no covariates
- staggered panel with covariates
- clustered example
- bootstrap example with fixed seed
- survey-weighted example if survey support is in scope

Fixture contents should include:

- group-time effects
- overall ATT
- event-study effects
- group effects
- efficient weights for selected cells
- condition numbers for selected cells
- stored EIF vectors for selected cells where feasible

Define tolerance bands:

- exact or near-exact for no-covariate point estimates
- small numerical tolerance for SEs
- slightly wider tolerance for DR and kernel-smoothed paths

## 10. Risks & Pitfalls

### Subtle differences between Python and R semantics

- `Inf`, `NA`, and `NaN` propagation differ in edge cases
- factor handling in R can silently alter ordering if not explicitly coerced
- `data.table` grouping order must be controlled to preserve deterministic unit ordering

### Parts of the code that are hard to translate

- conditional `Omega*(X)` with kernel smoothing
- survey-design variance parity if exact Python behavior is required
- pseudoinverse-based Hausman covariance difference logic
- cache management for nuisance estimates across many target cells

### Hidden assumptions in the Python implementation

- every unit must be observed in every period
- covariates are time-invariant
- pre-treatment placebo cells are estimated, not dropped
- `period_1` is excluded as a valid `t_pre`
- overall ATT is a cohort-share-weighted average of post-treatment cells, not the paper's equal-weight event-time average
- bootstrap is multiplier-based, not nonparametric resampling

## 11. Agent Execution Notes

### Order of implementation to minimize errors

Recommended order:

1. data validation and panel preparation
2. valid-pair enumeration
3. no-covariate cell estimator
4. analytical inference
5. aggregation + WIF
6. bootstrap
7. DR covariate path
8. survey support
9. Hausman pretest

### Components that can be implemented independently

- input validation
- valid-pair enumeration
- no-covariate `Omega*` and generated outcomes
- inference helpers
- bootstrap helper functions
- S3 print/summary/tidy methods

### Suggested checkpoints for validation

- Checkpoint 1: `prepare_edid_panel()` passes all validation tests
- Checkpoint 2: no-covariate `ATT(g, t)` matches Python on simple fixtures
- Checkpoint 3: overall/event-study/group aggregation matches Python
- Checkpoint 4: clustered analytical SEs and bootstrap results are numerically close
- Checkpoint 5: DR covariate path matches Python fixtures within tolerance
- Checkpoint 6: survey and Hausman features pass dedicated regression tests

### Implementation guidance for the agent

- centralize all internal constants: condition-number threshold, clipping bounds, minimum bootstrap count warnings
- keep a single source of truth for panel ordering
- avoid silent coercions or row drops
- return rich diagnostics, but gate large objects like EIF matrices behind `store_eif`
- do not optimize prematurely; establish mathematical parity first, then profile
