# Efficient DiD Implementation Guide

This document explains how `diff-diff` implements the Efficient Difference-in-Differences (EDiD) estimator from Chen, Sant'Anna, and Xie (2025), with the goal of making the method implementable in another package, especially in R.

The emphasis here is on methodology and implementation logic, not on Python-specific API details. Where the repository makes a concrete engineering choice, this guide distinguishes between:

- what is fundamental to the estimator, and
- what is a practical implementation choice that may be replaced in another language.

## 1. What the estimator is trying to do

The object of interest is a group-time average treatment effect on the treated:

`ATT(g, t) = E[Y_t(g) - Y_t(infinity) | G = g]`

where:

- `G = g` means a unit first receives treatment in period `g`
- `G = infinity` means never treated
- `t` is the calendar period

In staggered adoption settings, for a fixed target cell `(g, t)`, there are often many valid Difference-in-Differences style moments that identify the same ATT:

- different baseline periods `t_pre`
- different comparison cohorts `g'`

Under the stronger `PT-All` assumption, the model is overidentified. EDiD uses that overidentification to build a lower-variance estimator by optimally combining all valid moments instead of picking just one.

That is the central first principle:

1. Define all valid unbiased moments for a target `(g, t)`.
2. Estimate their covariance structure.
3. Combine them with minimum-variance weights subject to the weights summing to 1.

Without covariates, the combination weights are global for each `(g, t)`.
With covariates, the efficient weights become unit-specific functions `w(X_i)`.

## 2. First principles you should preserve in any reimplementation

If you implement EDiD in another language, these are the core ideas to preserve:

### 2.1 Balanced panel, absorbing treatment, cohort structure

The implementation assumes:

- a short balanced panel
- one row per unit-period
- absorbing treatment
- a well-defined first treatment period per unit

If these fail, the estimator should stop early rather than try to proceed.

### 2.2 Estimate each `ATT(g, t)` independently

The implementation loops over treatment cohorts `g` and time periods `t`, and treats each target cell as its own estimation problem. This is important because:

- the valid comparison set changes by target cell
- the efficient covariance matrix changes by target cell
- nuisance functions in the covariate path depend on the target cell

This also makes the method naturally parallelizable.

The repository estimates both:

- post-treatment cells, which are the main causal targets, and
- pre-treatment cells, which act as placebo or pre-trend diagnostics

It skips only the earliest observed period, which serves as the universal reference period in the PT-All path.

### 2.3 EDiD is a weighted combination of valid DiD moments

For each `(g, t)`, define an index set `H_gt` of valid pairs `(g', t_pre)`. Each pair generates one moment. The estimator is:

`ATT_hat(g, t) = sum_{j in H_gt} w_j * Y_hat_j`

in the no-covariate case, or

`ATT_hat(g, t) = E_n[ sum_{j in H_gt} w_j(X_i) * Y_hat_{ij} ]`

with covariates.

The efficient weights are always of the form:

`w = 1' Omega^{-1} / (1' Omega^{-1} 1)`

or, conditionally on covariates,

`w(X_i) = 1' Omega(X_i)^{-1} / (1' Omega(X_i)^{-1} 1)`

### 2.4 PT-All and PT-Post are different regimes

The repository supports two identification regimes:

- `PT-All`: stronger assumption, overidentified, uses all valid comparison pairs
- `PT-Post`: weaker assumption, just-identified, reduces to a single comparison moment

Under `PT-Post`, EDiD collapses to the standard single-baseline estimator (the implementation is designed to match Callaway-Sant'Anna for post-treatment effects).

### 2.5 Inference is influence-function based

The estimator is not just point estimation plus a generic bootstrap. The core analytical inference is based on the efficient influence function (EIF), computed for every target `(g, t)`, and then reused for:

- cell-level standard errors
- overall aggregation
- event-study aggregation
- group aggregation
- multiplier bootstrap

That architecture is worth preserving in an R implementation.

## 3. Files and conceptual roles in `diff-diff`

The repository splits the logic into five conceptual layers:

- `diff_diff/efficient_did.py`
  Main orchestration: validation, panel reshaping, target-cell loop, aggregation, analytical SEs, bootstrap, Hausman pretest.
- `diff_diff/efficient_did_weights.py`
  No-covariate path: valid pair enumeration, unconditional `Omega*`, efficient weights, generated outcomes, EIF.
- `diff_diff/efficient_did_covariates.py`
  Covariate path: outcome regression, sieve propensity ratios, inverse propensities, conditional `Omega*(X)`, per-unit weights, EIF.
- `diff_diff/efficient_did_bootstrap.py`
  Multiplier bootstrap based on stored EIF values.
- `diff_diff/efficient_did_results.py`
  Result container and reporting.

For an R implementation, the same conceptual split is useful even if the file structure differs.

## 4. Data representation used by the implementation

The implementation converts the long panel into unit-level structures:

- `outcome_wide`: matrix of shape `(n_units, n_periods)`
- `unit_cohorts`: vector of first treatment dates by unit
- `cohort_masks[g]`: logical vector for membership in cohort `g`
- `never_treated_mask`: logical vector for never-treated units
- `period_to_col`: map from calendar time to column index

If covariates are used, it also constructs:

- `covariate_matrix`: matrix of time-invariant unit covariates, one row per unit

This is the right representation for efficient numerical implementation in R as well.

## 5. Assumptions and validation rules

Before estimation, `diff-diff` checks:

1. Required columns exist.
2. Time and first-treatment variables are numeric.
3. The panel is balanced.
4. Outcomes are finite.
5. There are no duplicate `(unit, time)` rows.
6. Treatment is absorbing, meaning `first_treat` does not vary within unit.
7. A never-treated group exists, unless the implementation explicitly switches to a last-cohort-as-control design.
8. If covariates are supplied, they must:
   - exist,
   - be finite,
   - be constant within unit.

These checks are not incidental. They encode the data regime the estimator needs.

## 6. Control group handling

The repository supports two control group definitions:

- `never_treated`
- `last_cohort`

In the `last_cohort` option, the last treated cohort is reclassified as pseudo-never-treated, and time periods at or after that cohort's treatment start are dropped. This is a design choice for settings with no true never-treated units.

This is not the same as a "not-yet-treated" control definition. EDiD here is built around a never-treated comparison group, possibly created from the last cohort by trimming.

## 7. Step 1: enumerate valid comparison pairs

For each target `(g, t)`, the first key object is the set of valid pairs:

`H_gt = {(g', t_pre)}`

### 7.1 PT-Post

Under `PT-Post`, there is only one valid pair:

- comparison group: never treated
- baseline period: `g - 1 - anticipation`

So the model is just-identified.

### 7.2 PT-All

Under `PT-All`, the implementation includes all pairs such that:

- `g'` is either:
  - the never-treated group, or
  - any treatment cohort, including `g` itself
- `t_pre` is pre-treatment for the comparison cohort `g'`
- `t_pre` is not the universal baseline period `period_1`

Two implementation details matter:

1. Same-cohort comparisons `g' = g` are allowed.
   These are valid overidentifying moments.
2. Never-treated pairs `(infinity, t_pre)` are all retained under `PT-All`.
   Some are redundant, but the covariance-based weighting handles that redundancy.

In code this is done by `enumerate_valid_triples()`.

## 8. Step 2A: no-covariate path

When no covariates are supplied, the estimator uses the closed-form path.

### 8.1 Generated outcomes

For each valid pair `j = (g', t_pre)`, compute a scalar generated outcome:

`Y_hat_j = mean(Y_t - Y_1 | G=g)
         - mean(Y_t - Y_tpre | G=infinity)
         - mean(Y_tpre - Y_1 | G=g')`

Under `PT-Post`, `Y_1` is replaced by the cohort-specific baseline `Y_{g-1-anticipation}`.

In the repository this is `compute_generated_outcomes_nocov()`.

### 8.2 Unconditional covariance matrix `Omega*`

Then build the `|H_gt| x |H_gt|` covariance matrix `Omega*`. Its entries are the sample analog of the paper's covariance expression. The implementation constructs each matrix entry from up to five within-group covariance terms:

1. treated-group variance term
2. never-treated covariance term
3. treated/comparison cross term for the `j`th moment
4. treated/comparison cross term for the `k`th moment
5. comparison-group covariance term when `g'_j = g'_k`

The function is `compute_omega_star_nocov()`.

The important conceptual point is:

- `Omega*` is the covariance matrix of the available identifying moments
- efficient weighting is entirely driven by this covariance structure

### 8.3 Efficient weights

Once `Omega*` is available, compute:

`w = 1' Omega^{-1} / (1' Omega^{-1} 1)`

This is implemented in `compute_efficient_weights()`.

Practical safeguards in the repository:

- if `Omega*` is all zeros, use uniform weights
- if the condition number is too large, use the Moore-Penrose pseudoinverse
- if the weight denominator is nearly zero, use uniform weights

These are numerical fallback rules, not changes to the estimand.

### 8.4 Point estimate

The cell estimate is:

`ATT_hat(g, t) = w' Y_hat`

### 8.5 Efficient influence function

The no-covariate EIF is computed directly from:

- treated-group demeaned outcome change
- never-treated demeaned comparison change
- comparison-cohort demeaned baseline change

combined with the efficient weights.

In the repository this is `compute_eif_nocov()`.

That EIF is then used for analytical standard errors and all later aggregation.

## 9. Step 2B: covariate path

When covariates are supplied, `diff-diff` moves to a doubly robust implementation.

The theory requires nuisance objects for:

- outcome regressions
- generalized propensity-score ratios
- inverse generalized propensities
- conditional covariance matrices

### 9.1 What is essential versus what is a repository choice

Essential to the method:

- estimate `m_{g', a, b}(X) = E[Y_a - Y_b | G=g', X]`
- estimate ratios `r_{g, g'}(X) = p_g(X) / p_{g'}(X)`
- estimate inverse propensities `s_{g'}(X) = 1 / p_{g'}(X)`
- estimate conditional `Omega*(X)`
- form unit-specific efficient weights

Repository-specific implementation choices:

- linear OLS/WLS for outcome regression
- polynomial sieve basis for ratio and inverse-propensity estimation
- AIC/BIC to choose sieve degree
- Gaussian kernel with Silverman's rule for conditional covariance smoothing

An R implementation can change these choices if it preserves the identifying structure.

### 9.2 Outcome regression

For each needed group and time pair, the implementation estimates:

`m_hat_{g', a, b}(X) = E[Y_a - Y_b | G=g', X]`

using linear regression of `(Y_a - Y_b)` on:

- intercept
- covariates

within the selected group.

This is `estimate_outcome_regression()`.

The predictions are generated for all units, not just the training group.

### 9.3 Propensity-score ratio estimation

The ratio `r_{g, g'}(X)` is estimated by sieve minimum distance, using the paper's convex criterion. For each degree `K`:

`beta_K = arg min E_n[ G_g' * (psi_K(X)' beta)^2 - 2 * G_g * psi_K(X)' beta ]`

The first-order condition yields a linear system:

`(Psi_gprime' Psi_gprime) beta = sum_{i: G_i = g} psi_K(X_i)`

with weighted analogs when survey weights are present.

The repository then:

- tries `K = 1, ..., K_max`
- skips singular systems
- chooses the best `K` by AIC or BIC
- clips the fitted ratios to `[1 / ratio_clip, ratio_clip]`

This is `estimate_propensity_ratio_sieve()`.

### 9.4 Inverse propensity estimation

The implementation also estimates:

`s_{g'}(X) = 1 / p_{g'}(X)`

using an analogous sieve problem:

`beta_K = arg min E_n[ G_g' * (psi_K(X)' beta)^2 - 2 * psi_K(X)' beta ]`

Again, model selection is by AIC/BIC and fitted values are clipped to a plausible range.

This is `estimate_inverse_propensity_sieve()`.

### 9.5 Doubly robust generated outcomes

For each valid pair `j = (g', t_pre)` and each unit `i`, the implementation computes:

`gen_out[i, j]`

using the three-term expression from the paper:

1. treated contribution
2. never-treated contribution weighted by `r_hat_{g, infinity}(X_i)`
3. comparison-cohort contribution weighted by `r_hat_{g, g'}(X_i)`

This is the sample analog of the DR moment in Equation 4.4 and is implemented in `compute_generated_outcomes_cov()`.

The resulting object is a matrix of shape `(n_units, |H_gt|)`.

### 9.6 Conditional covariance `Omega*(X)`

The estimator then constructs a conditional covariance matrix for every unit:

`Omega_i = Omega*(X_i)`

This is a three-dimensional array of shape `(n_units, |H_gt|, |H_gt|)`.

In the repository:

- local covariances are estimated by Nadaraya-Watson kernel smoothing
- the kernel is Gaussian
- the bandwidth is Silverman's rule if not supplied
- covariance terms are scaled by the estimated inverse propensities `s_hat_g(X_i)`

This is implemented in `compute_omega_star_conditional()`.

Conceptually, this is the conditional analog of the unconditional `Omega*` from the no-covariate case.

### 9.7 Per-unit efficient weights

For each unit `i`, compute:

`w(X_i) = 1' Omega_i^{-1} / (1' Omega_i^{-1} 1)`

This is `compute_per_unit_weights()`.

As in the no-covariate path, the implementation uses:

- pseudoinverse fallback for ill-conditioned matrices
- uniform weights if the matrix is degenerate or the denominator is too small

### 9.8 Point estimate

The cell estimate is the sample mean of the per-unit efficient score:

`ATT_hat(g, t) = E_n[ sum_j w_j(X_i) * gen_out[i, j] ]`

In the repository this is implemented as:

1. compute `per_unit_scores[i] = sum_j w_j(X_i) * gen_out[i, j]`
2. average `per_unit_scores` over units

### 9.9 Efficient influence function

The covariate-path EIF is computed as:

`EIF_i = per_unit_score_i - ATT_hat(g, t)`

This is `compute_eif_cov()`.

The implementation uses the plug-in EIF and relies on the Neyman orthogonality logic discussed in the paper. In other words, once the nuisance estimates are plugged in, the score is centered at the final scalar ATT estimate.

## 10. Standard errors

### 10.1 Analytical SEs

For each cell `(g, t)`, the default analytical standard error is based on the EIF:

`SE_hat = sqrt( mean(EIF_i^2) / n )`

equivalently

`SE_hat = sqrt( sum(EIF_i^2) / n^2 )`

### 10.2 Cluster-robust analytical SEs

If clustering is requested, the repository:

1. aggregates EIF values within cluster
2. centers the cluster sums
3. computes a Liang-Zeger style sandwich variance
4. applies the small-sample factor `G / (G - 1)`

This is done through `_cluster_aggregate()` and `_compute_se_from_eif()`.

### 10.3 Survey-weighted SEs

The repository also supports survey design based inference, but this is not essential for a first R port unless survey support is in scope.

If survey support is included, the design principle is:

- point estimation uses unit-level survey weights
- EIF scores are passed into design-based variance routines

## 11. Aggregation after cell estimation

After all `ATT(g, t)` cells are computed, the implementation aggregates them in three ways.

### 11.1 Overall ATT

The overall ATT is a cohort-size weighted average of post-treatment cells:

`overall = sum_k q_k * ATT(g_k, t_k)`

where `q_k` is proportional to the cohort share `pi_g`.

Important note:

- this matches the repository's Callaway-Sant'Anna style "simple" aggregation
- it is not the same as the paper's equal-weight-over-event-times `ES_avg`

### 11.2 Event-study aggregation

For each relative time `e = t - g`, the implementation averages all finite `ATT(g, g+e)` values using cohort-size weights.

### 11.3 Group aggregation

For each cohort `g`, the implementation averages all post-treatment `ATT(g, t)` values equally over `t`.

## 12. Weight uncertainty in aggregation

A very important implementation detail is that aggregated standard errors are not obtained by simply averaging the cell EIFs.

Because the aggregation weights depend on estimated cohort shares, the repository adds a weight influence function (WIF) correction. This is handled in:

- `_compute_wif_contribution()`
- `_aggregate_overall()`
- `_aggregate_event_study()`

If you want inference for aggregated parameters to match the repository, you should preserve this correction.

This matters most for:

- overall ATT
- event-study effects

Group aggregation uses equal weights over time and does not need the same cohort-share WIF adjustment.

## 13. Multiplier bootstrap

The repository's bootstrap is not a nonparametric resampling bootstrap. It is a multiplier bootstrap built from the stored EIF vectors.

For each bootstrap draw:

`ATT_b(g, t) = ATT_hat(g, t) + (1/n) * sum_i xi_i * EIF_i(g, t)`

where the multipliers may be:

- Rademacher
- Mammen
- Webb

If clustering is used, the same multiplier is assigned to all units in a cluster.

The bootstrap then:

1. perturbs all cell estimates
2. recomputes overall, event-study, and group aggregates from the perturbed cells
3. computes bootstrap SEs, p-values, and confidence intervals

This is implemented in `efficient_did_bootstrap.py`.

For an R implementation, this is a practical and efficient inference approach because it avoids recomputing nuisance functions in every bootstrap repetition.

## 14. Numerical safeguards and edge cases

An implementation should handle these cases deliberately:

### 14.1 No valid comparison pairs

If `H_gt` is empty, return `NA` for that cell and mark its inference as undefined.

### 14.2 Single valid pair

If `|H_gt| = 1`, the efficient weight is trivially 1 and the estimator reduces to the unique valid DiD moment.

### 14.3 Singular or nearly singular covariance matrices

If `Omega*` or `Omega*(X_i)` is singular or poorly conditioned:

- use a pseudoinverse
- if needed, fall back to uniform weights

### 14.4 Extreme ratio estimates

If estimated propensity ratios are too large or too small, clip them and warn the user. This is an overlap diagnostic.

### 14.5 Failed sieve estimation

The repository falls back to simple defaults:

- ratio estimation failure -> constant ratio 1
- inverse-propensity estimation failure -> unconditional scaling

These fallbacks keep the estimator running and preserve DR consistency if the outcome model is correct.

### 14.6 Small or zero-weight cohorts

If a cohort is tiny, warn the user. If a cohort has zero effective weight, skip or fail cleanly depending on the design.

## 15. What is fundamental and what can change in an R implementation

### 15.1 Fundamental pieces to keep

- balanced-panel staggered adoption setup
- valid comparison pair construction
- efficient weighting by inverse covariance
- EIF-based analytical inference
- WIF correction for aggregated parameters
- DR generated outcomes in the covariate path
- conditional covariance weighting when covariates are used

### 15.2 Pieces you can reasonably replace

- OLS outcome regression with another learner
- polynomial sieve basis with splines or other sieves
- AIC/BIC model selection with cross-validation or another criterion
- Gaussian kernel with another smooth kernel
- Silverman bandwidth with cross-validated bandwidth
- exact object model and API

If you replace these pieces, the estimator remains methodologically faithful as long as the underlying nuisance objects and efficient weighting structure are preserved.

## 16. Recommended implementation order for an R package

If the goal is to get a reliable R version into production, the best order is:

1. Implement the no-covariate path first.
   This gives a stable, closed-form baseline and lets you validate pair enumeration, `Omega*`, weights, cell estimates, and EIFs.
2. Add analytical EIF-based SEs.
   This makes the estimator usable before bootstrap is added.
3. Add overall/event-study/group aggregation with WIF correction.
4. Add multiplier bootstrap from stored EIFs.
5. Add the covariate path.
   Start with simple nuisance estimators before attempting more flexible machine-learning options.

This order matches the estimator's internal dependency graph and greatly simplifies debugging.

## 17. Suggested pseudocode

Below is the implementation skeleton that most closely matches `diff-diff` while remaining language-agnostic.

```text
Input: balanced panel with unit id, time, outcome, first_treat, optional X

Validate data
Build sorted unit list and sorted time list
Pivot outcome to outcome_wide[n_units, n_periods]
Build unit_cohorts, cohort masks, never-treated mask
Build cohort fractions
If X present, build unit-level covariate matrix

For each treatment cohort g:
    Choose baseline reference:
        if PT-Post: base = g - 1 - anticipation
        else:       base = first observed period

    For each time t except the earliest observed period:
        Enumerate valid pairs H_gt
        If H_gt empty:
            store NA result and continue

        If no covariates:
            Compute generated outcomes Y_hat[j] for j in H_gt
            Compute Omega*[j, k]
            Compute efficient weights w
            ATT_hat(g, t) = sum_j w[j] * Y_hat[j]
            Compute EIF_i(g, t)

        If covariates:
            Estimate needed outcome regressions m_hat
            Estimate needed ratios r_hat
            Estimate needed inverse propensities s_hat
            Compute generated outcomes gen_out[i, j]
            Compute conditional Omega_i[j, k] for each unit i
            Compute per-unit weights w_i[j]
            score_i = sum_j w_i[j] * gen_out[i, j]
            ATT_hat(g, t) = mean(score_i)
            EIF_i(g, t) = score_i - ATT_hat(g, t)

        Compute analytical SE from EIF
        Store cell result and EIF

Aggregate post-treatment cells to:
    overall ATT
    event-study effects
    group effects
Add WIF correction to aggregated EIFs
Compute aggregate SEs and inference

If bootstrap requested:
    generate multiplier weights
    perturb each cell using stored EIFs
    re-aggregate perturbed cells
    replace inference with bootstrap inference

Return results object
```

## 18. What `diff-diff` is actually implementing, in one sentence

The repository implements EDiD as a cell-by-cell semiparametric minimum-variance combination of all valid DiD moments, with EIF-based inference and a doubly robust covariate extension that replaces global efficient weights with conditional efficient weights `w(X)`.

## 19. Practical guidance for an R port

If the goal is compatibility with this repository's behavior, try to match:

- valid pair construction exactly
- PT-Post baseline logic exactly
- cohort-size weighting in aggregations
- WIF correction for overall and event-study inference
- pseudoinverse fallback rules
- multiplier bootstrap based on stored EIFs

If the goal is methodological faithfulness rather than exact parity, the highest-priority pieces are:

- correct moment construction
- correct covariance-based efficient weighting
- correct EIF construction
- correct aggregation logic

## 20. Bottom line

The cleanest way to think about EDiD is:

- every `(g, t)` cell has several valid identifying moments
- the estimator computes all of them
- it learns their covariance structure
- it combines them in the most efficient way
- with covariates, that efficient combination becomes conditional on `X`
- inference is built from the resulting efficient influence function

That is the implementation core you should preserve in R.
