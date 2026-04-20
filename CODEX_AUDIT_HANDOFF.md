# CODEX audit handoff for Claude Code
## Goal
Audit and repair the `edid()` covariate path so that the implementation is defensible against Chen, Sant'Anna & Xie (2025), especially Eq. (3.9), Eq. (3.10), Eq. (3.12), and Eq. (4.1)-(4.4). Use the statsclawn plugin before changing code to derive the target formulas in code-ready notation and to document each mapping from paper to implementation.

Do not rewrite the package. Fix confirmed statistical and implementation bugs in the covariate path, add missing guards, and add executable tests and simulation evidence for bias, EIF correctness, variance calibration, coverage, and reproducibility.

## Context
The current code path is primarily in:
- `R/edid.R`
- `R/edid-fit.R`
- `R/edid-data.R`
- `R/edid-validate.R`
- `R/edid-cov.R`
- `R/edid-cov-eif.R`
- `R/edid-inference.R`

The current checked-in benchmark is `benchmark/compare_author_vs_edid.R`, but it is explicitly a no-covariate benchmark. The stronger covariate claims currently live in prose artifacts such as `audit_cov.md`, `spec.md`, and `test-spec.md`; those are not substitutes for executable tests.

## Highest-priority issues
### Must fix before merge
1. Verify and correct the generated-outcome formula in `R/edid-cov-eif.R`.
The current code uses `m_t - m_tp` keyed by `gp_lookup`, but Eq. (3.9)/(4.4) in the paper uses two distinct nuisance components:
- `m_{8,t,tpre}(X)` for the never-treated comparison trend
- `m_{g',tpre,1}(X)` for the comparison-cohort pretrend
Do not assume the current implementation is algebraically equivalent. Derive it with statsclawn and repair the code if it is not.

2. Verify and correct the covariate EIF in `R/edid-cov-eif.R`.
The current implementation is `drop(gen_out_mat %*% weights) - att_gt`, then centered. That is not a direct implementation of Eq. (3.10) / Theorem 3.2, which uses weighted pair-level influence functions, not just weighted generated outcomes minus a constant. Re-derive the EIF exactly and update the analytical SE path accordingly.

3. Verify and correct `Omega*` estimation in `R/edid-cov-eif.R`.
The paper’s Eq. (3.12) and Section 4 describe plug-in estimation of conditional covariance terms via kernel smoothing. The current code instead residualizes each generated-outcome column on a B-spline basis and smooths residual outer products. Either prove this is algebraically equivalent for the implemented estimator or replace it with a faithful plug-in estimator.

4. Fix `xformla` semantics in `R/edid-data.R` and `R/edid-validate.R`.
The current path uses `all.vars(xformla)` plus first-period raw columns. This silently ignores transformed terms like `I(x1^2)`, silently ignores all post-period values of time-varying covariates, and allows factor covariates in validation even though they fail later. Make formula handling coherent.

5. Make analytical covariate estimation reproducible.
The current cross-fitting folds are random and are not controlled by `seed`. Reuse the existing user-facing `seed` argument to control fold assignment for the analytical covariate path as well, without changing no-covariate behavior.

6. Add hard validation for unsupported covariate inputs.
At minimum, reject:
- any NA in covariates used by `xformla`
- any time-varying covariate within unit
- any formula feature you do not explicitly support
Do not silently coerce to partial behavior.

### Recommended improvements
1. Decide whether to support factor covariates properly via `model.matrix()` on one-row-per-unit invariant data, or to reject factors explicitly with an informative message. Do not keep the current “validation passes, runtime crashes” behavior.
2. Add diagnostics or stabilization for near-positivity and kernel/bandwidth pathologies. The constants in `R/edid-utils.R` suggest clipping/stabilization was planned but is currently unused.
3. Reassess clustered SE centering in `R/edid-inference.R`; `cluster_aggregate_edid()` centers cluster sums, but `safe_inference_edid()` currently does not use that helper.

## Required code inspections
Inspect these functions line by line against the paper with statsclawn:
- `compute_generated_outcomes_cov_edid()`
- `compute_omega_star_cov_edid()`
- `compute_eif_cov_edid()`
- `estimate_propensity_ratio_edid()`
- `estimate_all_conditional_means()`
- `fit_edid_cells()`
- `prepare_edid_panel()`
- `validate_edid_inputs()`

Also inspect the paper mappings in:
- `Efficient_DiD.pdf`
- `CSX2026.pdf`

Specific formula checks:
1. Map Eq. (4.1)-(4.2) to `estimate_propensity_ratio_edid()`.
2. Map Eq. (3.9)/(4.4) to `compute_generated_outcomes_cov_edid()`.
3. Map Eq. (3.10) and Theorem 3.2 to `compute_eif_cov_edid()`.
4. Map Eq. (3.12) and the kernel description in Section 4 to `compute_omega_star_cov_edid()`.
5. Confirm what must be cross-fitted and what need not be cross-fitted.

## Required code changes
### Must fix before merge
1. In `R/edid-cov.R`, restructure nuisance estimation so the conditional means required by the paper are estimable with unambiguous keys.
At minimum support the exact nuisance objects the paper requires, not just `(gp, period)` relative to period 1 if that is insufficient.

2. In `R/edid-cov-eif.R`, rewrite `compute_generated_outcomes_cov_edid()` to match the paper’s generated outcome exactly.
If the current self-comparison remapping `gp == g -> Inf` remains correct, keep it, but only after re-deriving the formula with statsclawn.

3. In `R/edid-cov-eif.R`, rewrite `compute_eif_cov_edid()` from the paper, not from the current “weighted generated outcomes minus constant” shortcut.
Make the returned EIF the actual influence function used by `safe_inference_edid()` and by the bootstrap path.

4. In `R/edid-cov-eif.R`, either:
- replace `compute_omega_star_cov_edid()` with a faithful plug-in estimator for Eq. (3.12), or
- add a derivation note in code comments and in documentation showing why the implemented estimator is equivalent
Do not leave this as an undocumented heuristic.

5. In `R/edid-data.R`, stop using `all.vars(xformla)` as the effective formula semantics.
Recommended implementation:
- build a one-row-per-unit covariate frame after enforcing time invariance
- use `model.frame()` / `model.matrix()` on that frame
- store the resulting numeric design matrix in `panel_obj$covariate_matrix`
If you intentionally restrict the allowed formula language, reject unsupported formulas early in `validate_edid_inputs()`.

6. In `R/edid-validate.R`, add:
- NA checks for all covariates referenced by `xformla`
- time-invariance checks for all covariates referenced by `xformla`
- explicit handling for factors and transformed terms
- informative failures for unsupported formula features

7. In `R/edid-fit.R` and `R/edid.R`, make fold assignment reproducible.
Use the existing `seed` argument to control analytical cross-fitting folds in the covariate path. Keep `xformla = NULL` and `xformla = ~1` behavior unchanged on the no-covariate path.

### Recommended improvements
1. Add optional storage of fold assignments in the fit object for debugging, or at least expose an internal deterministic fold builder used in tests.
2. Add explicit checks or warnings when a training fold has too few observations for the requested spline basis.
3. Add a clearer warning or failure mode for severe extrapolation from `splines::bs()` boundary knots.

## Required new tests
Create dedicated executable covariate-path tests. Suggested files:
- `tests/testthat/test-edid-cov-basic.R`
- `tests/testthat/test-edid-cov-formula.R`
- `tests/testthat/test-edid-cov-eif.R`
- `tests/testthat/test-edid-cov-variance.R`
- `tests/testthat/test-edid-cov-validation.R`

Add these tests exactly or equivalently:

1. Formula semantics tests.
- `xformla = ~x1` and `xformla = ~x1 + I(x1^2)` must not be identical on a nonlinear DGP when the same fold assignment is used.
- `xformla = ~x1 * x2` must either be supported correctly or fail early with an explicit unsupported-formula error.
- Factor covariates must either work through `model.matrix()` or fail in validation with an informative error.

2. Time-varying covariate rejection tests.
- Modify only post-period values of a covariate while keeping period 1 fixed.
- `edid(..., xformla = ~x1)` must error before estimation.
- Do not allow the current silent “use first row only” behavior.

3. NA covariate rejection tests.
- Add NA in period 1 for one unit and expect an informative error.
- Add NA only in a later period and still expect an informative error.

4. Reproducibility tests.
- Two analytical covariate-path runs with the same `seed` must be exactly identical.
- Two runs without resetting RNG but with the same user `seed` must still be identical.
- A test should confirm that fold assignment is deterministic when `seed` is fixed.

5. No-covariate regression tests.
- Keep the existing `xformla = NULL` vs `xformla = ~1` equality test.
- Add a test asserting that `~1` equality is due to deliberate routing to the no-covariate path, not because the covariate path was exercised.
- If you add any internal constant-design covariate-path test, keep it internal; do not change user-facing no-covariate behavior unless fixing a confirmed bug.

6. Generated-outcome formula tests.
- Construct a small synthetic example with known nuisance functions and verify `compute_generated_outcomes_cov_edid()` against a direct translation of Eq. (3.9)/(4.4).
- Include a self-comparison pair case `gp = g`.
- Include a cross-cohort pair `gp != g`.

7. EIF tests.
- Reconstruct the paper EIF directly from Eq. (3.10) / Theorem 3.2 and compare it to `compute_eif_cov_edid()` on the same data.
- Test more than zero mean. Zero mean alone is not informative here because the current implementation can be mechanically centered.
- Check that the returned EIF matches the direct formula elementwise up to numerical tolerance.

8. Variance tests.
- For each post-treatment cell, compare:
  - empirical Monte Carlo SD of the ATT estimator
  - mean reported analytical SE
  - mean EIF plug-in variance
- Fail if the average SE / empirical SD ratio is systematically far from 1 for large `n`.

9. Coverage tests.
- Estimate per-cell 95% coverage, not only pooled summaries.
- Include separate tests for `pt_assumption = "all"` and `pt_assumption = "post"`.

10. Near-positivity and bandwidth-pathology tests.
- Create a DGP with highly imbalanced treatment assignment by `X`.
- Verify that the code warns or fails in a controlled way rather than silently returning unstable results.
- Create a constant-covariate bandwidth case and verify the fallback is explicit and tested.

## Required simulation studies
Create an executable covariate benchmark, suggested file:
- `benchmark/edid_cov_sim.R`

Use the statsclawn plugin to restate the target estimand and the EIF before writing the simulation summary.

Run these simulation designs:

1. Paper-style covariate DGP suite.
Replicate the author’s covariate DGPs as closely as possible from the paper/slides or the author code, especially the DGP corresponding to the reported `(3,4)` SE issue.
Use sample sizes:
- `n = 200`
- `n = 500`
- `n = 1000`
- `n = 2000`

2. Per-cell bias study.
For each `(g,t)` post-treatment cell, report:
- Monte Carlo mean ATT
- true ATT
- bias
- RMSE

3. Per-cell variance calibration study.
For each `(g,t)` post-treatment cell, report:
- empirical SD across replications
- mean reported SE
- SE ratio = mean reported SE / empirical SD

4. Per-cell coverage study.
For each `(g,t)` post-treatment cell, report:
- nominal 95% coverage
- average CI length

5. EIF diagnostics.
For each post-treatment cell, report:
- Monte Carlo mean of the estimated EIF
- Monte Carlo variance implied by the EIF
- empirical ATT variance

6. Misspecification stress tests.
Run at least two designs:
- propensity ratio relatively well specified, outcome regression harder
- outcome regression relatively well specified, propensity ratio harder
The goal is not to prove exact double robustness in finite samples, but to detect whether the repaired implementation behaves in the expected direction.

7. Cross-fitting stability study.
For one moderate DGP, compare results across several fixed fold seeds.
Report variability induced by fold assignment itself.

## Acceptance criteria
### Must satisfy before merge
1. A written derivation note is added to the PR or documentation showing the exact mapping from paper Eq. (3.9), Eq. (3.10), Eq. (3.12), Eq. (4.1)-(4.4) to the repaired code.
2. `compute_generated_outcomes_cov_edid()`, `compute_eif_cov_edid()`, and `compute_omega_star_cov_edid()` are either shown equivalent to the paper or are rewritten to be faithful to it.
3. The analytical covariate path is deterministic when `seed` is fixed.
4. NA covariates and time-varying covariates fail early with informative errors.
5. Formula handling is coherent: transformed terms are either supported correctly or rejected explicitly.
6. Factor behavior is coherent: supported correctly or rejected explicitly in validation.
7. The no-covariate path remains unchanged for `xformla = NULL` and `xformla = ~1`, except for a confirmed bug fix that also affects the no-covariate estimator.
8. Dedicated executable covariate-path tests are present in `tests/testthat/`.
9. The benchmark script for the covariate path is checked in and runnable.

### Statistical acceptance targets
1. Bias must decrease with `n` on the main paper-style DGPs.
2. For the largest sample sizes, per-cell mean SE / empirical SD should be close to 1 and clearly not show a systematic severe inflation/deflation pattern.
3. Per-cell 95% coverage should be reported and should be reasonably near nominal for the main post-treatment cells.
4. A subtly wrong EIF sign/scaling convention must be detectable by the new tests.

## Guardrails while editing
1. Do not change behavior on the no-covariate path unless fixing a confirmed bug.
2. Do not silently broaden or narrow the `edid()` API. If you change what `xformla` supports, document the decision and add validation/tests for it.
3. Preserve user-facing argument names and object structure whenever possible.
4. Do not use `audit_cov.md`, `spec.md`, or `test-spec.md` as proof that the implementation is correct. Treat them as background only.
5. Prefer executable tests and simulation tables over prose claims.
6. Keep all statistical fixes isolated and documented. Do not mix formula repairs, validation changes, and simulation additions into one undocumented blob.

## Notes for documentation
1. Document every statistical change, not just code changes.
2. Document every new test and every new simulation design.
3. State clearly whether `xformla` now has full `model.matrix()` semantics or a restricted supported subset.
4. Document how analytical fold assignment is seeded and how users can reproduce covariate-path results.
5. If `Omega*` estimation remains computationally expensive, document the complexity and any practical limits.
