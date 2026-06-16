# Effective-n ridge/LW fix — full-adaptation audit (round complete)

**Branch:** `overnight-aggfix-robust` (repo `/Users/pcostag/Documents/GitHub/did`)
**Mandate (Pedro):** "the ridge penalty must make sense under weights … rock solid."
**Date:** 2026-06-15

## The fix

The vanishing ridge and Ledoit–Wolf intensities that stabilize the *weights*
(not the SE) divided by the raw unit count `panel_obj$n`. Under dispersed
`weightsname` observation weights the scale-invariant weighted `Omega*` is
backed by far fewer than `n` independent contributions (the heavily-weighted
units dominate), so the raw count **under-regularized the weighted efficient
leg by the factor `n / n_eff >= 1`** (the Kish design factor). The denominator
is now the per-cell Kish effective sample size

```
n_eff = if (is.null(unit_weights)) panel_obj$n              # unweighted: full n, byte-identical
        else  (sum w)^2 / sum(w^2)   over the cell's ACTIVE units   # weighted: Kish ESS
```

`n_eff` is implemented once in `R/edid-utils.R` (`n_eff_edid()`), with the
active-unit mask `active_mask_nocov_edid()` = treated cohort ∪ never-treated ∪
comparison cohorts (the nonzero rows of Ψ).

### Verified weights field (Task B)
`panel_obj$unit_weights` — `NULL` on the unweighted default; otherwise one
mean-1-normalized weight per unit (length `n`, aligned to `panel_obj$all_units`),
built in `prepare_edid_panel()` (`R/edid-data.R:109–116`). A cell's active units
are indexed by the cohort masks `panel_obj$cohort_masks[[as.character(g)]]` and
`panel_obj$never_treated_mask`.

## The three call sites (Task C)

| # | Site | Before | After |
|---|------|--------|-------|
| 1 | no-cov ridge, `R/edid-fit.R:949–960` | `lam_r <- Hh / panel_obj$n` | `n_eff <- n_eff_edid(unit_weights, active_mask_nocov_edid(g,pairs,panel_obj), panel_obj$n); lam_r <- Hh / n_eff` |
| 2a | no-cov LW intensity, `R/edid-nocov.R:375–380` | `b2 <- (q4/n^2 - n*sum(omega*omega))/n^2` | `b2_legacy <- (…)/n^2; n_eff <- n_eff_edid(…, n); b2 <- b2_legacy * (n / n_eff)` |
| 2b | no-cov LW EE chain rule, `R/edid-nocov.R:562–567` | `kappa_i <- -(2/(n*d2))*beta_i - …` | `n_eff <- n_eff_edid(…, n); kappa_i <- -(2/(n_eff*d2))*beta_i - …` |
| 3 | cov-path lift, `R/edid-cov-eif.R` `.edid_cov_ridge_lift_array/_pooled` (+ the EE `tr(C)/n` term) | `(H/n_full)`, `tr(C)/n` | `(H/n_eff)`, `tr(C)/.ridge_neff`; `n_eff = n_eff_edid(unit_weights, all-TRUE, n_full)` |

Call sites for the cov lift updated in `edid-cov-eif.R` (×2), `edid-cov-kernfast.R`
(×2), `edid-cov-sieve.R` (×2) to forward `panel_obj$unit_weights`.

**Byte-identity lever.** Site 2a keeps the *verbatim* legacy `b2_legacy`
expression and multiplies by `n / n_eff`, which is **exactly 1.0** unweighted
(`n_eff_edid` returns the same full `n`, identical doubles → no FP reassociation).
The first draft used `pi_hat/n_eff` directly and lost one ULP on one cell
(`1.776e-15`); the correction-factor form restores exact byte-identity.

## Properties that must hold — all verified

| Property | Result |
|---|---|
| **Byte-identical unweighted** (no-op at w-equal) | **`max|delta| = 0.000e+00`** over att/se/lambda/cond across 15 unweighted configs (no-cov + cov; none/ledoit_wolf/ridge; 3 staggered designs) |
| **`n_eff == n` unweighted** | unit-checked: `n_eff_edid(NULL, …, 10) == 10`; constant weight column → `n_eff == n_act` (4 for equal-w) |
| **λ → 0 asymptotically** | ridge `H/n_eff`: unweighted 0.050→0.003 (n 40→640); weighted 0.145→0.0068; `n_eff` grows ∝ n |
| **Sign-correct / well-conditioned weighted leg** | weighted ridge intensity LARGER than raw-n by `n/n_eff` (≈2.9× at n=40) → more shrinkage, as intended |
| **EE chain rule consistent (FD-oracled)** | weighted FD oracle (interior λ=0.475, n/n_eff=1.71): analytic vs FD Jacobian **9.4e-10** (tol 1e-5) |

## Full adaptation & audit battery

1. **Impact enumeration.** ratio_method {exp,direct}, weight_scheme
   {efficient,averaged,gmm,uniform}, omega smoother {ridge,ledoit_wolf,none},
   estimation_effect, misspec_robust, clustervars, both bootstraps, all
   aggregates, toolkit {weights,sargan,hausman,frontier,print,summary} —
   **option-matrix smoke sweep: ALL PASS (0 bad combinations)** on weighted
   no-cov data + unweighted cov. Unaffected (reason): covariate weighted path
   (scoped out — errors by design, so site 3 is a structural no-op today);
   PT-Post / uniform / H=1 (no weights estimated); `nocov_shrink="none"`
   (no intensity); thin-cohort guard / trim masks (operate upstream of the
   intensity).
2. **Estimation effects.** `n_eff` is a function of the FIXED weights only,
   constant in Ω̂ → the ridge term's `d/dΩ = I` is unchanged (plain-map EE branch
   unchanged). The LW `b²` chain rule derivative `d(b²) = -(2/n_eff)⟨Ω,dE⟩` was
   re-derived and **FD-oracled** (above). No channel skipped.
3. **Invariants.** Unweighted no-cov PT-All/PT-Post + covariate paths
   byte-identical (`max|delta| = 0`). FD oracle (`test-edid-nocov-estimation-effect.R`)
   passes with `NOT_CRAN=true`.
4. **Audit battery.** FD oracles (pass); **full testthat = 3977 PASS, 0 FAIL,
   89 pre-existing diagnostic WARN, 11 CRAN-skip**; targeted MC calibration
   (below); option-matrix sweep (pass).
5. **Docs same round.** roxygen regenerated (`n_eff_edid.Rd`,
   `active_mask_nocov_edid.Rd`, updated `shrink_omega_nocov_edid.Rd` /
   `compute_nocov_ee_correction_edid.Rd`); NEWS bullet added.
6. **No silent number changes — before/after table.** Only WEIGHTED ridge/LW
   fits move (recorded below); unweighted is byte-identical.

### Targeted MC (weighted ridge, disp=1.0, n=140, 400 reps, anchor ATT=1.7868)

| | MC SD | mean(SE) | mean(SE)/MC SD | 95% cover |
|---|---|---|---|---|
| BEFORE (raw n) | 0.1566 | 0.1389 | 0.887 | 0.895 |
| AFTER (n_eff)  | 0.1563 | 0.1401 | 0.896 | 0.900 |

The fix moves calibration in the correct direction (larger SE, better ratio and
coverage). The effect is modest because the ridge SE is the empirical variance
of the realized weighted IF (Ω does not enter it directly); the fix's primary
role is to make the *weight regularization* scale correctly under weights and be
a no-op unweighted — both achieved. The residual ~0.90 coverage at n=140 under
heavy dispersion is present BEFORE too (an analytic-plug-in-SE property, the
job of the multiplier bootstrap / estimation_effect), not introduced here.

## Artifacts
- `run_baseline.R`, `baseline.rds` (pre), `postfix.rds` (post), `compare.R`
- `lambda_vanish.R`, `fd_oracle_weighted.R`, `option_matrix_sweep.R`
- `mc_coverage.R`, `mc_beforeafter.R` (+ `.rds`/`.log`)
- full testthat log: `/tmp/gate_runs/ridgefix/full_suite.log`
