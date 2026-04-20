# Spec: edid() Covariate Path Audit and Repair

**Run**: REQ-20260417-EDID-COV
**Date**: 2026-04-17
**Status**: READY FOR BUILDER

---

## 1. Paper Formula Reference (R-ready notation)

### Eq. (3.9) — Generated Outcome for pair (g', tpre) targeting ATT(g,t)

```
Y_tilde[g',tpre] = (G_g / pi_g) * (Y_t - Y_1 - m_inf_t_tpre(X) - m_g'_tpre_1(X))
                  - r[g,inf](X) * (G_inf / pi_g) * (Y_t - Y_tpre - m_inf_t_tpre(X))
                  - r[g,g'](X)  * (G_g' / pi_g) * (Y_tpre - Y_1 - m_g'_tpre_1(X))
```

Where:
- `G_g = 1{G_i = g}`, `G_inf = 1{G_i = Inf}`, `G_g' = 1{G_i = g'}`
- `pi_g = E[G_g] = n_g/n`
- `m_inf_t_tpre(X) = E[Y_t - Y_tpre | G=Inf, X]` (never-treated trend from tpre to t)
- `m_g'_tpre_1(X) = E[Y_tpre - Y_1 | G=g', X]` (comparison cohort trend from 1 to tpre)
- `r[g,inf](X) = p_g(X) / p_inf(X)` (propensity ratio g vs never-treated)
- `r[g,g'](X) = p_g(X) / p_g'(X)` (propensity ratio g vs comparison cohort g')

**Self-comparison case (g' = g)**: The third term vanishes because there is no separate comparison cohort (G_g' = G_g is already used in the first term). The formula reduces to Eq. (3.2):

```
Y_tilde[g,tpre] = (G_g / pi_g) * (Y_t - Y_tpre - m_inf_t_tpre(X))
                 - r[g,inf](X) * (G_inf / pi_g) * (Y_t - Y_tpre - m_inf_t_tpre(X))
               = ((G_g / pi_g) - r[g,inf](X) * G_inf / pi_g) * (Y_t - Y_tpre - m_inf_t_tpre(X))
```

### Eq. (3.10) — Influence Function for pair (g', tpre)

```
IF[g',tpre] = Y_tilde[g',tpre] + (G_g / pi_g) * ATT(g,t)
```

Note the PLUS sign. The generated outcome has mean = -ATT(g,t) * pi_g, so IF has mean zero.

### Eq. (3.12) — Omega* conditional covariance matrix

The (j,k)-th element of Omega*(X) for pairs j=(g'_j, t'_j) and k=(g'_k, t'_k):

```
Omega*[j,k](X) = (1/p_g(X)) * Cov(Y_t - Y_1, Y_t - Y_1 | G=g, X)
               + (1/p_inf(X)) * Cov(Y_t - Y_{t'_j}, Y_t - Y_{t'_k} | G=Inf, X)
               - 1{g = g'_j} / p_g(X) * Cov(Y_t - Y_1, Y_{t'_j} - Y_1 | G=g, X)
               - 1{g = g'_k} / p_g(X) * Cov(Y_t - Y_1, Y_{t'_k} - Y_1 | G=g, X)
               + 1{g'_j = g'_k} / p_{g'_j}(X) * Cov(Y_{t'_j} - Y_1, Y_{t'_k} - Y_1 | G=g'_j, X)
```

### Eq. (4.3)-(4.4) — Plug-in Estimator

```
ATT_hat(g,t) = E_n[ w(X_i)' * Y_tilde_hat(X_i) ]
w(X) = Omega_hat*(X)^{-1} * 1 / (1' * Omega_hat*(X)^{-1} * 1)
```

### EIF for the estimator

```
EIF_i = w(X_i)' * IF_hat_i  [i.e., weighted IF across pairs]
      = w(X_i)' * Y_tilde_hat_i + (G_g_i / pi_g) * ATT_hat(g,t)
```

Note: This is NOT `w' * gen_out - ATT`. It is `w' * (gen_out + G_g/pi_g * ATT)`.
The current code `drop(gen_out_mat %*% weights) - att_gt` is WRONG because:
- It subtracts `att_gt` from ALL units instead of adding `G_g/pi_g * att_gt`
- The sign is wrong (should be +, not -)
- The weighting by G_g/pi_g is missing (only treated units get the ATT correction)

---

## 2. Bug-by-Bug Analysis and Repair Plan

### Bug 1: Generated Outcome Formula (edid-cov-eif.R:27-106)

**Current code**: Lines 84-99 compute:
```r
m_diff <- m_t - m_tp                     # m_{gp,t} - m_{gp,tpre}
y_diff <- ow[,col_t] - ow[,col_tp]      # Y_t - Y_tpre
phi_j  <- (Ig/pi_g - Igp_j * r_j) * (y_diff - m_diff)
```

**Paper Eq. (4.4)**: For a cross-cohort pair (g' != g mapped to Inf):
```
Y_tilde = (G_g/pi_g) * (Y_t - Y_1 - m_inf_t_1(X) - m_g'_tpre_1(X))
        - r[g,inf](X) * (G_inf/pi_g) * (Y_t - Y_tpre - m_inf_t_tpre(X))
        - r[g,g'](X) * (G_g'/pi_g) * (Y_tpre - Y_1 - m_g'_tpre_1(X))
```

**Diagnosis**: The current code computes `(G_g/pi_g - G_gp*r) * (Y_t - Y_tpre - (m_t - m_tpre))`. This is only correct for self-comparison pairs (g'=g), where the third term vanishes and `m_diff = m_inf_t_tpre(X)`. For cross-cohort pairs, the paper has THREE separate terms with different outcome differencing structures (Y_t - Y_1, Y_t - Y_tpre, Y_tpre - Y_1), not a single `(Y_t - Y_tpre)` difference.

**However**: Looking more carefully at Eq. (3.9), the structure for cross-cohort pairs with the current pair enumeration needs careful analysis. The `gp_lookup` remapping to `Inf` means the code treats ALL pairs as self-comparison against never-treated. This is a fundamental design choice that requires re-examination.

**Looking at pair enumeration**: `enumerate_valid_pairs_edid()` returns pairs where `gp` can be ANY treated cohort (including `target_g`). But the current generated-outcome code remaps `gp == g` to `Inf`. For `gp != g`, the code does NOT remap, meaning it should use `r[g,gp]` and `m_{gp,...}`.

**Required fix**: Implement the full three-term formula from Eq. (4.4) for cross-cohort pairs. For self-comparison pairs (gp==g, remapped to Inf), the simplified two-term formula is correct.

**Additional nuisance keys needed**:
- For cross-cohort pairs: need `m_inf_t_tpre(X)` (never-treated trend), `m_g'_tpre_1(X)` (comparison cohort pretrend), AND `r[g,g'](X)` AND `r[g,inf](X)`
- Currently only estimating one propensity ratio per gp and one conditional mean per (gp, period)

### Bug 2: EIF Formula (edid-cov-eif.R:243-248)

**Current code**:
```r
eif <- drop(gen_out_mat %*% weights) - att_gt
eif <- eif - mean(eif)
```

**Paper Eq. (3.10) + Theorem 3.2**: The EIF is:
```
EIF_i = sum_j w_j(X_i) * IF_{j,i}
      = sum_j w_j(X_i) * (Y_tilde_{j,i} + (G_g_i / pi_g) * ATT(g,t))
      = sum_j w_j(X_i) * Y_tilde_{j,i} + (G_g_i / pi_g) * ATT(g,t)
```

**Required fix**: Replace `- att_gt` with `+ (Ig / pi_g) * att_gt` where `Ig = 1{G_i == g}`.

The mechanical centering `eif - mean(eif)` was masking this bug.

### Bug 3: Omega* Estimation (edid-cov-eif.R:133-225)

**Current code**: Residualizes generated outcomes on B-splines, then smooths residual outer products with Nadaraya-Watson kernel.

**Paper Eq. (3.12) + Section 4**: The paper estimates each covariance term `Cov(Y_t - Y_s, Y_t - Y_s' | G=g', X=x)` directly using kernel smoothing on raw residuals from conditional means, NOT on generated-outcome residuals.

**Diagnosis**: The current approach estimates `Cov(phi_j, phi_k | X)` (covariance of generated outcomes) rather than `Omega*(X)` from Eq. (3.12). These are NOT the same. `Omega*(X)` is defined in terms of outcome change covariances, scaled by propensity scores. The generated outcomes include propensity weights, cohort indicators, and nuisance corrections — their covariance includes many cross-terms that don't appear in `Omega*(X)`.

**Required fix**: Replace with faithful plug-in of Eq. (3.12). Each entry needs specific conditional covariances of outcome changes within specific cohorts, NOT generated-outcome residual covariances.

### Bug 4: xformla Semantics (edid-data.R:145-159)

**Current code**: Uses `all.vars(xformla)` to extract variable names, then `as.matrix(cov_df)`.

**Problems**:
- `all.vars(~I(x1^2))` returns `"x1"`, not the transformed term — so `I(x1^2)` is silently treated as `x1`
- `all.vars(~x1*x2)` returns `c("x1","x2")` — interaction is lost
- Factor columns become numeric via `as.matrix()` without proper dummy coding
- Time-varying covariates: only first-period values used, no check that values are constant

**Required fix**:
1. Build one-row-per-unit data frame
2. Use `model.matrix(xformla, data=unit_df)` to expand formula (handles I(), interactions, factors)
3. Remove intercept column (handled by the estimator)
4. Store the expanded numeric matrix

### Bug 5: Seed Not Threaded (edid-fit.R:38)

**Current code**: `fold_id <- build_crossfit_folds_edid(n = panel_obj$n, K = 5L, seed = NULL)`

The `seed = NULL` is hardcoded. The user's `seed` argument from `edid()` is never passed through.

**Required fix**: Thread `seed` from `edid()` → `fit_edid_cells()` → `build_crossfit_folds_edid()`.

### Bug 6: Missing Validation (edid-validate.R)

**Missing checks**:
- No NA check on covariate columns referenced by `xformla`
- No time-invariance check for covariates
- No rejection of unsupported formula features
- Factors pass validation but crash later

**Required fix**: Add to `validate_edid_inputs()`:
1. NA check: `anyNA(data[[v]])` for each covariate variable
2. Time-invariance check: verify each covariate is constant within unit
3. Factor handling: either support via model.matrix or reject explicitly
4. Formula feature support: document which features are supported

---

## 3. Implementation Decisions

### Factor support: USE model.matrix()
- `model.matrix(xformla, data=unit_df)[, -1, drop=FALSE]` (drop intercept)
- This handles factors, I(), interactions, poly() automatically
- Validate that the result is numeric with no NA/Inf columns
- Reject formulas that produce zero columns after expansion

### Omega* estimation: REPLACE with faithful plug-in
- Cannot prove equivalence of residualized generated-outcome covariance to Eq. (3.12)
- Implement direct kernel estimation of each covariance term in (3.12)
- Use same Nadaraya-Watson kernel machinery already in place
- Accept O(n^2 * H^2) complexity with existing warning

### Generated outcome: REWRITE for cross-cohort pairs
- Self-comparison (gp==g): keep current two-term formula (correct after gp→Inf remap)
- Cross-cohort (gp!=g): implement full three-term Eq. (4.4)
- Need additional nuisance estimates: `r[g,inf]` for all cells, `m_inf_t_tpre` and `m_g'_tpre_1` as separate conditional means

---

## 4. Required Test Files

### test-edid-cov-basic.R
- Smoke test: `edid(data, ..., xformla = ~x1)` runs without error
- ATT estimates are finite and non-NA for post-treatment cells
- SEs are finite and positive
- `xformla = ~1` routes to no-covariate path (regression test)

### test-edid-cov-formula.R
- `~x1` vs `~x1 + I(x1^2)` produce different results on nonlinear DGP
- `~x1 * x2` either works or fails with explicit error
- Factor covariates either work via model.matrix or fail in validation
- Time-varying covariate rejection
- NA covariate rejection

### test-edid-cov-eif.R
- Construct small example with known nuisances
- Verify generated outcomes match Eq. (3.9)/(4.4) analytically
- Verify EIF matches Eq. (3.10) elementwise
- Include self-comparison and cross-cohort pair cases
- EIF has mean approximately zero
- Wrong-sign detection: deliberately flip a sign and check EIF changes

### test-edid-cov-variance.R
- Reproducibility: same seed → identical results
- Fold assignment is deterministic with fixed seed
- Two independent runs with same seed match exactly

### test-edid-cov-validation.R
- NA in covariates → error
- Time-varying covariates → error
- Factor covariates → handled (either works or clean error)
- Unsupported formula → error with message

---

## 5. Simulation Design (benchmark/edid_cov_sim.R)

### DGP
- n = {200, 500, 1000, 2000}
- 3 cohorts: g=3, g=5, never-treated
- 6 time periods
- X ~ N(0,1), treatment effect heterogeneous in X
- Parallel trends satisfied by construction
- True ATT(g,t) known analytically

### Studies
1. Bias study: per-cell Monte Carlo bias, should decrease with n
2. Variance calibration: SE/empirical SD ratio, should approach 1
3. Coverage: per-cell 95% CI coverage, should approach 0.95
4. Cross-fitting stability: same DGP, different fold seeds
5. Misspecification stress: one spec correct, other harder

### Replications
- 200 replications minimum per design
- Report: bias, RMSE, mean SE, SE ratio, coverage, CI length

