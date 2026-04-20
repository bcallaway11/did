# Test Specification — EDiD Covariate Extension

**Feature**: `xformla` argument for doubly-robust covariate adjustment in `edid()`
**Branch**: `feature-efficient-DiD-estimator`

---

## 1. Behavioral Contract

The covariate path must satisfy all of the following observable behaviors:

1. **Backward compatibility**: Calling `edid()` with `xformla = NULL` or `xformla = ~1` (or without the argument) produces results numerically identical to the current no-covariate implementation (tolerance 1e-10 on ATT estimates).

2. **Output structure invariance**: The `edid_fit` object returned with non-trivial `xformla` has identical S3 class and slot names to the no-covariate case.

3. **Correct doubly-robust estimator**: When either the propensity ratio or the conditional mean is correctly specified, the ATT estimate is consistent. Tested via: (a) correct propensity, misspecified outcome; (b) misspecified propensity, correct outcome.

4. **Dimension correctness**: For a dataset with n units, T periods, and K_cells = n_cohorts × (T-1) cells, the EIF matrix (when `store_eif = TRUE`) is n × K_cells.

5. **Nuisance function bounds**: Propensity ratio estimates are non-negative everywhere. Conditional mean estimates are finite (no NA, no Inf) when covariates are well-behaved.

6. **Validity flag**: Cell-level `inference_valid` flag follows the same logic as the no-covariate path (depends on n_pairs and SE computation).

---

## 2. Test Files to Create

### 2.1 `tests/testthat/test-edid-cov-basic.R`

Tests that verify backward compatibility, output structure, and the trivial-formula path.

### 2.2 `tests/testthat/test-edid-cov-nuisance.R`

Unit tests for individual nuisance functions: basis matrix, propensity ratio, conditional mean.

### 2.3 `benchmark/edid_cov_sim.R`

Simulation benchmark comparing edid() covariate path to author's reference implementation for DGPs 5–8 from Sim_10DGP.R.

---

## 3. Test Scenarios

### 3.1 Trivial xformla → Identical to no-covariate path

**File**: `test-edid-cov-basic.R`

**Setup**:
```r
set.seed(42)
n_units  <- 100
n_periods <- 6
# Balanced panel: cohorts 3, 5, Inf
# Generate standard DGP (no covariate-dependent treatment)
# [See DGP helper at end of this spec, Section 6]
df <- make_simple_panel_edid(n = 100, seed = 42)

fit_nocov <- edid(data = df, yname = "y", idname = "id",
                  tname = "period", gname = "first_treat",
                  pt_assumption = "all")

fit_trivial_null <- edid(data = df, yname = "y", idname = "id",
                         tname = "period", gname = "first_treat",
                         xformla = NULL, pt_assumption = "all")

fit_trivial_tilde1 <- edid(data = df, yname = "y", idname = "id",
                           tname = "period", gname = "first_treat",
                           xformla = ~1, pt_assumption = "all")
```

**Assertions**:
- `expect_equal(fit_nocov$att_gt$att, fit_trivial_null$att_gt$att, tolerance = 1e-10)`
- `expect_equal(fit_nocov$att_gt$att, fit_trivial_tilde1$att_gt$att, tolerance = 1e-10)`
- `expect_equal(fit_nocov$att_gt$se,  fit_trivial_tilde1$att_gt$se,  tolerance = 1e-10)`

### 3.2 Output structure invariance

**File**: `test-edid-cov-basic.R`

**Setup**:
```r
df_cov <- make_covariate_panel_edid(n = 150, seed = 123)
# df_cov has columns: id, period, y, first_treat, x1 (continuous covariate)

fit_cov <- edid(data = df_cov, yname = "y", idname = "id",
                tname = "period", gname = "first_treat",
                xformla = ~x1, pt_assumption = "all", store_eif = TRUE)
```

**Assertions**:
- `expect_s3_class(fit_cov, "edid_fit")`
- `expect_true(all(c("att_gt","overall","event_study","group","eif","bstrap","call") %in% names(fit_cov)))`
- `expect_true(is.data.frame(fit_cov$att_gt))`
- `expect_true(all(c("group","time","att","se","ci_lower","ci_upper","t_stat","p_value","n_pairs","is_pre") %in% names(fit_cov$att_gt)))`
- `expect_false(anyNA(fit_cov$att_gt$att))`  [all cells have valid estimates]
- `expect_equal(nrow(fit_cov$eif), 150L)`

### 3.3 Covariate path runs without error, 1D covariate

**File**: `test-edid-cov-basic.R`

**Setup**: Simple balanced panel with cohorts g=3, g=5, g=Inf; n=200; one continuous covariate `x1` that predicts treatment cohort assignment.

**Assertions**:
- `expect_no_error(edid(...))`
- ATT estimates are finite (no NA, no Inf)
- SE estimates are positive and finite
- EIF matrix has exactly n rows and n_cells columns

### 3.4 Covariate path, 2D covariate, PT-Post

**File**: `test-edid-cov-basic.R`

**Setup**: 2 cohorts (g=3, g=Inf), T=5 periods, n=200, covariates x1 and x2. Use `pt_assumption = "post"`.

**Assertions**:
- Runs without error
- H=1 pair per cell; n_pairs column in att_gt is 1 for all post cells
- `fit$att_gt$att` for pre-periods is NA (no valid pairs)

### 3.5 xformla validation errors

**File**: `test-edid-cov-basic.R`

| Input | Expected error message contains |
|-------|--------------------------------|
| `xformla = "x1"` (character, not formula) | "formula" |
| `xformla = ~nonexistent_col` | "nonexistent_col" |
| `xformla = ~x1` when x1 has NA values | "NA" |
| `xformla = ~x1` when x1 varies within unit | "time-invariant" or "constant" |

Use `expect_error(edid(...), regexp = "...")`.

### 3.6 xformla near-zero variance warning

**File**: `test-edid-cov-basic.R`

- Add a constant column to df. Call `edid(..., xformla = ~const_col)`.
- `expect_warning(..., regexp = "variance")` (or the near-zero variance message)

---

## 4. Nuisance Function Unit Tests

**File**: `tests/testthat/test-edid-cov-nuisance.R`

### 4.1 `build_crossfit_folds_edid()`

```r
folds <- build_crossfit_folds_edid(n = 100, K = 5L, seed = 1)
expect_length(folds, 100)
expect_true(all(folds %in% 1:5))
expect_equal(build_crossfit_folds_edid(100, K=5, seed=1),
             build_crossfit_folds_edid(100, K=5, seed=1))  # reproducible
```

### 4.2 `build_basis_matrix_edid()` — 1D input

```r
set.seed(1)
X <- matrix(rnorm(50), ncol = 1)
B <- build_basis_matrix_edid(X, bs_df = 4L)
expect_equal(nrow(B), 50L)
expect_equal(ncol(B), 4L)   # intercept = TRUE gives df columns for 1D
expect_true(all(is.finite(B)))
```

### 4.3 `build_basis_matrix_edid()` — 2D input

```r
X2 <- matrix(rnorm(100), ncol = 2)
B2 <- build_basis_matrix_edid(X2, bs_df = 4L)
expect_equal(nrow(B2), 50L)
expect_equal(ncol(B2), 4L + 3L)   # 4 + (4-1) = 7 columns
expect_true(all(is.finite(B2)))
```

### 4.4 `estimate_propensity_ratio_edid()` — dimension and non-negativity

```r
set.seed(42)
n_tr <- 80; n_te <- 20
X_tr <- matrix(rnorm(n_tr), ncol = 1)
X_te <- matrix(rnorm(n_te), ncol = 1)
G_tr <- sample(c(rep(3, 40), rep(0, 40)))
r_hat <- estimate_propensity_ratio_edid(X_tr, G_tr, X_te, g=3, gp=0, bs_df=4)
expect_length(r_hat, n_te)
expect_true(all(is.finite(r_hat)))
expect_true(all(r_hat >= 0))
```

### 4.5 `estimate_propensity_ratio_edid()` — constant-X special case

When all X values are the same, the ratio should be approximately p_g / p_{gp}:

```r
n <- 100
X_const  <- matrix(rep(1, n), ncol=1)   # all units have X=1
G_const  <- c(rep(3, 30), rep(0, 70))
X_te     <- matrix(rep(1, 10), ncol=1)
r_hat    <- estimate_propensity_ratio_edid(X_const, G_const, X_te, g=3, gp=0, bs_df=4)
# Expected: approximately 30/70 for all units (constant ratio)
expect_true(all(abs(r_hat - 30/70) < 0.15))
```

### 4.6 `estimate_conditional_mean_edid()` — dimension check

```r
n_tr <- 80; n_te <- 20
X_tr <- matrix(rnorm(n_tr), ncol = 1)
X_te <- matrix(rnorm(n_te), ncol = 1)
G_tr <- c(rep(0, 40), rep(3, 40))
Y_tr <- 2 + 0.5 * X_tr[,1] + rnorm(n_tr)
m_hat <- estimate_conditional_mean_edid(X_tr, Y_tr, G_tr, X_te, gp=0, bs_df=4)
expect_length(m_hat, n_te)
expect_true(all(is.finite(m_hat)))
```

### 4.7 `estimate_conditional_mean_edid()` — recovers linear mean

When the true conditional mean is linear and we use sufficient df:

```r
set.seed(7)
n <- 500
X  <- matrix(runif(n, -2, 2), ncol=1)
G  <- rep(0, n)   # all one cohort
Y  <- 1 + 2*X[,1] + rnorm(n, sd=0.1)
X_te <- matrix(seq(-1.5, 1.5, length.out=20), ncol=1)
m_hat <- estimate_conditional_mean_edid(X, Y, G, X_te, gp=0, bs_df=6)
true_m <- 1 + 2*X_te[,1]
expect_true(max(abs(m_hat - true_m)) < 0.5)   # loose tolerance; sieve approx
```

### 4.8 `compute_eif_cov_edid()` — zero mean

```r
set.seed(99)
n <- 100; H <- 3
gen_out <- matrix(rnorm(n*H), n, H)
w <- c(0.5, 0.3, 0.2)
att <- sum(w * colMeans(gen_out))
eif <- compute_eif_cov_edid(panel_obj = NULL, gen_out_mat = gen_out,
                              weights = w, att_gt = att)
# After centering, EIF should be zero-mean
expect_equal(mean(eif), 0, tolerance = 1e-12)
expect_length(eif, n)
```

(Note: `compute_eif_cov_edid` should not use `panel_obj` for the centering step; it operates on `gen_out_mat` directly. Passing `panel_obj = NULL` should work if panel_obj is not used in the body.)

---

## 5. Regression Test: xformla = ~1 → Identical to No-Covariate

This is the most critical regression test.

**File**: `tests/testthat/test-edid-cov-basic.R`

**Setup**: Use the existing test DGP from `test-edid-basic.R` or create a fresh one. Both PT-All and PT-Post must be tested.

**Expected**: All ATT values, SE values, and EIF columns agree to tolerance 1e-10.

```r
for (pt in c("all", "post")) {
  fit_base  <- edid(df, yname="y", idname="id", tname="period",
                    gname="first_treat", pt_assumption=pt, store_eif=TRUE)
  fit_tilde <- edid(df, yname="y", idname="id", tname="period",
                    gname="first_treat", xformla=~1, pt_assumption=pt, store_eif=TRUE)
  expect_equal(fit_base$att_gt$att, fit_tilde$att_gt$att, tolerance=1e-10,
               label=paste("ATT PT", pt))
  expect_equal(fit_base$att_gt$se,  fit_tilde$att_gt$se,  tolerance=1e-10,
               label=paste("SE PT", pt))
  if (!is.null(fit_base$eif) && !is.null(fit_tilde$eif)) {
    expect_equal(fit_base$eif, fit_tilde$eif, tolerance=1e-10,
                 label=paste("EIF PT", pt))
  }
}
```

---

## 6. Simulation Benchmark: Author's Reference Implementation

**File**: `benchmark/edid_cov_sim.R`

This script is NOT run as part of `R CMD check`. It is a manual validation tool run by the developer. It should produce a comparison table showing edid() ATT estimates vs the author's `simu()` estimates for DGPs 5–8.

### 6.1 DGP Definitions

The four DGPs replicate the covariate DGPs from `Sim_10DGP.R`:

**DGP 5**: Two cohorts (g=3, g=Inf), T=4 periods. Covariates: X1, X2 ~ Truncated-Normal(-3,3). Treatment assignment: logistic with `P(G=3|X) = logistic(-0.5 + 0.5*X1 - 0.25*X2)`. Conditional parallel trends. Index: `X1 + X2`. AR(1) errors with rho=0.5, heteroskedastic variance `exp(0.15*(X1-X2))`. True ATT(3,3) = 3, ATT(3,4) = 4.5.

**DGP 6** (author's DGP 7): Four cohorts (g=2, g=3, g=4, g=Inf), T=4. Same covariate structure. True ATTs: ATT(2,2)=2, ATT(2,3)=2.2, ATT(2,4)=2.4, ATT(3,3)=3, ATT(3,4)=4.5, ATT(4,4)=4.

**DGP 7** (author's DGP 8): Same as DGP 6 but with AR(1) errors with group-specific heteroskedastic variance `sd=1.5*g/3`.

**DGP 8** (author's DGP 9): Same cohort structure, non-linear trend `X1 + X2 + X1^2 + X2^2`.

### 6.2 Simulation design

```r
R       <- 500      # replications
n       <- 300      # sample size
seed    <- 20240417
set.seed(seed)

results_table <- data.frame(
  dgp = integer(),
  group = numeric(), time = numeric(),
  true_att = numeric(),
  edid_att_mean = numeric(), edid_att_bias = numeric(),
  edid_se_mean = numeric(),
  edid_cover_95 = numeric()
)
```

For each DGP and each replication:
1. Generate panel data from the DGP function.
2. Reshape to long format (id, period, y, first_treat, x1, x2).
3. Call `edid(..., xformla = ~x1 + x2, pt_assumption = "all")`.
4. Collect ATT estimates and SEs for each (g, t) cell.

### 6.3 Acceptance criteria for the benchmark

These are the numerical tolerances for declaring the covariate implementation correct:

| Metric | Threshold | Notes |
|--------|-----------|-------|
| Absolute bias of ATT | < 0.10 | At n=300, R=500; some simulation variability expected |
| SE ratio (est SE / empirical SD of ATT) | [0.7, 1.5] | Loose because Omega* smoothing adds variance |
| 95% CI coverage | [0.85, 0.99] | Wide band acceptable for a draft implementation |
| Failure rate (NA ATT) | < 0.05 | < 5% of replications |

These thresholds apply per (g, t) cell. If any cell fails, the output table should flag it.

### 6.4 Comparison column: author's implementation

Optionally call `simu(n=n, dgp=k)` from `Sim_10DGP.R` to get the author's ATT estimate for the same data realization and compute the difference. These should agree to within simulation noise (absolute difference < 0.20 for most cells in most replications).

---

## 7. DGP Helpers for testthat Files

Define in a `tests/testthat/helper-edid-cov.R` file:

```r
make_simple_panel_edid <- function(n = 100, seed = 42) {
  set.seed(seed)
  n_periods <- 6
  unit_ids  <- rep(1:n, each = n_periods)
  time_ids  <- rep(1:n_periods, times = n)
  cohort_assign <- rep(c(3, 5, Inf),
    times = c(ceiling(n/3), ceiling(n/3),
              n - 2*ceiling(n/3)))[1:n]
  first_treat_vec <- cohort_assign[unit_ids]
  treat_effect    <- as.numeric(time_ids >= first_treat_vec)
  y_vals <- 0.5*time_ids + treat_effect + rnorm(n*n_periods, sd=0.5)
  data.frame(id=unit_ids, period=time_ids, y=y_vals,
             first_treat=first_treat_vec)
}

make_covariate_panel_edid <- function(n = 150, seed = 123, n_periods = 6) {
  set.seed(seed)
  unit_ids <- rep(1:n, each=n_periods)
  time_ids <- rep(1:n_periods, times=n)
  # x1 is time-invariant, predicts cohort
  x1_units <- rnorm(n)
  # Cohort assignment depends on x1 (conditional PT setting)
  p_treat3 <- plogis(0.5 + 0.5*x1_units)
  u        <- runif(n)
  cohort_assign <- ifelse(u < p_treat3/2, 3,
                   ifelse(u < p_treat3, 5, Inf))
  first_treat_vec <- cohort_assign[unit_ids]
  x1_long   <- x1_units[unit_ids]
  # Outcome: covariate-dependent trend (conditional PT)
  treat_effect <- as.numeric(time_ids >= first_treat_vec)
  y_vals <- 0.5*time_ids + 0.3*x1_long*time_ids + treat_effect +
            rnorm(n*n_periods, sd=0.5)
  data.frame(id=unit_ids, period=time_ids, y=y_vals,
             first_treat=first_treat_vec, x1=x1_long)
}
```

---

## 8. Validation Commands

Run these after builder implements the feature:

```r
# Install and load
devtools::load_all("/path/to/did")

# Unit tests
testthat::test_file("tests/testthat/test-edid-cov-basic.R")
testthat::test_file("tests/testthat/test-edid-cov-nuisance.R")

# Full test suite (regression guard)
devtools::test()

# R CMD check (must pass with 0 ERRORs, 0 WARNINGs)
devtools::check()

# Benchmark (manual, not part of check)
source("benchmark/edid_cov_sim.R")
```

---

## 9. Property-Based Invariants

The following mathematical properties must hold for any valid input to `edid()` with covariates:

1. **EIF zero-mean**: `abs(mean(fit$eif[,k])) < 1e-8` for all cells k where `fit$eif` is non-NULL.

2. **SE positivity**: All `fit$att_gt$se` values are either positive or NA (NA only for cells with no valid pairs).

3. **CI containment**: `fit$att_gt$ci_lower < fit$att_gt$att` and `fit$att_gt$att < fit$att_gt$ci_upper` for all non-NA cells.

4. **Overall ATT in convex hull**: For post-treatment cells, `fit$overall$att` is between the minimum and maximum of finite cell-level ATTs (not a strict requirement, but a sanity check).

5. **Propensity ratio non-negativity**: The propensity ratio vector returned by `estimate_propensity_ratio_edid()` must have `all(r_hat >= 0)`.

6. **Weights sum to 1**: `abs(sum(weights) - 1) < 1e-12` for the output of `compute_efficient_weights_edid()`.

---

## 10. Edge Case Scenarios

| Scenario | Expected behavior |
|----------|------------------|
| Only 2 units per cohort | Warning emitted; falls back to nocov for affected cells |
| Single covariate is binary (0/1) | Runs; basis function may warn but does not error |
| xformla has 3+ covariates | Runs; basis grows to `bs_df + (d-1)*(bs_df-1)` columns |
| All units have identical X | Propensity ratio returns constant (p_g/p_gp); no error |
| n < 50 | Runs with warnings about nuisance estimation; does not error |
| PT-Post with xformla | One pair per cell; same formula applies; no error |
| clustervars with xformla | Cluster-robust SE computed via existing `safe_inference_edid()` |
| bstrap=TRUE with xformla | Bootstrap uses EIF from cov path; no changes to bootstrap code needed |
| control_group="notyettreated" with xformla | Runs; last cohort relabeled as never-treated; covariates handled correctly |

---

## 11. Cross-Reference with No-Covariate Path

The covariate path MUST reduce to the no-covariate path when `xformla = ~1`. The mechanism:

- `build_basis_matrix_edid` with intercept-only formula returns a column of 1s.
- `estimate_propensity_ratio_edid` with a constant basis returns `n_g / n_gp` for all units.
- `estimate_conditional_mean_edid` with a constant basis returns `mean(Y_delta[mask_gp])` for all units.
- The generated outcome formula then simplifies to the unconditional DiD moment.
- Omega*(X) with constant residuals collapses to the sample variance Omega (no X-dependence).

This reduction is verified by the regression test in Section 5.
