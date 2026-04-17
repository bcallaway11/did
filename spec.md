# Implementation Spec — PT-All Pair Enumeration Fix

**Target branch**: `feature-efficient-DiD-estimator`
**Files modified**: `R/edid-pairs.R`, `R/edid-nocov.R`
**Files read-only (no changes)**: `R/edid-data.R`, `R/edid-utils.R`, `R/edid-linalg.R`, `R/edid-fit.R`

---

## 1. Notation

| Symbol | Type | Description |
|--------|------|-------------|
| `g` | scalar | Target treatment cohort (`target_g`) |
| `t` | scalar | Target time period (`target_t`) |
| `g'` (`gp`) | scalar | Comparison cohort — always finite under PT-All after this fix |
| `t'` (`tpre`) | scalar | Pre-period used in the identifying moment |
| `period_1` | scalar | `min(time_periods)` — universal first period |
| `eff_start(g')` | scalar | `g' - anticipation` — effective treatment start of comparison cohort |
| `G_treated` | vector | `treatment_groups` — sorted unique finite cohort values |
| `n_g` | scalar | Number of units in cohort g |
| `n_inf` | scalar | Number of never-treated units |
| `n_{g'}` | scalar | Number of units in comparison cohort g' |
| `pi_g` | scalar | `n_g / n` — cohort fraction |
| `pi_inf` | scalar | `n_inf / n` |
| `delta_g_t_1` | vector length `n_g` | `Y_g(t) - Y_g(period_1)` — treated group change |
| `delta_inf_j` | vector length `n_inf` | `Y_inf(t) - Y_inf(tpre_j)` — never-treated time change for pair j |
| `delta_gp_j` | vector length `n_{gp_j}` | `Y_{gp_j}(tpre_j) - Y_{gp_j}(period_1)` — comparison cohort pre-trend |
| `H` | scalar integer | Number of valid pairs for cell (g, t) |
| `Omega` | H x H matrix | Moment covariance matrix |
| `w` | vector length H | Efficient weights, sum to 1 |

---

## 2. File 1: `R/edid-pairs.R` — Replace PT-All Loop

### 2.1 What to Replace

Replace the entire PT-All section (lines 50–76) — everything from the comment `# PT-All: loop over all candidate comparison cohorts` through the final `data.frame(...)` return — with the new logic below.

Do NOT touch lines 38–48 (the empty sentinel and the PT-Post branch). Those are correct and must be preserved exactly.

### 2.2 New PT-All Logic

```r
  # -------------------------------------------------------------------------
  # PT-All: loop over treated cohorts only (never-treated is NOT a comparison
  # cohort; it appears only as the time control E[Y_inf(t)-Y_inf(tpre)] inside
  # each moment).
  #
  # For g' == target_g: valid tpre = {s : s < eff_start(g')}
  #   -- INCLUDES period_1 (degenerate CS DiD moment; comparison EIF = 0)
  # For g' != target_g: valid tpre = {s : period_1 < s < eff_start(g')}
  #   -- EXCLUDES period_1 (non-degenerate moments only)
  # -------------------------------------------------------------------------
  out_gp   <- numeric(0L)
  out_tpre <- numeric(0L)

  for (gp in treatment_groups) {
    eff_start <- gp - anticipation
    if (gp == target_g) {
      # Self-pair: include period_1
      valid_tpre <- time_periods[time_periods < eff_start]
    } else {
      # Cross-pair: exclude period_1
      valid_tpre <- time_periods[
        time_periods > period_1 & time_periods < eff_start
      ]
    }
    if (length(valid_tpre) > 0L) {
      out_gp   <- c(out_gp,   rep(gp, length(valid_tpre)))
      out_tpre <- c(out_tpre, valid_tpre)
    }
  }

  if (length(out_gp) == 0L) return(empty)
  data.frame(gp = out_gp, tpre = out_tpre, stringsAsFactors = FALSE)
```

### 2.3 Roxygen `@details` Update

Replace the current `@details` paragraph (lines 13–17 in the existing file, contained within the roxygen block) with:

```
#' Under \strong{PT-Post}: returns exactly one pair \code{(Inf, target_g - 1 - anticipation)},
#' or a 0-row data.frame if that pre-period does not exist in \code{time_periods} or
#' equals \code{period_1}.
#'
#' Under \strong{PT-All}: iterates over treated cohorts \code{g'} only (the
#' never-treated group is the time control inside every moment, not a comparison
#' cohort). For \code{g' == target_g}: valid \code{tpre} are all periods strictly
#' less than \code{g' - anticipation}, including \code{period_1} (this is the
#' degenerate CS DiD moment whose comparison-cohort EIF term is identically zero).
#' For \code{g' != target_g}: valid \code{tpre} are periods strictly between
#' \code{period_1} and \code{g' - anticipation} (exclusive on both ends).
#' Returns a 0-row data.frame if no valid pairs exist (e.g., single cohort with
#' only one pre-period equal to \code{period_1}).
```

The `@param` lines, `@return` line, and `@keywords internal` line are unchanged.

---

## 3. File 2: `R/edid-nocov.R` — Remove Dead gp=Inf Branches (PT-All Only)

With the new pair rule, all `gp` values in PT-All pairs are finite. Three `else` branches that handle `gp = Inf` in PT-All context are now dead code. Remove them as specified below. **Do NOT touch the PT-Post branch (lines 47–60) in `compute_omega_star_nocov_edid` or the PT-Post branch (lines 221–227) in `compute_generated_outcomes_nocov_edid` or the PT-Post branch (lines 287–300) in `compute_eif_nocov_edid`.**

### 3.1 `compute_omega_star_nocov_edid` — Remove gp=Inf Branch in `delta_gp_cache` Loop

Current code (lines 82–95):

```r
  for (rr in seq_len(nrow(unique_gp_tpre))) {
    gp_val  <- unique_gp_tpre$gp[rr]
    tp_val  <- unique_gp_tpre$tpre[rr]
    key     <- paste0(gp_val, "_", tp_val)
    if (is.finite(gp_val)) {
      mask_gp <- panel_obj$cohort_masks[[as.character(gp_val)]]
      col_pre <- .col(panel_obj, tp_val)
      delta_gp_cache[[key]] <- ow[mask_gp, col_pre] - ow[mask_gp, col_1]
    } else {
      # Never-treated: delta_{inf, tpre, 1} = Y_{tpre} - Y_{period_1}
      col_pre <- .col(panel_obj, tp_val)
      delta_gp_cache[[key]] <- ow[mask_inf, col_pre] - ow[mask_inf, col_1]
    }
  }
```

Replace with (remove the `if`/`else` wrapper — only the finite branch remains):

```r
  for (rr in seq_len(nrow(unique_gp_tpre))) {
    gp_val  <- unique_gp_tpre$gp[rr]
    tp_val  <- unique_gp_tpre$tpre[rr]
    key     <- paste0(gp_val, "_", tp_val)
    mask_gp <- panel_obj$cohort_masks[[as.character(gp_val)]]
    col_pre <- .col(panel_obj, tp_val)
    delta_gp_cache[[key]] <- ow[mask_gp, col_pre] - ow[mask_gp, col_1]
  }
```

Note: `gp_val` is always finite here because the new pair enumeration never includes `Inf` as a comparison cohort in PT-All. No `is.finite()` guard needed.

### 3.2 `compute_omega_star_nocov_edid` — Simplify Inner Loop `n_gp` Lookups

In the inner `for (j in seq_len(H))` and `for (k in seq_len(H))` loops, the condition number lookups use:

```r
    n_gp_j <- if (is.finite(gp_j)) {
      sum(panel_obj$cohort_masks[[as.character(gp_j)]])
    } else {
      n_inf
    }
```

Replace both occurrences (one for `gp_j`, one for `gp_k`) with the direct finite-only form:

```r
    n_gp_j <- sum(panel_obj$cohort_masks[[as.character(gp_j)]])
```

```r
    n_gp_k <- sum(panel_obj$cohort_masks[[as.character(gp_k)]])
```

These lookups are already guarded by the `if (gp_j == gp_k)` check for Term D; the simplification is safe because `gp_j` and `gp_k` are always finite.

### 3.3 `compute_omega_star_nocov_edid` — Remove gp=Inf `gp_j == gp_k` Guard Adjustment

The Term D check `if (gp_j == gp_k)` at line 147 uses plain `==` which works for both finite and Inf. It continues to work correctly for finite values only. No change needed here.

### 3.4 `compute_generated_outcomes_nocov_edid` — Remove gp=Inf Else Branch

Current code (lines 242–248):

```r
    if (is.finite(gp_j)) {
      mask_gp  <- panel_obj$cohort_masks[[as.character(gp_j)]]
      term_gp  <- mean(ow[mask_gp, col_pre] - ow[mask_gp, col_1])
    } else {
      # gp = Inf: comparison is never-treated baseline change
      term_gp <- mean(ow[mask_inf, col_pre] - ow[mask_inf, col_1])
    }
```

Replace with:

```r
    mask_gp  <- panel_obj$cohort_masks[[as.character(gp_j)]]
    term_gp  <- mean(ow[mask_gp, col_pre] - ow[mask_gp, col_1])
```

### 3.5 `compute_eif_nocov_edid` — Remove gp=Inf Else Branch

Current code (lines 333–351):

```r
      # Comparison cohort contribution (subtract)
      if (is.finite(gp_j)) {
        mask_gp          <- panel_obj$cohort_masks[[as.character(gp_j)]]
        pi_gp            <- panel_obj$cohort_fractions[[as.character(gp_j)]]
        delta_gp_pre_base <- ow[mask_gp, col_pre] - ow[mask_gp, col_base]
        mean_gp_pre_base  <- mean(delta_gp_pre_base)
        eif[mask_gp] <- eif[mask_gp] -
          w_j * (delta_gp_pre_base - mean_gp_pre_base) / pi_gp
      } else {
        # gp_j == Inf: the comparison cohort is the never-treated group itself.
        # The moment has a baseline shift component E[Y_inf(tpre) - Y_inf(1)]
        # whose EIF contribution is:
        #   -(1/pi_inf) * demeaned(Y_inf(tpre) - Y_inf(1))
        # This is separate from the -delta_inf term (Y_inf(t) - Y_inf(tpre))
        # already accumulated above and must be added explicitly.
        delta_inf_pre_base <- ow[mask_inf, col_pre] - ow[mask_inf, col_base]
        mean_inf_pre_base  <- mean(delta_inf_pre_base)
        eif[mask_inf] <- eif[mask_inf] -
          w_j * (delta_inf_pre_base - mean_inf_pre_base) / pi_inf
      }
```

Replace with:

```r
      # Comparison cohort contribution (subtract)
      mask_gp          <- panel_obj$cohort_masks[[as.character(gp_j)]]
      pi_gp            <- panel_obj$cohort_fractions[[as.character(gp_j)]]
      delta_gp_pre_base <- ow[mask_gp, col_pre] - ow[mask_gp, col_base]
      mean_gp_pre_base  <- mean(delta_gp_pre_base)
      eif[mask_gp] <- eif[mask_gp] -
        w_j * (delta_gp_pre_base - mean_gp_pre_base) / pi_gp
```

**Important**: `col_base` is defined as `col_1` (the column index for `period_1`) at line 305 and is unchanged. The EIF formula `delta_gp_pre_base = Y_{gp}(tpre) - Y_{gp}(period_1)` remains correct.

---

## 4. Degenerate Pair Analysis — Omega Entries When `g' == target_g`, `tpre == period_1`

When `gp == target_g` and `tpre == period_1`:

- `delta_gp_j = Y_g(period_1) - Y_g(period_1) = 0` (zero vector)
- In Omega: Term C_j = `Cov(delta_g_t_1, 0) / n_g = 0`; Term D = `Cov(0, ...) = 0` or `Cov(..., 0) = 0`
- The EIF contribution from the comparison cohort term: `w_j * (0 - 0) / pi_g = 0`
- The generated outcome: `term_gp = mean(0) = 0`, so `y_hat_j = term_g - term_inf`

No special case handling is needed. The existing code in `compute_omega_star_nocov_edid` and `compute_eif_nocov_edid` handles the zero vector correctly via the general formula. The `cov_nn_edid(delta_gp_j, ...)` call on a zero vector returns 0 cleanly. No divide-by-zero risk.

---

## 5. Unchanged Paths — Explicit Confirmation

These must remain byte-for-byte identical after the edit:

1. **`pt_assumption == "post"` branch in `enumerate_valid_pairs_edid`** (lines 40–48): unchanged.
2. **`pt_assumption == "post"` branch in `compute_omega_star_nocov_edid`** (lines 47–61): unchanged.
3. **`pt_assumption == "post"` branch in `compute_generated_outcomes_nocov_edid`** (lines 221–227): unchanged.
4. **`pt_assumption == "post"` branch in `compute_eif_nocov_edid`** (lines 287–300): unchanged.
5. **All of `R/edid-data.R`**: unchanged.
6. **All of `R/edid-utils.R`**: unchanged.
7. **All of `R/edid-linalg.R`**: unchanged.
8. **All of `R/edid-fit.R`**: unchanged.

---

## 6. Edge Case Handling

### 6.1 Single Treated Cohort

If `length(treatment_groups) == 1`, then `g' == target_g` always. The only valid pairs are `{(target_g, s) : s < target_g - anticipation}`. If `target_g - anticipation <= period_1`, no periods satisfy `s < target_g - anticipation`, so `pairs` is empty and the cell is set to NA (existing fallback in `edid-fit.R` line 73–96).

If `target_g - anticipation == period_1 + 1` (one period before effective start), then `valid_tpre = {period_1}` — exactly one pair, the degenerate CS DiD moment. This is correct behavior.

### 6.2 Anticipation > 0

`eff_start <- gp - anticipation`. All tpre conditions use `eff_start` as the upper bound. The self-pair condition `tpre < eff_start` means `tpre < gp - anticipation`. The cross-pair condition `period_1 < tpre < eff_start` means `period_1 < tpre < gp - anticipation`. This correctly restricts comparison cohorts to use only their clean pre-treatment periods.

### 6.3 `control_group = "last_cohort"`

`prepare_edid_panel()` in `edid-data.R` relabels the last finite cohort as `Inf` in `unit_cohorts` and removes it from `treatment_groups` before any estimation. Therefore `enumerate_valid_pairs_edid` receives a `treatment_groups` vector that already excludes the last cohort. No changes needed to `edid-pairs.R` for this case.

### 6.4 Empty Cross-Pairs (`g' != target_g` with no valid tpre)

If `gp != target_g` and `gp - anticipation <= period_1 + 1` (i.e., `eff_start - period_1 <= 1`), there are no valid cross-pair tpre values and the `length(valid_tpre) == 0` guard skips that `gp`. This is correct: cohorts that are treated too early to have interior pre-periods contribute no cross-pair moments.

---

## 7. Input Validation — No Changes Needed

`enumerate_valid_pairs_edid` performs no input validation beyond the `pt_assumption` branch dispatch and length checks. The existing interface contract (`treatment_groups` is a sorted numeric vector of finite values, `time_periods` is sorted, `period_1 == time_periods[1]`) is preserved.

---

## 8. Summary of All Line-Level Changes

### `R/edid-pairs.R`

| Location | Action |
|----------|--------|
| Lines 13–17 (roxygen `@details`) | Replace with new description (Section 2.3) |
| Lines 50–76 (PT-All block) | Replace entirely with new loop (Section 2.2) |

### `R/edid-nocov.R`

| Function | Location | Action |
|----------|----------|--------|
| `compute_omega_star_nocov_edid` | Lines 82–95 (`delta_gp_cache` loop) | Remove `if/else`; keep only finite branch (Section 3.1) |
| `compute_omega_star_nocov_edid` | Lines 103–107 (`n_gp_j` lookup) | Remove `if/else`; direct lookup (Section 3.2) |
| `compute_omega_star_nocov_edid` | Lines 119–123 (`n_gp_k` lookup) | Remove `if/else`; direct lookup (Section 3.2) |
| `compute_generated_outcomes_nocov_edid` | Lines 242–248 (`gp_j` branch) | Remove `if/else`; keep only finite branch (Section 3.4) |
| `compute_eif_nocov_edid` | Lines 333–351 (comparison cohort branch) | Remove `if/else`; keep only finite branch (Section 3.5) |
