# Validation: FULL testthat suite — over-id eigen-ridge fix

**Repo:** `/Users/pcostag/Documents/GitHub/did` · **branch:** `overnight-aggfix-robust`
**Date:** 2026-06-15 · **Task:** run the FULL testthat suite (all 64 files), compare to known-good **3977/0**, FAIL on any new failure.
**Load:** `pkgload::load_all(".")` (read-only, no source edits) · R 4.6.0 · testthat 3.3.2

---

## VERDICT: PASS

The full testthat suite is **byte-for-byte at the known-good baseline** under the standard harness:

| metric | known-good (ridgefix round) | this run (v3) | match |
|---|---|---|---|
| files | 64 | 64 | yes |
| **PASS** | **3977** | **3977** | yes |
| **FAIL** | **0** | **0** | yes |
| **ERROR** | 0 | **0** | yes |
| WARN | 89 | 89 | yes |
| SKIP | 11 | 11 | yes |

No new failure. No new error. WARN/SKIP counts identical to baseline. The over-id eigen-ridge fix
(`.edid_if_diff_quadform` weight-dispersion noise floor, `EDID_OVERID_DISP_C = 1.5`) introduces **zero
test regressions**.

The fix is present and live in this tree (`R/edid-hausman.R:225` `EDID_OVERID_DISP_C <- 1.5`;
`floor_rel`/`lift` logic at `R/edid-hausman.R:248-270`; threaded into `R/edid-sargan.R:289`).

---

## Fix-relevant CRAN-gated tests executed and PASSED (the trust question)

This run set `NOT_CRAN=true`, so the directly-fix-relevant `skip_on_cran()` tests RAN (they were
skipped in the no-env run, see "Harness reconciliation"):

- `test-edid-toolkit.R:44` — **"Hausman test has approximately correct size under PT-All AND power
  under violation"** → PASS. The test-enforced power property (the fix still REJECTS genuine PT
  violations) holds under the eigen-ridge lift.
- `test-edid-toolkit.R:614` — **"edid_sargan detects the violated moment under a PT-All violation"**
  → PASS. Sargan retains detection power.
- Plus the full thin-cohort MC (`test-edid-thin-cohort.R` n=1500 spillover/calibration), cov-variance
  coverage MC (`n=200, R=50`), JEL replication, mboot cluster-sum, adaptive-inference quadrature MC —
  all PASS.

(Note: the size-vs-power *calibration* of `c = 1.5` — the MC sweep that justifies the constant — is a
separate gate, not the testthat task. This run confirms the in-suite size/power assertion passes; it
does not re-derive the calibration table.)

---

## Harness reconciliation (why the first run looked broken — important)

Three runs were performed; only the third matches the known-good harness. The first two are recorded
to document the harness sensitivity, NOT as evidence against the fix.

| run | `export_all` | `NOT_CRAN` | PASS | FAIL | ERROR | WARN | SKIP | status |
|---|---|---|---|---|---|---|---|---|
| v1 | **FALSE** | unset | 2630 | 24 | 117 | 67 | 54 | **harness artifact — invalid** |
| v2 | TRUE | unset | 3394 | 0 | 0 | 67 | 55 | clean, but CRAN-gated tests skipped |
| **v3** | **TRUE** | **true** | **3977** | **0** | **0** | **89** | **11** | **MATCHES known-good — verdict** |

- **v1's 24 FAIL + 117 ERROR were 100% a harness artifact, not the fix.** With `export_all = FALSE`,
  internal (non-exported) functions are invisible to the test environment. Many edid tests call
  internal functions *unqualified* (e.g. `prepare_edid_panel(...)`, `.edid_aks_core(...)`,
  `validate_edid_inputs(...)`, `get_wide_data(...)` — none in NAMESPACE), so they errored with
  `could not find function "..."` / `threw an error with unexpected message`. Every captured v1
  failure/error traced to this single cause. The standard testthat harness uses `export_all = TRUE`
  (the `load_all` default; what `devtools::test` does), under which all of these resolve.
- **v2 vs v3 (3394 → 3977):** the gap is the 55 `skip_on_cran()` assertions. `NOT_CRAN` was empty in
  the shell, so v2 skipped all CRAN-gated tests (SKIP 55). v3 set `NOT_CRAN=true`, running 44 of them
  (the remaining 11 skips are genuine environment skips), which is exactly the baseline's SKIP 11 and
  lifts PASS to 3977. This reproduces the ridgefix verdict's own note: "Pass count rose 3394→3977
  because this round added/modified tests" — i.e. 3977 is the NOT_CRAN-on count, 3394 the off count.

---

## WARN / SKIP composition (all pre-existing, none from the fix)

**89 warnings** — two benign, expected diagnostic families, identical to baseline:
1. *Extreme propensity ratios (max > 100) ... thin PAIRWISE overlap ... trimmed* — the trim/keep-mask
   diagnostic on deliberately thin-overlap test designs.
2. *higher_order: could not recover the overall-aggregate weights ... the overall IF is not in the
   column span ... increment skipped* — the documented, expected `group`-aggregation message.

**11 skips** — genuine environment skips, none related to over-id:
- did v2.1.2 not on CRAN (7, `test-inference.R`); fork-unsafe macOS Accelerate BLAS (2,
  `test-edid-parallel.R`); known DRDID crash bug (1, `test-user_bug_fixes.R`); orthogonal-ACH sign
  undefined under uniform weights (1, `test-edid-ach-correction.R`).

---

## Artifacts

- `quality_reports/drafts/gate_runs/overid/full_testthat_v3.log` — the verdict run (NOT_CRAN=true), full skip/warn listing + summary.
- `quality_reports/drafts/gate_runs/overid/full_testthat_df_v3.rds` — per-test result frame (v3).
- `quality_reports/drafts/gate_runs/overid/full_testthat_v2.log` / `..._df_v2.rds` — export_all=TRUE, CRAN off (3394/0 clean).
- `quality_reports/drafts/gate_runs/overid/full_testthat_raw.log` / `..._df.rds` — v1 (export_all=FALSE) — harness-artifact run, retained for the diagnosis.

**Bottom line:** full testthat = **3977 PASS / 0 FAIL / 0 ERROR**, an exact match to known-good. No new
failure introduced by the over-id eigen-ridge fix.
