# Option-matrix smoke sweep — effective-n ridge/LW fix (independent re-run)

**Branch:** `overnight-aggfix-robust` · **Date:** 2026-06-15
**Script:** `option_matrix_sweep_full.R` · **Log:** `/tmp/gate_runs/ridgefix/option_matrix_sweep_full.log`
**Result:** **ALL PASS — 90 combinations, 0 bad, process exit 0 (deterministic over 2 runs).**

## Axis coverage (every axis the audit rule enumerates)

| Axis | Values exercised | Sections |
|---|---|---|
| `ratio_method` | exp, direct | C, H |
| `weight_scheme` | efficient, averaged, gmm, uniform | A, B, F, H |
| omega smoother | **kernel, kernel_orig, sieve** (cov path, via `edid_omega_method`) | G, H, I |
| `omega_cov_shrink` (nocov_shrink) | **none, ledoit_wolf, ridge** | A, B, C, G |
| `weightsname` | **NULL (unweighted) and weighted col "w"** | A vs B; all |
| `misspec_robust` / `estimation_effect` / `higher_order` | TRUE/FALSE; ho on cov path (no-cov guard verified) | D, I |
| `clustervars` | one-way cluster col | E, I |
| both bootstraps | `bstrap=TRUE` (mult.) + `cband_method` {analytic, multiplier} | E, I |
| `aggregate` | event_study, group, calendar, overall, none | F |
| `moment_set` | NULL + a valid `(g,gp,tpre)` data.frame | F |
| `bs_df` | 4 (default) + `"ic"` (cov path) | I |
| `pt_assumption` | all, post | F, K |
| toolkit | **edid_weights, edid_sargan, edid_hausman, edid_frontier, edid_adaptive** | K (weighted), L (cov) |
| print / summary | both, on weighted ridge + cov ridge | K, L |
| `$args` refit snapshot | present-and-non-null check | K |
| weighted-cov contract | `weightsname`+`xformla` errors **by design** (scoped out) | J |

## Rigor checks (fail-loud)

- **No garbage:** 74 fit-producing combos → att ∈ [1.3687, 1.7817] (true cohort effects 1.5/2.5; overall ~1.4–1.8), se ∈ [0.1325, 1.1101], **all finite, all se>0**.
- **kernel == kernel_orig:** `max|Δatt| = 0`, `max|Δse| = 0` (smoother axis consistent through the fix).
- **Weighted leg live & distinct:** weighted-ridge-eff att 1.5400 vs unweighted-ridge-eff att 1.6928 (Δ=0.153) — n_eff branch active; unweighted leg unaffected.
- **Determinism:** identical summary line + exit 0 on two independent invocations.

## Triage of the 4 initial flags (all HARNESS bugs, NOT package regressions)

The first draft of the harness fed 4 invalid argument combinations; the package
responded with **correct documented behavior** in every case (clean errors / correct
output shape, never garbage), so none is a regression:

1. `higher_order=TRUE` on a no-cov fit → documented `stop()` (`R/edid.R:760`):
   no covariates ⇒ zero higher-order variance by construction. (×2 flags.)
   Harness fix: ho exercised on the cov path (section I); a no-cov `ho=FALSE`
   explicit pins the guard.
2. `aggregate="none"` returns only `att_gt` (no aggregate object) — harness picker
   yielded NA. Verified `att_gt` has 12 finite-att/se cells. Harness now checks `att_gt`.
3. `moment_set="own"` is invalid — `moment_set` must be a data.frame with columns
   `g, gp, tpre` (`R/edid.R:827`). Harness now passes a valid data.frame (att 1.5749).

After the harness fixes, the corrected sweep is **90/90 PASS, 0 bad**.

## Scope note

Per the fix's design, the **weighted-cov** path is scoped out and errors by design
(section J confirms), so the cov-path ridge lift (site 3) is a structural
byte-identical no-op today and is exercised on the **unweighted** cov path
(kernel/kernel_orig/sieve), where it must be invariant — confirmed.
