# Effective-n ridge/LW fix — AUDIT VERDICT + Schmitt re-certification

**Branch:** `overnight-aggfix-robust` (`/Users/pcostag/Documents/GitHub/did`, fix present as
uncommitted working-tree diff)
**Date:** 2026-06-15
**Mandate (Pedro):** "the ridge penalty must make sense under weights … rock solid."
**Verdict:** **AUDIT-CLEAN — ready to commit.** No FAIL across the full-adaptation battery.

---

## The fix (one line)

The vanishing ridge and Ledoit–Wolf intensities that stabilize the *weights* (not the SE)
divided by the raw unit count `panel_obj$n`; under dispersed observation weights the
scale-invariant weighted `Ω*` is backed by far fewer than `n` independent contributions, so
the raw count **under-regularized the weighted efficient leg by the Kish design factor
`n/n_eff ≥ 1`.** The denominator is now the per-cell Kish effective sample size
`n_eff = (Σw)²/Σw²` over active units (and **exactly the full `n` when unweighted →
byte-identical**), via the single helper `n_eff_edid()` (`R/edid-utils.R:256`). Three no-cov
sites + one cov-path lift, all wired through it.

---

## 2) AUDIT VERDICT — did the full-adaptation battery pass?

| Check | Result | Note |
|---|:--:|---|
| **Byte-identical unweighted** (no-X + cov; none/LW/ridge; no-X==trivial-X; PT-Post) | **PASS** | `max\|delta\| = 0.000e+00` on att/se/lambda/cond across all 15 unweighted configs, verified on LIVE patched code. Structural: `n_eff_edid()` early-returns full `n` when `w=NULL`. |
| **FD-oracle** (weighted ridge plain-map under n_eff; cov-path lift incl. trace) | **PASS** | no-cov plain-map D vs FD Jacobian rel 3.2e-10 (wtd) / 3.9e-10 (unwtd); cov-path lift 17/18 cells rel ~1e-10–1e-9 (the one 4.5e-8 is a central-difference truncation artifact, drops to 5.2e-9 at h=1e-5). LW chain rule 9.4e-10. |
| **testthat full suite** | **PASS** | 3977 PASS / **0 FAIL** / 0 ERROR / 89 pre-existing diagnostic WARN / 11 CRAN-skip (64 files). Pass count *rose* (3394→3977) because this round added/modified tests; none removed. |
| **Option-matrix smoke sweep** | **PASS** | **90/90** combinations pass, 0 bad, deterministic over 2 runs; all axes the rule enumerates covered; weighted leg confirmed LIVE (wtd-ridge att 1.540 ≠ unwtd 1.693) and unwtd leg unaffected; kernel vs kernel_orig max\|Δ\|=0. |
| **MC calibration** (95% CI coverage, ES_avg & ES(0), Schmitt~54% / Bailey~48% Kish, 500 reps) | **PASS** | All dCoverage within ±0.006 (inside ±0.019 MC band); mean(SE) **rises** in every cell (correct direction); mean(SE)/MC-SD moves toward 1.0. No regression. |
| **Weighted apps before/after** (Schmitt / Gadenne / Bailey-GB) | **PASS** | Ridge cond# improves by exactly the per-app Kish factor n/n_eff (Bailey ×0.13, Schmitt ×0.54, Gadenne ×0.97); all signs correct; 0 Inf cells; unweighted byte-identical to 12 digits. |

**FAIL list: NONE.** Two reported, non-blocking caveats (neither a regression, neither
introduced by this fix):
1. **Doc-wording precision (asymptotically-negligible).** The no-cov plain-map EE D does not
   differentiate the ridge's `mean(diag(Ω))` data-dependence — a *pre-existing*
   `O(H/n_eff)` ridge-trace channel that vanishes as λ→0 (identical in character unweighted,
   where `n_eff==n`). The cov path **does** include this channel (FD-oracled). Action: soften
   "the EXACT weight-estimation correction for w(Ω̂+cI)" to "exact up to the
   asymptotically-negligible ridge-trace term" in `effective_n_ridge_fix_audit.md`.
2. **Tooling:** `devtools` absent → `testthat::test_local` used (loads via `pkgload`, same 64
   files). Cosmetic.

---

## 1) SCHMITT RECOVERY — weighted (`w_base`, cluster `h`), FIXED n_eff ridge

Re-run on the LIVE patched tree (`schmitt_recert.R` → `schmitt_recert.{log,rds}`). Panel
`/tmp/gate_runs/rangel-60/panel_bal.rds`, 1293 hospitals, cohorts 2000–2010 + never, years
1996–2014, cluster = hospital `h` (1293 clusters).

**Schmitt n_eff (Kish ESS).** Global Kish ESS of the mean-1 `w_base` (1996 baseline
discharges) = **701.8 of 1293 units (n/n_eff = 1.843, CV(w) = 0.918)** — moderate dispersion.

### The gatekeeper: PLUG-IN (none) over-id — NEVER ridged

Contiguous joint Hausman over-id (U = PT-Post weighted, R = PT-All weighted), plug-in `none`,
grown from e=0. This is the honest window certifier and the n_eff fix correctly does **not**
touch it (it is on the un-ridged `none` path):

| e\[0,emax\] | df | chi² | p | decision |
|---|---|---|---|---|
| **e\[0,0\]** | 1 | 1.96 | **0.1613** | **PASS** |
| e\[0,1\] | 2 | 14.23 | 0.0008 | fail |
| e\[0,2\]…e\[0,14\] | … | 14.6→176 | ≤0.0022 | fail |

**→ PLUG-IN (none) CERTIFIED WINDOW (from e=0) = `e[0,0]`** (only the contemporaneous effect).

**Effective-n over-id-size caveat.** The over-id χ² over-rejects at small `n_eff` (its
quadratic form is backed by ESS, not `n`); at Schmitt's **moderate** n_eff≈702 the rejection
is unlikely to be pure size distortion — χ²=14.23 (df=2) at e\[0,1\] is far past nominal, and
the statistic grows monotonically (→176 by e=14), so the narrow window reflects a **genuine**
PT-Post-vs-PT-All weighted-moment mismatch past e=0, not a small-sample artifact. (Contrast:
the **unweighted** Schmitt certifies cleanly to e\[0,5\]; weighting toward high-volume
hospitals genuinely shortens the honest window.)

### Recovered weighted gain (certified window, EE-corrected, ARE = (CS-SE/eff-SE)²)

| Version | ES_avg (SE) | ARE_avg vs CS-never | ES(0) (SE) | ARE_0 vs CS-never |
|---|---:|---:|---:|---:|
| A none (plug-in) | 0.0135 (0.0123) | 1.77 (uncorrected) | 0.0684 (0.0149) | 1.00 |
| A+ none + EE | 0.0135 (0.0138) | 1.41 | 0.0684 (0.0156) | 0.90 |
| **C ridge + EE (default)** | **0.0629 (0.0150)** | **1.18** | **0.0395 (0.0136)** | **1.19** |
| B lw + EE | 0.0745 (0.0186) | 0.78 | 0.0593 (0.0176) | 0.72 |
| CS-never (anchor) | 0.0581 (0.0164) | 1.00 | 0.0438 (0.0149) | 1.00 |

Post-fix ridge ES_avg 0.0629 (SE 0.0150) and cond# 649/392 reproduce the audited weighted-apps
POST exactly. The fix **raises** the ridge SE (pre-fix ridge was ES_avg 0.0572, SE 0.0143,
ARE_avg 1.31) — i.e. it gives back the over-claimed gain (ARE 1.31→**1.18**), the honest
direction, while halving the conditioning (cond# ×0.54).

**Verdict (one line):** *Modest, well-conditioned weighted win — efficient ridge buys ≈1.18×
(ES_avg) to 1.19× (ES(0)) the CS estimator with a well-conditioned Ω (cond# ≈ 650, 0 Inf
cells) — but the honest over-id window is narrow (plug-in certifies only e\[0,0\]), so the
weighted Schmitt remains fragile relative to the cleanly-certified unweighted e\[0,5\].*

---

## 3) Bailey-GB — correct-signed diagnostic, no win claim

Bailey-GB (Construction B, `w_pop` 1960 county population, CV=2.57, global Kish ESS 402/3059,
n/n_eff ≈ 7.6×, the fix's stress case). Post-fix ridge+EE: **ES_avg −8.96 (SE 2.66), ES(0)
−6.57** — sign **negative and correct** (CHC establishment lowers age-adjusted mortality),
cond# falls ×0.13 (= the Kish factor, 2724→362), 0 Inf cells. The SE (2.66) is large relative
to the point, the over-id window is not a clean efficiency story, and the heavy dispersion is
exactly where under-regularization bit hardest — so this is reported as a **correct-signed
diagnostic, NOT a win claim.** The strongest demonstration that the fix engages: cond# drops
by precisely n/n_eff. (LW byte-identical pre→post because λ saturates at the 1.0 cap in every
cell — correct hard-invariant behaviour.)

---

## Bottom line

**The fix is AUDIT-CLEAN and ready to commit.** Byte-identical unweighted (structural,
`max|Δ|=0`), FD-oracled, full testthat 3977/0, option-matrix 90/90, MC no-regression,
weighted apps sign-correct with Kish-scaled conditioning. The only open item is a one-line
doc-wording softening (ridge-trace channel "exact" → "exact up to the asymptotically-
negligible ridge-trace term") — cosmetic, not a code change. **Schmitt:** modest
well-conditioned weighted win (ARE ≈ 1.18–1.19, cond# ≈ 650) on the plug-in-certified window
`e[0,0]`; honest window is narrow, so still-fragile relative to the unweighted e[0,5].
**Bailey-GB:** correct-signed diagnostic only.

*Artifacts: `schmitt_recert.{R,log,rds}`, `effective_n_ridge_fix_audit.md`,
`weighted_apps_before_after.md`, `mc_es_coverage_table.md`, `option_matrix_sweep_RESULT.md`,
`testthat_full.log`, plus the FD-oracle scripts in this directory.*
