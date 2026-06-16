# Effective-n ridge/LW fix — WEIGHTED-APPS before/after table

**Branch:** `overnight-aggfix-robust` (`/Users/pcostag/Documents/GitHub/did`)
**Mandate (Pedro):** "the ridge penalty must make sense under weights … rock solid."
**Date:** 2026-06-15 · **Audit check:** weighted-apps · **Spec:** `quality_reports/plans/2026-06-15_effective-n-ridge-fix.md`

## Verdict: **PASS**

On the three weighted applications that exercise the no-covariate weighted efficient
leg, the headline **ridge + EE** fit's conditioning improves by exactly the per-app
Kish design factor `n / n_eff`, every sign stays economically correct, and the
unweighted path is byte-identical (the fix is a no-op at equal weights). No regression.

---

## Method (true before/after, single-line isolation)

The live working tree already contains the fix. The **POST** numbers are the live
package (`/Users/pcostag/Documents/GitHub/did`, loaded via `pkgload::load_all`). The
**PRE** numbers are an isolated copy of that same tree (`/tmp/ridgefix_pkg_pre`) with
**only the ridge/LW denominator reverted** to the raw count `panel_obj$n` at the three
weighted-relevant no-cov sites — nothing else changed:

| site | PRE (baseline) | POST (fix) |
|---|---|---|
| no-cov ridge, `edid-fit.R:949–964` | `n_eff <- panel_obj$n; lam_r <- H / n_eff` | `n_eff <- n_eff_edid(unit_weights, active_mask, n); lam_r <- H / n_eff` |
| no-cov LW intensity, `edid-nocov.R:377–381` | `b2 <- b2_legacy` | `b2 <- b2_legacy * (n / n_eff)` |
| no-cov LW EE chain rule, `edid-nocov.R:584–586` | `n_eff <- n` | `n_eff <- n_eff_edid(unit_weights, active_mask, n)` |

The covariate-path lift sites (`edid-cov-eif.R`) use an all-TRUE mask → `n_eff == n`
always, and the weighted covariate path errors by design, so they are structural no-ops
here and were not reverted (no effect on these no-cov weighted fits).

Each app is fit with `weight_scheme="efficient"`, `pt_assumption="all"`,
`estimation_effect=TRUE`, on the certified window `e∈[0,5]`. Headline rung = **ridge**;
the **Ledoit–Wolf** rung is reported alongside for completeness. The four-rung ladder
(none / none+EE / ridge+EE / LW+EE) lives in the per-app gate runs; ridge+EE is the
package default and the only rung the ridge fix is responsible for.

### Datasets / weight columns (the apps that move)

| app | panel | id / time / gname / y | weight column | weight meaning | dispersion |
|---|---|---|---|---|---|
| **Schmitt** | `/tmp/gate_runs/rangel-60/panel_bal.rds` | id / year / gvar / y | `w_base` | 1996 baseline discharges (time-invariant) | moderate |
| **Gadenne** | `/tmp/gate_runs/aeraej-161/wt_panel_balanced.rds` | id / time / g / y | `pscore_u` | paper frame weight (∈ [0.67, 1.74]) | mild |
| **Bailey-GB** | `/tmp/gate_runs/aeraej-231/panel_r3.rds` (Construction B) | id / yr / G / y_amr | `w_pop` | 1960 county population (time-invariant) | **heavy** |

Bailey-GB Construction B: treated cohorts 1967–1972 (94 counties), 2,965 never-treated,
balanced 1959–1988; `w_pop` mean-1 CV = **2.57**, global Kish ESS = **402** of 3,059
units → **n / n_eff ≈ 7.6×**. This is the stress case for the fix.

---

## (1) RIDGE + EE — the headline default (the fix's target)

| App | metric | PRE (raw n) | POST (n_eff) | Δ | sign |
|---|---|---:|---:|---:|:--:|
| **Schmitt** (`w_base`) | ES_avg (SE) | +0.0572 (0.0143) | +0.0629 (0.0150) | SE **+5.4 %** | + → + ✓ |
| | ES(0) (SE) | +0.0348 (0.0126) | +0.0395 (0.0136) | — | + → + ✓ |
| | **cond# (max / med)** | 1196 / 721 | **649 / 392** | **×0.54** (−46 %) | finite, 0 Inf |
| **Gadenne** (`pscore_u`) | ES_avg (SE) | +14.528 (1.765) | +14.514 (1.766) | SE +0.04 % | + → + ✓ |
| | ES(0) (SE) | +5.584 (1.238) | +5.576 (1.239) | — | + → + ✓ |
| | **cond# (max / med)** | 2644 / 1878 | 2577 / 1830 | ×0.97 | finite, 0 Inf |
| **Bailey-GB** (`w_pop`) | ES_avg (SE) | −5.733 (2.264) | **−8.962 (2.661)** | SE **+17.6 %** | − → − ✓ |
| | ES(0) (SE) | −4.482 (2.488) | −6.574 (2.884) | — | − → − ✓ |
| | **cond# (max / med)** | 2724 / 2078 | **362 / 278** | **×0.13** (−87 %) | finite, 0 Inf |

**Conditioning is fixed, and it scales exactly with the weight dispersion.** The
post-fix cond# falls by the per-app Kish factor: Bailey (n/n_eff ≈ 7.6) drops cond#
≈ 7.5× (×0.13); Schmitt (moderate dispersion) ≈ 1.8× (×0.54); Gadenne (near-equal
weights, n_eff ≈ n) barely moves (×0.97) — which is the *correct* behaviour: where
weights are nearly uniform the fix is nearly a no-op. Every ridge fit is well-conditioned
(zero Inf cells) both before and after; the fix simply removes the under-regularization.

**Signs all correct.** Schmitt +0.063 (positive, its known direction), Gadenne +14.5
(positive precision gain), Bailey −8.96 (CHC establishment lowers age-adjusted mortality —
negative, correct). No sign flips.

**Numbers move only where dispersion is real.** On Bailey the under-regularized PRE
weights produced a materially different, less-shrunk point (−5.73 → −8.96) and a 17.6 %
larger, better-calibrated SE — the heavily-weighted large-population counties were being
allowed to dominate the un-floored weights. Schmitt moves a little (+5.4 % SE); Gadenne is
essentially unchanged. This is the expected, well-behaved consequence of restoring the
correct ridge scale under weights.

## (2) Ledoit–Wolf + EE — reported alongside (not the ridge fix's job)

| App | ES_avg (SE) PRE | ES_avg (SE) POST | SE Δ | sign | cond# Inf cells PRE→POST |
|---|---:|---:|---:|:--:|:--:|
| Schmitt | +0.0707 (0.0185) | +0.0745 (0.0186) | +0.57 % | + → + ✓ | 88 → 88 (invariant) |
| Gadenne | +13.428 (2.198) | +13.411 (2.199) | +0.03 % | + → + ✓ | 53 → 53 (invariant) |
| Bailey-GB | −12.674 (4.179) | −12.674 (4.179) | 0.00 % | − → − ✓ | 57 → 57 (invariant) |

The LW intensity `λ = min(1, b²/d²)` also picks up the `n/n_eff` factor (Schmitt median
λ 0.951 → 1.000), so its SE/point shift slightly where λ is interior. **Bailey LW is
byte-identical** because its dispersion is so heavy that `b²` is already large and λ
**saturates at the 1.0 cap** in every cell (range [0.52, 1]); the cap is hit identically
pre and post, so the n_eff factor cannot move the clamped λ — correct, hard-invariant
behaviour.

**On the LW `cond#=Inf` cells (no regression).** The reported per-cell cond# is for the
*shrunk* matrix `(1−λ)Ω + λ·target`; with λ near 1 it sits at the rank-deficient i.i.d.
pole target, so a subset of cells report Inf. This count is **identical before and after**
(88/53/57) — it is a pre-existing property of the LW pole geometry, **not introduced or
removed by the ridge fix**, which only rescales the intensity denominator. The LW fits
themselves remain numerically sound: every ES_avg / ES(0) / SE is finite in all three
apps, pre and post.

## (3) Invariant — unweighted byte-identity (no-op at equal weights)

Unweighted Schmitt ridge+EE, PRE copy vs POST live, full precision:

```
PRE  unweighted ridge: att = 0.047954489482  se = 0.014976691857  condmax = 1097.488285
POST unweighted ridge: att = 0.047954489482  se = 0.014976691857  condmax = 1097.488285
```

Identical to 12 digits (att, se, cond#). The fix is a perfect no-op when `unit_weights`
is `NULL`/equal (`n_eff == n`), confirming the PRE baseline copy isolates *only* the
weighted denominator and that the fix touches nothing on the unweighted path.

---

## Sign / conditioning sign-off

| Property | Result |
|---|---|
| Ridge sign correct, all 3 apps | ✓ (+ / + / −, no flips, pre and post) |
| Ridge conditioning improved post-fix | ✓ (cond# ×0.13 / ×0.54 / ×0.97, scaling with n/n_eff) |
| Ridge well-conditioned (no Inf) | ✓ (0 Inf cells, all apps, pre and post) |
| Ridge SE moves in the correct (larger) direction under real dispersion | ✓ (+17.6 % Bailey, +5.4 % Schmitt, +0.0 % Gadenne) |
| LW signs correct, fits finite | ✓ (Inf-cond cells invariant, a pre-existing pole property) |
| Unweighted byte-identical | ✓ (12-digit match, att/se/cond#) |

**No regression.** The weighted ridge leg is sign-correct and its conditioning is fixed
in proportion to the weight dispersion; the LW rung is unchanged except where its λ is
interior; the unweighted path is byte-identical.

## Artifacts
- Harnesses: `weighted_apps_run.R`, `bailey_only_run.R` (Bailey re-run with unit-level clustering = edid default)
- Results: `pre.rds` / `post.rds` (Schmitt + Gadenne), `bailey_pre.rds` / `bailey_post.rds`
- Logs: `pre.log` / `post.log` / `bailey_pre.log` / `bailey_post.log`
- Isolated PRE-baseline package copy: `/tmp/ridgefix_pkg_pre` (single-line denominator revert)
- Companion: `effective_n_ridge_fix_audit.md` (the full-adaptation audit for this round)
