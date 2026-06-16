# Effective-n ridge fix — targeted MC calibration (event-study estimands)

**Audit item 4d.** Coverage of the 95% CI for the dynamic event-study average
(`ES_avg` = `fit$overall$overall.*`) and the on-impact effect (`ES(0)` =
`fit$event_study$att.egt/se.egt` at event time 0), under dispersed time-invariant
(per-unit) weights, comparing the ridge/LW intensity denominator BEFORE (raw n)
vs AFTER (per-cell Kish ESS `n_eff`).

- DGP: staggered, two treated cohorts (g3 effect 1.5, g5 effect 2.5), 60 never-treated, n=140, T=6, ridge shrinkage.
- Weights: log-normal `exp(N(0,s^2))`, calibrated to the per-cell Kish ESS targets.
- BEFORE = `n_eff_edid` overridden in the namespace to return the raw active count (reproduces pre-fix intensity); AFTER = the implemented fix. Same seeds for both arms.
- R = 500 reps per cell. Estimands anchored to the weighted plim per regime via a 60-fit large-sample average.
- Harness: `quality_reports/drafts/gate_runs/ridgefix/mc_es_coverage.R`; results `/tmp/gate_runs/ridgefix/mc_es_coverage.rds`.

## Realized dispersion

| Regime | s | target Kish | realized mean panel Kish | n/n_eff (g-cell, s) |
|---|---|---|---|---|
| Schmitt-like | 0.82 | ~54% | 0.536 | ~2.0 |
| Bailey-like  | 0.93 | ~48% | 0.457 | 2.544 (verified) |

n/n_eff = 2.544 under Bailey dispersion: the fix more than doubles the ridge
intensity, so this is a real test (not a no-op). Equal weights -> n_eff == n_active
(delta 0); NULL weights -> full n. Byte-identity preserved.

## Coverage table (500 reps)

### Schmitt-like (Kish ~54%)  — anchors ES_avg=1.7710, ES(0)=2.0175

| estimand | arm | MC SD | mean(SE) | mean(SE)/MC SD | 95% cover | bias |
|---|---|---|---|---|---|---|
| ES_avg | BEFORE (raw n) | 0.1404 | 0.1281 | 0.912 | 0.934 | -0.0234 |
| ES_avg | AFTER (n_eff)  | 0.1400 | 0.1285 | 0.918 | 0.932 | -0.0238 |
| ES(0)  | BEFORE (raw n) | 0.1644 | 0.1558 | 0.948 | 0.940 | -0.0237 |
| ES(0)  | AFTER (n_eff)  | 0.1634 | 0.1561 | 0.956 | 0.938 | -0.0241 |

### Bailey-like (Kish ~48%)  — anchors ES_avg=1.7703, ES(0)=2.0170

| estimand | arm | MC SD | mean(SE) | mean(SE)/MC SD | 95% cover | bias |
|---|---|---|---|---|---|---|
| ES_avg | BEFORE (raw n) | 0.1503 | 0.1336 | 0.889 | 0.924 | -0.0230 |
| ES_avg | AFTER (n_eff)  | 0.1497 | 0.1344 | 0.898 | 0.918 | -0.0237 |
| ES(0)  | BEFORE (raw n) | 0.1772 | 0.1637 | 0.924 | 0.928 | -0.0236 |
| ES(0)  | AFTER (n_eff)  | 0.1756 | 0.1644 | 0.936 | 0.930 | -0.0243 |

## Gate checks (tol = 2 x binomial SE on coverage = 0.019 at 500 reps)

| regime | estimand | dCover | dSE | dRatio | verdict |
|---|---|---|---|---|---|
| Schmitt | ES_avg | -0.002 | +0.00043 | +0.006 | OK |
| Schmitt | ES(0)  | -0.002 | +0.00035 | +0.008 | OK |
| Bailey  | ES_avg | -0.006 | +0.00079 | +0.009 | OK |
| Bailey  | ES(0)  | +0.002 | +0.00066 | +0.012 | OK |

## Verdict: PASS

- No coverage regression in any regime x estimand: all dCover within +/-0.006, far inside the +/-0.019 MC band (the small negative ES_avg moves are noise — mean(SE) rose, so coverage cannot truly drop from the fix).
- mean(SE) rises AFTER in EVERY cell (+0.0003 to +0.0008): the fix lifts the regularization in the correct (calibration-improving) direction; it never shrinks the reported SE.
- Gain retained and slightly improved: mean(SE)/MC-SD moves toward 1.0 in all four cells (+0.006 to +0.012).
- Effect is modest by construction (the ridge SE is the empirical variance of the realized weighted IF; the fix's primary role is correct weight regularization, which is a 2.5x change here yet leaves the IF — and thus the SE — only marginally moved). The residual ~0.92-0.94 coverage at n=140 under heavy dispersion is present BEFORE too (a plug-in-SE / finite-sample property, the job of the bootstrap / estimation-effect leg), not introduced by this fix.
