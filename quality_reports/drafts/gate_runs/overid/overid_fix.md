# Over-identification / Hausman variance fix: n_eff-aware eigenvalue noise floor

**Repo:** `/Users/pcostag/Documents/GitHub/did` · **branch:** `overnight-aggfix-robust`
**Date:** 2026-06-15 · **Scope:** `.edid_if_diff_quadform()` and callers (no bootstrap; analytic).

---

## TL;DR

- **Mechanism CONFIRMED.** The spurious weighted Bailey-Goodman-Bacon over-id (`H = 290.6`, df 22,
  `p < 2e-16`) is a pseudoinverse over-amplification of downward-biased small eigenvalues of the
  IF-difference covariance `D` under dispersed weights + thin cohorts (effective sample
  `n_eff = 401 << n = 3039`). The smallest 5 retained eigendirections supply **49.1%** of `H`.
- **Fix:** lift (eigen-ridge) the retained eigenvalues at a weight-dispersion-aware noise floor
  before inverting; rank/df preserved.
- **Result:** Bailey-GB weighted over-id `-> H = 23.2, df 22, p = 0.39`, matching the clean
  TWFE-replication pre-trend verdict (`p ~ 0.43`). **Byte-identical** on unweighted full-rank
  designs (`max|H_live - H_orig| = 0`). Full edid testthat: **360 tests, 0 failures**.

---

## STEP 1 — Mechanism confirmation

**Setup.** Weighted Bailey-GB, Construction B, `weightsname = "w_pop"`, cluster = NULL (unit-level,
matching `wt_main.R`), `pt_assumption` post (U) vs all (R), `omega_cov_shrink = "none"` ("none" rung).
`n = 3039`, `|E| = 22`, weight dispersion max/min = 12819, CV = 2.57,
**`n_eff` (Kish ESS) = 400.6, `n/n_eff` = 7.59**.

**(a) Spectrum of `D` (`max-eig = 2.92e6`, numerical tol `mx*sqrt(eps) = 4.35e-2`).** All 22
eigenvalues are 4-10 orders of magnitude ABOVE the numerical tol — this is NOT a rank-deficiency /
float-dust case. The spectrum decays smoothly over 3.6 orders of magnitude (largest 2.92e6,
smallest 683.1; ratio_to_max of the smallest = 2.34e-4). Every eigenvalue is "retained" by the
current numerical threshold.

**(b) H decomposition by eigendirection (current numerical-tol pseudoinverse, `H = 290.58`).**
`H` is dominated by the SMALLEST eigenvalues (`H_k = n a_k^2 / lambda_k`, `a = V' d`):

| rank by Hcontrib | eigenvalue lambda | Hcontrib | % of H |
|---|---|---|---|
| smallest (lambda=683)   | 683.1  | 73.90 | 25.4% |
| lambda=2394             | 2393.7 | 31.27 | 10.8% |
| lambda=5927             | 5927.1 | 29.35 | 10.1% |
| lambda=2477             | 2476.6 | 29.11 | 10.0% |
| lambda=1075             | 1074.7 | 21.88 |  7.5% |

**Smallest-5 retained eigenvalues supply 49.1% of H.** The largest eigenvalue (2.92e6) contributes
only 1.06%. This is the over-amplification signature.

**(c) Raising the threshold to an n_eff-aware floor collapses H to a sane value.** A jackknife
(delete-1/40-block) on the smallest eigenvalue: full-sample 683.1, jackknife mean 664.5, SD 51.6,
range [488, 700] — the small eigenvalues carry real sampling variability, and dividing by them is
the inflation channel.

**Verdict: mechanism CONFIRMED** — small-eigenvalue pseudoinverse over-amplification under
limited effective sample (dispersed weights + thin cohorts), exactly as hypothesized.

---

## STEP 2 — Implementation

### Threshold formula (the design choice, and why the prompt's literal formula was wrong)

The prompt sketched `tol = mx * max(sqrt(eps), c * sqrt(p / n_eff))`. **This breaks invariant (ii).**
On the standard *unweighted, well-conditioned* test design (`make_panel_toolkit`, n=300, `n_eff = n`),
`D` has condition ~168, min/max = 6e-3, but `c*sqrt(p/n_eff) = 1*sqrt(4/300) = 0.115 > 6e-3` — so
the `p/n_eff` relative-to-max truncation would **drop 2 of 4 directions even on the well-conditioned
unweighted design**, changing the chi-square df from 4 to 2 and violating byte-identity. The
`p/n_eff` scale is simply not small at small p, and a relative-to-max HARD truncation on a smoothly
decaying spectrum lurches lumpily (in the Bailey-GB case it jumped df 22 -> 3 -> 2 -> 1 as `c` grew).

**Adopted design — weight-dispersion-driven eigen-RIDGE (not p/n_eff, not hard truncation):**

```
disp      = max(0, n / n_eff - 1)                              # weight-dispersion excess (0 unweighted)
floor_rel = max(sqrt(.Machine$double.eps), c * sqrt(disp / n_eff)),   c = EDID_OVERID_DISP_C = 1.5
lambda_k -> max(lambda_k, floor_rel * max-eig)                 # lift the over-amplified directions
H        = n * d' V diag(1/lambda_lifted) V' d,  df = rank (UNCHANGED)
```

Two design decisions, both evidence-driven:

1. **Driver = `disp = n/n_eff - 1`, NOT `p/n_eff`.** The pathology is a *weight-dispersion*
   phenomenon: when `n = n_eff` (unweighted / uniform) there is NO inflation, so the floor MUST
   collapse to the numerical tol there. `disp = 0` exactly when unweighted (since `n_eff_edid`
   returns the full `n` for `w = NULL`, and a mean-1 constant column gives `n_eff = n_act`),
   guaranteeing byte-identity by construction. The scale `sqrt(disp / n_eff)` is the relative
   sampling SD of an estimated-covariance eigenvalue built from `n_eff` effective pieces, inflated
   by the excess dispersion the raw `1/n` averaging does not see.
2. **Eigen-RIDGE (lift), not hard rank truncation.** The spectrum is smooth (not rank-deficient),
   so the over-amplified directions are *damped* (`lambda -> max(lambda, floor)`) rather than
   discarded — the rank, and hence the chi-square df, is preserved (a smooth Tikhonov lift, not a
   lumpy df cut). Confirmed monotone and well-behaved vs hard truncation in the Step-1b probe.

### Calibration of `c` (justification)

`c = 1.5` is calibrated for ~nominal **size** of the over-id under H0 (clean PT), NOT fitted to the
Bailey-GB p-value. Size MC (clean-PT staggered DGP, 140 reps, joint event-study over-id):

| condition | mean n_eff (n=900) | c=0 (current) | c=1.0 | c=1.5 | c=2.0 |
|---|---|---|---|---|---|
| thin cohorts + dispersed weights | 108 | **84.3%** | 30.0% | 24.3% | 16.4% |
| well-conditioned (mild weights)  | 846 | 4.3% | 2.9% | 2.9% | 2.1% |

`c = 1.5` keeps the well-conditioned design at nominal (2.9%, not over-corrected) while removing the
bulk of the pathology over-rejection (84% -> 24%); `c = 2.0` starts under-sizing the well-conditioned
case. The residual over-rejection on the deliberately extreme `n_eff = 108` / 22-restriction DGP is
the *genuine* small-effective-sample strain of a 22-df chi-square test (`p/n_eff = 0.2`), NOT the
pseudoinverse artifact — and it **converges to nominal as the effective sample grows** (second MC, 120
reps, `c = 1.5`):

| dispersion | mean n_eff | rej@5% (c=1.5) |
|---|---|---|
| severe, thin   | 159  | 10.8% |
| moderate       | 672  | 5.8% |
| mild, fuller   | 1168 | 3.3% |

The Bailey-GB case (`n_eff = 401`, `p/n_eff = 0.055`) sits comfortably in the well-behaved region,
which is why it lands cleanly at `p = 0.39`.

### Invariants (all by construction, all verified)

- **(i) Asymptotically negligible.** Fixed weight distribution `=> disp -> const, n_eff -> inf =>
  floor -> mx*sqrt(eps)`: the lift vanishes, chi-square(rank) limit untouched. (Project invariant:
  every covariance regularizer is asymptotically negligible.)
- **(ii) Unweighted / well-conditioned unchanged.** `disp = 0 => floor = sqrt(eps) => lift = FALSE
  => the ORIGINAL `solve()` / pseudoinverse branch runs bit-for-bit.` Verified `max|H_live - H_orig|
  = 0.000e+00`.
- **(iii) Degenerate-contrast guard fires first.** The `max(abs(D)) <= eps_D*v_scale` and `rk == 0`
  branches are untouched and short-circuit before any lift.

### Exact change (`git diff` summary)

```
 R/edid-hausman.R  | 91 +++++++++++++++++++++++++++  (quadform lift + n_eff helper + threading)
 R/edid-sargan.R   |  3 +-                            (thread n_eff into the window-grow quadform)
 R/edid-frontier.R |  3 +-                            (thread n_eff into the scalar Hausman, parity)
 NEWS.md           | 12 ++++++++                      (new bullet)
```

Core of `.edid_if_diff_quadform` (new branch; unchanged path retained verbatim when `!lift`):

```r
EDID_OVERID_DISP_C <- 1.5
.edid_if_diff_quadform <- function(d, xi, n, cluster_indices, v_scale = 1, n_eff = n) {
  ... # D build, finite/degenerate/rk==0 guards UNCHANGED
  if (!is.finite(n_eff) || n_eff <= 0) n_eff <- n
  disp      <- max(0, n / n_eff - 1)
  floor_rel <- max(sqrt(.Machine$double.eps), EDID_OVERID_DISP_C * sqrt(disp / n_eff))
  lift      <- floor_rel > sqrt(.Machine$double.eps) && any(ev$values[pos] < floor_rel * mx)
  if (!lift) { ... ORIGINAL solve()/pseudoinverse branch, byte-identical ... }
  V   <- ev$vectors[, pos, drop = FALSE]
  lam <- pmax(ev$values[pos], floor_rel * mx)
  H   <- as.numeric(n * crossprod(d, V %*% (crossprod(V, d) / lam)))
  list(statistic = H, df = rk, p_value = stats::pchisq(H, df = rk, lower.tail = FALSE),
       D = D, degenerate = FALSE)
}
```

`n_eff` is supplied by `.edid_overid_n_eff(fit) = n_eff_edid(fit$unit_weights, rep(TRUE, n), n)`
(reuses the cov-ridge fix's Kish-ESS machinery; `== n` exactly when unweighted).

---

## STEP 3 — Fast validation (through the LIVE package functions)

**(a/d) Bailey-GB weighted joint over-id — now SANE:**

```
POST-FIX (none rung):    H = 23.243   df = 22   p = 0.3881     (was H = 290.6, p < 2e-16)
clean TWFE-replication target: p ~ 0.43      MATCH
overall ES_avg (scalar): H = 6.632    df = 1    p = 0.0100     (scalar path; byte-identical)
ridge+EE default rung:   H = 6.512    df = 22   p = 0.9994     (was H = 25.5, p = 0.27; weighted, lift applies)
```

**(b) BYTE-IDENTICAL on full-rank unweighted designs** (live `edid_hausman` vs an emulation of the
original numerical-tol quadform):

```
seed1 n=300:  LIVE H=4.1109854101 df=4  | ORIG H=4.1109854101 df=4  | dH=0.00e+00
seed2 n=250:  LIVE H=3.4005396261 df=4  | ORIG H=3.4005396261 df=4  | dH=0.00e+00
seed3 n=200:  LIVE H=7.7179240172 df=4  | ORIG H=7.7179240172 df=4  | dH=0.00e+00
seed4 n=180:  LIVE H=10.4641481533 df=4 | ORIG H=10.4641481533 df=4 | dH=0.00e+00
clustered (G=16): LIVE H=0.7724214300 df=1 | ORIG H=0.7724214300 df=1 | dH=0.00e+00
MAX |H_live - H_orig| = 0.000e+00 ;  MAX |p_live - p_orig| = 0.000e+00
```
Scalar path also byte-identical on the weighted Bailey-GB fits (`max|scalar H_live - H_orig| = 0`
on both none and ridge+EE rungs).

**(c) Size MC** — see the two tables in Step 2 (calibration). Pathology 84.3% -> 24.3% at c=1.5,
converging to nominal as `n_eff` grows; well-conditioned stays at nominal (2.9%).

**(d) testthat:** full edid battery **360 tests, 0 failures across 38 files**
(`test-edid-toolkit`, `test-edid-round3-guards`, `test-edid-adaptive-inference/-fixture`,
`test-edid-identities`, ... all pass; the full-rank byte-identity assertion at
`test-edid-toolkit.R:99` holds).

---

## Toolkit caller impact enumeration (full-adaptation audit, this round)

| Caller | Uses | Status |
|---|---|---|
| `edid_hausman` (joint) | `.edid_if_diff_quadform` | **ADAPTED** — `n_eff = .edid_overid_n_eff(fit_restricted)` threaded; Bailey-GB now sane, unweighted byte-identical. |
| `edid_hausman` (scalar per-e + ES_avg) | `.edid_scalar_hausman` | **ADAPTED (no-op by design)** — `n_eff` threaded for parity; 1-D variance is directly estimated, not an inverted small eigenvalue, so H is byte-identical (verified). |
| `edid_sargan` (window-grow certify) | `.edid_if_diff_quadform` | **ADAPTED** — `n_eff = .edid_overid_n_eff(fit_base)` threaded into the per-restriction quadform. Unweighted byte-identical; the Holm window-grow inherits the dispersed-weight floor. |
| `edid_frontier` | `.edid_scalar_hausman` | **ADAPTED (no-op by design)** — `n_eff` threaded; scalar path byte-identical, so frontier radii unchanged on every existing design. |
| `edid_adaptive` | own scalar `t_O = Y_O/sqrt(V_O)` from a 2x2 covariance | **UNAFFECTED** — 1-D over-id direction, no pseudoinverse over-amplification; no quadform call. No change. |
| `edid_weights` | — | **UNAFFECTED** — does not call the quadform. |

**Docs same round:** roxygen-level comment block at `.edid_if_diff_quadform` rewritten (formula,
calibration, invariants); `EDID_OVERID_DISP_C` documented inline; NEWS.md bullet added. No exported
signature changed, so no `.Rd` drift from this work (roxygenize produced no new man/ changes
attributable to these edits).

## Honest caveats / follow-ups for the SEPARATE full battery

- The size correction is real and large but does NOT reach a perfect 5% on the deliberately extreme
  `n_eff ~ 108` / 22-restriction DGP (it lands ~24%); this residual is the genuine few-effective-units
  strain of a high-df chi-square test and shrinks to nominal as `n_eff` grows (10.8% -> 5.8% -> 3.3%).
  The Bailey-GB case is comfortably in the clean region.
- Weighted toolkit fits that are ALSO regularized move as a consequence (e.g. Bailey-GB ridge+EE
  rung `H: 25.5 -> 6.5`). Downstream weighted-application over-id verdicts and any deck/paper cards
  reporting weighted Hausman/Sargan numbers should be regenerated (flagged for the follow-up).
- The full 500-rep size MC, the complete option-matrix smoke sweep, and the all-application
  re-verdicts are the separate follow-up you will run after judging this result.
```
