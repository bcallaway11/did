# Weighted-covariate path rollout — audit trail

Target: enable a correct, fully-audited observation-weighted (`weightsname`) covariate
(`xformla`) path in `edid()`. Plan: `Efficient_DiD_Claude/quality_reports/plans/crystalline-herding-rose.md`.
Branch: `overnight-aggfix-robust`.

Standing gate after every phase (must stay green):
```
Rscript quality_reports/wcov/invariant_harness.R check   # -> ALL BYTE-IDENTICAL (PASS), 24 configs
```

---

## Phase 1 — Weighted nuisance WLS  ✅ COMPLETE (2026-06-16)

**Scope.** Gave the six sieve / Riesz nuisance fitters an obs-weight (`weights=` / `obsw=`)
slot, threaded `panel_obj$unit_weights` through the three `estimate_all_*` aggregators.

**Files changed.** `R/edid-cov.R` only:
- `estimate_propensity_ratio_edid` — WLS Gram `B_gp'W_gp B_gp`, obs-weighted col-sums, weighted
  score (`w * B(G_gp r - G_g)`), weighted IC loss. `H_inv = n*pinv(BtB_gp)` unchanged (Gram carries W).
- `estimate_inverse_propensity_edid` — WLS Gram, obs-weighted `colSums(w*B)`, weighted score, IC loss.
- `estimate_conditional_mean_edid` — `solve_ols_edid(..., weights=w_gp)`, weighted score/Hessian,
  weighted-mean degenerate fallback, weighted IC RSS.
- `fit_exp_riesz_edid` — `obsw` arg; weighted comparison indicator `cw = obsw*comp` folds into the
  loss / grad / Hessian (FOC `E_n[obsw·psi·comp·e^{eta}] = tcol`). `live` (support) stays weight-invariant.
- `exp_riesz_warmstart_edid` — `obsw` arg; weighted seed loss + WLS LS candidate.
- `estimate_propensity_ratio_exp_edid`, `estimate_inverse_propensity_exp_edid` — obs-weighted `tcol`,
  weighted `const_level`, pass `obsw` to solver + warm start, weighted M-estimator aux
  (`ow_s * score`, `crossprod(B, ow_s·diag·B)`; ridge-rescue `pen` terms stay unweighted), weighted IC loss.
- 3 aggregators (`estimate_all_propensity_ratios`, `estimate_all_inverse_propensities`,
  `estimate_all_conditional_means`) — extract `w_vec <- panel_obj$unit_weights`, pass
  `weights = w_vec[train_idx]` to each per-fold fit. **The 4 call sites in edid-fit.R / edid-boot.R
  need no edit** — they already pass `panel_obj`, which carries the weights.

**NULL-guard discipline.** Every weighted line is an explicit `if (is.null(weights))` branch whose
NULL arm is the *verbatim* original expression ⇒ byte-identical by construction. Where a scalar
multiplier was cleaner (`ow_s <- if (is.null(weights)) 1 else weights`), `1*x == x` in IEEE754 keeps
the unweighted result bit-identical.

**Gates (all green).**
1. New unit test `tests/testthat/test-edid-weighted-cov-wls.R` (6 tests, all pass):
   - WLS(w≡1) == OLS **byte-identical incl. aux** (`pred`/`s_hat`, `beta`, `score_mat`, `H_inv`) for
     all six fitters.
   - WLS(w≡c), c≠1 — fitted nuisance **prediction-invariant** (loss scales by c ⇒ same minimizer).
   - Dispersed mean-1 weight genuinely **moves** the conditional-mean fit (weights not inert).
2. Byte-identity invariant harness: **24 configs, worst |diff| = 0.000e+00, ALL BYTE-IDENTICAL (PASS)**.
3. Unweighted cov-path regression files (test-edid-cov-basic/eif/formula/ridge/validation/variance,
   ach-correction, ptpost-cov): **116 pass / 0 fail** (warnings pre-existing — deliberate edge cases).

**Docs.** `@param weights` added to the three Rd-generating fitters; `@param obsw` to the exp solver
+ warm start; exp wrappers inherit via `@inheritParams`. NEWS + the comparison note are batched into
Phase 6 (consolidated docs), per plan.

**Invariants still holding.** No-cov path and unweighted-cov path byte-identical (gate 2). The
end-to-end weighted-cov path is still blocked by the `edid.R:732` guard (removed in Phase 4), so the
new weighted branches are reachable only via direct fitter calls (gate 1) until then.

**Open follow-through into later phases.** The weighted aux (`score_mat`, `H_inv`) produced here is
consumed by the ACH correction — Phase 3 must verify the weighted estimation-effect against the FD
oracle. Phase 2 (weighted Omega*(X)) consumes the weighted nuisance *predictions* only.

---

## Phase 2 — Weighted Omega*(X) (3 smoothers)  ✅ COMPLETE (2026-06-16)

**Implemented** the convention below across all three builders. **Files:** `R/edid-cov-kernfast.R`,
`R/edid-cov-sieve.R`, `R/edid-cov-eif.R` (`compute_omega_star_cov_edid`).
- Each builder gains `uw <- panel_obj$unit_weights`, `.nrm <- if NULL n else sum(uw)`, and a pooling
  helper `wmean_o` (= `mean` when NULL).
- **Conditional moment**: kernel builders fold `w` into the kernel columns inside `get_kp`
  (`Kg <- Kg * rep(uw[idx], each = n)`; `Ks = rowSums` of the weighted `Kg`) → `cmean`/`ccov`/
  `kernel_cond_cov_kp`/`term_psi` all inherit the weighted NW. Sieve folds `w` into the WLS group Gram
  (`.sieve_group_pieces`: `crossprod(B_grp, w_grp*B_grp)`, stores `w_grp`) and the projection targets
  in `cmean`/`ccov` (`B'W v`, `B'W (AB)`).
- **Pooling**: per-unit builds `mean(o) → wmean_o(o)`; averaged builds weight the `avg_block` (row
  weight `uw*prefac`, normalize `/.nrm`; kernel keeps the column weight in `Kg`, sieve carries it via
  `w_grp*h` and `BtB_inv = (B'WB)^{-1}`), and `wmean_o` on the term1/term3-4 pooled blocks.
- **`pi_inf`** → `sum(uw[mask_inf])/.nrm` (NULL-guard). `pi_g` untouched (already weighted via
  `cohort_fractions`). Shrinkage-λ internals (`shape_var`/`samp_var`/`m_eff`) left unweighted
  (vanishing regularizers; the leading Omega-bar pooling IS weighted).
- The kernel-orig psi/EE channel inherits the weighted `Kg` (correct-direction for Phase 3); it is
  unreachable end-to-end until the Phase-4 guard drop, and stays byte-identical at `uw=NULL`.

**Gates (all green).**
1. New `tests/testthat/test-edid-weighted-cov-omega.R` (6 tests, all pass): weighted branch `uw≡1` ==
   unweighted **byte-identical** for all 3 builders × {averaged, per-unit}; a dispersed mean-1 weight
   **moves** Omega (not inert); a constant weight column normalizes to all-ones → Omega == unweighted
   (end-to-end panel build). This is STRONGER than the plan's byte-identity-only Phase-2 gate.
2. Byte-identity invariant harness (uw=NULL across kernel/kernel_orig/sieve × exp/direct):
   **24 configs, worst |diff| = 0.000e+00, ALL BYTE-IDENTICAL (PASS)**.
3. Broader cov-path regression files (test-edid-cov-*, ach-correction, ptpost-cov, weighted-cov-wls,
   weighted-cov-omega): **145 pass / 0 fail**.

**Correctness boundary (per plan).** Weighted-Omega *correctness* (not just w≡1 collapse) is validated
end-to-end by Phase 3's weighted FD oracle + Phase 5's MC calibration. The byte-identity gate + the
w≡1/dispersed/normalization tests confirm the unweighted path is untouched and the weighted code is
live and non-inert.

---

### (superseded) Phase 2 design notes

**Convention (decided, principled; correctness to be confirmed by the Phase-3 FD oracle + Phase-5 MC).**
Omega*(X) is the *conditional* covariance of the per-unit moment given X (a DGP feature scaled by 1/p
prefactors), NOT a variance-of-a-mean — so obs weights enter **linearly** here (the design-Bessel
`Σw²/(Σw)²` belongs to the EE/SE step, Phase 3), in two independent places:
  1. **Weighted Nadaraya-Watson / WLS conditional moment** — fold `w` into the smoother:
     kernel `K_iℓ → w_ℓ K_iℓ`; sieve `B'B → B'WB`, `B'v → B'Wv`. The shift-stabilizing pre-centering
     by the (unweighted) group mean can stay — Cov_K is shift-invariant; weighting enters only via the
     weighted E_K/projection.
  2. **Weighted pooling over the marginal X** — `mean(omega_jk_i) → weighted.mean(omega_jk_i, uw)`.
  3. **`pi_inf`** — `sum(mask_inf)/n → sum(uw[mask_inf])/n` (NULL-guard). `pi_g` already weighted via
     `cohort_fractions` (edid-data.R:123) — DO NOT touch.
Shrinkage λ internals (`shape_var`, `samp_var`, `m_eff`) are vanishing regularizers — leave unweighted
(byte-identical when NULL); only the leading Omega-bar pooling is weighted.

**Exact edit points (verified line numbers).**
- **Kernel** `compute_omega_star_cov_edid` (R/edid-cov-eif.R): add `uw <- panel_obj$unit_weights` near
  top; `pi_inf` at L472; fold `w` into `Kg`/`Ks` inside `get_kp` (L502-514) → ALL downstream
  (`cmean_psi`, `.cov_psi`, `kernel_cond_cov_kp`, `term_psi`) inherit weighting; weighted pooling at
  L797 (`omega_jk <- mean(omega_jk_i)`).
- **Kernel-fast** `compute_omega_star_kernel_fast_edid` (R/edid-cov-kernfast.R): `pi_inf` L29; fold `w`
  into `Kg`/`Ks` at L49 and into the averaged block `wrow`/`cWp` (L134-137: `wrow = prefac/Ks`, needs
  the column weights too); per-unit average L124/164. psi channel is blocked here (deferred to the
  kernel builder), so no EE work in this file.
- **Sieve** `compute_omega_star_sieve_edid` (R/edid-cov-sieve.R): `pi_inf` L45; weight the group Gram
  in `.sieve_group_pieces` (`BtB <- crossprod(B_grp)` L22-24 → `crossprod(B_grp, w[idx]*B_grp)`) and the
  projection target (`crossprod(B_grp, vc[idx])` → `crossprod(B_grp, w[idx]*vc[idx])`, L69-70); weighted
  pooling L223; the averaged block `cB`/`h` (L230-233). **Sieve has an INLINE psi/EE channel
  (L101-190)** — its OLS-projection IF (`aAB = BtB_inv B'Ws`, residuals) must be weighted in the SAME
  round as the value (full-adaptation), and FD-oracled in Phase 3.
- Group-size guards (`length(idx) < 2`) stay raw-count (structural support, weight-invariant), like the
  exp-Riesz `live` mask.

**Gate (per approved plan): w=NULL byte-identical (1e-12) on all three smoothers** via the invariant
harness (it already covers kernel/kernel_orig/sieve × exp/direct). Weighted-Omega *correctness* is NOT
gated at Phase 2 — it is validated end-to-end by Phase 3's weighted FD oracle and Phase 5's MC.
## Phase 3 — Weighted plug-in + EE/ACH + FD oracle  ✅ 3a/3b COMPLETE (2026-06-16); 3c PENDING

**Discovery (plan understated this):** the covariate path did NOT obs-weight its plug-in moment at all
(weights reached only the cov-ridge `n_eff`). So Phase 3 has a foundational prerequisite — the
obs-weighted (Hajek) plug-in moment/EIF — that precedes the EE.

**3a — obs-weighted plug-in (Hajek).** `R/edid-fit.R`: efficient `att_gt <- weighted.mean(wY_i, uw)`
(L727), averaged `att_gt` uses obs-weighted column means (L813). `R/edid-cov-eif.R`
`compute_eif_cov_edid`: the EIF carries the per-unit `uw` in both the no-trim and trim branches
(`pi_g` already weighted; per-pair `att_j` obs-weighted). SE is unchanged — the eif carries `uw` and
`sum(w)=n`, so `sqrt(sum(eif^2)/n^2)` is the Hajek design variance (mirrors the no-cov convention).

**3b — obs-weighted ACH (estimation_effect channel).** `score_mat`/`H_inv` already carry `uw` (Phase 1);
the only change is folding `uw` into the moment derivative `Γ = (1/n)B'(uw·s)` in BOTH
`compute_ach_correction_analytic_cov_edid` (analytic) and the `edid_ach="fd"` oracle (`m0`/`Γ` via
`weighted.mean(·, uw)`), so the oracle stays a valid ground truth for the obs-weighted moment.

**Guard.** `R/edid.R:732` hard stop removed → weighted-covariate path runs end-to-end (required to FD-
oracle it). Comment documents the validation status; nothing is committed until Phase 5 is green.

**Gates (all green).**
1. New `tests/testthat/test-edid-weighted-cov-e2e.R` (8 assertions, all pass):
   - **(B) weighted FD oracle**: dispersed mean-1 weights, `misspec_robust=FALSE`, `estimation_effect=TRUE`
     — analytic ACH SE == FD-oracle SE **to 1e-5** (the EE channel is correct under weights), and the
     ACH moves the SE by >1e-3 (a real correction).
   - **(A) constant-weight collapse**: a constant weight column (→ `unit_weights≡1`) makes the entire
     weighted-cov fit (att + SE) == the unweighted-cov fit to 1e-8, for `estimation_effect ∈ {F,T}` AND
     under `misspec_robust=TRUE`.
   - **(D)** dispersed weights move the covariate-path estimate (not inert).
2. Byte-identity harness (uw=NULL): **24 configs, 0 diff, ALL BYTE-IDENTICAL** — the unweighted path is
   untouched by the Hajek/EIF/ACH edits (all NULL-guarded).
3. Broader cov+nocov+weighted regression (test-edid-cov-*, ach-correction, ptpost-cov, nocov*,
   weighted-cov-wls/omega/e2e): **324 pass / 0 fail**.

**3c — EMPIRICALLY VALID (conservative); design-Bessel precision refinement REMAINING (research-grade).**
The `misspec_robust` Sigma_Omega weight-estimation channel under DISPERSED weights: (a) constant-weight
collapse byte-identical; (b) Phase 2 weights the kernel `Kg` feeding `term_psi`; (c) **MC coverage ~nominal
(0.958–0.971)** under dispersed weights — VALID inference. The weighted SE is CONSERVATIVE (variance ratio
analytic/MC ≈ 2.76 wt vs 1.45 unw → the EXTRA factor ≈ 1.9 ≈ Kish `n/n_eff` for the test weights).

**Diagnosis (traced `term_psi`, edid-cov-eif.R:654–676):** `psi_omega[ℓ]` already carries `uw_ℓ` (via the
Phase-2 weighted `Kg`: `S0[ℓ] = Σ_i sc_i Kg[i,ℓ] = uw_ℓ Σ_i sc_i K_iℓ`). The marginal `Σ_i` over eval units
does NOT carry `uw_i`, but adding it is ~NEUTRAL here because the weights are ⊥ X (the local kernel average
of `uw_i` ≈ 1), so that is NOT the fix. The real over-statement is structural: the weight-estimation
variance is a DEGENERATE second-order U-statistic whose true magnitude is `O(1/n_eff)`, but the first-order
`psi_omega` fold gives a naive `Σuw²/n²` (first-order) variance — exactly the term the no-cov path corrects
with the design-Bessel `s2_g = Σw²/(Σw)²` and the cov path does not yet have.

**Weight-channel FD oracle BUILT + a genuine bug FIXED (2026-06-16, session 2).**
`quality_reports/wcov/psiomega_fd_oracle.R` compares the analytic `psi_omega[ℓ]` to a case-weight FD
(perturb unit ℓ's weight ONLY in Ω̂, φ + outer-Hajek frozen). It revealed the analytic psi channel was
MISSING the outer Hajek marginal weight `uw_i` over eval units (Phase 2 weighted the value-pooling
`wmean_o` but NOT the psi channel's `Σ_i`): weighted cor(analytic,FD) = **0.42**. Fix — fold `uw_i` into
(a) `term_psi`'s `sc` (kernel, edid-cov-eif.R) and (b) the inv-p `coupled_C`; mirror in the sieve psi
channel (edid-cov-sieve.R `add_term`: marginal `uw_i` on the `Σ_i` accumulator, the perturbing-unit
WLS weight `w_ℓ` on the coefficient-IF residuals, and the WLS projection target `B'Wv` in `rawfit` to
match the value cmean). Result: weighted cor **0.42 → 0.93** and slope ratio wt/unw **0.276 → 1.04** —
the per-unit weight-channel IF now matches the FD oracle, same as unweighted. **All NULL-guarded →
byte-identity 0 diff; 37 weighted-cov tests pass; sieve byte-identity 0 diff.** This was a real
correctness bug (matters for clustered / aggregated SEs, where per-unit IF shape — not just its norm —
drives the variance).

**Residual aggregate conservativeness** (adversarial-DGP MC: wT cell-SE ratio still ~1.65, ≈ unchanged
by the structural fix because the per-unit NORM was preserved). This is the misspec channel's robust-SE
nature: it is present UNWEIGHTED too (uT ratio 1.20) and is amplified on the adversarial DGP (extreme
propensity ratios → heavy trimming; the plug-in itself severely under-covers there, wF ratio 0.62). The
first-order `psi_omega` variance vs the degenerate-U truth gap.

**VERDICT — Phase 3c SOLVED (healthy-overlap MC `mc_healthy.R`, 800 reps, n=500, gentle propensity).**
With the marginal-`uw` fix, the weighted `misspec_robust` SE is **nominal** on healthy overlap, BOTH
smoothers (true ATT=1, cell (2,4); ratio = mean_SE/MC_SD): kernel unw (1.04, cover 0.966) | **kernel wt
(1.00, 0.955)** | sieve unw (1.05, 0.964) | **sieve wt (1.01, 0.956)**. The weighted SE covers at nominal,
as well as unweighted. So the ~1.65 adversarial-DGP conservativeness was the robust misspec channel's
intended behavior (present unweighted too; plug-in itself under-covers there) — NOT a weighting bug. No
degenerate-U design-Bessel surgery needed: with the marginal Hajek weight correct, the weighted channel
tracks the unweighted one (nominal on healthy designs, conservative-but-valid on adversarial).

**TOOLKIT verified under weighted-cov** (`test-edid-weighted-cov-toolkit.R`, 25/25): a constant weight
column reproduces the UNWEIGHTED toolkit output to machine precision for ALL five —
`edid_weights` 6e-17, `edid_sargan` 1e-13, `edid_hausman` 9e-14, `edid_frontier` 3e-14,
`edid_adaptive` 2e-14 — and each runs + returns finite output under dispersed weights. The toolkit
consumes the (now weight-correct) eif/att consistently; the prior round's Kish-`n_eff` eigen-ridge
hardening covers the design factor.

**Cell-Hessian (`higher_order`) under weights:** the option-matrix sweep runs `higher_order=TRUE` with
weights without error; a dedicated weighted FD oracle for the cell Hessian is the one remaining
nice-to-have (it shares the ACH machinery already FD-validated in 3b).

---

## Phase 4 — Unblock guard + disclosure  ✅ guard done; disclosure polish PENDING
Guard relaxed (see Phase 3). `print`/`summary` already render the weighted-covariate fit correctly and
the `$args` refit snapshot works (the toolkit refits run on weighted-cov fits — sweep below). A one-line
print note of the supported `weightsname × xformla` combo is a remaining nicety (folds into Phase 6).
Updated `test-edid-weightsname.R` (the old "covariate + weights => error" assertion → now asserts the
path RUNS and a constant weight column reproduces the unweighted covariate fit to 1e-8).

## Phase 5 — Audit battery  ⏳ IN PROGRESS
- **Option-matrix smoke sweep** (`quality_reports/wcov/option_matrix_sweep.R`): **28/28 weighted-cov fits
  OK, 7/7 toolkit OK — ALL CLEAN.** Covers {kernel, kernel_orig, sieve} × {direct, exp} × {plugin, ee,
  misspec} + {averaged, gmm} + bs_df="ic" + higher_order + multiplier bootstrap + clustering + {group,
  event_study, calendar, overall} aggregation, all WITH `weightsname`; toolkit = weights / sargan /
  hausman / frontier / adaptive / summary / print on weighted-cov fits. **Flag:** `sieve|exp|misspec`
  SE (1.08 vs plug-in 0.56) is inflated — the un-validated 3c Sigma_Omega channel under dispersed weights.
- **Full testthat**: **3457 pass / 0 fail** (warn 67, skip 55) — ALL PASS (was 3394/0; the new
  weighted-cov tests — wls/omega/e2e/toolkit — add coverage). One guard-assertion test updated
  (`test-edid-weightsname.R`); no unweighted numbers move (byte-identity harness + cov/nocov regression).
  Option-matrix sweep re-run after the 3c fix: still **28/28 fits + 7/7 toolkit, ALL CLEAN**.
- **MC calibration** (`mc_calibration.R`, 2000 reps, config: efficient, **misspec_robust=FALSE**,
  estimation_effect=TRUE; true ATT=1): cell (2,4) bias +0.099, mean_SE/MC_SD = 0.69, coverage 0.83;
  cell (3,4) bias +0.034, ratio 0.80, coverage 0.88. **Interpretation:** this config OMITS the Omega
  weight-estimation variance (that is exactly what `misspec_robust=TRUE` / the 3c Sigma_Omega channel
  supplies), so under-coverage is the EXPECTED plug-in-weights behavior, not by itself a weighting bug.
  Decisive check is `mc_diagnostic.R` (weighted vs unweighted × misspec F/T on one DGP): weighting is
  correct iff weighted ≈ unweighted within each misspec level.
- **MC diagnostic** (`mc_diagnostic.R`, 800 reps, kernel) — RESULT (bias / ratio mean_SE÷MC_SD / coverage):
  - `misspec=FALSE`: unw c24 (+0.071, 0.62, 0.805) ≈ wt c24 (+0.069, 0.70, 0.839); unw c34 (0.80, 0.882)
    ≈ wt c34 (0.83, 0.875). **Bias identical wt vs unw; weighted tracks (even slightly beats) the
    unweighted under-coverage** → the FALSE under-coverage is the plug-in-weights property, NOT a
    weighting bug. **Phase 3a/3b weighting is calibration-correct.**
  - `misspec=TRUE`: unw c24 (1.20, **0.964**), wt c24 (1.66, **0.971**); unw c34 (1.06, 0.954), wt c34
    (1.42, **0.958**). **The Sigma_Omega channel restores ~nominal coverage under weights too.** The
    weighted SE is slightly CONSERVATIVE (ratio 1.4–1.7 vs unw 1.1–1.2) — valid inference in the safe
    direction; tightening it toward nominal is exactly the 3c design-Bessel refinement (a precision
    improvement, not a correctness fix).

## Phase 6 — Docs (NEWS + comparison note)  ⏳ PENDING
roxygen `@param weights`/`@param obsw` already added (Phases 1-2). Remaining: NEWS entry + a print/summary
note of the supported combo + this audit doc as the comparison record. No before/after table needed on
existing paths (proven byte-identical); the weighted-cov SE is genuinely new (its "before" was an error).

---

## Phase 3d — FULL weight-propagation audit ("weights must propagate fully into every option")  ✅ COMPLETE (2026-06-16, session 2)

Triggered by the user mandate to ensure observation weights propagate through EVERY option. Audited every
unweighted aggregation over per-unit quantities in the cov path. Two further gaps found + fixed (both
NULL-guarded => byte-identical unweighted):

1. **`m_common` (overlap treated-mass) was UNWEIGHTED** while `pi_g` (cohort_fractions) is weighted, so the
   overlap renorm `renorm_fac = pi_g / m_common = mean(uw | g) != 1` even with NO trimming. Fix
   (`edid_cell_trim_structure`, edid-cov-eif.R): `m_common = E_n[uw * G_g * keep]` (weighted Hajek mass);
   the dead-pair mass check likewise weighted. Now `m_common == pi_g` under no-trim => `renorm == 1`,
   consistent with the EIF/Hessian `m_kept = pi_g`. This was the residual that made the higher-order cell
   Hessian 1.0585x off (= `mean(uw|g)`); after it, the weighted analytic cell Hessian matches a clean
   central-FD of the weighted att, the higher-order FD-oracle test passes, the att stays UNBIASED, and MC
   coverage stays nominal (kernel wt 0.958, sieve wt 0.958).

2. **`weight_scheme = "gmm"` ignored weights** in three places: the gmm weight inverted an UNWEIGHTED
   `cov(gen_out_mat)`; the gmm data-channel `Cmat`/`mbar`/`psi_plug` were unweighted (and with the weighted
   gmm weight the `w_chk` gate would have silently SKIPPED the channel); the gmm ACH correction
   (`compute_gmm_weight_correction_cov_edid`) used unweighted `colMeans`/`mean`. Fix: weighted covariance
   (`cov.wt`), weighted `mbar`, obs-weighted `psi_plug`, weighted correction means.

**Already-correct (verified):** the aggregation cohort shares (edid-mp.R passes the per-unit `.w` to
`compute.aggte()` => ES/overall/group/calendar shares are obs-weighted); the "averaged" psi channel (shares
the now-weighted `term_psi`); clustering (the cluster-robust SE = `rowsum(eif, cluster)` of the uw-carrying
EIF, and `sigma_quad`'s clustered `V` uses the uw-carrying scores -- so clustering inherits weighting
through the EIF).

**Gates (all green).**
- **Option-matrix compatibility** (`test-edid-weighted-cov-options.R`, 88/0): a CONSTANT weight column ==
  the unweighted fit (att + SE + aggregate, 1e-7) across the FULL matrix -- kernel/kernel_orig/sieve x
  direct/exp x efficient/averaged/gmm/uniform x plugin/ee/misspec/higher_order x shrink{none,LW,ridge} x
  bs_df="ic" x PT-Post x trim x {group,event_study,calendar,overall} x clustering x multiplier bootstrap.
  => weighting is fully compatible with every option.
- `test-edid-weighted-cov-gmm.R` (6/0): gmm/averaged constant-weight == unweighted; dispersed runs finite.
- Byte-identity harness 0 diff; FULL testthat green; healthy-overlap MC nominal coverage after both fixes.
