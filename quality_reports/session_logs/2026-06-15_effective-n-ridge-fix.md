# Session Log: 2026-06-15 - Effective-n ridge/LW fix

## Goal
Make the ridge / Ledoit-Wolf moment-covariance intensity "make sense under weights":
replace raw panel_obj$n with per-cell Kish ESS (n_eff) at the 3 intensity sites, with
byte-identical unweighted behavior. Full-adaptation & audit, no shortcuts.

## Key Decisions
- Weights field confirmed: panel_obj$unit_weights (NULL unweighted; else mean-1 per-unit).
- n_eff_edid() returns FULL panel_obj$n when unweighted (NOT n_act) -- matches the legacy
  denominator at every site, guaranteeing byte-identity in EVERY design (active-set Kish ESS
  only on the weighted branch). This is the mandate's own definition.
- LW b2: kept verbatim legacy expression * (n/n_eff) so the unweighted factor is EXACTLY 1.0
  (first draft used pi_hat/n_eff and lost one ULP -> fixed).
- LW EE chain rule kappa_i: 2/(n*d2) -> 2/(n_eff*d2), consistent with the new b2; FD-oracled
  under dispersed weights (interior lambda) to 9.4e-10.
- Cov-path lift + its EE tr(C)/n term wired to n_eff; weighted-cov path is scoped out (errors),
  so site 3 is a structural BYTE-IDENTICAL no-op today, correct-by-construction later.

## Work Done
- Sites: edid-fit.R (no-cov ridge), edid-nocov.R (LW b2 + kappa_i + roxygen),
  edid-cov-eif.R (2 lift helpers + EE tr term + 2 call sites), edid-cov-kernfast.R / sieve.R
  (2 call sites each). New helpers n_eff_edid + active_mask_nocov_edid in edid-utils.R.
- Audit: byte-identity max|delta|=0 (15 unweighted configs); weighted ridge/LW move (recorded);
  lambda-vanish (ridge H/n_eff 0.05->0.003 unw, 0.145->0.0068 wtd); weighted FD oracle 9.4e-10;
  full testthat 3977 PASS 0 FAIL; option-matrix sweep ALL PASS; MC before/after (cover 0.895->0.900).
- Docs: roxygen regenerated (n_eff_edid.Rd, active_mask_nocov_edid.Rd, 2 updated); NEWS bullet;
  audit doc quality_reports/drafts/gate_runs/ridgefix/effective_n_ridge_fix_audit.md.

## Open Questions
- Residual ~0.90 analytic-SE coverage at n=140 heavy dispersion (present before too; bootstrap/EE job).

## Next Steps
- None for this fix. Round complete; ready for review/commit per user.
