#!/usr/bin/env Rscript
# Lambda-vanish check: the ridge intensity (H/n_eff) and the LW lambda must shrink
# as the sample grows (O(1/n_eff)), under BOTH unweighted and dispersed weights.
suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE)))

mk <- function(n_per_grp, seed = 11L, dispersion = NULL) {
  set.seed(seed)
  n_treat <- n_per_grp; n_never <- n_per_grp; n_periods <- 5L
  n <- n_treat + n_never
  uid <- rep(seq_len(n), each = n_periods); tid <- rep(seq_len(n_periods), times = n)
  ft  <- c(rep(3L, n_treat * n_periods), rep(Inf, n_never * n_periods))
  y <- rep(rnorm(n, 0, 1), each = n_periods) +
       rep(seq(0, 0.5, length.out = n_periods), times = n) + rnorm(n * n_periods, 0, 0.5) +
       2 * ((uid <= n_treat) & (tid >= 3L))
  df <- data.frame(unit = uid, time = tid, outcome = y, first_treat = ft)
  if (!is.null(dispersion)) {
    w_unit <- exp(rnorm(n, 0, dispersion))
    df$w <- w_unit[df$unit]
  }
  df
}

# Direct ridge intensity probe: lam_r = H / n_eff, mean(diag) scale.
probe_lambda <- function(df, weightsname = NULL) {
  po <- prepare_edid_panel(df, yname = "outcome", idname = "unit",
                           tname = "time", gname = "first_treat",
                           weightsname = weightsname)
  pr <- enumerate_valid_pairs_edid(3, po$treatment_groups, po$time_periods, po$period_1, "all")
  am <- active_mask_nocov_edid(3, pr, po)
  n_eff <- n_eff_edid(po$unit_weights, am, po$n)
  H <- nrow(pr)
  # LW lambda for a representative cell (g=3, t=last)
  t_last <- max(po$time_periods)
  om <- compute_omega_star_nocov_edid(3, t_last, pr, po, "all")
  sh <- shrink_omega_nocov_edid(om, 3, t_last, pr, po)
  list(n = po$n, n_eff = n_eff, H = H, ridge_int = H / n_eff, lw_lambda = sh$lambda)
}

cat("=== UNWEIGHTED: ridge intensity H/n_eff and LW lambda vs n ===\n")
for (m in c(20L, 40L, 80L, 160L, 320L)) {
  p <- probe_lambda(mk(m, seed = 11L))
  cat(sprintf("  n=%4d  n_eff=%7.1f  H=%d  ridge_int=H/n_eff=%.5f  LW_lambda=%.5f\n",
              p$n, p$n_eff, p$H, p$ridge_int, p$lw_lambda))
}
cat("\n=== WEIGHTED (dispersion=1.0): ridge intensity and LW lambda vs n ===\n")
for (m in c(20L, 40L, 80L, 160L, 320L)) {
  p <- probe_lambda(mk(m, seed = 11L, dispersion = 1.0), weightsname = "w")
  cat(sprintf("  n=%4d  n_eff=%7.1f  H=%d  ridge_int=H/n_eff=%.5f  LW_lambda=%.5f\n",
              p$n, p$n_eff, p$H, p$ridge_int, p$lw_lambda))
}
cat("\n(both columns must DECREASE toward 0 as n grows; n_eff grows ~proportionally to n)\n")
