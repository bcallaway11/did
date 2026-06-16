#!/usr/bin/env Rscript
# FD-oracle of the WEIGHTED Ledoit-Wolf chain-rule EE term (n_eff != n).
# Verifies that compute_nocov_ee_correction_edid's analytic per-unit directions D
# (which now use kappa_i with 2/(n_eff*d2)) reproduce brute-force finite differences
# of the weight map w(Omega_sh(Omega)) when the LW b2 carries the n/n_eff correction.
suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE)))

# --- weighted panel with dispersed time-invariant weights ---------------------
mk <- function(n = 360L, seed = 4L, Tn = 6L, disp = 0.8, rho = 0.6) {
  set.seed(seed)
  G <- sample(c(3L, 5L, Inf), n, TRUE, c(.35, .30, .35))
  # heteroskedastic, serially-correlated shocks -> Omega genuinely OFF the i.i.d. pole
  e <- matrix(rnorm(n * Tn), n, Tn)
  for (s in 2:Tn) e[, s] <- rho * e[, s - 1L] + sqrt(1 - rho^2) * e[, s]
  scale_i <- exp(rnorm(n, 0, 0.5))                 # cross-unit variance heterogeneity
  e <- e * scale_i
  Y <- rnorm(n) + matrix(0.3 * seq_len(Tn), n, Tn, byrow = TRUE) + e
  for (g in c(3L, 5L)) for (t in g:Tn) Y[G == g, t] <- Y[G == g, t] + 1 + 0.3 * (t - g)
  w_unit <- exp(rnorm(n, 0, disp))
  data.frame(id = rep(seq_len(n), each = Tn), time = rep(seq_len(Tn), n),
             y = as.vector(t(Y)), g = rep(G, each = Tn),
             w = rep(w_unit, each = Tn))
}

df <- mk(seed = 33L, rho = 0.70)               # lands LW lambda ~ 0.475 (interior) under weights
pn <- prepare_edid_panel(df, "y", "id", "time", "g", weightsname = "w")
stopifnot(!is.null(pn$unit_weights))

gg <- 5L; tt <- 5L
prs <- enumerate_valid_pairs_edid(gg, pn$treatment_groups, pn$time_periods, pn$period_1, "all", 0L)
H   <- nrow(prs)
Om  <- compute_omega_star_nocov_edid(gg, tt, prs, pn, "all")
psi <- compute_psi_moments_nocov_edid(gg, tt, prs, pn)
q4  <- sum(rowSums(psi * psi)^2)
S   <- compute_pole_structure_nocov_edid(gg, tt, prs, pn)
n   <- pn$n
am  <- active_mask_nocov_edid(gg, prs, pn)
n_eff <- n_eff_edid(pn$unit_weights, am, n)
cat(sprintf("H=%d  n=%d  n_eff=%.3f  (n/n_eff = %.4f)\n", H, n, n_eff, n / n_eff))
stopifnot(n_eff < n - 1)   # genuinely dispersed: the correction is exercised

# Reference weight map with the WEIGHTED LW b2 (carries n/n_eff), S/n/n_eff/q4 fixed.
w_of_omega <- function(Omx) {
  ss <- sum(S * S); sigma2 <- sum(Omx * S) / ss; target <- sigma2 * S
  d2 <- sum((Omx - target)^2)
  b2_legacy <- (q4 / n^2 - n * sum(Omx * Omx)) / n^2
  b2  <- b2_legacy * (n / n_eff)
  lam <- min(1, max(0, b2) / d2)
  Os  <- (1 - lam) * Omx + lam * target
  A   <- solve(Os); u <- drop(A %*% rep(1, nrow(Os))); u / sum(u)
}

# analytic EE (returns per-unit directions D)
shr <- shrink_omega_nocov_edid(Om, gg, tt, prs, pn)
lam <- shr$lambda; Omu <- shr$omega
w   <- compute_efficient_weights_edid(Omu)
ee  <- compute_nocov_ee_correction_edid(gg, tt, prs, pn, omega_raw = Om, omega_used = Omu,
                                        weights = w, shrink_lambda = lam, return_D = TRUE)
cat(sprintf("LW lambda = %.5f  (interior in (0,1): %s)\n", lam, lam > 0 && lam < 1))
stopifnot(is.finite(lam), lam > 0, lam < 1, isTRUE(ee$applied))

# FD: directional derivative of w along each per-unit v_i = psi_i psi_i'/n - Om.
# Analytic claims D[i,] = J[v_i] = d/dh w(Om + h v_i) |_{h=0}. Compare to central FD.
Dana <- ee$D
hh <- 1e-6
max_rel <- 0
for (i in sample(seq_len(n), 40)) {        # 40 random units (incl. heavily-weighted)
  vi <- outer(psi[i, ], psi[i, ]) / n - Om
  num <- (w_of_omega(Om + hh * vi) - w_of_omega(Om - hh * vi)) / (2 * hh)
  den <- max(1, max(abs(num)))
  max_rel <- max(max_rel, max(abs(num - Dana[i, ])) / den)
}
cat(sprintf("max relative |D_analytic - D_FD| over 40 units = %.3e\n", max_rel))
cat(sprintf(">>> WEIGHTED FD-ORACLE: %s (tol 1e-5)\n", if (max_rel < 1e-5) "PASS" else "FAIL"))
