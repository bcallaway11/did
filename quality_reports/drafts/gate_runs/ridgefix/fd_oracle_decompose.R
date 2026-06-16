#!/usr/bin/env Rscript
# DECOMPOSITION oracle: isolate WHAT the plain-map analytic D actually differentiates.
# Three candidate realized maps, all evaluated at the SAME ridged omega_used = M:
#   (A) TRUE ridge map:   w(Om_raw + lam_r*mean(diag(Om_raw))*I)        [ridge moves with Om]
#   (B) frozen-c map:     w(Om_raw + c0*I), c0 = lam_r*mean(diag(Om_raw0)) FIXED   [additive const]
#   (C) pure map at M:    w(M + t*V_i), i.e. perturb M directly by V_i (no ridge wrapper)
# The package plain-map D uses omega_used = M, direction v_i = psi psi'/n - Om_raw.
# Hypothesis: plain-map D == FD of map (C)  (it differentiates w at the point M in the
# RAW direction v_i), and is NOT the derivative of the realized data map (A).
# This tells us whether the EE correction is the EXACT weight-estimation correction for
# the realized ridge estimator, or only for a 'frozen-ridge' approximation.

suppressMessages(pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE))
options(digits = 12)

mk_wt_panel <- function(n, seed, disp = 1.0, Tn = 6L, rho = 0.5) {
  set.seed(seed)
  G <- sample(c(3L, 5L, 0L), n, TRUE, c(.35, .30, .35))
  e <- matrix(rnorm(n * Tn), n, Tn)
  for (s in 2:Tn) e[, s] <- rho * e[, s - 1L] + sqrt(1 - rho^2) * e[, s]
  Y <- rnorm(n) + matrix(0.3 * seq_len(Tn), n, Tn, byrow = TRUE) + e
  for (g in c(3L, 5L)) for (t in g:Tn) Y[G == g, t] <- Y[G == g, t] + 1 + 0.3 * (t - g)
  raw_w <- exp(disp * rnorm(n))
  data.frame(id = rep(seq_len(n), each = Tn), time = rep(seq_len(Tn), n),
             y = as.vector(t(Y)), g = rep(G, each = Tn), w = rep(raw_w, each = Tn))
}

run_cell <- function(df, gg, tt, weighted, tag) {
  df$g <- ifelse(df$g == 0L, Inf, df$g)
  pn <- if (weighted) prepare_edid_panel(df, "y","id","time","g", anticipation=0L, weightsname="w")
        else          prepare_edid_panel(df, "y","id","time","g", anticipation=0L)
  prs <- enumerate_valid_pairs_edid(gg, pn$treatment_groups, pn$time_periods, pn$period_1,"all",0L)
  H <- nrow(prs); if (H < 2L) return(invisible())
  Om <- compute_omega_star_nocov_edid(gg, tt, prs, pn, "all")
  psi <- compute_psi_moments_nocov_edid(gg, tt, prs, pn); n_u <- pn$n
  amask <- active_mask_nocov_edid(gg, prs, pn)
  n_eff <- n_eff_edid(pn$unit_weights, amask, pn$n)
  lam_r <- H / n_eff
  c0 <- lam_r * mean(diag(Om))
  M  <- Om + c0 * diag(H)
  w  <- compute_efficient_weights_edid(M)
  ee <- compute_nocov_ee_correction_edid(gg, tt, prs, pn, omega_raw = Om, omega_used = M,
                                         weights = w, shrink_lambda = NA_real_, return_D = TRUE)
  if (!isTRUE(ee$applied)) { cat(sprintf("  %s: n/a (%s)\n", tag, ee$reason)); return(invisible()) }
  Dan <- ee$D
  pure_w <- function(MM) { A <- solve(MM); u <- drop(A %*% rep(1,H)); u/sum(u) }
  h <- 1e-6 * max(abs(Om))
  DA <- DB <- DC <- matrix(NA_real_, n_u, H)
  for (i in seq_len(n_u)) {
    Vi <- tcrossprod(psi[i, ]) / n_u - Om
    # (A) true ridge map: ridge recomputed each perturbation
    rA <- (pure_w((Om+h*Vi) + (lam_r*mean(diag(Om+h*Vi)))*diag(H)) -
           pure_w((Om-h*Vi) + (lam_r*mean(diag(Om-h*Vi)))*diag(H))) / (2*h)
    # (B) frozen additive constant c0
    rB <- (pure_w((Om+h*Vi) + c0*diag(H)) - pure_w((Om-h*Vi) + c0*diag(H))) / (2*h)
    # (C) perturb M directly by Vi
    rC <- (pure_w(M + h*Vi) - pure_w(M - h*Vi)) / (2*h)
    DA[i,] <- rA; DB[i,] <- rB; DC[i,] <- rC
  }
  rel <- function(Dfd) max(abs(Dan - Dfd)) / max(abs(Dfd), 1e-300)
  cat(sprintf("  %-22s | n/n_eff=%.3f lam_r=%.4f\n", tag, n_u/n_eff, lam_r))
  cat(sprintf("     D vs (A) TRUE ridge map        rel = %.3e\n", rel(DA)))
  cat(sprintf("     D vs (B) frozen-const ridge    rel = %.3e\n", rel(DB)))
  cat(sprintf("     D vs (C) perturb M directly    rel = %.3e\n", rel(DC)))
  # also: does (B) == (C)? (frozen-const additive vs direct M perturbation)
  cat(sprintf("     (B) vs (C) FD-FD               rel = %.3e\n",
              max(abs(DB - DC))/max(abs(DC),1e-300)))
}

cat("=== DECOMPOSITION: what does plain-map D differentiate? ===\n\n")
run_cell(mk_wt_panel(120,101,1.0), 3L, 4L, TRUE,  "WEIGHTED disp=1.0")
run_cell(mk_wt_panel(200,404,2.0), 5L, 5L, TRUE,  "WEIGHTED disp=2.0")
run_cell(mk_wt_panel(120,101,1.0), 3L, 4L, FALSE, "UNWEIGHTED")
