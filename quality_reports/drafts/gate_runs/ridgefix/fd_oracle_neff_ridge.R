#!/usr/bin/env Rscript
# FD-ORACLE: weighted ridge (and LW) corrected SE under the n_eff (Kish ESS) lift.
# Claim under audit (plan/audit doc): with the ridge intensity lam_r = H/n_eff (n_eff a
# function of the FIXED weights only, CONSTANT in Omega-hat), the EE "plain map" branch
# (use_chain = FALSE), evaluated at omega_used = the n_eff-RIDGED Omega, is the EXACT
# weight-estimation correction -- i.e. the analytic per-unit Jacobian directions D match
# the finite-difference Jacobian of the ACTUAL realized ridge weight map.
#
# The subtlety this oracle is built to CATCH: the realized ridge matrix
#     M(Omega) = Omega + lam_r * mean(diag(Omega)) * I
# depends on Omega-hat through mean(diag(Omega)) too. So when Omega-hat moves by the
# per-unit direction v_i = psi_i psi_i'/n - Omega_raw, the TRUE perturbation to the
# inverted matrix is  dM = v_i + lam_r * mean(diag(v_i)) * I, not just v_i. The FD map
# below uses the TRUE M(Omega) (recomputing the ridge each perturbation), so if the
# plain-map analytic D omitted the isotropic mean(diag) channel the oracle would FAIL.
#
# Tolerance: 1e-8 (relative, max over all entries).  Reported both relative and absolute.

suppressMessages(pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE))
set.seed(20260615)
options(digits = 12)

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- weighted staggered panel builder -------------------------------------
mk_wt_panel <- function(n, seed, disp = 1.0, Tn = 6L, rho = 0.5) {
  set.seed(seed)
  G <- sample(c(3L, 5L, 0L), n, TRUE, c(.35, .30, .35))
  e <- matrix(rnorm(n * Tn), n, Tn)
  for (s in 2:Tn) e[, s] <- rho * e[, s - 1L] + sqrt(1 - rho^2) * e[, s]
  Y <- rnorm(n) + matrix(0.3 * seq_len(Tn), n, Tn, byrow = TRUE) + e
  for (g in c(3L, 5L)) for (t in g:Tn) Y[G == g, t] <- Y[G == g, t] + 1 + 0.3 * (t - g)
  # dispersed positive unit weights (lognormal); replicated across a unit's rows
  raw_w <- exp(disp * rnorm(n))
  data.frame(id = rep(seq_len(n), each = Tn), time = rep(seq_len(Tn), n),
             y = as.vector(t(Y)), g = rep(G, each = Tn),
             w = rep(raw_w, each = Tn))
}

## ---- the EXACT realized ridge weight map (mirrors edid-fit.R:956-960) -------
# Given a RAW Omega-hat, apply the n_eff ridge and return the efficient weights.
# lam_r is passed in as a FIXED scalar (= H / n_eff), exactly as the fit computes it
# (n_eff is a function of the fixed weights, NOT of Omega-hat).
w_of_omega_ridge <- function(Om_raw, lam_r) {
  H <- nrow(Om_raw)
  M <- Om_raw + (lam_r * mean(diag(Om_raw))) * diag(H)   # ridge term DEPENDS on Om_raw
  A <- solve(M); u <- drop(A %*% rep(1, H)); u / sum(u)
}

## ---- one cell's analytic-vs-FD comparison ---------------------------------
oracle_cell <- function(df, gg, tt, weighted = TRUE, tol = 1e-8, h_scale = 1e-6) {
  df$g <- ifelse(df$g == 0L, Inf, df$g)
  pn <- if (weighted) prepare_edid_panel(df, "y", "id", "time", "g", anticipation = 0L,
                                         weightsname = "w")
        else          prepare_edid_panel(df, "y", "id", "time", "g", anticipation = 0L)
  prs <- enumerate_valid_pairs_edid(gg, pn$treatment_groups, pn$time_periods,
                                    pn$period_1, "all", 0L)
  H <- nrow(prs)
  if (H < 2L) return(NULL)
  Om_raw <- compute_omega_star_nocov_edid(gg, tt, prs, pn, "all")   # exact Psi'Psi/n^2 (weighted)
  psi    <- compute_psi_moments_nocov_edid(gg, tt, prs, pn)         # n x H
  n_u    <- pn$n

  # n_eff exactly as the fit does it (active-unit Kish ESS)
  amask <- active_mask_nocov_edid(gg, prs, pn)
  n_eff <- n_eff_edid(pn$unit_weights, amask, pn$n)
  lam_r <- H / n_eff
  M     <- Om_raw + (lam_r * mean(diag(Om_raw))) * diag(H)          # the ridged Omega (omega_used)
  w     <- compute_efficient_weights_edid(M)

  # analytic plain-map D from the package (use_chain = FALSE: shrink_lambda = NA)
  ee <- compute_nocov_ee_correction_edid(gg, tt, prs, pn, omega_raw = Om_raw,
                                         omega_used = M, weights = w,
                                         shrink_lambda = NA_real_, return_D = TRUE)
  if (!isTRUE(ee$applied)) return(list(applied = FALSE, reason = ee$reason,
                                       g = gg, t = tt, weighted = weighted))
  Dan <- ee$D

  # FD Jacobian of the TRUE realized ridge map, direction v_i (raw-Omega perturbation)
  h   <- h_scale * max(abs(Om_raw))
  Dfd <- matrix(NA_real_, n_u, H)
  for (i in seq_len(n_u)) {
    Vi <- tcrossprod(psi[i, ]) / n_u - Om_raw
    wp <- w_of_omega_ridge(Om_raw + h * Vi, lam_r)
    wm <- w_of_omega_ridge(Om_raw - h * Vi, lam_r)
    Dfd[i, ] <- (wp - wm) / (2 * h)
  }

  scaleD <- max(abs(Dfd), 1e-300)
  rel_D  <- max(abs(Dan - Dfd)) / scaleD
  abs_D  <- max(abs(Dan - Dfd))

  # assembled var_add via the SAME assembly but with FD directions
  a    <- drop(psi %*% w)
  qfd  <- -sum(a * rowSums(Dfd * psi)) / n_u^3
  va_fd <- ee$delta_df + 2 * qfd
  rel_va <- abs(ee$var_add - va_fd) / max(abs(ee$var_add), 1e-300)

  # also report the raw-n (legacy) intensity for context: how big is n/n_eff?
  list(applied = TRUE, g = gg, t = tt, weighted = weighted, H = H, n = n_u,
       n_eff = n_eff, ridge_factor = n_u / n_eff, lam_r = lam_r,
       rel_D = rel_D, abs_D = abs_D, rel_va = rel_va,
       maxw = max(abs(w)), pass = (rel_D < tol && rel_va < tol))
}

## ===========================================================================
## RUN: weighted staggered designs + an unweighted control + a sieve/cov touch
## ===========================================================================
cat("=== FD-ORACLE: n_eff weighted RIDGE plain-map EE (tol 1e-8) ===\n\n")

configs <- list(
  list(tag = "WEIGHTED disp=1.0", n = 120, seed = 101, disp = 1.0, g = 3L, t = 4L, wtd = TRUE),
  list(tag = "WEIGHTED disp=1.0", n = 120, seed = 101, disp = 1.0, g = 3L, t = 5L, wtd = TRUE),
  list(tag = "WEIGHTED disp=1.5", n = 160, seed = 202, disp = 1.5, g = 5L, t = 6L, wtd = TRUE),
  list(tag = "WEIGHTED disp=0.7", n = 140, seed = 303, disp = 0.7, g = 3L, t = 6L, wtd = TRUE),
  list(tag = "WEIGHTED disp=2.0", n = 200, seed = 404, disp = 2.0, g = 5L, t = 5L, wtd = TRUE),
  list(tag = "UNWEIGHTED ctrl",   n = 120, seed = 101, disp = 1.0, g = 3L, t = 4L, wtd = FALSE)
)

rows <- list()
allpass <- TRUE
for (cf in configs) {
  df <- mk_wt_panel(cf$n, cf$seed, cf$disp)
  res <- oracle_cell(df, cf$g, cf$t, weighted = cf$wtd, tol = 1e-8)
  if (is.null(res)) { cat(sprintf("  [skip] %s g=%d t=%d : H<2\n", cf$tag, cf$g, cf$t)); next }
  if (!isTRUE(res$applied)) {
    cat(sprintf("  [n/a]  %s g=%d t=%d : %s\n", cf$tag, cf$g, cf$t, res$reason)); next
  }
  ok <- res$pass; allpass <- allpass && ok
  cat(sprintf("  [%s] %-18s g=%d t=%d | n=%d n_eff=%7.2f  n/n_eff=%.3f  lam_r=%.4f\n",
              if (ok) "PASS" else "FAIL", res$tag <- cf$tag, cf$g, cf$t,
              res$n, res$n_eff, res$ridge_factor, res$lam_r))
  cat(sprintf("           rel|D-Dfd|=%.3e  abs|D-Dfd|=%.3e  rel|var_add-FD|=%.3e  max|w|=%.3f\n",
              res$rel_D, res$abs_D, res$rel_va, res$maxw))
  rows[[length(rows) + 1L]] <- res
}

cat(sprintf("\nworst rel|D-Dfd| over applied cells = %.3e\n",
            max(vapply(rows, function(r) r$rel_D, 0))))
cat(sprintf("worst rel|var_add-FD| over applied cells = %.3e\n",
            max(vapply(rows, function(r) r$rel_va, 0))))
cat(sprintf("\nOVERALL: %s\n", if (allpass) "PASS (all cells < 1e-8)" else "FAIL"))
saveRDS(rows, "/tmp/gate_runs/ridgefix/fd_oracle_neff_ridge_rows.rds")
