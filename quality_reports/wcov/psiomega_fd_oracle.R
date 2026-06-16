## =====================================================================
## quality_reports/wcov/psiomega_fd_oracle.R
## Weight-estimation (Sigma_Omega) channel FD oracle, UNDER OBS WEIGHTS.
##
## The analytic misspec_robust SE folds psi_omega (the IF of theta_hat through
## the estimated efficient weights W = Omega_hat^{-1}). This oracle computes the
## SAME IF by a case-weight finite difference: perturb unit l's weight ONLY in
## the Omega_hat estimation (phi and the OUTER Hajek average frozen), recompute
## W, and read theta_hat's response. The empirical IF_l^FD is compared to the
## analytic psi_omega[l].
##
## Calibrate on the UNWEIGHTED fit (analytic trusted there), then read the
## weighted discrepancy: if analytic/FD variance ratio jumps by ~n/n_eff under
## dispersed weights, the analytic over-states the weight-channel variance by the
## Kish factor (the degenerate-U / design-Bessel gap).
## Usage: Rscript quality_reports/wcov/psiomega_fd_oracle.R
## =====================================================================
suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE)
options(edid_mc_cores = 1L, edid_omega_method = "kernel_orig")   # psi-capable kernel builder

make_panel <- function(n = 400, seed = 21, disp = 0.8) {
  set.seed(seed)
  Tn  <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P   <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  wu   <- exp(disp * rnorm(n)); wu <- wu / mean(wu)
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = gc,            # never-treated coded as Inf (direct prepare_edid_panel)
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

# Build the per-cell weight-channel pieces for cell (g, t), at a given unit_weights vector.
cell_pieces <- function(panel, g, t, pairs) {
  pfn <- pairs
  self <- is.finite(pfn$gp) & pfn$gp == g; if (any(self)) pfn$gp[self] <- Inf
  crossp <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(crossp) > 0L) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(crossp$tpre))))
  pr <- suppressWarnings(estimate_all_propensity_ratios(panel, g, pfn, bs_df = 4L, K_folds = 1L,
            fold_id = rep(1L, panel$n), return_aux = TRUE, ratio_method = "exp"))
  ip <- suppressWarnings(estimate_all_inverse_propensities(panel, g, pairs, bs_df = 4L, K_folds = 1L,
            fold_id = rep(1L, panel$n), return_aux = TRUE, ratio_method = "exp"))
  cm_combos <- unique(rbind(data.frame(gp = pfn$gp, period = t), data.frame(gp = pfn$gp, period = pfn$tpre)))
  cm <- suppressWarnings(estimate_all_conditional_means(panel, pfn, t_val = t, bs_df = 4L, K_folds = 1L,
            fold_id = rep(1L, panel$n), return_aux = FALSE))
  list(pairs = pairs, pfn = pfn, pr = pr$predictions, ip = ip, cm = cm)
}

# theta_hat through the WEIGHT channel only: Omega from `uw_omega`, outer Hajek avg with `uw_outer`, phi frozen.
theta_wc <- function(panel, g, t, pc, phi_frozen, uw_omega, uw_outer) {
  p2 <- panel; p2$unit_weights <- uw_omega
  oa <- suppressWarnings(compute_omega_star_cov_edid(p2, g, t, pc$pairs, pc$pr, pc$cm, pc$ip,
            return_pointwise = TRUE))
  W  <- compute_pointwise_weights_edid(oa, d = ncol(panel$covariate_matrix))
  if (is.list(W)) W <- W$W
  wY <- rowSums(phi_frozen * W)
  if (is.null(uw_outer)) mean(wY) else stats::weighted.mean(wY, uw_outer)
}

run_oracle <- function(weighted, seed = 21, nL = 70L) {
  df    <- make_panel(seed = seed)
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1,
                              weightsname = if (weighted) "w" else NULL)
  n  <- panel$n
  uw <- panel$unit_weights                                   # NULL (unw) or mean-1 (wt)
  g <- 2; t <- 4
  pairs <- enumerate_valid_pairs_edid(g, panel$treatment_groups,
             as.numeric(names(panel$period_to_col)), panel$period_1, "all", 0L)
  pc  <- cell_pieces(panel, g, t, pairs)
  phi <- suppressWarnings(compute_generated_outcomes_cov_edid(panel, g, t, pc$pairs, pc$pr, pc$cm, "all"))
  if (is.list(phi)) phi <- phi$gen_out

  uw1  <- if (is.null(uw)) rep(1, n) else uw                 # working weight vector for Omega
  base <- theta_wc(panel, g, t, pc, phi, uw1, uw)

  # ---- analytic psi_omega for this cell (replicate edid-fit.R's assembly) ----
  oa <- suppressWarnings(compute_omega_star_cov_edid(panel, g, t, pc$pairs, pc$pr, pc$cm, pc$ip,
            return_pointwise = TRUE))
  pw <- compute_pointwise_weights_edid(oa, d = ncol(panel$covariate_matrix),
            gen_out_mat = phi, need_coup = TRUE)
  po <- suppressWarnings(compute_omega_star_cov_edid(panel, g, t, pc$pairs, pc$pr, pc$cm, pc$ip,
            return_pointwise = TRUE,
            psi_qw = list(pointwise = TRUE, Q = pw$Q, W = pw$W, lambda = attr(oa, "shrink_lambda"),
                          C = pw$C, ridge = !is.null(attr(oa, "ridge_lift")))))
  corr_an <- compute_invp_correction_analytic_cov_edid(n, attr(pc$ip, "aux"), po$coupled_C)
  psi_an  <- po$psi - corr_an                                # analytic weight-channel IF (length n)

  # ---- FD case-weight IF for a sample of units l ----
  set.seed(99); Ls <- sort(sample.int(n, min(nL, n)))
  delta <- 1e-4
  iffd <- numeric(length(Ls))
  for (k in seq_along(Ls)) {
    l <- Ls[k]
    up <- uw1; up[l] <- uw1[l] + delta
    um <- uw1; um[l] <- uw1[l] - delta
    iffd[k] <- (theta_wc(panel, g, t, pc, phi, up, uw) -
                theta_wc(panel, g, t, pc, phi, um, uw)) / (2 * delta)
  }
  # case-weight IF convention: the eif entry carries the perturbing unit's own weight, so
  # eif_l = uw_l * n * d(theta)/d(uw_l)  (the structural influence n*dtheta/duw_l, times uw_l, matching
  # the analytic psi_omega which carries uw_l via the weighted kernel). uw_l = 1 on the unweighted path.
  uw_l   <- if (is.null(uw)) rep(1, length(Ls)) else uw[Ls]
  eif_fd <- uw_l * n * iffd
  an_s   <- psi_an[Ls]
  # robust slope (analytic ~ c * FD); compare c across weighted/unweighted to detect a scale gap
  ok <- is.finite(eif_fd) & is.finite(an_s) & abs(eif_fd) > 1e-8
  slope <- sum(an_s[ok] * eif_fd[ok]) / sum(eif_fd[ok]^2)
  corr  <- suppressWarnings(stats::cor(an_s[ok], eif_fd[ok]))
  neff  <- if (is.null(uw)) n else (sum(uw)^2) / sum(uw^2)
  cat(sprintf("[%s] cell(2,4) nL=%d  cor(an,FD)=%.3f  slope(an/FD)=%.3f  n/n_eff=%.3f  sqrt=%.3f\n",
              if (weighted) "WT " else "UNW", sum(ok), corr, slope, n / neff, sqrt(n / neff)))
  invisible(list(slope = slope, corr = corr, neff = neff, n = n))
}

cat("=== weight-channel (psi_omega) FD oracle: analytic vs case-weight FD ===\n")
u <- run_oracle(FALSE)
w <- run_oracle(TRUE)
cat(sprintf("\n[RATIO] slope_wt / slope_unw = %.3f  (if ~sqrt(n/n_eff)=%.3f the analytic over-states the\n",
            w$slope / u$slope, sqrt(w$n / w$neff)))
cat("        weighted psi_omega by the Kish factor; if ~1 the analytic IF is correctly scaled.)\n")
