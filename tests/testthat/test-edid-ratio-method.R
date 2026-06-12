library(testthat)

# ---------------------------------------------------------------------------
# ratio_method: coherent (multinomial-logit sieve) cross-cohort nuisances.
#
# Regression guards for the audited with-X PT-All degeneracy: the legacy per-pair
# LS sieve for r_{g,g'} (and the per-cohort LS sieve for 1/p_{g'}) solves a system
# whose Gram matrix uses ONLY the n_{g'} comparison-cohort observations, so basis
# directions thin on g' explode -- producing large NEGATIVE fitted "ratios" on
# sizable shares of the comparison cohort, |r| > 1e4 tails, and fitted inverse
# propensities orders of magnitude above the 1/pi_{g'} scale, even under healthy
# cohort-vs-never-treated overlap. The coherent system (per-cohort ridge-logistic
# h_c = log(p_c/p_NT) on a shared basis; r = exp(h_g - h_g'); 1/p_c analytic from
# the implied multinomial shares) removes the thin-denominator Grams.
# ---------------------------------------------------------------------------

# Staggered DGP with one THIN comparison cohort (the audited failure mode):
# cohorts 3 (large), 4 (thin), never-treated; 2 covariates shifting cohort
# membership so the sieve has real signal to fit.
make_thin_cohort_panel <- function(n = 600L, seed = 42L, p_thin = 0.05) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- runif(n, -1, 1)
  sc <- exp(cbind(0, 0.8 * x1 - 0.3 * x2, log(p_thin / (1 - p_thin)) + 0.5 * x2))
  P  <- sc / rowSums(sc)
  g  <- vapply(seq_len(n), function(i) sample(c(Inf, 3, 4), 1L, prob = P[i, ]), numeric(1))
  df <- do.call(rbind, lapply(1:6, function(tt) {
    tau <- 1 * (is.finite(g) & tt >= g)
    # Inf-coded never-treated (accepted by edid() and required by direct prepare_edid_panel calls)
    data.frame(id = 1:n, t = tt, g = g, x1 = x1, x2 = x2,
               y = 0.5 * x1 + 0.2 * x2 + 0.3 * tt + tau + rnorm(n, 0, 0.7))
  }))
  df
}

test_that("direct per-pair LS ratio sieve degenerates on a thin comparison cohort; coherent does not", {
  df <- make_thin_cohort_panel()
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts
  pfn <- data.frame(gp = c(Inf, 4), tpre = c(2, 2))     # target g = 3; cross comparison g' = 4 (thin)
  fid <- rep(1L, panel$n)

  pr_dir <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g = 3, pairs = pfn, bs_df = 4L, K_folds = 1L, fold_id = fid, ratio_method = "direct"))
  pr_coh <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g = 3, pairs = pfn, bs_df = 4L, K_folds = 1L, fold_id = fid, ratio_method = "coherent"))

  i4 <- which(G == 4)                                   # the ONLY units where r_{3,4} enters the moment
  expect_gt(length(i4), 5L)
  r_dir <- pr_dir[["4"]][i4]
  r_coh <- pr_coh[["4"]][i4]
  # the documented failure: the direct fit assigns NEGATIVE "ratios" to part of the
  # comparison cohort and/or explodes far beyond the plausible scale
  expect_true(any(r_dir < 0) || max(abs(r_dir)) > 50 * max(r_coh))
  # the coherent ratio is positive by construction and bounded at the consumed units
  expect_true(all(r_coh > 0))
  expect_lt(max(r_coh), 1e3)
  # the never-treated ratio is the SAME LS fit under both methods (bitwise)
  expect_identical(pr_dir[["Inf"]], pr_coh[["Inf"]])
})

test_that("coherent ratios are algebraically coherent and the implied masses are calibrated", {
  df <- make_thin_cohort_panel(n = 800L, seed = 7L, p_thin = 0.15)
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts
  hs <- fit_coherent_logodds_edid(panel$covariate_matrix, G, cohorts = c(3, 4), bs_df = 4L)
  r34 <- exp(hs[["3"]]$h - hs[["4"]]$h)
  r43 <- exp(hs[["4"]]$h - hs[["3"]]$h)
  expect_equal(r34 * r43, rep(1, panel$n), tolerance = 1e-12)   # exact coherence r_{g,g'} r_{g',g} = 1
  # Hajek mass calibration: E_n[r_{g,g'} G_{g'}] ~ pi_g (the LS fit enforces it exactly;
  # the coherent fit should be close on a well-specified logistic DGP)
  expect_equal(mean(r34 * (G == 4)), mean(G == 3), tolerance = 0.35 * mean(G == 3))
})

test_that("direct LS inverse propensities explode on a thin cohort; coherent ones sit at the right scale", {
  df <- make_thin_cohort_panel()
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts
  pairs_g <- data.frame(gp = c(3, 4), tpre = c(2, 2))
  fid <- rep(1L, panel$n)
  ip_dir <- suppressWarnings(estimate_all_inverse_propensities(
    panel, g = 3, pairs = pairs_g, bs_df = 4L, K_folds = 1L, fold_id = fid, ratio_method = "direct"))
  ip_coh <- suppressWarnings(estimate_all_inverse_propensities(
    panel, g = 3, pairs = pairs_g, bs_df = 4L, K_folds = 1L, fold_id = fid, ratio_method = "coherent"))
  s4_scale <- 1 / mean(G == 4)                          # the unconditional 1/pi_4 anchor
  # coherent: mean within an order of magnitude of 1/pi_4, strictly positive
  expect_true(all(ip_coh[["4"]] > 0))
  expect_lt(mean(ip_coh[["4"]]), 20 * s4_scale)
  expect_gt(mean(ip_coh[["4"]]), s4_scale / 20)
  # direct: the audited s-channel pathology -- wild scale and/or a large clamped-to-zero share
  expect_true(mean(ip_dir[["4"]]) > 50 * s4_scale || mean(ip_dir[["4"]] == 0) > 0.2)
  # the never-treated inverse propensity is the SAME LS fit under both methods (bitwise)
  expect_identical(ip_dir[["Inf"]], ip_coh[["Inf"]])
})

test_that("with-X PT-All on the thin-cohort design: coherent default has sane SEs, direct inflates", {
  skip_on_cran()
  df <- make_thin_cohort_panel()
  f_nox <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "event_study",
                                 cband = FALSE, seed = 1))
  f_coh <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                                 aggregate = "event_study", cband = FALSE, seed = 1))
  f_dir <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                                 aggregate = "event_study", cband = FALSE, seed = 1,
                                 ratio_method = "direct"))
  se_nox <- f_nox$overall$overall.se
  se_coh <- f_coh$overall$overall.se
  se_dir <- f_dir$overall$overall.se
  expect_lt(se_coh, 6 * se_nox)              # restored: a sane multiple of the no-X SE
  expect_lt(se_coh, se_dir)                  # and strictly better than the legacy construction
  # the point estimate recovers the homogeneous true effect (1.0) within sampling noise
  expect_lt(abs(f_coh$overall$overall.att - 1), 4 * se_coh + 0.25)
})

test_that("PT-Post-X, no-X, and own-cohort moment sets are invariant to ratio_method", {
  df <- make_thin_cohort_panel(n = 400L, seed = 3L, p_thin = 0.10)
  # PT-Post with covariates: single never-treated moment -> bitwise identical
  fp_c <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2, pt_assumption = "post",
                                aggregate = "none", cband = FALSE))
  fp_d <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2, pt_assumption = "post",
                                aggregate = "none", cband = FALSE, ratio_method = "direct"))
  expect_identical(fp_c$att_gt$att, fp_d$att_gt$att)
  expect_identical(fp_c$att_gt$se,  fp_d$att_gt$se)
  # no-X: the argument is a covariate-path object -> bitwise identical
  fn_c <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none", cband = FALSE))
  fn_d <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none", cband = FALSE,
                                ratio_method = "direct"))
  expect_identical(fn_c$att_gt$att, fn_d$att_gt$att)
  expect_identical(fn_c$att_gt$se,  fn_d$att_gt$se)
  # own-cohort moment set: no cross pairs -> the cross-cohort RATIO construction never enters;
  # with fixed (uniform) weights the variance-channel prefactors 1/p_g do not enter either, so
  # the fit is bitwise invariant. (Under the efficient/averaged schemes own-set MOMENTS are
  # still identical, but the Omega prefactor 1/p_g legitimately moves the estimated weights:
  # the coherent system replaces the audited LS s-channel there too.)
  tg <- c(3, 4); tp <- 1:6
  ms <- do.call(rbind, lapply(tg, function(g) data.frame(g = g, gp = g, tpre = tp[tp < g])))
  fo_c <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2, aggregate = "none",
                                cband = FALSE, moment_set = ms, weight_scheme = "uniform"))
  fo_d <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2, aggregate = "none",
                                cband = FALSE, moment_set = ms, weight_scheme = "uniform",
                                ratio_method = "direct"))
  expect_identical(fo_c$att_gt$att, fo_d$att_gt$att)
  expect_identical(fo_c$att_gt$se,  fo_d$att_gt$se)
})

test_that("ratio_method is validated, stored on the fit, and snapshotted for refits", {
  df <- make_thin_cohort_panel(n = 300L, seed = 9L, p_thin = 0.2)
  expect_error(edid(df, "y", "id", "t", "g", ratio_method = "bogus"), "coherent")
  f <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                             cband = FALSE, ratio_method = "direct"))
  expect_identical(f$ratio_method, "direct")
  expect_identical(f$args$ratio_method, "direct")
  f2 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                              cband = FALSE))
  expect_identical(f2$ratio_method, "coherent")
})

# ---------------------------------------------------------------------------
# Ratio-targeted trim masks + the trim mask reaching the Omega/psi channel
# ---------------------------------------------------------------------------

test_that("finite-cohort trim masks key on the pair's ratio only; the never-treated mask keeps r AND s", {
  n <- 8L
  pr <- list("Inf" = c(1, 1, 500, 1, 1, 1, 1, 1),  "4" = c(1, 500, 1, 1, 1, 1, 1, 1))
  ip <- list("Inf" = c(1, 1, 1, 500, 1, 1, 1, 1),  "4" = c(500, 1, 1, 1, 1, 1, 1, 1),
             "3"   = rep(500, n))
  tk <- build_trim_keep_edid(pr, ip, trim_level = 200, n = n)
  expect_false(tk[["Inf"]][3])                     # r_{g,Inf} extreme -> trimmed
  expect_false(tk[["Inf"]][4])                     # 1/p_NT extreme -> trimmed (legacy never-treated mask)
  expect_false(tk[["4"]][2])                       # r_{g,4} extreme -> trimmed
  expect_true(tk[["4"]][1])                        # 1/p_4 extreme does NOT trim (variance channel only)
  expect_true(all(tk[["3"]]))                      # s-only finite key: never trimmed by s
  expect_null(build_trim_keep_edid(pr, ip, trim_level = Inf, n = n))
})

test_that("the cell-common keep mask reaches the Omega builders (trimmed units' prefactors zeroed)", {
  df <- make_thin_cohort_panel(n = 300L, seed = 11L, p_thin = 0.3)
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  g <- 3; t <- 4
  pairs <- enumerate_valid_pairs_edid(g, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  pfn <- pairs; sc <- is.finite(pfn$gp) & pfn$gp == g; pfn$gp[sc] <- Inf
  cr <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(cr)) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cr$tpre))))
  fid <- rep(1L, panel$n)
  pr <- suppressWarnings(estimate_all_propensity_ratios(panel, g, pfn, 4L, 1L, fid, ratio_method = "coherent"))
  cm <- suppressWarnings(estimate_all_conditional_means(panel, pfn, t, 4L, 1L, fid))
  ip <- suppressWarnings(estimate_all_inverse_propensities(panel, g, pairs, 4L, 1L, fid, ratio_method = "coherent"))

  keep <- rep(1, panel$n); keep[1:25] <- 0          # a hand-built cell-common trim mask
  old <- options(edid_shrink_lambda = 0)            # disable the pooled blend so zero rows stay zero
  on.exit(options(old), add = TRUE)
  arr0 <- suppressWarnings(compute_omega_star_kernel_fast_edid(panel, g, t, pairs, pr, cm, ip,
                                                               return_pointwise = TRUE))
  arr1 <- suppressWarnings(compute_omega_star_kernel_fast_edid(panel, g, t, pairs, pr, cm, ip,
                                                               return_pointwise = TRUE, keep = keep))
  expect_equal(max(abs(arr1[1:25, , ])), 0)         # trimmed units: every Omega entry zeroed
  expect_equal(arr1[26:50, , ], arr0[26:50, , ], tolerance = 1e-12)  # kept units: unchanged slices
  # keep = NULL and keep = all-ones agree
  arr2 <- suppressWarnings(compute_omega_star_kernel_fast_edid(panel, g, t, pairs, pr, cm, ip,
                                                               return_pointwise = TRUE, keep = rep(1, panel$n)))
  expect_equal(arr2, arr0, tolerance = 1e-15)
  # pooled (averaged) Omega: keep changes the average; the two kernel builders stay build-invariant
  om_fast <- suppressWarnings(compute_omega_star_kernel_fast_edid(panel, g, t, pairs, pr, cm, ip, keep = keep))
  om_orig <- suppressWarnings(compute_omega_star_cov_edid(panel, g, t, pairs, pr, cm, ip, keep = keep))
  expect_equal(unclass(om_fast), unclass(om_orig), tolerance = 1e-9, check.attributes = FALSE)
})

test_that("legacy floor remains reachable via options(edid_legacy_floor = TRUE)", {
  df <- make_thin_cohort_panel(n = 300L, seed = 13L, p_thin = 0.3)
  f_new <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                 cband = FALSE, weight_scheme = "averaged"))
  old <- options(edid_legacy_floor = TRUE)
  f_leg <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                 cband = FALSE, weight_scheme = "averaged"))
  options(old)
  expect_false(isTRUE(all.equal(f_new$att_gt$att, f_leg$att_gt$att, tolerance = 1e-12)))
})
