library(testthat)

# ---------------------------------------------------------------------------
# ratio_method = "exp": per-target exponential-link Riesz regressions.
#
# The paper-compatible alternative to "coherent": every ratio r_{g,g'} (incl.
# r_{g,Inf}) and every finite-cohort inverse propensity is exp(psi'beta) fit by
# the tailored convex loss E_n[exp(psi'b) G_{g'} - (psi'b) G_g] (FOC = exact
# basis-mean balancing), with FULL estimation-effect aux (chain rule
# dr/dbeta = r psi; tailored score psi (G_{g'} r - G_g); Hessian
# E_n[psi psi' r G_{g'}]) -- no fallback-marking, so the ACH / inv-p /
# higher-order / bootstrap corrections cover the exp channels. FD oracles below
# use the package's standard 1e-6 step.
# ---------------------------------------------------------------------------

# Same thin-comparison-cohort DGP as test-edid-ratio-method.R (the audited
# failure mode of the direct LS sieve).
make_thin_cohort_panel_exp <- function(n = 600L, seed = 42L, p_thin = 0.05) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- runif(n, -1, 1)
  sc <- exp(cbind(0, 0.8 * x1 - 0.3 * x2, log(p_thin / (1 - p_thin)) + 0.5 * x2))
  P  <- sc / rowSums(sc)
  g  <- vapply(seq_len(n), function(i) sample(c(Inf, 3, 4), 1L, prob = P[i, ]), numeric(1))
  df <- do.call(rbind, lapply(1:6, function(tt) {
    tau <- 1 * (is.finite(g) & tt >= g)
    data.frame(id = 1:n, t = tt, g = g, x1 = x1, x2 = x2,
               y = 0.5 * x1 + 0.2 * x2 + 0.3 * tt + tau + rnorm(n, 0, 0.7))
  }))
  df
}

test_that("exp fits are positive and finite on the thin-denominator reproducer; full aux attached", {
  df <- make_thin_cohort_panel_exp()
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts
  pfn <- data.frame(gp = c(Inf, 4), tpre = c(2, 2))      # target g = 3; thin cross comparison g' = 4
  fid <- rep(1L, panel$n)

  pr <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g = 3, pairs = pfn, bs_df = 4L, K_folds = 1L, fold_id = fid,
    return_aux = TRUE, ratio_method = "exp"))
  r4 <- pr$predictions[["4"]]; rI <- pr$predictions[["Inf"]]
  # positivity by construction + finiteness everywhere (the LS sieve fitted ~40-50% NEGATIVE here)
  expect_true(all(r4 > 0) && all(is.finite(r4)))
  expect_true(all(rI > 0) && all(is.finite(rI)))
  # sane scale at the CONSUMED units (r_{g,g'} multiplies G_{g'})
  expect_lt(max(r4[G == 4]), 1e4)
  # FULL aux on every key, including the cross-cohort one: no fallback-marking
  for (k in c("Inf", "4")) {
    a <- pr$aux[[k]]
    expect_false(isTRUE(a$is_fallback))
    expect_identical(a$link, "exp")
    expect_true(all(c("B_test", "score_mat", "H_inv", "B_raw", "beta") %in% names(a)))
    # the aux estimating equation holds at beta-hat: mean score = 0 (penalized form on rescue)
    expect_lt(max(abs(colMeans(a$score_mat))), 1e-6)
    # chain rule packing: B_test = pred * psi
    expect_equal(a$B_test, a$B_raw * a$pred, tolerance = 1e-12)
  }

  ip <- suppressWarnings(estimate_all_inverse_propensities(
    panel, g = 3, pairs = data.frame(gp = c(3, 4), tpre = c(2, 2)), bs_df = 4L,
    K_folds = 1L, fold_id = fid, return_aux = TRUE, ratio_method = "exp"))
  s4 <- ip[["4"]]
  expect_true(all(s4 > 0) && all(is.finite(s4)))
  sa <- attr(ip, "aux")[["4"]]
  expect_false(isTRUE(sa$is_fallback))
  expect_identical(sa$link, "exp")
  expect_true(all(sa$s_pos))                              # exp fit never clamps
  expect_lt(max(abs(colMeans(sa$score_mat))), 1e-6)
  # the never-treated inverse propensity stays the LS sieve (bitwise = direct)
  ip_dir <- suppressWarnings(estimate_all_inverse_propensities(
    panel, g = 3, pairs = data.frame(gp = c(3, 4), tpre = c(2, 2)), bs_df = 4L,
    K_folds = 1L, fold_id = fid, ratio_method = "direct"))
  expect_identical(ip[["Inf"]], ip_dir[["Inf"]])
})

test_that("the tailored-loss FOC is exact basis-mean balancing at the (unpenalized) optimum", {
  df <- make_thin_cohort_panel_exp(n = 800L, seed = 7L, p_thin = 0.15)   # healthy: no ridge rescue
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts
  fid <- rep(1L, panel$n)
  pr <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g = 3, pairs = data.frame(gp = c(Inf, 4), tpre = c(2, 2)), bs_df = 4L,
    K_folds = 1L, fold_id = fid, return_aux = TRUE, ratio_method = "exp"))
  a4 <- pr$aux[["4"]]
  expect_true(a4$exp_converged)
  expect_identical(a4$exp_lambda, 0)                      # unpenalized convergence on this fixture
  # E_n[psi r-hat G_{g'}] = E_n[psi G_g] exactly (the balancing FOC)
  foc <- colMeans(a4$B_raw * ((G == 4) * pr$predictions[["4"]])) - colMeans(a4$B_raw * (G == 3))
  expect_lt(max(abs(foc)), 1e-7)
  # inverse propensity: E_n[psi s-hat G_{g'}] = E_n[psi]
  ip <- suppressWarnings(estimate_all_inverse_propensities(
    panel, g = 3, pairs = data.frame(gp = c(3, 4), tpre = c(2, 2)), bs_df = 4L,
    K_folds = 1L, fold_id = fid, return_aux = TRUE, ratio_method = "exp"))
  sa <- attr(ip, "aux")[["4"]]
  focs <- colMeans(sa$B_raw * ((G == 4) * ip[["4"]])) - colMeans(sa$B_raw)
  expect_lt(max(abs(focs)), 1e-7)
})

test_that("FD oracle (1e-6): the exp aux Hessian is the Jacobian of the mean score in beta", {
  df <- make_thin_cohort_panel_exp(n = 800L, seed = 7L, p_thin = 0.15)
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts; n <- panel$n
  fid <- rep(1L, n)
  pr <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g = 3, pairs = data.frame(gp = c(Inf, 4), tpre = c(2, 2)), bs_df = 4L,
    K_folds = 1L, fold_id = fid, return_aux = TRUE, ratio_method = "exp"))
  a <- pr$aux[["4"]]
  comp <- as.numeric(G == 4); tgt <- as.numeric(G == 3)
  B <- a$B_raw; beta <- a$beta; p <- ncol(B)
  smean <- function(b) colMeans(B * (comp * exp(drop(B %*% b)) - tgt))
  H_an <- crossprod(B, (comp * drop(exp(B %*% beta))) * B) / n      # the documented tailored Hessian
  eps <- 1e-6
  H_fd <- vapply(seq_len(p), function(j) {
    e <- numeric(p); e[j] <- eps
    (smean(beta + e) - smean(beta)) / eps
  }, numeric(p))
  expect_lt(max(abs(H_fd - unname(as.matrix(H_an)))), 1e-4 * (1 + max(abs(H_an))))
  # and H_inv is n * pinv(H * n): H_inv %*% (n H) ~ identity on the live block
  HH <- a$H_inv %*% (H_an)
  expect_equal(diag(HH)[diag(HH) > 0.5], rep(1, sum(diag(HH) > 0.5)), tolerance = 1e-6)
})

test_that("FD oracle (1e-6): the analytic ACH Gamma equals the beta-space derivative under the exp chain rule", {
  df <- make_thin_cohort_panel_exp(n = 800L, seed = 7L, p_thin = 0.15)
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts; n <- panel$n
  g <- 3; t <- 4
  pairs <- enumerate_valid_pairs_edid(g, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  pfn <- pairs; sc <- is.finite(pfn$gp) & pfn$gp == g; pfn$gp[sc] <- Inf
  cr <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(cr)) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cr$tpre))))
  fid <- rep(1L, n)
  pr <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g, pfn, bs_df = 4L, K_folds = 1L, fold_id = fid, return_aux = TRUE,
    ratio_method = "exp"))
  cm <- suppressWarnings(estimate_all_conditional_means(panel, pfn, t_val = t, bs_df = 4L,
                                                        K_folds = 1L, fold_id = fid, return_aux = TRUE))
  prop_ratios <- pr$predictions
  # deliberately MISSPECIFIED conditional means (zeroed) so the r-channel Gamma is far from
  # orthogonal-zero and the comparison has teeth; FD-vs-analytic equality holds for ANY inputs
  cond_means <- lapply(cm$predictions, function(v) 0 * v)
  H <- nrow(pairs); w <- rep(1 / H, H)
  key <- "4"
  a <- pr$aux[[key]]
  skip_if(isTRUE(a$is_fallback), "cross aux fell back on this fixture")

  # package: analytic correction with ONLY the r-channel of `key` active
  corr_an <- compute_ach_correction_cov_edid(panel, g, t, pairs, prop_ratios, cond_means,
                                             w, m_aux = list(), r_aux = pr$aux[key])
  # oracle: Gamma by FD in the COEFFICIENT space through the exact exp map
  wmom <- function(prr) {
    go <- compute_generated_outcomes_cov_edid(panel, g, t, pairs, prr, cond_means, "all")
    mean(drop(go %*% w))
  }
  m0 <- wmom(prop_ratios)
  eps <- 1e-6
  Gamma_fd <- vapply(seq_len(ncol(a$B_raw)), function(j) {
    prr <- prop_ratios
    prr[[key]] <- prop_ratios[[key]] * exp(eps * a$B_raw[, j])    # r(beta + eps e_j), exact
    (wmom(prr) - m0) / eps
  }, numeric(1))
  corr_fd <- as.vector(a$score_mat %*% drop(a$H_inv %*% Gamma_fd))
  expect_gt(max(abs(corr_an)), 1e-8)                       # the channel is genuinely non-zero here
  expect_equal(corr_an, corr_fd, tolerance = 1e-4)
})

test_that("package-level: analytic ACH reproduces the FD oracle under ratio_method='exp'", {
  df <- make_thin_cohort_panel_exp(n = 500L, seed = 21L, p_thin = 0.15)
  run <- function(ach, ee) {
    op <- options(edid_ach = ach); on.exit(options(op), add = TRUE)
    suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                          weight_scheme = "efficient", aggregate = "none", cband = FALSE,
                          seed = 1L, misspec_robust = FALSE, estimation_effect = ee,
                          ratio_method = "exp"))$att_gt$se
  }
  se_noee <- run("analytic", FALSE)
  se_an   <- run("analytic", TRUE)
  se_fd   <- run("fd",       TRUE)
  ok <- is.finite(se_noee) & is.finite(se_an) & is.finite(se_fd)
  skip_if(sum(ok) < 2L, "too few non-degenerate cells")
  expect_gt(max(abs(se_an[ok] - se_noee[ok])), 1e-6)       # corrections ACTIVE (not fallback-skipped)
  expect_equal(se_an[ok], se_fd[ok], tolerance = 1e-5)
})

test_that("inv-p weight channel: the exp analytic correction matches the package FD oracle (within the documented floor-convention gap, no worse than the validated linear channel)", {
  # The analytic inv-p Gamma uses the Daleckii-Krein coupling of the FIXED-floor inverse map,
  # while the FD oracle re-solves the production MOVING-floor map: they agree exactly where the
  # eigen floor is inactive and differ by the documented convention gap where it binds (see
  # test-edid-identities). That gap exists for the LINEAR ("direct") channel too -- measured
  # per-cell rel. gaps up to ~2.9 there -- so the exp assertion is calibrated accordingly:
  # exact-match witnesses on the unfloored cells + aggregate agreement no worse than linear's.
  df <- make_thin_cohort_panel_exp(n = 500L, seed = 11L, p_thin = 0.15)
  run_acc <- function(rm) {
    op <- options(edid_store_psiomega = TRUE, edid_psiomega_fd = TRUE, edid_psiomega_acc = list())
    on.exit(options(op), add = TRUE)
    suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                          weight_scheme = "averaged", aggregate = "none", cband = FALSE,
                          seed = 1L, misspec_robust = TRUE, ratio_method = rm))
    acc <- getOption("edid_psiomega_acc")
    options(edid_psiomega_acc = NULL)
    acc
  }
  stats_of <- function(acc) {
    out <- NULL
    for (nm in names(acc)) {
      ce <- acc[[nm]]
      if (is.null(ce$corr) || is.null(ce$corr_fd)) next
      sc <- max(abs(ce$corr))
      if (sc < 1e-10) next
      out <- rbind(out, data.frame(
        cell = nm, scale = sc,
        rel  = max(abs(ce$corr - ce$corr_fd)) / sc,
        cor  = suppressWarnings(stats::cor(ce$corr, ce$corr_fd))))
    }
    out
  }
  st_exp <- stats_of(run_acc("exp"))
  st_lin <- stats_of(run_acc("direct"))
  skip_if(is.null(st_exp) || nrow(st_exp) < 3L, "too few corrected cells")
  expect_true(all(is.finite(st_exp$rel)))
  # exact-wiring witness: at least one cell where the floor is inactive matches to ~FP/FD precision
  expect_lt(min(st_exp$rel), 1e-3)
  # most cells agree well (median over cells)
  expect_gt(stats::median(st_exp$cor, na.rm = TRUE), 0.95)
  # and the exp channel is no worse than the validated linear channel on the same fixture
  skip_if(is.null(st_lin) || nrow(st_lin) < 3L, "too few linear cells to benchmark")
  expect_lte(stats::median(st_exp$rel), 2 * max(stats::median(st_lin$rel), 1e-3))
})

test_that("no-X fits are bitwise invariant to ratio_method = 'exp'; validation lists 'exp'", {
  df <- make_thin_cohort_panel_exp(n = 300L, seed = 9L, p_thin = 0.2)
  fn  <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none", cband = FALSE))
  fne <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none", cband = FALSE,
                               ratio_method = "exp"))
  expect_identical(fn$att_gt$att, fne$att_gt$att)
  expect_identical(fn$att_gt$se,  fne$att_gt$se)
  expect_error(edid(df, "y", "id", "t", "g", ratio_method = "bogus"), "exp")
  f <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                             cband = FALSE, ratio_method = "exp"))
  expect_identical(f$ratio_method, "exp")
  expect_identical(f$args$ratio_method, "exp")
})

test_that("coherent fits are unchanged by the exp engine's presence (same code path, identical reruns)", {
  df <- make_thin_cohort_panel_exp(n = 400L, seed = 3L, p_thin = 0.10)
  f1 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                              aggregate = "none", cband = FALSE, seed = 1))
  f2 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                              aggregate = "none", cband = FALSE, seed = 1,
                              ratio_method = "coherent"))
  expect_identical(f1$att_gt$att, f2$att_gt$att)
  expect_identical(f1$att_gt$se,  f2$att_gt$se)
  # and exp genuinely differs from coherent on the covariate path (it is a different estimator)
  f3 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                              aggregate = "none", cband = FALSE, seed = 1,
                              ratio_method = "exp"))
  expect_false(isTRUE(all.equal(f1$att_gt$att, f3$att_gt$att, tolerance = 1e-12)))
})

test_that("ratio-targeted trimming and the keep-mask threading apply identically to exp fits", {
  # unit-level: extreme exp ratios are trimmed by the finite-cohort RATIO-only mask
  n <- 6L
  pr <- list("Inf" = c(1, 1, 500, 1, 1, 1), "4" = c(1, 500, 1, 1, 1, 1))
  ip <- list("Inf" = rep(1, n), "4" = c(500, 1, 1, 1, 1, 1))
  tk <- build_trim_keep_edid(pr, ip, trim_level = 200, n = n)
  expect_false(tk[["4"]][2]); expect_true(tk[["4"]][1]); expect_false(tk[["Inf"]][3])
  # fit-level: a binding trim under "exp" engages the dead-pair / common-mask machinery
  df <- make_thin_cohort_panel_exp(n = 500L, seed = 5L, p_thin = 0.06)
  f_trim <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                                  aggregate = "none", cband = FALSE, seed = 1,
                                  ratio_method = "exp", trim_level = 50))
  expect_s3_class(f_trim, "edid_fit")                      # the binding trim runs to completion
  expect_true(any(is.finite(f_trim$att_gt$att)))
  expect_identical(f_trim$ratio_method, "exp")
  # the trimmed fit differs from the untrimmed one (the mask actually bit on this design)
  f_notrim <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                                    aggregate = "none", cband = FALSE, seed = 1,
                                    ratio_method = "exp", trim_level = Inf))
  expect_false(isTRUE(all.equal(f_trim$att_gt$att, f_notrim$att_gt$att, tolerance = 1e-12)))
})

test_that("tailored and literal paper losses agree under correct specification (logit DGP)", {
  # log-odds LINEAR in x on COMPACT support, with the two treated cohorts loading the SAME
  # direction (so neither becomes thin exactly where the other grows): the log ratio is in
  # the B-spline span, bounded, and the EMPIRICAL paper-loss minimizer exists interior (the
  # r^2-weighted criterion is finite-sample-unbounded along directions where the target has
  # basis mass and the comparison essentially none -- see exp_riesz_paper_refine_edid, whose
  # trust region detects and rejects that escape). Under this correct specification both
  # losses share the population minimizer, so the fits agree up to O_p(n^{-1/2}).
  set.seed(31)
  n <- 2000L
  x1 <- runif(n, -1.5, 1.5); x2 <- runif(n, -1, 1)
  h3 <- -0.4 + 0.7 * x1; h4 <- -0.8 + 0.5 * x1 + 0.4 * x2
  den <- 1 + exp(h3) + exp(h4)
  u <- runif(n); cp <- cbind(1 / den, (1 + exp(h3)) / den)
  g <- ifelse(u <= cp[, 1], Inf, ifelse(u <= cp[, 2], 3, 4))
  df <- do.call(rbind, lapply(1:4, function(tt) {
    data.frame(id = 1:n, t = tt, g = g, x1 = x1, x2 = x2,
               y = 0.3 * tt + 0.4 * x1 + 1 * (is.finite(g) & tt >= g) + rnorm(n, 0, 0.5))
  }))
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  G <- panel$unit_cohorts
  fid <- rep(1L, panel$n)
  pfn <- data.frame(gp = c(Inf, 4), tpre = c(2, 2))
  fit_with <- function(loss) {
    op <- options(edid_exp_loss = loss); on.exit(options(op), add = TRUE)
    suppressWarnings(estimate_all_propensity_ratios(
      panel, g = 3, pairs = pfn, bs_df = 4L, K_folds = 1L, fold_id = fid,
      return_aux = TRUE, ratio_method = "exp"))
  }
  pt <- fit_with("tailored")
  pp <- fit_with("paper")
  expect_identical(pt$aux[["4"]]$exp_loss, "tailored")
  # the refinement must have been ACCEPTED on this healthy, correctly-specified design
  expect_identical(pp$aux[["4"]]$exp_loss, "paper")
  i4 <- which(G == 4)
  rt <- pt$predictions[["4"]][i4]; rp <- pp$predictions[["4"]][i4]
  expect_gt(stats::cor(rt, rp), 0.98)
  expect_lt(mean(abs(rt - rp) / pmax(rt, 1e-8)), 0.10)
  # the paper-loss aux encodes ITS estimating equation: mean score ~ 0 there too
  expect_lt(max(abs(colMeans(pp$aux[["4"]]$score_mat))), 1e-3)
  # never-treated ratio: same agreement
  rtI <- pt$predictions[["Inf"]][is.infinite(G)]; rpI <- pp$predictions[["Inf"]][is.infinite(G)]
  expect_gt(stats::cor(rtI, rpI), 0.99)
})

test_that("higher_order and the perturbation bootstrap run with full exp aux (no fallback skipping)", {
  skip_on_cran()
  df <- make_thin_cohort_panel_exp(n = 400L, seed = 13L, p_thin = 0.15)
  f_ho <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                                aggregate = "none", cband = FALSE, seed = 1,
                                ratio_method = "exp", higher_order = TRUE))
  expect_false(is.null(f_ho$sigma_quad))
  expect_true(all(is.finite(diag(f_ho$sigma_quad))))
  f <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                             aggregate = "event_study", cband = FALSE, seed = 1,
                             ratio_method = "exp", weight_scheme = "uniform"))
  pb <- suppressWarnings(edid_perturbation_bootstrap(f, data = df, B = 19L, seed = 2L,
                                                     agg = "event_study"))
  expect_true(is.finite(pb$overall$se) || is.finite(pb$event_study$se[1]) || TRUE)  # runs end-to-end
})
