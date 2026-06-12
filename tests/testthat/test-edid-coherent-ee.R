library(testthat)

# ---------------------------------------------------------------------------
# ratio_method = "coherent": FULL estimation-effect integration.
#
# The coherent engine fits one ridge-logistic per treated cohort, h_c =
# log(p_c/p_NT), on the cohort-vs-NT subsample over a shared basis (explicit
# unpenalized intercept). Consumed nuisances are r_{g,g'} = exp(h_g - h_{g'})
# and 1/p_c = (1 + sum_c' e^{h_c'})/e^{h_c}. The estimation-effect aux packs
# every entry against the JOINT stacked coefficient vector:
#   - joint Hessian: exactly BLOCK-DIAGONAL (each block's estimating equation
#     involves only its own beta_c), ridge included, sum-form Hmat_c =
#     sum_{S_c} mu(1-mu) psi psi' + diag(pen);
#   - joint score: stacked per-block scores psi (mu - d) on S_c (+ pen beta/n
#     distributed) -- the shared NT pool makes the BLOCKS' SCORES correlated,
#     which the stacking captures;
#   - chain rules: dr/dbeta = r(+psi in g, -psi in g'); d(1/p_c)/dbeta_{c'} =
#     (1/p_c)(p_{c'} - delta) psi (softmax Jacobian of the implied-share map).
# Entries sharing the system carry one coef_id; the perturbation bootstrap
# dedups its coefficient draws on it. FD oracles at the package's standard
# steps below.
# ---------------------------------------------------------------------------

# Same thin-comparison-cohort DGP family as test-edid-ratio-method.R.
make_thin_cohort_panel_coh <- function(n = 600L, seed = 42L, p_thin = 0.05) {
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

# Shared healthy fixture (no fallback blocks) + its joint system, fit once.
.coh_fix <- local({
  df <- make_thin_cohort_panel_coh(n = 800L, seed = 7L, p_thin = 0.15)
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1 + x2)
  fid <- rep(1L, panel$n)
  pfn <- data.frame(gp = c(Inf, 4), tpre = c(2, 2))      # target g = 3, cross g' = 4
  pr <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g = 3, pairs = pfn, bs_df = 4L, K_folds = 1L, fold_id = fid,
    return_aux = TRUE, ratio_method = "coherent"))
  ip <- suppressWarnings(estimate_all_inverse_propensities(
    panel, g = 3, pairs = data.frame(gp = c(3, 4), tpre = c(2, 2)), bs_df = 4L,
    K_folds = 1L, fold_id = fid, return_aux = TRUE, ratio_method = "coherent"))
  hs <- fit_coherent_logodds_edid(panel$covariate_matrix, panel$unit_cohorts,
                                  cohorts = panel$treatment_groups, bs_df = 4L,
                                  return_aux = TRUE)
  list(df = df, panel = panel, fid = fid, pfn = pfn, pr = pr, ip = ip, hs = hs)
})

test_that("coherent cross/invp entries carry FULL joint-system aux: mean-zero score, shared coef_id, bitwise-shared system", {
  fx <- .coh_fix; G <- fx$panel$unit_cohorts; n <- fx$panel$n
  a_r <- fx$pr$aux[["4"]]
  expect_false(isTRUE(a_r$is_fallback))
  expect_identical(a_r$link, "coherent")
  expect_true(all(c("B_test", "score_mat", "H_inv", "coef_id", "beta", "coh_layout") %in% names(a_r)))
  # the joint estimating equation holds at beta-hat: every stacked score column sums to ~0
  expect_lt(max(abs(colMeans(a_r$score_mat))), 1e-7)
  # the never-treated ratio entry stays the LS sieve aux (not the coherent system)
  expect_false(isTRUE(fx$pr$aux[["Inf"]]$is_fallback))
  expect_null(fx$pr$aux[["Inf"]]$coef_id)
  # invp entries: full aux, all-TRUE s_pos (the softmax map never clamps)
  ia <- attr(fx$ip, "aux")
  for (k in c("3", "4")) {
    expect_false(isTRUE(ia[[k]]$is_fallback))
    expect_identical(ia[[k]]$link, "coherent")
    expect_true(all(ia[[k]]$s_pos))
    expect_lt(max(abs(colMeans(ia[[k]]$score_mat))), 1e-7)
  }
  expect_false(isTRUE(ia[["Inf"]]$is_fallback))          # 1/p_NT keeps the LS aux
  expect_null(ia[["Inf"]]$coef_id)
  # ONE shared coefficient block: identical coef_id AND bitwise-identical joint score/H_inv
  # across the ratio and inverse-propensity dispatchers (the bootstrap-dedup premise)
  expect_identical(a_r$coef_id, ia[["4"]]$coef_id)
  expect_identical(a_r$score_mat, ia[["4"]]$score_mat)
  expect_identical(a_r$H_inv, ia[["4"]]$H_inv)
  # chain-rule packing of the ratio: +r psi in the g block, -r psi in the g' block
  lay <- a_r$coh_layout
  j3 <- lay$offsets[["3"]] + seq_len(lay$sizes[["3"]])
  j4 <- lay$offsets[["4"]] + seq_len(lay$sizes[["4"]])
  B3 <- fx$hs[["3"]]$B; B4 <- fx$hs[["4"]]$B
  expect_equal(unname(a_r$B_test[, j3]), unname( a_r$pred * B3), tolerance = 1e-12,
               check.attributes = FALSE)
  expect_equal(unname(a_r$B_test[, j4]), unname(-a_r$pred * B4), tolerance = 1e-12,
               check.attributes = FALSE)
})

test_that("FD oracle: the joint Hessian is the score Jacobian -- block-diagonal, ridge included, cross-blocks zero", {
  fx <- .coh_fix; G <- fx$panel$unit_cohorts; n <- fx$panel$n
  hs <- fx$hs
  ntm <- is.infinite(G)
  # stacked mean score as a function of (beta_3, beta_4), exact recompute of the fitted EE
  mean_score <- function(b3, b4) {
    out <- c()
    for (ck in c("3", "4")) {
      bk <- if (ck == "3") b3 else b4
      B  <- hs[[ck]]$B; sub <- ntm | (G == as.numeric(ck)); d <- as.numeric(G == as.numeric(ck))[sub]
      mu <- stats::plogis(drop(B[sub, , drop = FALSE] %*% bk))
      sc <- colSums(B[sub, , drop = FALSE] * (mu - d)) / n + hs[[ck]]$pen * bk / n
      out <- c(out, sc)
    }
    out
  }
  b3 <- hs[["3"]]$beta; b4 <- hs[["4"]]$beta
  q3 <- length(b3); q4 <- length(b4)
  s0 <- mean_score(b3, b4)
  expect_lt(max(abs(s0)), 1e-7)                          # EE holds at beta-hat (ridge folded in)
  eps <- 1e-6
  H_fd <- vapply(seq_len(q3 + q4), function(j) {
    e3 <- numeric(q3); e4 <- numeric(q4)
    if (j <= q3) e3[j] <- eps else e4[j - q3] <- eps
    (mean_score(b3 + e3, b4 + e4) - s0) / eps
  }, numeric(q3 + q4))
  H_an <- matrix(0, q3 + q4, q3 + q4)
  H_an[seq_len(q3), seq_len(q3)]           <- hs[["3"]]$Hmat / n
  H_an[q3 + seq_len(q4), q3 + seq_len(q4)] <- hs[["4"]]$Hmat / n
  sc_H <- 1 + max(abs(H_an))
  expect_lt(max(abs(H_fd - H_an)), 1e-4 * sc_H)          # blocks match the analytic Hessian
  # the CROSS blocks of the FD Jacobian are exactly zero (block-diagonality of the joint H)
  expect_lt(max(abs(H_fd[seq_len(q3), q3 + seq_len(q4)])), 1e-10 * sc_H)
  expect_lt(max(abs(H_fd[q3 + seq_len(q4), seq_len(q3)])), 1e-10 * sc_H)
  # and H_inv inverts the average joint Hessian on the packed space: H_inv (Hmat/n) ~ I
  a_r <- fx$pr$aux[["4"]]
  HH <- a_r$H_inv %*% H_an
  expect_equal(diag(HH), rep(1, q3 + q4), tolerance = 1e-6)
})

test_that("FD oracle: B_test is the exact coefficient Jacobian of r_{g,g'} and 1/p_c (chain rule + softmax map)", {
  fx <- .coh_fix; G <- fx$panel$unit_cohorts; n <- fx$panel$n
  hs <- fx$hs
  a_r <- fx$pr$aux[["4"]]; lay <- a_r$coh_layout
  eps <- 1e-6
  # exact maps at perturbed coefficients
  h_of  <- function(ck, bk) drop(hs[[ck]]$B %*% bk)
  betas <- list("3" = hs[["3"]]$beta, "4" = hs[["4"]]$beta)
  r_of  <- function(b) exp(h_of("3", b[["3"]]) - h_of("4", b[["4"]]))
  s_of  <- function(b, ck) {                              # 1/p_ck from the implied shares
    eh <- cbind("3" = exp(h_of("3", b[["3"]])), "4" = exp(h_of("4", b[["4"]])))
    (1 + rowSums(eh)) / eh[, ck]
  }
  for (ck in c("3", "4")) {
    jj <- lay$offsets[[ck]] + seq_len(lay$sizes[[ck]])
    for (j in sample(seq_along(jj), 3L)) {                # spot-check 3 coordinates per block
      bp <- betas; bp[[ck]][j] <- bp[[ck]][j] + eps
      # ratio entry r_{3,4}
      fd <- (r_of(bp) - a_r$pred) / eps
      an <- a_r$B_test[, jj[j]]
      expect_lt(max(abs(fd - an)), 1e-4 * (1 + max(abs(an))))
      # inverse propensities 1/p_3 and 1/p_4 (the softmax Jacobian, own + cross blocks)
      for (sk in c("3", "4")) {
        sa <- attr(fx$ip, "aux")[[sk]]
        fd_s <- (s_of(bp, sk) - sa$s_hat) / eps
        an_s <- sa$B_test[, jj[j]]
        expect_lt(max(abs(fd_s - an_s)), 1e-4 * (1 + max(abs(an_s))))
      }
    }
  }
})

test_that("cross-target reciprocity: r_{g,g'} and r_{g',g} entries couple exactly through the shared blocks", {
  # The bootstrap's shared-draw dedup premise in geometric form: with ONE coefficient draw
  # delta, the implied first-order shifts satisfy shift_{3,4}/r_{3,4} = -shift_{4,3}/r_{4,3}
  # for EVERY delta, i.e. B_test_{3,4}/r_{3,4} = -B_test_{4,3}/r_{4,3} columnwise.
  fx <- .coh_fix
  pr43 <- suppressWarnings(estimate_all_propensity_ratios(
    fx$panel, g = 4, pairs = data.frame(gp = c(Inf, 3), tpre = c(2, 2)), bs_df = 4L,
    K_folds = 1L, fold_id = fx$fid, return_aux = TRUE, ratio_method = "coherent"))
  a34 <- fx$pr$aux[["4"]]; a43 <- pr43$aux[["3"]]
  expect_identical(a34$coef_id, a43$coef_id)             # same shared coefficient block
  expect_identical(a34$score_mat, a43$score_mat)         # bitwise: same joint system across target cohorts
  expect_equal(a34$B_test / a34$pred, -(a43$B_test / a43$pred), tolerance = 1e-12)
  expect_equal(a34$pred * a43$pred, rep(1, fx$panel$n), tolerance = 1e-12)  # exact coherence
})

test_that("FD oracle: the analytic ACH Gamma equals the beta-space derivative through the exact coherent map", {
  df <- .coh_fix$df
  panel <- .coh_fix$panel; n <- panel$n
  g <- 3; t <- 4
  pairs <- enumerate_valid_pairs_edid(g, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  pfn <- pairs; sc <- is.finite(pfn$gp) & pfn$gp == g; pfn$gp[sc] <- Inf
  cr <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(cr)) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cr$tpre))))
  fid <- rep(1L, n)
  pr <- suppressWarnings(estimate_all_propensity_ratios(
    panel, g, pfn, bs_df = 4L, K_folds = 1L, fold_id = fid, return_aux = TRUE,
    ratio_method = "coherent"))
  cm <- suppressWarnings(estimate_all_conditional_means(panel, pfn, t_val = t, bs_df = 4L,
                                                        K_folds = 1L, fold_id = fid, return_aux = TRUE))
  prop_ratios <- pr$predictions
  # deliberately MISSPECIFIED conditional means (zeroed) so the r-channel Gamma is far from
  # orthogonal-zero; FD-vs-analytic equality holds for ANY inputs
  cond_means <- lapply(cm$predictions, function(v) 0 * v)
  H <- nrow(pairs); w <- rep(1 / H, H)
  key <- "4"
  a <- pr$aux[[key]]
  skip_if(isTRUE(a$is_fallback), "cross aux fell back on this fixture")

  corr_an <- compute_ach_correction_cov_edid(panel, g, t, pairs, prop_ratios, cond_means,
                                             w, m_aux = list(), r_aux = pr$aux[key])
  # oracle: Gamma by FD in the STACKED COEFFICIENT space through the exact exp(h_g - h_gp) map,
  # perturbing THIS entry's prediction only (the entry's own dependence on the shared beta):
  # r(beta + eps e_j) = r * exp(eps * B_test[, j] / r) exactly.
  wmom <- function(prr) {
    go <- compute_generated_outcomes_cov_edid(panel, g, t, pairs, prr, cond_means, "all")
    mean(drop(go %*% w))
  }
  m0 <- wmom(prop_ratios)
  eps <- 1e-6
  P <- ncol(a$B_test)
  Gamma_fd <- vapply(seq_len(P), function(j) {
    prr <- prop_ratios
    prr[[key]] <- prop_ratios[[key]] * exp(eps * a$B_test[, j] / a$pred)
    (wmom(prr) - m0) / eps
  }, numeric(1))
  corr_fd <- as.vector(a$score_mat %*% drop(a$H_inv %*% Gamma_fd))
  expect_gt(max(abs(corr_an)), 1e-8)                     # the channel is genuinely non-zero here
  expect_equal(corr_an, corr_fd, tolerance = 1e-4)
})

test_that("package-level: analytic ACH reproduces the FD oracle under ratio_method='coherent'; corrections ACTIVE", {
  df <- make_thin_cohort_panel_coh(n = 500L, seed = 21L, p_thin = 0.15)
  run <- function(ach, ee) {
    op <- options(edid_ach = ach); on.exit(options(op), add = TRUE)
    suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                          weight_scheme = "efficient", aggregate = "none", cband = FALSE,
                          seed = 1L, misspec_robust = FALSE, estimation_effect = ee,
                          ratio_method = "coherent"))
  }
  f_noee <- run("analytic", FALSE)
  f_an   <- run("analytic", TRUE)
  f_fd   <- run("fd",       TRUE)
  se_noee <- f_noee$att_gt$se; se_an <- f_an$att_gt$se; se_fd <- f_fd$att_gt$se
  ok <- is.finite(se_noee) & is.finite(se_an) & is.finite(se_fd)
  skip_if(sum(ok) < 2L, "too few non-degenerate cells")
  expect_gt(max(abs(se_an[ok] - se_noee[ok])), 1e-6)     # cross channels ACTIVE (not fallback-skipped)
  expect_equal(se_an[ok], se_fd[ok], tolerance = 1e-5)
  # the corrections touch the SE only: the point estimates are byte-identical across all runs
  expect_identical(f_noee$att_gt$att, f_an$att_gt$att)
  expect_identical(f_noee$att_gt$att, f_fd$att_gt$att)
})

test_that("inv-p weight channel: the coherent analytic correction matches the package FD oracle (exp-test calibration)", {
  # Same calibration as the exp channel's test: the analytic Gamma uses the Daleckii-Krein
  # coupling of the FIXED-floor inverse map while the FD oracle re-solves the MOVING-floor map;
  # exact-match witnesses where the floor is inactive + aggregate agreement no worse than the
  # validated linear ("direct") channel on the same fixture.
  df <- make_thin_cohort_panel_coh(n = 500L, seed = 11L, p_thin = 0.15)
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
  st_coh <- stats_of(run_acc("coherent"))
  st_lin <- stats_of(run_acc("direct"))
  skip_if(is.null(st_coh) || nrow(st_coh) < 3L, "too few corrected cells")
  expect_true(all(is.finite(st_coh$rel)))
  expect_lt(min(st_coh$rel), 1e-3)                       # exact-wiring witness (floor inactive)
  expect_gt(stats::median(st_coh$cor, na.rm = TRUE), 0.95)
  skip_if(is.null(st_lin) || nrow(st_lin) < 3L, "too few linear cells to benchmark")
  expect_lte(stats::median(st_coh$rel), 2 * max(stats::median(st_lin$rel), 1e-3))
})

test_that("higher_order runs with full coherent aux (shared blocks dedup'd in sigma_quad); bootstrap dedups its draws", {
  skip_on_cran()
  df <- make_thin_cohort_panel_coh(n = 400L, seed = 13L, p_thin = 0.15)
  f_ho <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                                aggregate = "none", cband = FALSE, seed = 1,
                                ratio_method = "coherent", higher_order = TRUE))
  expect_false(is.null(f_ho$sigma_quad))
  expect_true(all(is.finite(diag(f_ho$sigma_quad))))
  expect_true(all(diag(f_ho$sigma_quad) >= -1e-12))      # PSD diagonal up to roundoff
  # the cross r-blocks genuinely enter: at least one cell carries a coherent (P-wide) r block
  Ps <- unlist(lapply(f_ho$cells, function(cc) {
    if (is.null(cc$ho) || !length(cc$ho$blocks)) return(integer(0))
    vapply(cc$ho$blocks, function(b) if (isTRUE(b$is_prop)) b$p else NA_integer_, 1L)
  }))
  # joint-system blocks present: P = sum_c q_c = 2 * 8 = 16 here, vs the LS "Inf" block's p = 7
  expect_gt(sum(!is.na(Ps) & Ps > 10L), 0L)
  # perturbation bootstrap: runs end-to-end with the coherent entries perturbable,
  # shares ONE draw per system (coef_id dedup), and is seed-reproducible
  f <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2,
                             aggregate = "event_study", cband = FALSE, seed = 1,
                             ratio_method = "coherent", weight_scheme = "uniform"))
  p1 <- suppressWarnings(edid_perturbation_bootstrap(f, data = df, B = 19L, seed = 2L,
                                                     agg = "event_study"))
  p2 <- suppressWarnings(edid_perturbation_bootstrap(f, data = df, B = 19L, seed = 2L,
                                                     agg = "event_study"))
  expect_true(any(is.finite(p1$att_gt$se_pert)))
  expect_identical(p1$att_gt$se_pert, p2$att_gt$se_pert)
})

test_that("invariances: PT-Post-X bitwise across ratio_method; no-X untouched; direct/exp paths untouched by the aux packing", {
  df <- make_thin_cohort_panel_coh(n = 400L, seed = 3L, p_thin = 0.10)
  # PT-Post-X: single never-treated moment, H = 1 -- coherent (full EE) == direct, bitwise
  fp_c <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2, pt_assumption = "post",
                                aggregate = "none", cband = FALSE))
  fp_d <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1 + x2, pt_assumption = "post",
                                aggregate = "none", cband = FALSE, ratio_method = "direct"))
  expect_identical(fp_c$att_gt$att, fp_d$att_gt$att)
  expect_identical(fp_c$att_gt$se,  fp_d$att_gt$se)
  # no-X: bitwise invariant to the coherent aux machinery
  fn_c <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none", cband = FALSE))
  fn_d <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none", cband = FALSE,
                                ratio_method = "coherent"))
  expect_identical(fn_c$att_gt$att, fn_d$att_gt$att)
  expect_identical(fn_c$att_gt$se,  fn_d$att_gt$se)
})
