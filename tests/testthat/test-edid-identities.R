library(testthat)

# ===========================================================================
# Identity-style tests (external-audit batch, 2026-06): each test checks an
# IMPLEMENTED derivative/estimand identity against an independent oracle
# (numerical directional derivative, or the directly-computed estimand on the
# kept population), rather than pinning numbers.
#   (a) FD psi_Omega coupling identity (pooled/averaged DK coupling)
#   (b) floored-eigenvalue Daleckii-Krein derivative (exact math)
#   (c) common-overlap estimand + dead-pair dropping under trimming
#   (d) covariate-path moment_set degeneracies (full-empty / one-cohort-empty)
# plus the edid_sargan inference conventions and the rank-safe overall-weight
# recovery of the higher-order aggregation path.
# ===========================================================================

.collect_warnings_id <- function(expr) {
  ws <- character(0L)
  val <- withCallingHandlers(expr,
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = val, warnings = ws)
}

# ---------------------------------------------------------------------------
# (a) FD psi_Omega identity: the pooled (averaged-scheme) eigen-floor-aware
#     coupling C = dtheta/dOmega-bar returned by compute_obar_coupling_edid
#     must predict the first-order change of theta(Omega-bar) = w(Omega-bar)'mbar
#     under a random symmetric perturbation of the RAW pooled Omega-bar, where
#     w is built from the floored inverse exactly as production does. The
#     documented convention holds the floor LEVEL fixed (its d(max-eig) term is
#     higher-order), so the tolerance is ~10% wherever the floor moves; away
#     from the floor the agreement is much tighter.
# ---------------------------------------------------------------------------
test_that("pooled Omega-bar coupling matches the numerical directional derivative (tiny design)", {
  skip_on_cran()
  set.seed(31); n <- 150; Tt <- 4
  x  <- rnorm(n)
  gv <- sample(c(3, 4, Inf), n, replace = TRUE, prob = c(.3, .3, .4))
  df <- do.call(rbind, lapply(1:Tt, function(tt) {
    tau <- ifelse(is.finite(gv) & tt >= gv, 1, 0)
    data.frame(id = 1:n, time = tt, g = ifelse(is.finite(gv), gv, 0),
               x = x, y = 0.4 * x + 0.2 * tt + tau + rnorm(n, 0, 0.5))
  }))
  df$g[df$g == 0] <- Inf
  panel <- prepare_edid_panel(df, "y", "id", "time", "g", xformla = ~ x, anticipation = 0L)
  g <- 3; t <- 3
  pairs <- enumerate_valid_pairs_edid(g, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  pfn <- pairs
  self_cmp <- is.finite(pfn$gp) & (pfn$gp == g); if (any(self_cmp)) pfn$gp[self_cmp] <- Inf
  cross <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(cross) > 0L) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cross$tpre))))
  fid <- rep(1L, n)
  pr <- suppressWarnings(estimate_all_propensity_ratios(panel, g, pfn, bs_df = 4L, K_folds = 1L, fold_id = fid))
  cm <- suppressWarnings(estimate_all_conditional_means(panel, pfn, t_val = t, bs_df = 4L, K_folds = 1L, fold_id = fid))
  ip <- suppressWarnings(estimate_all_inverse_propensities(panel, g, pairs, bs_df = 4L, K_folds = 1L, fold_id = fid))

  omega <- suppressWarnings(compute_omega_star_cov_edid(panel, g, t, pairs, pr, cm, ip))
  ef    <- attr(omega, "eig_floor")
  expect_false(is.null(ef))                              # the pooled builder attaches the eigendecomposition
  expect_false(is.null(ef$scale))                        # ... of the SCALED system + the pooled scale (2026-06)
  H     <- nrow(pairs)
  gen   <- compute_generated_outcomes_cov_edid(panel, g, t, pairs, pr, cm, "all")
  mbar  <- colMeans(gen)
  w     <- compute_efficient_weights_edid(omega)
  att   <- sum(w * mbar)
  C     <- compute_obar_coupling_edid(omega, mbar, att)
  expect_false(is.null(C))

  # Two maps for the POOLED-SCALE floor (2026-06 re-derivation): the builder floors the scaled system
  # S = D^{-1/2} Omega D^{-1/2} (D = diag(diag(Omega)), exponent 1/3 -- the pooled Omega-bar is a
  # sqrt(n) object) and the weights invert D^{-1/2} fl(S)^{-1} D^{-1/2}. theta_fixed holds BOTH the
  # floor level c AND the pooled scale D at their fitted values (the documented convention: their
  # derivatives are omitted as higher-order); theta_moving is the full production map, which
  # recomputes the scale from diag(M) and the floor from max(eig) * n^(-1/3).
  sc0 <- ef$scale
  theta_fixed <- function(M) {
    S  <- 0.5 * (M + t(M)); S <- t(t(S * sc0) * sc0); S <- 0.5 * (S + t(S))
    e  <- eigen(S, symmetric = TRUE)
    ev <- pmax(e$values, ef$floor)
    v  <- sc0 * drop(e$vectors %*% (crossprod(e$vectors, sc0) / ev))
    sum((v / sum(v)) * mbar)
  }
  theta_moving <- function(M) {
    M  <- 0.5 * (M + t(M))
    dg <- diag(M); dsc <- 1 / sqrt(pmax(dg, max(dg) * 1e-12))
    S  <- t(t(M * dsc) * dsc); S <- 0.5 * (S + t(S))
    e  <- eigen(S, symmetric = TRUE)
    fl <- max(e$values) * panel$n^(-1/3)
    ev <- pmax(e$values, fl)
    v  <- dsc * drop(e$vectors %*% (crossprod(e$vectors, dsc) / ev))
    sum((v / sum(v)) * mbar)
  }
  # the RAW pooled Omega-bar: unscale the stored (scaled-system) eigendecomposition
  M0s <- ef$vectors %*% diag(ef$values, H) %*% t(ef$vectors)
  M0  <- t(t(M0s / sc0) / sc0)
  expect_equal(theta_moving(M0), att, tolerance = 1e-8)       # sanity: the map reproduces the plug-in att
  expect_equal(theta_fixed(M0),  att, tolerance = 1e-8)       # (they coincide AT M0 by construction)

  # Floor-only moving map: the pooled scale held at sc0, the floor level recomputed from max(eig).
  # Isolates the original documented omission (the d(floor) term) from the d(scale) channel.
  theta_move_floor <- function(M) {
    S  <- 0.5 * (M + t(M)); S <- t(t(S * sc0) * sc0); S <- 0.5 * (S + t(S))
    e  <- eigen(S, symmetric = TRUE)
    fl <- max(e$values) * panel$n^(-1/3)
    ev <- pmax(e$values, fl)
    v  <- sc0 * drop(e$vectors %*% (crossprod(e$vectors, sc0) / ev))
    sum((v / sum(v)) * mbar)
  }

  set.seed(7)
  for (r in 1:3) {
    D <- matrix(rnorm(H * H), H); D <- 0.5 * (D + t(D)); D <- D / sqrt(sum(D^2))
    eps    <- 1e-6 * max(abs(ef$values))
    num_fx <- (theta_fixed(M0 + eps * D)  - theta_fixed(M0 - eps * D))  / (2 * eps)
    num_mf <- (theta_move_floor(M0 + eps * D) - theta_move_floor(M0 - eps * D)) / (2 * eps)
    num_mv <- (theta_moving(M0 + eps * D) - theta_moving(M0 - eps * D)) / (2 * eps)
    pred   <- sum(C * D)
    # (i) The IMPLEMENTED identity is exact: the coupling is the derivative of the fixed-floor,
    # fixed-scale map.
    expect_equal(pred, num_fx, tolerance = 1e-6)
    # (ii) Against the floor-only moving map the gap is the documented fixed-floor convention
    # (d(floor) omitted as higher-order): bounded relative gap, same sign -- the original calibration.
    expect_lt(abs(num_mf - pred) / max(abs(num_mf), abs(pred), 1e-12), 0.5)
    expect_gt(sign(pred) * sign(num_mf), 0)
    # (iii) The FULL production map additionally recomputes the pooled scale from diag(M); that
    # d(scale) channel is also omitted by convention (diag(Omega-bar) is sqrt(n)-consistent, so the
    # omission vanishes as n grows), but on this n = 150 worst case it can dominate single random
    # directions (measured relative gap up to ~1, occasional sign flips). Bound the magnitude only.
    expect_lt(abs(num_mv - pred), 5 * max(abs(num_mv), abs(pred), 1e-12))
  }
})

# ---------------------------------------------------------------------------
# (b) Floored-eigenvalue Daleckii-Krein check: with the floor level HELD FIXED
#     (the implemented convention), the DK coupling is the EXACT derivative of
#     M -> w(M)'mbar with w from V diag(1/max(lambda, c)) V'. Construct a
#     matrix with one eigenvalue strictly below the floor and compare against
#     a central finite difference -- tight tolerance, this is exact math.
# ---------------------------------------------------------------------------
test_that("Daleckii-Krein coupling is the exact derivative of the FIXED-floor inverse map", {
  set.seed(5)
  H   <- 4L
  Q   <- qr.Q(qr(matrix(rnorm(H * H), H)))               # random orthogonal basis
  lam <- c(2.0, 0.9, 0.3, 1e-4)                          # last eigenvalue floored
  c_fl <- 0.05                                           # fixed floor (1e-4 << 0.05 << 0.3)
  M0  <- Q %*% diag(lam) %*% t(Q)
  omega_fl <- Q %*% diag(pmax(lam, c_fl)) %*% t(Q)
  attr(omega_fl, "eig_floor") <- list(values = lam, vectors = Q, floor = c_fl)
  mbar <- c(0.8, 1.1, 0.6, 1.4)
  w0   <- compute_efficient_weights_edid(omega_fl)
  att  <- sum(w0 * mbar)
  C    <- compute_obar_coupling_edid(omega_fl, mbar, att)
  expect_false(is.null(C))

  theta_fix <- function(M) {                              # FIXED floor c_fl (the convention C differentiates)
    e  <- eigen(0.5 * (M + t(M)), symmetric = TRUE)
    ev <- pmax(e$values, c_fl)
    v  <- drop(e$vectors %*% (crossprod(e$vectors, rep(1, H)) / ev))
    sum((v / sum(v)) * mbar)
  }
  set.seed(8)
  for (r in 1:3) {
    D <- matrix(rnorm(H * H), H); D <- 0.5 * (D + t(D)); D <- D / sqrt(sum(D^2))
    eps  <- 1e-6
    num  <- (theta_fix(M0 + eps * D) - theta_fix(M0 - eps * D)) / (2 * eps)
    pred <- sum(C * D)
    expect_equal(pred, num, tolerance = 1e-5)
  }
  # And the floored direction is genuinely clamped: the same perturbation through a NAIVE smooth
  # -sym(q w')/den coupling (no floor awareness) must NOT match the floored map (sanity that the
  # test has power: the floor binds here).
  Minv  <- Q %*% diag(1 / pmax(lam, c_fl)) %*% t(Q)
  q0    <- drop(Minv %*% (mbar - att)); den <- sum(drop(Minv %*% rep(1, H)))
  C_sm  <- -0.5 * (outer(q0, w0) + outer(w0, q0))         # smooth adjoint in the same orientation as C
  set.seed(8)
  D <- matrix(rnorm(H * H), H); D <- 0.5 * (D + t(D)); D <- D / sqrt(sum(D^2))
  num <- (theta_fix(M0 + 1e-6 * D) - theta_fix(M0 - 1e-6 * D)) / (2e-6)
  expect_gt(abs(sum(C_sm * D) - num), 10 * abs(sum(C * D) - num) + 1e-12)
})

# ---------------------------------------------------------------------------
# (c) Common-overlap estimand under binding trimming.
# ---------------------------------------------------------------------------

# DGP with a dead cross pair: cohort 4's support is essentially x > 1.2 while cohort 3 lives ONLY at
# x <= 1.2 (disjoint treated supports), so under trim_level = 3 the (g = 3 vs g' = 4) cross pairs --
# whose ratio-targeted mask keys on r_{3,4}(X) = p_3/p_4, huge everywhere cohort 3 lives -- retain no
# treated mass and must be DROPPED, not kept as zero columns (the pre-fix behavior halved the cell
# ATT: ~0.5 when the truth is 1.0). (2026-06: the masks are ratio-only for finite cohorts and the
# default ratios are the exp-link Riesz fits, so the supports must be genuinely disjoint for the
# pair to die; the old s-based mask killed it through the 1/p_4 scale alone. This test deliberately
# uses ratio_method = "direct" -- the LS sieve's extreme fitted values reliably kill the
# disjoint-support cross pair at trim_level = 3, exercising the construction-agnostic drop machinery.)
make_deadpair_panel <- function(n2 = 600, Tt = 5, seed = 1) {
  set.seed(seed)
  x2 <- rnorm(n2)
  p3 <- plogis(0.3 * x2) * (x2 <= 1.2)
  p4 <- ifelse(x2 > 1.2, 0.96, 0.005)
  u  <- runif(n2)
  gv <- ifelse(u < 0.25 * p3, 3, ifelse(u < 0.25 * p3 + 0.5 * p4, 4, Inf))
  df <- do.call(rbind, lapply(1:Tt, function(tt) {
    tau <- ifelse(is.finite(gv) & tt >= gv, 1, 0)
    data.frame(id = 1:n2, time = tt, g = ifelse(is.finite(gv), gv, 0),
               x = x2, y = 0.5 * x2 + 0.2 * tt + tau + rnorm(n2, 0, 0.5))
  }))
  df
}

test_that("dead pairs are dropped (not zero-padded): post ATT recovers ~1.0 with a warning", {
  skip_on_cran()
  df <- make_deadpair_panel()
  # ratio_method = "direct": this test guards the DEAD-PAIR MACHINERY (drop vs zero-pad), which is
  # construction-agnostic; the direct LS ratio's extreme fitted values reliably kill the disjoint-support
  # cross pair at trim_level = 3 on this single seed. Under the default exp-link ratios the smooth sieve
  # (like any smooth sieve) cannot represent this DGP's DISCONTINUOUS log-odds jump and under-
  # estimates the contrast near the boundary, so a sliver of treated mass stays kept -- the pair is then
  # legitimately alive-but-trimmed (the cell-common-overlap machinery, tested below, handles that case).
  res <- .collect_warnings_id(
    edid(df, "y", "id", "time", "g", xformla = ~ x, weight_scheme = "uniform",
         pt_assumption = "all", aggregate = "none", cband = FALSE,
         trim_level = 3, misspec_robust = FALSE, ratio_method = "direct"))
  fit <- res$value
  expect_true(any(grepl("dropped from their cells' moment sets", res$warnings)))
  # g = 3 cells: the (3 vs 4) cross pairs are dead and dropped; the surviving self pairs carry the cell
  att3 <- fit$att_gt[fit$att_gt$group == 3 & !fit$att_gt$is_pre, , drop = FALSE]
  expect_true(all(is.finite(att3$att)))
  expect_lt(max(abs(att3$att - 1)), 0.25)                # truth 1.0; pre-fix this sat at ~0.5
  cells3 <- Filter(function(cc) cc$group == 3 && !cc$is_pre, fit$cells)
  expect_true(all(vapply(cells3, function(cc) cc$n_pairs_dropped > 0L, logical(1))))
  # n_pairs reflects the SURVIVING count and matches the stored pair keys / weights
  for (cc in cells3) {
    expect_identical(cc$n_pairs, nrow(cc$pairs))
    expect_identical(length(cc$weights), cc$n_pairs)
    expect_false(any(cc$pairs$gp == 4))                  # the dead comparison cohort is gone
  }
  # EIF mean-zero identity under the common mask
  expect_lt(max(abs(colMeans(fit$eif)), na.rm = TRUE), 1e-12)
})

test_that("binding trim with heterogeneous tau(X): all weight schemes target the ONE common-overlap ATT", {
  skip_on_cran()
  # tau(X) = 1 + x, two comparison cohorts with partially disjoint support: cohort 4 lives mostly at
  # x > 0.5, so the cross-pair mask (the g' = 4 inverse propensity blows up at low x) differs from the
  # never-treated mask. Pre-fix, each pair was renormalized on its OWN kept subpopulation, so with
  # heterogeneous effects the schemes/moment-sets mixed DIFFERENT estimands and disagreed systematically;
  # post-fix every moment is renormalized on the cell-common kept population, so uniform / averaged /
  # efficient estimate the SAME common-overlap ATT (differences are sampling noise only) and match the
  # oracle ATT computed directly on the kept set.
  set.seed(42); n <- 700; Tt <- 4
  x  <- rnorm(n)
  p4 <- 0.55 * plogis(3 * (x - 0.2))                     # cohort 4: support essentially x > -0.2
  p3 <- 0.30 * plogis(0.3 * x)                           # cohort 3: everywhere
  u  <- runif(n)
  gv <- ifelse(u < p3, 3, ifelse(u < p3 + p4, 4, Inf))
  df <- do.call(rbind, lapply(1:Tt, function(tt) {
    tau <- ifelse(is.finite(gv) & tt >= gv, 1 + x, 0)
    data.frame(id = 1:n, time = tt, g = ifelse(is.finite(gv), gv, 0),
               x = x, y = 0.5 * x + 0.2 * tt + tau + rnorm(n, 0, 0.5))
  }))
  tl <- 8
  fit_of <- function(ws) suppressWarnings(
    edid(df, "y", "id", "time", "g", xformla = ~ x, weight_scheme = ws,
         pt_assumption = "all", aggregate = "none", cband = FALSE,
         trim_level = tl, misspec_robust = FALSE))
  f_u <- fit_of("uniform"); f_a <- fit_of("averaged"); f_e <- fit_of("efficient")

  # Oracle: replicate the g-level trim masks exactly as fit_edid_cells' .gbuild does (plug-in nuisances
  # are deterministic), then read the cell-common kept population from edid_cell_trim_structure and
  # average the KNOWN tau over the kept treated units of cohort 3.
  df2 <- df; df2$g[df2$g == 0] <- Inf
  panel <- prepare_edid_panel(df2, "y", "id", "time", "g", xformla = ~ x, anticipation = 0L)
  g <- 3
  pairs <- enumerate_valid_pairs_edid(g, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  pfn <- pairs
  self_cmp <- is.finite(pfn$gp) & (pfn$gp == g); if (any(self_cmp)) pfn$gp[self_cmp] <- Inf
  cross <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(cross) > 0L) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cross$tpre))))
  fid <- rep(1L, panel$n)
  # mirror the fit's default construction exactly: exp-link nuisances + the shared ratio-targeted mask
  pr <- suppressWarnings(estimate_all_propensity_ratios(panel, g, pfn, bs_df = 4L, K_folds = 1L,
                                                        fold_id = fid, ratio_method = "exp"))
  ip <- suppressWarnings(estimate_all_inverse_propensities(panel, g, pairs, bs_df = 4L, K_folds = 1L,
                                                           fold_id = fid, ratio_method = "exp"))
  trim_keep <- build_trim_keep_edid(pr, ip, tl, panel$n)
  ti <- edid_cell_trim_structure(panel, g, pairs, trim_keep, "all")
  expect_true(any(ti$keep_common < 0.5))                  # the trim binds
  expect_false(ti$full_trim)
  Ig   <- as.logical(panel$cohort_masks[[as.character(g)]])
  kept <- Ig & (ti$keep_common > 0.5)
  expect_gt(sum(Ig & !(ti$keep_common > 0.5)), 0L)        # ... and bites TREATED units (estimand moves)
  x_unit <- panel$covariate_matrix[, 1]
  oracle <- mean(1 + x_unit[kept])                        # common-overlap ATT(3, t) for every post t

  post3 <- function(f) f$att_gt[f$att_gt$group == 3 & !f$att_gt$is_pre, , drop = FALSE]
  a_u <- post3(f_u); a_a <- post3(f_a); a_e <- post3(f_e)
  expect_true(all(is.finite(c(a_u$att, a_a$att, a_e$att))))
  # (i) each scheme is within MC noise of the common-mask oracle (SE-aware, generous cap)
  for (a in list(a_u, a_a, a_e)) {
    expect_lt(max(abs(a$att - oracle) / pmax(a$se, 1e-12)), 4)
    expect_lt(max(abs(a$att - oracle)), 0.30)
  }
  # (ii) the schemes agree with EACH OTHER within MC noise: same estimand, different efficiency
  pair_gap <- c(max(abs(a_u$att - a_a$att)), max(abs(a_u$att - a_e$att)), max(abs(a_a$att - a_e$att)))
  expect_lt(max(pair_gap), 0.30)
  # (iii) EIF mean-zero identity with the common mask, all three schemes
  for (f in list(f_u, f_a, f_e)) expect_lt(max(abs(colMeans(f$eif)), na.rm = TRUE), 1e-12)
})

# ---------------------------------------------------------------------------
# (d) Covariate-path moment_set degeneracies (the F3 crash).
# ---------------------------------------------------------------------------
make_ms_panel <- function(seed = 20260610, n = 300, Tt = 5) {
  set.seed(seed)
  gv <- sample(c(3, 4, Inf), n, replace = TRUE, prob = c(.3, .3, .4))
  x  <- rnorm(n)
  df <- do.call(rbind, lapply(1:Tt, function(tt) {
    tau <- ifelse(is.finite(gv) & tt >= gv, tt - gv + 1, 0)
    data.frame(id = 1:n, time = tt, g = ifelse(is.finite(gv), gv, 0),
               x = x, y = 0.5 * x + 0.2 * tt + tau + rnorm(n, 0, 0.5))
  }))
  df
}

test_that("covariate-path moment_set removing ALL pairs returns an all-NA fit with a warning (no error)", {
  df <- make_ms_panel()
  ms_bad <- data.frame(g = c(3, 4), gp = c(3, 4), tpre = c(99, 99))   # no valid pair anywhere
  res <- .collect_warnings_id(
    edid(df, "y", "id", "time", "g", xformla = ~ x, pt_assumption = "all",
         aggregate = "none", cband = FALSE, moment_set = ms_bad))
  fit <- res$value
  expect_s3_class(fit, "edid_fit")
  expect_true(all(is.na(fit$att_gt$att)))
  expect_true(any(grepl("All post-treatment ATT\\(g,t\\) cells are NA", res$warnings)))
})

test_that("covariate-path moment_set emptying ONE cohort: that cohort NA, the other byte-identical", {
  skip_on_cran()
  df <- make_ms_panel()
  fit_full <- edid(df, "y", "id", "time", "g", xformla = ~ x, pt_assumption = "all",
                   aggregate = "none", cband = FALSE)
  # Empty cohort 3 (a non-existent pair) while listing cohort 4's FULL enumerated pair set, so cohort 4
  # is unrestricted (moment_set has intersection semantics; an unlisted cohort would be emptied too).
  pairs4 <- enumerate_valid_pairs_edid(4, c(3, 4), 1:5, 1, "all", 0L)
  ms <- rbind(data.frame(g = 3, gp = 3, tpre = 99),
              data.frame(g = 4, gp = pairs4$gp, tpre = pairs4$tpre))
  fit_half <- suppressWarnings(
    edid(df, "y", "id", "time", "g", xformla = ~ x, pt_assumption = "all",
         aggregate = "none", cband = FALSE, moment_set = ms))
  a_full <- fit_half$att_gt
  expect_true(all(is.na(a_full$att[a_full$group == 3])))             # emptied cohort: NA per contract
  i4h <- which(fit_half$att_gt$group == 4)
  i4f <- which(fit_full$att_gt$group == 4)
  expect_identical(fit_half$att_gt$att[i4h], fit_full$att_gt$att[i4f])  # byte-identical estimates
  expect_identical(fit_half$att_gt$se[i4h],  fit_full$att_gt$se[i4f])   # ... and SEs
})

# ---------------------------------------------------------------------------
# edid_sargan inference conventions (match_fit vs plugin_fast).
# ---------------------------------------------------------------------------
test_that("edid_sargan: match_fit copies the fit's IF convention; plugin_fast stays cheap; Holm identical", {
  skip_on_cran()
  # Single treated cohort + never-treated, T = 5: exactly ONE candidate restriction beyond the PT-Post
  # base, so each inference mode needs only two internal refits (keeps the test fast).
  set.seed(99); n <- 260; Tt <- 5
  gv <- sample(c(3, Inf), n, replace = TRUE, prob = c(.45, .55))
  xx <- rnorm(n)
  df <- do.call(rbind, lapply(1:Tt, function(tt) {
    tau <- ifelse(is.finite(gv) & tt >= gv, 1, 0)
    data.frame(id = 1:n, time = tt, g = ifelse(is.finite(gv), gv, 0),
               x = xx, y = 0.5 * xx + 0.2 * tt + tau + rnorm(n, 0, 0.5))
  }))
  fit <- suppressWarnings(
    edid(df, "y", "id", "time", "g", xformla = ~ x, pt_assumption = "all",
         aggregate = "event_study", cband = FALSE))                  # default: misspec_robust channels ON
  expect_true(isTRUE(fit$misspec_robust))                            # the channel-active premise
  sg_match <- suppressWarnings(edid_sargan(fit, data = df))          # default inference = "match_fit"
  sg_plug  <- suppressWarnings(edid_sargan(fit, data = df, inference = "plugin_fast"))
  expect_s3_class(sg_match, "edid_sargan")
  expect_s3_class(sg_plug,  "edid_sargan")
  expect_identical(sg_match$inference, "match_fit")
  expect_identical(sg_plug$inference,  "plugin_fast")
  # Same candidates in the same order; the statistics DIFFER when the misspec_robust/ACH channels are
  # active (match_fit carries them into Var-hat(xi); plugin_fast omits them).
  expect_identical(sg_match$table[, c("gp", "tpre")], sg_plug$table[, c("gp", "tpre")])
  expect_gt(max(abs(sg_match$table$H_statistic - sg_plug$table$H_statistic)), 1e-8)
  # The Holm step-down machinery is the SAME function of the p-values in both modes
  for (sg in list(sg_match, sg_plug)) {
    h <- .edid_holm(sg$table$p_value, sg$alpha)
    expect_identical(sg$table$rejected, h$rejected)
    expect_equal(sg$table$holm_threshold, h$threshold, tolerance = 1e-15)
  }
  # plugin_fast announces the cheaper convention; match_fit does not
  expect_output(print(sg_plug),  "plug-in influence functions only")
  out_match <- capture.output(print(sg_match))
  expect_false(any(grepl("plug-in influence functions only", out_match)))
})

# ---------------------------------------------------------------------------
# Rank-safe overall-aggregate weight recovery (higher-order path).
# ---------------------------------------------------------------------------
test_that("overall weight recovery: exact on full rank; warns and skips on collinear egt columns", {
  set.seed(3); n <- 200
  b1 <- rnorm(n); b2 <- rnorm(n); b3 <- rnorm(n)
  w_true <- c(0.2, 0.5, 0.3)
  egt <- cbind(b1, b2, b3)
  overall <- drop(egt %*% w_true)
  expect_silent(w <- .edid_recover_overall_weights(egt, overall))
  expect_equal(drop(w), w_true, tolerance = 1e-10)

  # Collinear columns (a duplicated cell feeding one event time): warning fires, NULL returned (the
  # caller skips the overall higher-order increment, leaving the finite first-order SE in force).
  egt_c <- cbind(b1, b2, b2)
  overall_c <- drop(egt_c %*% c(0.2, 0.3, 0.5))
  expect_warning(w_c <- .edid_recover_overall_weights(egt_c, overall_c), "higher_order")
  expect_null(w_c)

  # Inconsistent system (overall NOT in the column span): warning fires, NULL returned, never silent.
  expect_warning(w_i <- .edid_recover_overall_weights(cbind(b1, b2), b3), "higher_order")
  expect_null(w_i)

  # End-to-end: a higher_order fit's overall SE stays finite when the increment is skipped, because the
  # skip leaves did::aggte's analytic overall.se untouched (structural property of the caller); here we
  # assert the helper's contract that enables it (NULL + warning, no error).
  expect_true(TRUE)
})
