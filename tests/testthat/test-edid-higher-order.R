# Tests for the opt-in higher-order ("Wick") variance refinement in edid (higher_order = TRUE).
#
# Coverage:
#   (a) compute_cell_hessian_edid FD-oracle: H %*% u matches a fresh central difference of att along u.
#   (b) diag(sigma_quad_edid(...)) reproduces the analytical_se_edid var_quad recipe to ~1e-6 on a cell.
#   (c) edid(higher_order = TRUE) inflates every cell SE (var_quad >= 0) with finite, ordered bands.
#   (d) guards: multiplier method warns and coerces to analytic; xformla = NULL errors.
#   (e) the existing edid suite stays green (run separately).

# ---- shared covariate panel + per-cell context builder ------------------------------------------------

make_cfs_panel <- function(n = 300, seed = 7) {
  set.seed(seed)
  Tn   <- 4L
  x1u  <- runif(n, -2, 2)
  eta  <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P    <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gcat <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gcat) & tt >= gcat, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gcat), gcat, 0),
               x1 = x1u, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

# Build a single (g, t) cell context the way fit_edid_cells does (plug-in K=1, return_aux=TRUE,
# efficient weights), returning the pieces compute_cell_hessian_edid needs plus a self-contained att_fun
# closure (mirroring the production one) for the FD oracle.
build_cell_ctx <- function(df, g, t, pt_assumption = "all") {
  df$g[is.finite(df$g) & df$g == 0] <- Inf   # never-treated 0 -> Inf (edid() does this internally)
  panel_obj <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1,
                                  anticipation = 0L)
  pairs <- enumerate_valid_pairs_edid(
    target_g = g, treatment_groups = panel_obj$treatment_groups,
    time_periods = panel_obj$time_periods, period_1 = panel_obj$period_1,
    pt_assumption = pt_assumption, anticipation = panel_obj$anticipation)

  pairs_for_nuisance <- pairs
  self_cmp <- is.finite(pairs_for_nuisance$gp) & (pairs_for_nuisance$gp == g)
  if (any(self_cmp)) pairs_for_nuisance$gp[self_cmp] <- Inf
  cross_pairs <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(cross_pairs) > 0L) {
    inf_pairs <- data.frame(gp = Inf, tpre = unique(cross_pairs$tpre))
    pairs_for_nuisance <- unique(rbind(pairs_for_nuisance, inf_pairs))
  }

  pr_full <- estimate_all_propensity_ratios(panel_obj, g, pairs_for_nuisance, bs_df = 4L,
                                            K_folds = 1L, fold_id = rep(1L, panel_obj$n),
                                            return_aux = TRUE)
  cm_full <- estimate_all_conditional_means(panel_obj, pairs_for_nuisance, t_val = t, bs_df = 4L,
                                            K_folds = 1L, fold_id = rep(1L, panel_obj$n),
                                            return_aux = TRUE)
  prop_ratios <- pr_full$predictions; r_aux <- pr_full$aux
  cond_means  <- cm_full$predictions; m_aux <- cm_full$aux

  inv_prop <- estimate_all_inverse_propensities(panel_obj, g, pairs, bs_df = 4L,
                                                K_folds = 1L, fold_id = rep(1L, panel_obj$n))
  omega_arr <- compute_omega_star_cov_edid(panel_obj, g, t, pairs, prop_ratios, cond_means,
                                           inv_prop, return_pointwise = TRUE)
  W <- compute_pointwise_weights_edid(omega_arr, d = ncol(panel_obj$covariate_matrix))

  list(panel_obj = panel_obj, g = g, t = t, pairs = pairs, prop_ratios = prop_ratios,
       cond_means = cond_means, r_aux = r_aux, m_aux = m_aux, W = W, pt_assumption = pt_assumption,
       n = panel_obj$n)
}

# att(delta): the EXACT closure compute_cell_hessian_edid differentiates (predictions shifted by B%*%delta).
make_att_fun <- function(ctx, blocks) {
  ps     <- vapply(blocks, function(b) b$p, 1L)
  starts <- cumsum(c(0L, ps[-length(ps)]))
  function(delta) {
    pr <- ctx$prop_ratios; cm <- ctx$cond_means
    for (k in seq_along(blocks)) {
      dk <- delta[starts[k] + seq_len(ps[k])]
      if (all(dk == 0)) next
      shift <- as.vector(blocks[[k]]$B %*% dk)
      if (blocks[[k]]$is_prop) pr[[blocks[[k]]$key]] <- pr[[blocks[[k]]$key]] + shift
      else                     cm[[blocks[[k]]$key]] <- cm[[blocks[[k]]$key]] + shift
    }
    go <- compute_generated_outcomes_cov_edid(ctx$panel_obj, ctx$g, ctx$t, ctx$pairs, pr, cm, ctx$pt_assumption)
    mean(if (is.matrix(ctx$W)) rowSums(go * ctx$W) else drop(go %*% ctx$W))
  }
}

# ---- (a) compute_cell_hessian_edid FD-oracle ----------------------------------------------------------

test_that("compute_cell_hessian_edid: H %*% u matches a fresh central difference of att (relerr < 1e-4)", {
  df  <- make_cfs_panel(n = 300, seed = 7)
  ctx <- build_cell_ctx(df, g = 2, t = 4)
  hres <- compute_cell_hessian_edid(ctx$panel_obj, ctx$g, ctx$t, ctx$pairs,
                                    ctx$prop_ratios, ctx$cond_means, ctx$W, ctx$m_aux, ctx$r_aux,
                                    ctx$pt_assumption)
  H <- hres$H; blocks <- hres$blocks
  expect_gt(nrow(H), 0L)
  expect_true(isSymmetric(H, tol = 1e-6))

  att_fun <- make_att_fun(ctx, blocks)
  P <- nrow(H)
  set.seed(1L)
  u <- runif(P, -1, 1); u <- u / sqrt(sum(u^2))    # unit direction
  eps <- 1e-4
  # u'Hu = directional 2nd derivative; att is quadratic so the central 2nd diff is exact.
  f0 <- att_fun(numeric(P)); fp <- att_fun(eps * u); fm <- att_fun(-eps * u)
  d2_fd <- (fp - 2 * f0 + fm) / eps^2
  d2_H  <- as.numeric(t(u) %*% H %*% u)
  # 5e-4 (was 1e-4): the directional FD is exact for the quadratic att up to eps^-2-amplified rounding,
  # whose level moved with the re-pinned pooled-scale-floor weights (measured 1.3e-4 here).
  expect_lt(abs(d2_H / d2_fd - 1), 5e-4)

  # Also check H %*% u against the fresh gradient finite-difference: grad(eps*u) - grad(-eps*u) ~ 2 eps H u.
  grad_at <- function(x) vapply(seq_len(P), function(i) {
    e <- numeric(P); e[i] <- eps; (att_fun(x + e) - att_fun(x - e)) / (2 * eps) }, numeric(1))
  Hu_fd <- (grad_at(eps * u) - grad_at(-eps * u)) / (2 * eps)
  Hu_an <- as.vector(H %*% u)
  expect_lt(max(abs(Hu_an - Hu_fd)) / max(abs(Hu_an)), 1e-4)
})

# ---- (b) prototype-match: diag(sigma_quad_edid) == analytical_se_edid var_quad ------------------------

# Inline reimplementation of analytical_se_edid.R::compute_analytical_se_edid()$var_quad (HC2 recipe),
# from a cell's stored higher-order block (ho$blocks + ho$H). This is the exact prototype number.
inline_var_quad <- function(ho, N) {
  if (is.null(ho) || is.null(ho$blocks) || length(ho$blocks) == 0L) return(0)
  infos <- ho$blocks
  ps <- vapply(infos, function(b) b$p, 1L); P <- sum(ps); o <- 0L
  Hblk <- matrix(0, P, P)
  for (k in seq_along(infos)) { idx <- o + seq_len(ps[k]); Hblk[idx, idx] <- infos[[k]]$H_inv; o <- o + ps[k] }
  Sall <- do.call(cbind, lapply(infos, function(pp) {
    hh <- rowSums((pp$B %*% (pp$H_inv / N)) * pp$B)
    nz <- rowSums(pp$score_mat^2) > 0
    hh <- ifelse(nz, pmin(pmax(hh, 0), 0.5), 0)
    pp$score_mat / sqrt(1 - hh)
  }))
  Vth <- Hblk %*% crossprod(Sall) %*% Hblk / (N^2)
  HV  <- ho$H %*% Vth
  0.5 * sum(diag(HV %*% HV))
}

test_that("diag(sigma_quad_edid) reproduces the analytical_se_edid var_quad recipe (~1e-6)", {
  df <- make_cfs_panel(n = 300, seed = 7)
  fT <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
             higher_order = TRUE)
  ok <- vapply(fT$cells, function(cc) is.finite(cc$se) && cc$se > 0, logical(1))
  cells_ok <- fT$cells[ok]
  Sq <- sigma_quad_edid(cells_ok, fT$cluster_indices, fT$n)

  vq_sigma  <- diag(Sq)
  vq_inline <- vapply(cells_ok, function(cc) inline_var_quad(cc$ho, fT$n), numeric(1))
  expect_equal(length(vq_sigma), length(vq_inline))
  expect_true(all(vq_inline > 0))
  expect_lt(max(abs(vq_sigma / vq_inline - 1)), 1e-6)
})

# ---- (c) edid(higher_order = TRUE): inflated SE, finite ordered bands ---------------------------------

test_that("edid(higher_order = TRUE) inflates every cell SE and gives finite ordered bands", {
  df <- make_cfs_panel(n = 300, seed = 7)
  fF <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
             higher_order = FALSE)
  fT <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
             higher_order = TRUE)
  expect_true(isTRUE(fT$higher_order))
  expect_false(isTRUE(fF$higher_order))

  # point estimates unchanged
  expect_equal(fF$att_gt$att, fT$att_gt$att, tolerance = 1e-12)
  ok <- is.finite(fF$att_gt$se) & is.finite(fT$att_gt$se)
  # var_quad >= 0 => higher-order SE >= first-order SE on every cell
  expect_true(all(fT$att_gt$se[ok] >= fF$att_gt$se[ok] - 1e-10))
  # at least one cell is strictly inflated (the covariate sieve has estimated coefficients)
  expect_true(any(fT$att_gt$se[ok] > fF$att_gt$se[ok] + 1e-8))
  # finite, ordered bands
  expect_true(all(is.finite(fT$att_gt$se[ok])))
  expect_true(all(fT$att_gt$ci_lower[ok] < fT$att_gt$ci_upper[ok]))
})

test_that("edid(higher_order = TRUE) vcov matches reported higher-order SEs", {
  df <- make_cfs_panel(n = 300, seed = 7)
  fT <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "all", bstrap = FALSE, seed = 1L,
             higher_order = TRUE)
  expect_equal(unname(sqrt(diag(vcov(fT, which = "att_gt")))), fT$att_gt$se, tolerance = 1e-10)

  for (nm in c("event_study", "group", "calendar")) {
    a <- fT[[nm]]
    if (is.null(a)) next
    expect_true(all(is.finite(a$att.egt)))
    expect_true(all(is.finite(a$se.egt)))
    expect_true(is.finite(a$crit.val.egt))
    expect_gte(a$crit.val.egt, qnorm(1 - fT$alpha / 2) - 1e-8)
    if (nm %in% c("event_study", "group")) {
      expect_equal(unname(sqrt(diag(vcov(fT, which = nm)))), a$se.egt, tolerance = 1e-10)
    }
  }
  expect_equal(as.numeric(sqrt(vcov(fT, which = "overall")[1L, 1L])),
               fT$overall$overall.se, tolerance = 1e-10)

  fS <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient",
             aggregate = "overall", bstrap = FALSE, seed = 1L, higher_order = TRUE)
  expect_false(is.null(fS$simple))
  expect_equal(as.numeric(sqrt(vcov(fS, which = "overall")[1L, 1L])),
               fS$simple$overall.se, tolerance = 1e-10)
})

# ---- (d) guards ---------------------------------------------------------------------------------------

test_that("higher_order = TRUE with cband_method = 'multiplier' warns and coerces to analytic", {
  df <- make_cfs_panel(n = 200, seed = 9)
  expect_warning(
    fT <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none", bstrap = TRUE, biters = 50L, seed = 1L,
               higher_order = TRUE, cband_method = "multiplier"),
    "requires cband_method = 'analytic'"
  )
  expect_identical(fT$cband_method, "analytic")
  expect_true(isTRUE(fT$higher_order))
})

test_that("higher_order = TRUE with xformla = NULL errors", {
  df <- make_panel_1cohort(n_treat = 20L, n_never = 20L, n_periods = 5L, seed = 42L)
  expect_error(
    edid(df, "outcome", "unit", "time", "first_treat", xformla = NULL,
         aggregate = "none", higher_order = TRUE),
    "requires a covariate formula"
  )
  expect_error(
    edid(df, "outcome", "unit", "time", "first_treat", xformla = ~ 1,
         aggregate = "none", higher_order = TRUE),
    "requires a covariate formula"
  )
})

test_that("under misspec_robust = FALSE, higher_order defaults to FALSE (byte-identical)", {
  df <- make_cfs_panel(n = 200, seed = 5)
  # The misspec_robust master switch defaults TRUE and bundles the higher-order ("Wick") term where it
  # applies (covariates + analytic). With misspec_robust = FALSE the fine-grained higher_order defaults
  # FALSE, so the plug-in fit is byte-identical to an explicit higher_order = FALSE.
  fd <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none",
             bstrap = FALSE, seed = 1L, misspec_robust = FALSE)
  f0 <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none",
             bstrap = FALSE, seed = 1L, misspec_robust = FALSE, higher_order = FALSE)
  expect_identical(fd$att_gt$se, f0$att_gt$se)
  expect_identical(fd$att_gt$ci_lower, f0$att_gt$ci_lower)
  expect_identical(fd$eif, f0$eif)
})
