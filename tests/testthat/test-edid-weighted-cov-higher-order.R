library(testthat)

# ============================================================
# Phase 3c (higher_order under obs weights): the per-cell "Wick" Hessian
# of the OBS-WEIGHTED (Hajek) att(theta) must carry the marginal weight.
#
#   (A) FD oracle: u' H_analytic u matches a CLEAN central difference of the
#       weighted att(theta) along u (eps=1e-4) -- the weighted analytic cell
#       Hessian is correct (the production edid_hessian="fd" path is a coarser
#       eps=5e-2 oracle; this uses a clean eps for a tight check).
#   (B) constant weight column => higher_order SE == unweighted higher_order SE.
#   (C) higher_order MOVES the weighted SE vs the plug-in (real refinement).
# ============================================================

make_wh_panel <- function(n = 320, seed = 11) {
  set.seed(seed); Tn <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 0.5 * x1u
  P   <- exp(cbind(0, eta, 0.55 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.4 * x1u, 1)
  wu   <- exp(0.6 * rnorm(n)); wu <- wu / mean(wu)
  rows <- lapply(1:Tn, function(tt) {
    ht <- (tt - 1) * (0.4 * x1u); tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

# weighted per-cell context (mirrors fit_edid_cells; weightsname threaded into the panel)
build_wcell_ctx <- function(df, g, t, wt) {
  df$g[is.finite(df$g) & df$g == 0] <- Inf
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1, weightsname = wt, anticipation = 0L)
  pairs <- enumerate_valid_pairs_edid(g, panel$treatment_groups,
             panel$time_periods, panel$period_1, "all", panel$anticipation)
  pfn <- pairs; sc <- is.finite(pfn$gp) & pfn$gp == g; if (any(sc)) pfn$gp[sc] <- Inf
  cp <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
  if (nrow(cp)) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cp$tpre))))
  pr <- estimate_all_propensity_ratios(panel, g, pfn, bs_df = 4L, K_folds = 1L, fold_id = rep(1L, panel$n), return_aux = TRUE)
  cm <- estimate_all_conditional_means(panel, pfn, t_val = t, bs_df = 4L, K_folds = 1L, fold_id = rep(1L, panel$n), return_aux = TRUE)
  ip <- estimate_all_inverse_propensities(panel, g, pairs, bs_df = 4L, K_folds = 1L, fold_id = rep(1L, panel$n))
  oa <- compute_omega_star_cov_edid(panel, g, t, pairs, pr$predictions, cm$predictions, ip, return_pointwise = TRUE)
  W  <- compute_pointwise_weights_edid(oa, d = ncol(panel$covariate_matrix)); if (is.list(W)) W <- W$W
  list(panel = panel, g = g, t = t, pairs = pairs, pr = pr$predictions, cm = cm$predictions,
       r_aux = pr$aux, m_aux = cm$aux, W = W)
}

att_fun_w <- function(ctx, blocks) {
  ps <- vapply(blocks, function(b) b$p, 1L); starts <- cumsum(c(0L, ps[-length(ps)]))
  uw <- ctx$panel$unit_weights
  function(delta) {
    pr <- ctx$pr; cm <- ctx$cm
    for (k in seq_along(blocks)) {
      dk <- delta[starts[k] + seq_len(ps[k])]; if (all(dk == 0)) next
      shift <- as.vector(blocks[[k]]$B %*% dk)
      if (blocks[[k]]$is_prop) pr[[blocks[[k]]$key]] <- pr[[blocks[[k]]$key]] + shift
      else                     cm[[blocks[[k]]$key]] <- cm[[blocks[[k]]$key]] + shift
    }
    go <- compute_generated_outcomes_cov_edid(ctx$panel, ctx$g, ctx$t, ctx$pairs, pr, cm, "all")
    v  <- if (is.matrix(ctx$W)) rowSums(go * ctx$W) else drop(go %*% ctx$W)
    if (is.null(uw)) mean(v) else stats::weighted.mean(v, uw)
  }
}

test_that("(A) weighted analytic cell Hessian matches a clean central-FD of the weighted att", {
  ctx <- build_wcell_ctx(make_wh_panel(), g = 2, t = 4, wt = "w")
  hres <- compute_cell_hessian_edid(ctx$panel, ctx$g, ctx$t, ctx$pairs, ctx$pr, ctx$cm,
                                    ctx$W, ctx$m_aux, ctx$r_aux, "all")   # analytic default
  H <- hres$H; blocks <- hres$blocks
  expect_gt(nrow(H), 0L); expect_true(isSymmetric(H, tol = 1e-6))
  af <- att_fun_w(ctx, blocks); P <- nrow(H)
  set.seed(1L); u <- runif(P, -1, 1); u <- u / sqrt(sum(u^2)); eps <- 1e-4
  d2_fd <- (af(eps * u) - 2 * af(numeric(P)) + af(-eps * u)) / eps^2
  d2_H  <- as.numeric(t(u) %*% H %*% u)
  expect_lt(abs(d2_H / d2_fd - 1), 1e-3)                  # weighted analytic == clean weighted FD
})

test_that("(B) constant weight column => higher_order SE == unweighted higher_order SE", {
  df <- make_wh_panel(); df_c <- df; df_c$w <- 4
  ho <- function(d, wt) suppressWarnings(edid(d, "y","id","t","g", xformla = ~ x1, weightsname = wt,
    weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
    misspec_robust = FALSE, higher_order = TRUE))$att_gt$se
  su <- ho(df, NULL); sc <- ho(df_c, "w"); k <- is.finite(su) & is.finite(sc)
  expect_equal(sc[k], su[k], tolerance = 1e-8)
})

test_that("(C) higher_order moves the weighted-cov SE vs the plug-in", {
  df <- make_wh_panel()
  base <- function(ho) suppressWarnings(edid(df, "y","id","t","g", xformla = ~ x1, weightsname = "w",
    weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
    misspec_robust = FALSE, higher_order = ho))$att_gt$se
  s0 <- base(FALSE); s1 <- base(TRUE); k <- is.finite(s0) & is.finite(s1)
  expect_gt(max(abs(s1[k] - s0[k])), 1e-6)
})
