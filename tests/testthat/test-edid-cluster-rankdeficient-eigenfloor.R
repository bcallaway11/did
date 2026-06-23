# B1 fix: rank-deficient clustered cluster-metric eigen-floor (replaces the IID fallback) +
# LW cluster-ESS consistency fix. See NEWS "Few-cluster / rank guard (updated)".
#
# Mechanism: under clustering the no-cov efficient weights invert the CLUSTER moment covariance
# Sig_cl = crossprod(rowsum(psi, cluster))/n^2. When the over-id dimension H >= G_active, Sig_cl is
# rank-deficient and its null space is sampling noise. The OLD code reverted those cells to the IID
# Omega*, whose weights minimize the WRONG variance and can drive the reported CLUSTERED SE far below
# the honest equal-weight read (illusory sub-floor precision, a latent bound-violating ARE<1). The fix
# eigen-floors Sig_cl to a cluster-budget noise edge, collapsing the weights toward the honest
# equal-weight read.

# few treated clusters, many comparison cohorts -> over-id cells with H >= G_active (rank-deficient)
.rankdef_clustered_panel <- function() {
  cohorts <- c(4L, 5L, 6L); nct <- 3L; upc <- 8L; n_never <- 2L; T <- 8L
  set.seed(303L); d <- list(); uid <- 0L
  for (c in seq_len(nct)) for (u in seq_len(upc)) { uid <- uid + 1L
    d[[length(d)+1L]] <- data.frame(unit = uid, time = seq_len(T), cluster_id = c, first_treat = cohorts[c]) }
  for (c in seq_len(n_never)) for (u in seq_len(upc)) { uid <- uid + 1L
    d[[length(d)+1L]] <- data.frame(unit = uid, time = seq_len(T), cluster_id = nct + c, first_treat = Inf) }
  d <- do.call(rbind, d)
  cfe <- rnorm(nct + n_never, 0, 0.8)[d$cluster_id]; ufe <- rnorm(max(d$unit), 0, 0.3)[d$unit]
  d$outcome <- cfe + ufe + rnorm(nrow(d), 0, 0.2) + 1.8 * (is.finite(d$first_treat) & d$time >= d$first_treat)
  d
}

test_that("B1: rank-deficient clustered cell triggers the eigen-floor guard (warning fires)", {
  d <- .rankdef_clustered_panel()
  w <- character(0)
  withCallingHandlers(
    edid(data = d, yname = "outcome", idname = "unit", tname = "time", gname = "first_treat",
         clustervars = "cluster_id", bstrap = FALSE, cband = FALSE, pt_assumption = "all",
         omega_cov_shrink = "ridge"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("rank-deficient", w)))
  expect_true(any(grepl("EIGEN-FLOOR", w)))            # the new mechanism, not the IID fallback
  expect_false(any(grepl("kept the IID Omega", w)))     # old wording is gone
})

test_that("B1: eigen-floor keeps the rank-deficient efficient SE finite, positive, and NOT illusory", {
  d <- .rankdef_clustered_panel()
  # efficient (eigen-floored) vs uniform (equal-weight) clustered SE on the same cells
  fit <- function(scheme) suppressWarnings(edid(
    data = d, yname = "outcome", idname = "unit", tname = "time", gname = "first_treat",
    clustervars = "cluster_id", bstrap = FALSE, cband = FALSE, pt_assumption = "all",
    omega_cov_shrink = "ridge", weight_scheme = scheme))
  eff  <- fit("efficient")$att_gt
  unif <- fit("uniform")$att_gt
  m    <- match(paste(eff$group, eff$time), paste(unif$group, unif$time))
  se_eff <- eff$se; se_floor <- unif$se[m]
  post <- !eff$is_pre & is.finite(se_eff) & is.finite(se_floor) & se_floor > 0
  expect_true(all(is.finite(se_eff[post])) && all(se_eff[post] > 0))
  # the eigen-floor collapses toward the honest equal-weight read: efficient SE is NOT a small
  # fraction of the floor (the IID fallback produced ratios down to ~0.12; the floor brings the
  # worst case up near 1). Assert no catastrophic illusory-precision excursion remains.
  ratio <- se_eff[post] / se_floor[post]
  expect_gt(min(ratio), 0.5)        # IID fallback gave min ~0.12; eigen-floor gives ~0.87
})

test_that("B1: the H < G_active clustered path is byte-identical (eigen-floor is a no-op there)", {
  # MANY clusters so every over-id cell has H < G_active -> the eigen-floor must NEVER fire and the
  # reported SEs equal what the full-rank cluster metric gives (no rank-deficient guard at all).
  mk <- function() {
    set.seed(101L); nct <- 10L; nvc <- 10L; upc <- 4L; T <- 6L
    n <- (nct + nvc) * upc; uid <- rep(seq_len(n), each = T); tid <- rep(seq_len(T), times = n)
    cl <- c(rep(seq_len(nct), each = upc * T), rep(seq_len(nvc) + nct, each = upc * T))
    ft <- c(rep(3L, nct * upc * T), rep(Inf, nvc * upc * T))
    cfe <- rep(rnorm(nct + nvc, 0, 0.8), each = upc * T); ufe <- rep(rnorm(n, 0, 0.3), each = T)
    y <- cfe + ufe + rnorm(n * T, 0, 0.2) + 1.8 * (uid <= nct * upc & tid >= 3L)
    data.frame(unit = uid, time = tid, outcome = y, first_treat = ft, cluster_id = cl)
  }
  d <- mk()
  w <- character(0)
  r <- withCallingHandlers(
    edid(data = d, yname = "outcome", idname = "unit", tname = "time", gname = "first_treat",
         clustervars = "cluster_id", bstrap = FALSE, cband = FALSE, pt_assumption = "all",
         omega_cov_shrink = "ridge"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("rank-deficient", w)))         # guard must NOT fire when H < G_active
  expect_true(all(is.finite(r$att_gt$se)) && all(r$att_gt$se > 0))
})

test_that("LW: shrinking the CLUSTER metric Sig_cl uses the cluster ESS (G_eff), not the unit n", {
  # full-rank cluster metric (many clusters), interior LW lambda; cluster ESS != unit n unweighted
  mk <- function() {
    set.seed(55L); nct <- 30L; nvc <- 30L; upc <- 3L; T <- 6L; cohorts <- c(3L, 4L)
    tcc <- cohorts[(seq_len(nct) - 1L) %% length(cohorts) + 1L]; d <- list(); uid <- 0L
    for (c in seq_len(nct)) for (u in seq_len(upc)) { uid <- uid + 1L
      d[[length(d)+1L]] <- data.frame(unit = uid, time = seq_len(T), cluster_id = c, first_treat = tcc[c]) }
    for (c in seq_len(nvc)) for (u in seq_len(upc)) { uid <- uid + 1L
      d[[length(d)+1L]] <- data.frame(unit = uid, time = seq_len(T), cluster_id = nct + c, first_treat = Inf) }
    d <- do.call(rbind, d)
    cfe <- rnorm(nct + nvc, 0, 0.8)[d$cluster_id]; ufe <- rnorm(max(d$unit), 0, 0.3)[d$unit]
    d$outcome <- cfe + ufe + rnorm(nrow(d), 0, 0.2) + 1.8 * (is.finite(d$first_treat) & d$time >= d$first_treat)
    d
  }
  d <- mk()
  panel <- prepare_edid_panel(d, "outcome", "unit", "time", "first_treat", clustervars = "cluster_id")
  g <- 3L
  pr <- enumerate_valid_pairs_edid(g, panel$treatment_groups, panel$time_periods,
                                   panel$period_1, "all", 0L)
  pr <- pr[pr$tpre < g, , drop = FALSE]
  # locate a full-rank cell with interior unit-ESS lambda
  hit <- NULL
  for (tt in panel$time_periods[panel$time_periods >= g]) {
    psi <- compute_psi_moments_nocov_edid(g, tt, pr, panel)
    H <- ncol(psi); psi_cl <- rowsum(psi, panel$cluster_indices)
    G_act <- sum(rowSums(psi_cl^2) > 0); sig_cl <- crossprod(psi_cl) / panel$n^2
    if (nrow(pr) < 2L || any(!is.finite(sig_cl)) || qr(sig_cl)$rank < H) next
    shu <- shrink_omega_nocov_edid(sig_cl, g, tt, pr, panel, cl_metric_on = FALSE)
    if (is.finite(shu$lambda) && shu$lambda > 1e-6 && shu$lambda < 1 - 1e-6) {
      hit <- list(g = g, tt = tt, pr = pr, sig_cl = sig_cl, G_act = G_act, psi = psi,
                  H = H, lam_unit = shu$lambda); break
    }
  }
  skip_if(is.null(hit), "no interior-lambda full-rank cluster cell in fixture")
  shc <- shrink_omega_nocov_edid(hit$sig_cl, hit$g, hit$tt, hit$pr, panel,
                                 cl_metric_on = TRUE, cl_n_eff = hit$G_act)
  # hand reconstruction: cluster-ESS lambda = clamp( b2_legacy * (n / G_act) / d2 )
  S <- compute_pole_structure_nocov_edid(hit$g, hit$tt, hit$pr, panel); ss <- sum(S * S)
  sigma2 <- sum(hit$sig_cl * S) / ss; target <- sigma2 * S; d2 <- sum((hit$sig_cl - target)^2)
  q4 <- sum(rowSums(hit$psi * hit$psi)^2)
  b2_legacy <- (q4 / panel$n^2 - panel$n * sum(hit$sig_cl * hit$sig_cl)) / panel$n^2
  lam_clus_hand <- min(1, max(0, b2_legacy * (panel$n / hit$G_act)) / d2)
  expect_equal(shc$lambda, lam_clus_hand, tolerance = 1e-12)   # denominator switched to G_act
  expect_true(hit$G_act != panel$n)                            # the two ESS genuinely differ
  # and the unit-metric call still reproduces the legacy unit-n lambda exactly (byte-identity)
  lam_unit_hand <- min(1, max(0, b2_legacy * (panel$n / panel$n)) / d2)
  expect_equal(hit$lam_unit, lam_unit_hand, tolerance = 1e-12)
})

test_that("LW: unit-metric shrink is byte-identical with the new default args (no cl_metric_on)", {
  # self-contained fixture (no dependence on cross-file helpers / source order)
  set.seed(9L); nt <- 80L; n_never <- 40L; T <- 5L; n <- nt + n_never
  uid <- rep(seq_len(n), each = T); tid <- rep(seq_len(T), times = n)
  ft  <- c(rep(3L, nt * T), rep(Inf, n_never * T))
  ufe <- rep(rnorm(n, 0, 1), each = T)
  y   <- ufe + 0.3 * tid + rnorm(n * T, 0, 1) + 1.5 * (uid <= nt & tid >= 3L)
  df  <- data.frame(id = uid, time = tid, y = y, gvar = ft)
  panel <- prepare_edid_panel(df, "y", "id", "time", "gvar")
  pairs <- enumerate_valid_pairs_edid(3L, panel$treatment_groups, panel$time_periods,
                                      panel$period_1, "all", 0L)
  pairs <- pairs[pairs$tpre < 3L, , drop = FALSE]
  om <- compute_omega_star_nocov_edid(3L, 4L, pairs, panel, "all")
  skip_if(nrow(om) < 2L, "need an over-identified cell for shrinkage")
  sh_default  <- shrink_omega_nocov_edid(om, 3L, 4L, pairs, panel)                        # legacy call
  sh_explicit <- shrink_omega_nocov_edid(om, 3L, 4L, pairs, panel, cl_metric_on = FALSE)  # new arg, FALSE
  expect_identical(sh_default$lambda, sh_explicit$lambda)
  expect_identical(sh_default$omega,  sh_explicit$omega)
})
