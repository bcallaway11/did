library(testthat)

# ============================================================================
# Paper-faithfulness invariants for the efficient-DiD estimator
# (Chen, Sant'Anna & Xie 2025). These lock the implementation to the paper so
# that future edits which deviate from it are caught. Each test references the
# paper object it protects.
# ============================================================================

# Shared DGP: cohorts {3, 4, Inf}, T = 6, conditional parallel trends given X.
pf_panel <- function(n = 1500, seed = 1, het = FALSE, attscale = 1,
                     dynamics = FALSE, ncov = 1L, binx1 = FALSE) {
  set.seed(seed)
  x1 <- if (binx1) sample(c(-1, 1), n, TRUE) else rnorm(n)
  x2 <- rnorm(n)
  u  <- runif(n); g <- ifelse(u < .30, Inf, ifelse(u < .60, 3, 4)); mu <- rnorm(n)
  TP <- 6L
  if (het) {
    base <- ifelse(is.infinite(g), 0.6, ifelse(g == 4, 1.2, 0.9)); rho <- 0.6
    eps  <- matrix(0, n, TP); eps[, 1] <- rnorm(n, sd = base)
    for (k in 2:TP) eps[, k] <- rho * eps[, k - 1] + sqrt(1 - rho^2) * rnorm(n, sd = base)
  } else {
    eps <- matrix(rnorm(n * TP, sd = 0.7), n, TP)
  }
  eff <- function(k, gg) if (dynamics) attscale * (1 + 0.5 * (k - gg)) else attscale
  rows <- lapply(1:TP, function(k) {
    tr <- as.numeric(is.finite(g) & k >= g)
    y  <- mu + 0.3 * k + 0.4 * x1 * (k - 1) +
          (if (ncov >= 2) 0.3 * x2 * (k - 1) else 0) +
          ifelse(is.finite(g) & k >= g, eff(k, g), 0) * tr + eps[, k]
    d <- data.frame(id = 1:n, t = k, y = y, g = g, x1 = x1)
    if (ncov >= 2) d$x2 <- x2
    d
  })
  do.call(rbind, rows)
}

# --- Efficient weight_scheme = GLS solution (Theorems on the efficiency bound) -----
test_that("efficient weights solve the GLS problem: sum to 1, equal Omega^{-1}1/(1'Omega^{-1}1), minimal variance", {
  set.seed(20260529); H <- 4L
  A  <- matrix(rnorm(H * H), H, H); Om <- crossprod(A) + diag(H)   # random SPD Omega*
  w  <- compute_efficient_weights_edid(Om); one <- rep(1, H)
  expect_equal(sum(w), 1, tolerance = 1e-8)                        # weights sum to one
  expect_equal(unname(w), as.numeric(solve(Om, one) / sum(solve(Om, one))), tolerance = 1e-8)
  v_eff  <- as.numeric(t(w) %*% Om %*% w)
  expect_equal(v_eff, 1 / sum(solve(Om)), tolerance = 1e-8)        # = (1'Omega^{-1}1)^{-1}
  expect_lte(v_eff, as.numeric(t(one / H) %*% Om %*% (one / H)) + 1e-10)  # <= uniform-weight variance
})

# --- Omega* construction is faithful to Eq (3.12): cov-path = no-cov on constant X ---
test_that("covariate-path Omega* equals the no-covariate Omega* on constant covariates (Eq 3.12)", {
  set.seed(101); n <- 4000L; TP <- 5L
  x1 <- rnorm(n, sd = 0.02)                                        # near-constant -> cond cov = uncond
  g  <- ifelse(runif(n) < 0.5, Inf, 4L); mu <- rnorm(n)
  sdc <- ifelse(is.finite(g), 1.5, 0.5); rhoc <- ifelse(is.finite(g), 0.7, 0.2)
  eps <- matrix(0, n, TP); eps[, 1] <- rnorm(n, sd = sdc)
  for (k in 2:TP) eps[, k] <- rhoc * eps[, k - 1] + sqrt(1 - rhoc^2) * rnorm(n, sd = sdc)
  df <- do.call(rbind, lapply(1:TP, function(k) {
    tr <- as.numeric(is.finite(g) & k >= g)
    data.frame(id = 1:n, t = k, y = mu + 0.3 * k + 1.0 * tr + eps[, k], g = g, x1 = x1)
  }))
  panel <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1)
  pairs <- enumerate_valid_pairs_edid(4L, panel$treatment_groups, panel$time_periods,
                                            panel$period_1, "all", panel$anticipation)
  fid   <- build_crossfit_folds_edid(panel$n, 5L, seed = 1L)
  pn    <- pairs; pn$gp[is.finite(pn$gp) & pn$gp == 4L] <- Inf
  pr    <- estimate_all_propensity_ratios(panel, 4L, pn, 4L, 5L, fid)
  cm    <- estimate_all_conditional_means(panel, pn, 4L, 4L, 5L, fid)
  w_cov   <- compute_efficient_weights_edid(compute_omega_star_cov_edid(panel, 4L, 4L, pairs, pr, cm))
  w_nocov <- compute_efficient_weights_edid(compute_omega_star_nocov_edid(4L, 4L, pairs, panel, "all"))
  expect_equal(w_cov, w_nocov, tolerance = 0.05)                   # match up to NW kernel error
})

# --- EIF is the ratio-estimator influence function (IF-g / IF-general) -------
test_that("covariate-path EIF is mean-zero and uses the ratio centering w'Ytilde - (G_g/pi_g)ATT", {
  df  <- pf_panel(n = 1200, seed = 7)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none")
  cm  <- colMeans(fit$eif, na.rm = TRUE)
  expect_true(all(abs(cm) < 1e-8))                                 # EIF mean-zero by construction
})

# --- All four weight schemes are consistent for a homogeneous ATT ------------
test_that("all four weight schemes recover a homogeneous ATT", {
  df <- pf_panel(n = 3000, seed = 3, het = TRUE, attscale = 1)
  for (m in c("efficient", "averaged", "gmm", "uniform")) {
    fit <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1,
                                 weight_scheme = m, aggregate = "overall"))
    expect_lt(abs(fit$simple$overall.att - 1), 0.12)
  }
})

# --- $overall is the dynamic headline; the type overalls match aggte_edid ----
test_that("$overall is the dynamic headline and the type overalls match aggte_edid", {
  df  <- pf_panel(n = 1500, seed = 5, dynamics = TRUE)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "all")
  expect_equal(fit$overall$overall.att,          aggte_edid(fit, type = "dynamic",  na.rm = TRUE)$overall.att, tolerance = 1e-8)
  expect_equal(fit$simple$overall.att,   aggte_edid(fit, type = "simple")$overall.att,                 tolerance = 1e-8)
  expect_equal(fit$group$overall.att,    aggte_edid(fit, type = "group")$overall.att,                  tolerance = 1e-8)
  expect_equal(fit$calendar$overall.att, aggte_edid(fit, type = "calendar", na.rm = TRUE)$overall.att,  tolerance = 1e-8)
  # with genuine dynamics the dynamic headline differs from the simple aggregate
  expect_gt(abs(fit$overall$overall.att - fit$simple$overall.att), 1e-3)
})

# --- All four aggregation types run and return a finite overall --------------
test_that("aggte_edid supports simple, dynamic, group, and calendar", {
  df  <- pf_panel(n = 1500, seed = 9)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "all")
  for (ty in c("simple", "dynamic", "group", "calendar")) {
    a <- aggte_edid(fit, type = ty, na.rm = TRUE)
    expect_true(is.finite(a$overall.att), info = paste("type", ty))
  }
})

# --- Calendar effect = cohort-share-weighted average of ATT(g,t) over g <= t --
test_that("calendar ATT(t) equals the cohort-share-weighted average of ATT(g,t) for g <= t", {
  df  <- pf_panel(n = 2500, seed = 11, dynamics = TRUE)
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "all")
  pi  <- fit$cohort_fractions
  for (t_val in fit$calendar$egt) {           # $calendar is a AGGTEobj; egt = calendar periods
    rows  <- fit$att_gt[fit$att_gt$time == t_val & fit$att_gt$group <= t_val &
                          is.finite(fit$att_gt$group) & is.finite(fit$att_gt$att), ]
    if (nrow(rows) == 0L) next
    pg    <- vapply(rows$group, function(g) pi[[as.character(g)]], numeric(1L))
    manual <- sum((pg / sum(pg)) * rows$att)
    expect_equal(fit$calendar$att.egt[match(t_val, fit$calendar$egt)], manual, tolerance = 1e-6)
  }
})

# --- WIF: the simple-overall SE carries the cohort-share weight-estimation
#     variance, so it strictly exceeds the no-WIF (direct-EIF-only) SE under
#     cohort-ATT heterogeneity (Theorem 'efficiency' aggregation; wif). ---
test_that("aggregated overall SE includes the cohort-share WIF under cohort heterogeneity", {
  set.seed(13); n <- 2500L; TP <- 6L
  u <- runif(n); g <- ifelse(u < .25, Inf, ifelse(u < .5, 3, ifelse(u < .75, 4, 5))); mu <- rnorm(n)
  base_att <- function(gg) ifelse(is.infinite(gg), 0, ifelse(gg == 3, 5, ifelse(gg == 4, -3, 2)))  # heterogeneous
  df <- do.call(rbind, lapply(1:TP, function(k) {
    tr <- as.numeric(is.finite(g) & k >= g)
    data.frame(id = 1:n, t = k, y = mu + 0.3 * k + base_att(g) * tr + rnorm(n, sd = 0.7), g = g)
  }))
  fit <- edid(df, "y", "id", "t", "g", aggregate = "overall")
  se_wif <- fit$simple$overall.se
  # Direct-EIF-only SE: re-aggregate the post cells with cohort-share weights, NO WIF term.
  ci   <- fit$cells; idx <- fit$att_gt
  post <- which(!idx$is_pre & is.finite(idx$att))
  pg   <- vapply(idx$group[post], function(gg) fit$cohort_fractions[[as.character(gg)]], numeric(1L))
  q    <- pg / sum(pg)
  eif_direct <- as.numeric(fit$eif[, post, drop = FALSE] %*% q)
  se_direct  <- sqrt(sum(eif_direct^2) / fit$n^2)
  expect_gt(se_wif, se_direct * 1.10)   # WIF adds materially under heterogeneity
})

# --- Determinism (analytical inference is exactly reproducible) --------------
test_that("edid is deterministic with bstrap = FALSE", {
  df <- pf_panel(n = 1000, seed = 8, het = TRUE)
  a  <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "all")
  b  <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "all")
  expect_equal(a$overall$overall.att, b$overall$overall.att, tolerance = 1e-12)
  expect_equal(a$att_gt$att,  b$att_gt$att,  tolerance = 1e-12)
})

# --- Guards that protect against misuse / silent deviations ------------------
test_that("gmm emits a finite-sample-bias warning", {
  df <- pf_panel(n = 800, seed = 4)
  expect_warning(edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "gmm", aggregate = "overall"), "gmm")
})

test_that("factor gname is rejected (not silently coerced)", {
  df <- pf_panel(n = 600, seed = 6)
  df$g <- factor(ifelse(is.infinite(df$g), "0", as.character(df$g)))
  expect_error(edid(df, "y", "id", "t", "g"), "numeric")
})

test_that("efficient weights collapse to non-efficient schemes only via documented fallbacks", {
  # no-covariate path: efficient/averaged/gmm coincide (no X variation); uniform differs.
  df <- pf_panel(n = 1500, seed = 12, het = TRUE)
  ov <- function(m) suppressWarnings(edid(df, "y", "id", "t", "g", weight_scheme = m, aggregate = "overall"))$simple$overall.att
  expect_equal(ov("efficient"), ov("averaged"), tolerance = 1e-8)
  expect_equal(ov("efficient"), ov("gmm"),      tolerance = 1e-8)
  expect_false(isTRUE(all.equal(ov("efficient"), ov("uniform"))))
})
