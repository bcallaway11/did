# test-edid-weightsname.R
# Observation/sampling weights (`weightsname`) on the no-covariate path:
#   (a) weightsname = NULL / a constant column reproduces the unweighted fit exactly;
#   (b) the weighted Hajek estimand matches the standard weighted DR/CS oracle
#       (did::att_gt(weightsname=) on the PT-Post anchor; DRDID on a 2x2);
#   (c) the weighted no-covariate estimation-effect correction is structurally sound
#       and FD-consistent;
#   (d) input validation and the covariate-path scope-out gate.

# ---- small staggered fixtures ---------------------------------------------------
mk_panel <- function(n = 300, Tn = 5, seed = 1, never = c(0, Inf), tau = 1.2) {
  set.seed(seed)
  nt <- sample(never, 1L)
  g  <- sample(c(nt, 3, 4), n, replace = TRUE, prob = c(.5, .25, .25))
  ai <- stats::rnorm(n)
  w  <- stats::runif(n, 0.4, 4)
  cl <- sample(1:10, n, replace = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    t  <- seq_len(Tn)
    trt <- is.finite(g[i]) && g[i] != 0
    y  <- ai[i] + 0.3 * t + stats::rnorm(Tn) +
      ifelse(trt & t >= g[i], tau * (t - g[i] + 1), 0)
    data.frame(id = i, yr = t, g = g[i], y = y, w = w[i], cl = cl[i])
  }))
}

# (a) ------------------------------------------------------------------------------
test_that("weightsname = NULL reproduces the default fit byte-for-byte", {
  d  <- mk_panel(seed = 11)
  f0 <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "all", pt_assumption = "all")
  fN <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "all", pt_assumption = "all",
             weightsname = NULL)
  expect_identical(f0$att_gt, fN$att_gt)
  expect_identical(as.numeric(f0$eif), as.numeric(fN$eif))
  expect_identical(f0$overall$overall.att, fN$overall$overall.att)
})

test_that("a constant weight column reproduces the unweighted fit to ~1e-12 across the matrix", {
  d <- mk_panel(seed = 12); d$wc <- 7
  for (pt in c("all", "post")) for (ws in c("efficient", "averaged", "uniform"))
    for (sh in c("ridge", "ledoit_wolf", "none")) {
      fu <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "all",
                 pt_assumption = pt, weight_scheme = ws, omega_cov_shrink = sh)
      fw <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "all",
                 pt_assumption = pt, weight_scheme = ws, omega_cov_shrink = sh, weightsname = "wc")
      expect_lt(max(abs(fu$att_gt$att - fw$att_gt$att)), 1e-12)
      expect_lt(max(abs(fu$att_gt$se  - fw$att_gt$se), na.rm = TRUE), 1e-11)
      expect_lt(max(abs(as.numeric(fu$eif) - as.numeric(fw$eif))), 1e-11)
    }
})

test_that("the weighted Omega* == crossprod(weighted psi)/n^2 identity holds exactly", {
  d  <- mk_panel(seed = 13, never = Inf)
  pn <- prepare_edid_panel(d, "y", "id", "yr", "g", NULL, NULL, NULL, 0L, "w")
  prs <- enumerate_valid_pairs_edid(3, pn$treatment_groups, pn$time_periods, pn$period_1, "all", 0L)
  skip_if(nrow(prs) < 2L)
  om  <- compute_omega_star_nocov_edid(3, 5, prs, pn, "all")
  psi <- compute_psi_moments_nocov_edid(3, 5, prs, pn)
  expect_lt(max(abs(om - crossprod(psi) / pn$n^2)), 1e-12 * max(abs(om)))
})

# (b) oracle matches ---------------------------------------------------------------
test_that("weighted PT-Post anchor matches did::att_gt(weightsname=) to ~1e-8", {
  skip_if_not_installed("did")
  d  <- mk_panel(n = 400, seed = 21, never = 0)   # att_gt 0-never convention
  fp <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "none",
             pt_assumption = "post", weightsname = "w")
  ag <- did::att_gt(yname = "y", tname = "yr", idname = "id", gname = "g", xformla = ~1,
                    data = d, weightsname = "w", control_group = "nevertreated",
                    bstrap = FALSE, cband = FALSE, base_period = "universal", est_method = "reg")
  e <- fp$att_gt[fp$att_gt$time >= fp$att_gt$group, c("group", "time", "att", "se")]
  m <- merge(e, data.frame(group = ag$group, time = ag$t, att_did = ag$att, se_did = ag$se),
             by = c("group", "time"))
  expect_lt(max(abs(m$att - m$att_did)), 1e-8)
  expect_lt(max(abs(m$se  - m$se_did), na.rm = TRUE), 1e-8)
})

test_that("weighted 2x2 matches DRDID weighted to ~1e-8", {
  skip_if_not_installed("DRDID")
  set.seed(31); n <- 600
  ai <- stats::rnorm(n); g <- sample(c(0, 2), n, replace = TRUE); w <- stats::runif(n, .5, 3)
  d <- do.call(rbind, lapply(seq_len(n), function(i) {
    t <- 1:2
    y <- ai[i] + 0.4 * t + stats::rnorm(2) + ifelse(g[i] > 0 & t >= g[i], 2.0, 0)
    data.frame(id = i, yr = t, g = g[i], y = y, w = w[i])
  }))
  fp <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "none",
             pt_assumption = "post", weightsname = "w")
  e <- fp$att_gt[fp$att_gt$group == 2 & fp$att_gt$time == 2, ]
  dd <- d; dd$D <- as.integer(dd$g == 2)
  out <- DRDID::drdid(yname = "y", tname = "yr", idname = "id", dname = "D", xformla = ~1,
                      data = dd, panel = TRUE, weightsname = "w", estMethod = "trad")
  expect_lt(abs(e$att - out$ATT), 1e-8)
  expect_lt(abs(e$se  - out$se),  1e-8)
})

# (c) weighted estimation effect ---------------------------------------------------
test_that("weighted no-cov estimation-effect correction: structural identities + Bessel reduction", {
  d  <- mk_panel(n = 200, seed = 41, never = Inf)
  pn <- prepare_edid_panel(d, "y", "id", "yr", "g", NULL, NULL, NULL, 0L, "w")
  prs <- enumerate_valid_pairs_edid(3, pn$treatment_groups, pn$time_periods, pn$period_1, "all", 0L)
  skip_if(nrow(prs) < 2L)
  om  <- compute_omega_star_nocov_edid(3, 5, prs, pn, "all")
  w   <- compute_efficient_weights_edid(om)
  ee  <- compute_nocov_ee_correction_edid(3, 5, prs, pn, omega_raw = om, omega_used = om,
                                          weights = w, shrink_lambda = NA, return_D = TRUE)
  expect_true(isTRUE(ee$applied))
  expect_lt(max(abs(colSums(ee$D))), 1e-8 * max(abs(ee$D)))   # sum_i d_i = 0 exactly
  expect_equal(ee$q_opt, -ee$cov_lead, tolerance = 1e-10)
  expect_gt(ee$delta_df, 0)
  expect_gte(ee$q_opt, 0)
  # weighted Bessel fraction reduces to 1/(m-1) when weights are constant
  pn_u <- prepare_edid_panel(transform(d, w = 1), "y", "id", "yr", "g", NULL, NULL, NULL, 0L, "w")
  ee_u <- compute_nocov_ee_correction_edid(3, 5, prs, pn_u, omega_raw = om, omega_used = om,
                                           weights = w, shrink_lambda = NA)
  ee_n <- compute_nocov_ee_correction_edid(
    3, 5, prs,
    prepare_edid_panel(d, "y", "id", "yr", "g", NULL, NULL, NULL, 0L, NULL),
    omega_raw = om, omega_used = om, weights = w, shrink_lambda = NA)
  expect_equal(ee_u$delta_df, ee_n$delta_df, tolerance = 1e-12)
})

# (d) validation + gates -----------------------------------------------------------
test_that("weightsname validation rejects bad columns and the covariate path is scoped out", {
  d <- mk_panel(seed = 51)
  expect_error(edid(d, "y", "id", "yr", "g", weightsname = "nope"), "not a column")
  d2 <- d; d2$wneg <- -1
  expect_error(edid(d, "y", "id", "yr", "g", weightsname = "y"), "coincides")
  expect_error(edid(transform(d, wn = -1), "y", "id", "yr", "g", weightsname = "wn"), "nonnegative")
  expect_error(edid(transform(d, wz = 0), "y", "id", "yr", "g", weightsname = "wz"), "all zero")
  # time-varying weight rejected
  dv <- d; dv$wv <- dv$yr
  expect_error(edid(dv, "y", "id", "yr", "g", weightsname = "wv"), "time-invariant")
  # covariate path + weights: now SUPPORTED (weighted-covariate rollout) -- runs (no guard stop), and a
  # CONSTANT weight column reproduces the unweighted covariate fit (the FD-oracle validation of the
  # weighted estimation-effect lives in test-edid-weighted-cov-e2e.R).
  dx <- d; set.seed(1); xmap <- stats::rnorm(length(unique(d$id))); dx$x <- xmap[match(dx$id, sort(unique(dx$id)))]
  fwx <- suppressWarnings(edid(dx, "y", "id", "yr", "g", xformla = ~x, weightsname = "w",
              weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, misspec_robust = FALSE))
  expect_false(is.null(fwx$att_gt))
  dxc <- dx; dxc$w <- 7                                 # constant column -> unit_weights == 1
  fxc <- suppressWarnings(edid(dxc, "y", "id", "yr", "g", xformla = ~x, weightsname = "w",
              weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, misspec_robust = FALSE))
  fxu <- suppressWarnings(edid(dx,  "y", "id", "yr", "g", xformla = ~x,
              weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, misspec_robust = FALSE))
  ok <- is.finite(fxc$att_gt$se) & is.finite(fxu$att_gt$se)
  expect_equal(fxc$att_gt$se[ok], fxu$att_gt$se[ok], tolerance = 1e-8)
})

test_that("the weighted estimand differs from the unweighted one (sanity)", {
  d  <- mk_panel(n = 400, seed = 61)
  fu <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "all", pt_assumption = "all")
  fw <- edid(d, "y", "id", "yr", "g", xformla = ~1, aggregate = "all", pt_assumption = "all",
             weightsname = "w")
  expect_false(isTRUE(all.equal(fu$att_gt$att, fw$att_gt$att, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(fu$overall$overall.att, fw$overall$overall.att, tolerance = 1e-6)))
})
