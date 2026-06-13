library(testthat)

# ===========================================================================
# AKS B-FLCI inference layer of edid_adaptive() (Armstrong, Kline & Sun,
# Econometrica 2025, Section 4.2.2 / eqs. (7)-(8), arXiv-v6 numbering),
# validated against aks_inference_fixtures.rds: the recorded output of the
# AUTHORS' OWN R code (calculate_adaptive_estimates.R, calculate_simple_CI.R,
# calculate_B_FLCI.R of github.com/lsun20/MissAdapt, commit 98d823a, run
# verbatim against their lookup tables) on the README/dCdH vignette inputs
# and 5 synthetic input sets spanning |corr| in [0.07, 0.97] -- each with
# B in {0, 1, 9}.
#
# Reproduction contract (per the fixture file's own header):
#   - critical values and CI bounds: EXACT (deterministic spline lookups;
#     tolerance 1e-12 relative, observed difference is bitwise 0);
#   - min/max coverages: the fixtures record both the authors' Monte Carlo
#     (set.seed(1), 100000 draws over b = seq(-9, 9, 0.025)) and a
#     deterministic quadrature cross-check (fields quad_*). The package's
#     deterministic quadrature must match the quad_* fields to ~1e-10 and the
#     authors' MC within its noise (<= 5e-3).
#   - quad_*_cov_min_inB restricts to |b/sigma_O| <= B-row value, the region
#     where AKS eq. (8) guarantees >= 95%: the ADAPTIVE intervals meet it
#     everywhere (min 0.9508 across all fixtures); the SHIPPED soft-threshold
#     table fails it at large B*|rho| (0.868 at |rho| = 0.97, B = 9; 0.743 at
#     the on-grid corner rho = -0.995, B = 9) because that table is calibrated
#     to a different (extrapolated) threshold than the soft-threshold estimate
#     -- which is why st_cv = "exact" (the default) re-solves eq. (8) at the
#     correct lambda*(rho). Fixture md5: 26b811fce37b5037413d10df2dc31752.
# ===========================================================================

aks_inf_fixtures <- function() {
  path <- test_path("aks_inference_fixtures.rds")
  skip_if(!file.exists(path), "aks_inference_fixtures.rds not available")
  readRDS(path)
}

aks_core_for <- function(p, tables) {
  did:::.edid_aks_core(YR = p$YR, VR = p$VR, YU = p$YU, VU = p$VU, VUR = p$VUR,
                       tables = tables)
}

st_delta_for <- function(lambda) {
  function(t) (t > lambda) * (t - lambda) + (t < -lambda) * (t + lambda)
}

B_SET <- c(0, 1, 9)
B_KEY <- c("B0", "B1", "B9")

# ---------------------------------------------------------------------------
# 1. Exact reproduction of every fixture estimate, critical value, and bound
# ---------------------------------------------------------------------------

test_that("derived components reproduce the authors' code on every fixture set", {
  f   <- aks_inf_fixtures()
  tab <- did:::.edid_aks_lookup()
  for (nm in names(f)) {
    r    <- f[[nm]]
    core <- aks_core_for(r$inputs, tab)
    d    <- r$derived
    expect_equal(core$YO,             d$YO,             tolerance = 1e-12, info = nm)
    expect_equal(core$VO,             d$VO,             tolerance = 1e-12, info = nm)
    expect_equal(core$VUO,            d$VUO,            tolerance = 1e-12, info = nm)
    expect_equal(core$tO,             d$tO,             tolerance = 1e-12, info = nm)
    expect_equal(core$corr,           d$corr,           tolerance = 1e-12, info = nm)
    expect_equal(core$GMM,            d$GMM,            tolerance = 1e-12, info = nm)
    expect_equal(core$se_GMM,         d$se_GMM,         tolerance = 1e-12, info = nm)
    expect_equal(core$adaptive,       d$adaptive,       tolerance = 1e-12, info = nm)
    expect_equal(core$adaptive_st,    d$adaptive_st,    tolerance = 1e-12, info = nm)
    expect_equal(core$soft_threshold, d$soft_threshold, tolerance = 1e-12, info = nm)
  }
})

test_that("B-FLCI cvs and bounds reproduce the authors' code exactly (st_cv = 'missadapt')", {
  f   <- aks_inf_fixtures()
  tab <- did:::.edid_aks_lookup()
  for (nm in names(f)) {
    r    <- f[[nm]]
    core <- aks_core_for(r$inputs, tab)
    ci   <- did:::.edid_aks_ci(core, tab, B_SET, st_cv = "missadapt")
    ad   <- ci[ci$variant == "adaptive", ]
    st   <- ci[ci$variant == "soft_threshold", ]
    expect_identical(unique(st$cv_source), "missadapt_table")
    for (i in seq_along(B_SET)) {
      z <- r$b_flci[[B_KEY[i]]]
      lab <- paste(nm, B_KEY[i])
      # the effective table row (B = 0 -> the B-tilde = 0.01 row)
      expect_equal(ad$B_tilde[i], z$B_row_value, tolerance = 1e-12, info = lab)
      # critical values: exact spline-lookup reproduction
      expect_equal(ad$cv[i], z$cv_adaptive, tolerance = 1e-12, info = lab)
      expect_equal(st$cv[i], z$cv_st,       tolerance = 1e-12, info = lab)
      # bounds: center +- cv * sqrt(VU) (sigma_U scale, never se_GMM/sigma_R)
      expect_equal(ad$lower[i], z$adaptive_lo, tolerance = 1e-12, info = lab)
      expect_equal(ad$upper[i], z$adaptive_hi, tolerance = 1e-12, info = lab)
      expect_equal(st$lower[i], z$st_lo,       tolerance = 1e-12, info = lab)
      expect_equal(st$upper[i], z$st_hi,       tolerance = 1e-12, info = lab)
    }
    expect_equal(unique(ci$sigma_U), sqrt(r$inputs$VU), tolerance = 1e-12, info = nm)
  }
})

test_that("README/dCdH headline values are pinned at full precision", {
  f   <- aks_inf_fixtures()
  tab <- did:::.edid_aks_lookup()
  core <- aks_core_for(f$readme_dcdh$inputs, tab)
  expect_equal(core$corr,     -0.769619216097768,   tolerance = 1e-14)
  expect_equal(core$tO,       -1.74735919034965,    tolerance = 1e-14)
  expect_equal(core$adaptive,    0.00356524756118796, tolerance = 1e-14)
  expect_equal(core$adaptive_st, 0.00360684586877713, tolerance = 1e-14)
  ci <- did:::.edid_aks_ci(core, tab, B_SET, st_cv = "missadapt")
  expect_equal(ci$cv[ci$variant == "adaptive"],
               c(1.53826804222675, 1.74223004325498, 2.32252928068568),
               tolerance = 1e-14)
  expect_equal(ci$cv[ci$variant == "soft_threshold"],
               c(1.62214812951357, 1.76618645610213, 2.11419155201473),
               tolerance = 1e-14)
  # the B = 1 adaptive FLCI of AKS Table 2 Panel B (sqrt(VU) = 0.0014)
  ad1 <- ci[ci$variant == "adaptive" & ci$B == 1, ]
  expect_equal(ad1$lower, 0.00112612550063099, tolerance = 1e-14)
  expect_equal(ad1$upper, 0.00600436962174494, tolerance = 1e-14)
  # printed-precision regression against AKS Table 2 Panel B (1.54/1.74/2.32,
  # 1.62/1.77/2.11)
  expect_identical(round(ci$cv[ci$variant == "adaptive"], 2), c(1.54, 1.74, 2.32))
  expect_identical(round(ci$cv[ci$variant == "soft_threshold"], 2), c(1.62, 1.77, 2.11))
})

# ---------------------------------------------------------------------------
# 2. Deterministic quadrature coverage: equals the fixture cross-checks,
#    agrees with the authors' Monte Carlo, and verifies the eq.-(8) guarantee
# ---------------------------------------------------------------------------

test_that("quadrature coverage matches the fixture quad_* fields and the authors' MC", {
  skip_on_cran()
  f   <- aks_inf_fixtures()
  tab <- did:::.edid_aks_lookup()
  bg  <- did:::.edid_aks_b_grid()
  sg <- tanh(seq(-3, -0.05, 0.05))
  for (nm in names(f)) {
    r    <- f[[nm]]
    core <- aks_core_for(r$inputs, tab)
    st_d <- st_delta_for(core$soft_threshold)

    # simple CI (center +- 1.96 sigma_U): bounds exact, coverage band
    sc <- r$simple_ci
    sigU <- sqrt(r$inputs$VU)
    expect_equal(core$adaptive - 1.96 * sigU, sc$adaptive_lo, tolerance = 1e-12, info = nm)
    expect_equal(core$adaptive + 1.96 * sigU, sc$adaptive_hi, tolerance = 1e-12, info = nm)
    expect_equal(core$adaptive_st - 1.96 * sigU, sc$st_lo, tolerance = 1e-12, info = nm)
    expect_equal(core$adaptive_st + 1.96 * sigU, sc$st_hi, tolerance = 1e-12, info = nm)
    qa <- did:::.edid_aks_flci_coverage(core$psi_fun, core$corr, 1.96, bg)
    qs <- did:::.edid_aks_flci_coverage(st_d,         core$corr, 1.96, bg)
    expect_equal(min(qa), sc$quad_adaptive_cov_min, tolerance = 1e-10, info = nm)
    expect_equal(max(qa), sc$quad_adaptive_cov_max, tolerance = 1e-10, info = nm)
    expect_equal(min(qs), sc$quad_st_cov_min,       tolerance = 1e-10, info = nm)
    expect_equal(max(qs), sc$quad_st_cov_max,       tolerance = 1e-10, info = nm)
    expect_lt(abs(min(qa) - sc$adaptive_cov_min), 5e-3)
    expect_lt(abs(max(qa) - sc$adaptive_cov_max), 5e-3)
    expect_lt(abs(min(qs) - sc$st_cov_min), 5e-3)
    expect_lt(abs(max(qs) - sc$st_cov_max), 5e-3)

    # B-FLCIs: coverage of the table cvs over the full +-9 grid and within B
    ci <- did:::.edid_aks_ci(core, tab, B_SET, st_cv = "missadapt")
    # the soft threshold the authors' calculate_B_FLCI.R actually uses inside
    # its coverage simulation: thresholds.mat splined against the SIGNED grid
    # but evaluated at abs(corr) -- an off-grid extrapolation yielding
    # ~0.45-0.54 for every rho instead of lambda*(rho). Their b_flci MC
    # coverages were generated under THIS threshold (and the shipped st cv
    # table is calibrated to it), so reproducing them requires it; the
    # quad_st_* fixture fields use the CORRECT lambda*(rho) instead.
    lam_extrap <- stats::splinefun(sg, tab$st, method = "fmm",
                                   ties = mean)(abs(core$corr))
    st_d_extrap <- st_delta_for(lam_extrap)
    for (i in seq_along(B_SET)) {
      z   <- r$b_flci[[B_KEY[i]]]
      lab <- paste(nm, B_KEY[i])
      cva <- ci$cv[ci$variant == "adaptive"][i]
      cvs <- ci$cv[ci$variant == "soft_threshold"][i]
      qa  <- did:::.edid_aks_flci_coverage(core$psi_fun, core$corr, cva, bg)
      qs  <- did:::.edid_aks_flci_coverage(st_d,         core$corr, cvs, bg)
      inB <- abs(bg) <= z$B_row_value + 1e-12
      expect_equal(min(qa), z$quad_adaptive_cov_min, tolerance = 1e-10, info = lab)
      expect_equal(max(qa), z$quad_adaptive_cov_max, tolerance = 1e-10, info = lab)
      expect_equal(min(qs), z$quad_st_cov_min,       tolerance = 1e-10, info = lab)
      expect_equal(max(qs), z$quad_st_cov_max,       tolerance = 1e-10, info = lab)
      expect_equal(min(qa[inB]), z$quad_adaptive_cov_min_inB, tolerance = 1e-10, info = lab)
      expect_equal(min(qs[inB]), z$quad_st_cov_min_inB,       tolerance = 1e-10, info = lab)
      # the authors' MC agrees within its noise: directly for the adaptive
      # interval; for the soft-threshold interval only under the extrapolated
      # threshold their simulation uses (observed deviation <= 0.0022 across
      # all 18 combinations -- numerical proof of the calibration mispairing,
      # since the correct-lambda quadrature deviates by up to 0.16 from their
      # MC at |rho| = 0.97)
      expect_lt(abs(min(qa) - z$adaptive_cov_min), 5e-3)
      expect_lt(abs(max(qa) - z$adaptive_cov_max), 5e-3)
      qs_x <- did:::.edid_aks_flci_coverage(st_d_extrap, core$corr, cvs, bg)
      expect_lt(abs(min(qs_x) - z$st_cov_min), 5e-3)
      expect_lt(abs(max(qs_x) - z$st_cov_max), 5e-3)
      # eq.-(8) guarantee for the HEADLINE (adaptive) interval: >= 95% within
      # |b| <= B at the plug-in rho (the fixtures sit at 0.9508-0.9522;
      # slight conservativeness from the 2-decimal cv tabulation)
      expect_gte(min(qa[inB]), 0.9505)
    }
  }
})

test_that("quadrature reproduces the authors' MC for the naive-YR and pre-test foils", {
  skip_on_cran()
  f  <- aks_inf_fixtures()
  bg <- did:::.edid_aks_b_grid()
  z  <- seq(-8.5, 8.5, by = 0.01)
  wz <- stats::dnorm(z) * 0.01
  for (nm in names(f)) {
    r <- f[[nm]]; p <- r$inputs
    YO <- p$YR - p$YU; VO <- p$VR - 2 * p$VUR + p$VU; VUO <- p$VUR - p$VU
    corr <- VUO / sqrt(VO) / sqrt(p$VU); s <- sqrt(1 - corr^2)
    k <- 1 + VO / VUO; thr <- 1.96 * sqrt(p$VR / p$VU)
    # naive restricted CI YR +- 1.96 sigma_R: statistic corr*(k*t_b - b),
    # threshold 1.96*sqrt(VR/VU) on the sigma_U scale (their calculate_simple_CI)
    qYR <- vapply(bg, function(b) {
      m <- corr * (k * (z + b) - b)
      sum((stats::pnorm((thr - m) / s) - stats::pnorm((-thr - m) / s)) * wz)
    }, numeric(1L))
    expect_lt(abs(min(qYR) - r$simple_ci$YR_cov_min), 5e-3)
    expect_lt(abs(max(qYR) - r$simple_ci$YR_cov_max), 5e-3)
    # pre-test CI: switches between YU +- 1.96 sigma_U and YR +- 1.96 sigma_R
    # on |t_O| >< 1.96 (discontinuous in Z1, so the z-grid Riemann sum is
    # coarser here: tolerance 1e-2)
    qPT <- vapply(bg, function(b) {
      tb <- z + b
      m1 <- corr * (tb - b); m2 <- corr * (k * tb - b)
      rej1 <- stats::pnorm((-1.96 - m1) / s) + 1 - stats::pnorm((1.96 - m1) / s)
      rej2 <- stats::pnorm((-thr - m2) / s) + 1 - stats::pnorm((thr - m2) / s)
      1 - sum(((abs(tb) > 1.96) * rej1 + (abs(tb) < 1.96) * rej2) * wz)
    }, numeric(1L))
    expect_lt(abs(min(qPT) - r$simple_ci$pretest_cov_min), 1e-2)
    expect_lt(abs(max(qPT) - r$simple_ci$pretest_cov_max), 1e-2)
  }
})

# ---------------------------------------------------------------------------
# 3. The corrected soft-threshold cv (st_cv = "exact"): eq.-(8) coverage gate
# ---------------------------------------------------------------------------

test_that("corrected ST cvs restore the eq.-(8) guarantee where the shipped table fails", {
  skip_on_cran()
  f   <- aks_inf_fixtures()
  tab <- did:::.edid_aks_lookup()
  bg  <- did:::.edid_aks_b_grid()

  # (a) all six fixture correlations x B in {0, 1, 9}: solved cv covers
  # >= 0.949 within |b| <= B at the CORRECT lambda*(rho)
  for (nm in names(f)) {
    r    <- f[[nm]]
    core <- aks_core_for(r$inputs, tab)
    st_d <- st_delta_for(core$soft_threshold)
    ci   <- did:::.edid_aks_ci(core, tab, B_SET, st_cv = "exact")
    st   <- ci[ci$variant == "soft_threshold", ]
    expect_identical(unique(st$cv_source), "exact")
    for (i in seq_along(B_SET)) {
      qb  <- bg[abs(bg) <= st$B_tilde[i] + 1e-12]
      qs  <- did:::.edid_aks_flci_coverage(st_d, core$corr, st$cv[i], qb)
      expect_gte(min(qs), 0.949)
    }
  }

  # (b) the on-grid corner where the SHIPPED table fails: rho = tanh(-3)
  # (= -0.99505..., column 1 of the tables -- no corr interpolation at all),
  # B-tilde = 9 (row 91). With the correct lambda*(rho) = thresholds.mat
  # column 1, the shipped cv undercovers grossly; the exact solve restores
  # the guarantee.
  rho_corner <- tanh(-3)
  lam_corner <- tab$st[1L]                       # lambda* at |rho| = 0.99505
  st_d       <- st_delta_for(lam_corner)
  inB9       <- bg
  cv_shipped <- tab$flci_cv_st[91L, 1L]
  cov_shipped <- did:::.edid_aks_flci_coverage(st_d, rho_corner, cv_shipped, inB9)
  expect_lt(min(cov_shipped), 0.78)              # documented failure: ~0.743
  expect_gt(min(cov_shipped), 0.70)
  cv_exact <- did:::.edid_aks_st_cv_exact(rho_corner, lam_corner, 9)
  cov_exact <- did:::.edid_aks_flci_coverage(st_d, rho_corner, cv_exact, inB9)
  expect_gte(min(cov_exact), 0.949)
  expect_gt(cv_exact, cv_shipped)                # the correction enlarges the cv

  # neighboring on-grid column (rho = tanh(-2.95) = -0.99454...): shipped
  # min coverage ~0.751, exact >= 0.949
  rho2 <- tanh(-2.95); lam2 <- tab$st[2L]
  st_d2 <- st_delta_for(lam2)
  cov_shipped2 <- did:::.edid_aks_flci_coverage(st_d2, rho2, tab$flci_cv_st[91L, 2L], inB9)
  expect_lt(min(cov_shipped2), 0.90)
  cv_exact2 <- did:::.edid_aks_st_cv_exact(rho2, lam2, 9)
  expect_gte(min(did:::.edid_aks_flci_coverage(st_d2, rho2, cv_exact2, inB9)), 0.949)
})

test_that("exact ST solver: bisection invariants and symmetry in the sign of rho", {
  skip_on_cran()
  # coverage at the returned cv clears the level; a slightly smaller cv does not
  cv <- did:::.edid_aks_st_cv_exact(-0.7, 0.6, 1)
  bg <- did:::.edid_aks_b_grid(); qb <- bg[abs(bg) <= 1 + 1e-12]
  st_d <- st_delta_for(0.6)
  expect_gte(min(did:::.edid_aks_flci_coverage(st_d, -0.7, cv, qb)), 0.95 - 1e-6)
  expect_lt(min(did:::.edid_aks_flci_coverage(st_d, -0.7, cv - 1e-3, qb)), 0.95)
  # exact sign invariance (the pivot's coverage is even in rho)
  expect_equal(did:::.edid_aks_st_cv_exact(0.7, 0.6, 1), cv, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 4. Lookup conventions: B-row tolerance matching, |corr| mirror, clamping
# ---------------------------------------------------------------------------

test_that("B values are matched to the grid within tolerance, never by float equality", {
  tab <- did:::.edid_aks_lookup()
  # B = 0.3 and 0.7 crash MissAdapt's exact `==` match (the Matlab-written
  # grid doubles differ from the R literals in the last bit); here they work
  rows <- did:::.edid_aks_flci_rows(c(0, 0.3, 0.7, 1, 9, Inf), tab$flci_B_grid)
  expect_identical(rows$row, c(1L, 4L, 8L, 11L, 91L, 91L))
  expect_equal(rows$B_tilde, c(0.01, 0.3, 0.7, 1, 9, 9), tolerance = 1e-8)
  # the requested off-grid value is preserved as the label; the row value is
  # the tabulated double
  expect_identical(rows$B[6], Inf)
  # off-grid / invalid requests error informatively
  expect_error(did:::.edid_aks_flci_rows(0.05, tab$flci_B_grid), "tabulated")
  expect_error(did:::.edid_aks_flci_rows(9.5,  tab$flci_B_grid), "tabulated")
  expect_error(did:::.edid_aks_flci_rows(-1,   tab$flci_B_grid), "nonnegative")
  expect_error(did:::.edid_aks_flci_rows(NA_real_, tab$flci_B_grid), "missing")
})

test_that("cv lookup at |corr| equals the authors' signed-grid lookup (mirror identity)", {
  tab <- did:::.edid_aks_lookup()
  signed_grid <- tanh(seq(-3, -0.05, 0.05))
  for (corr in c(-0.769619216097768, -0.3, -0.97)) {
    for (r in c(1L, 11L, 91L)) {
      ours   <- did:::.edid_aks_flci_cv(tab$flci_cv_adaptive, r, tab$corr_grid, abs(corr))
      theirs <- stats::splinefun(signed_grid, tab$flci_cv_adaptive[r, ],
                                 method = "fmm", ties = mean)(corr)
      expect_equal(ours, theirs, tolerance = 1e-14)
      ours_st   <- did:::.edid_aks_flci_cv(tab$flci_cv_st, r, tab$corr_grid, abs(corr))
      theirs_st <- stats::splinefun(signed_grid, tab$flci_cv_st[r, ],
                                    method = "fmm", ties = mean)(corr)
      expect_equal(ours_st, theirs_st, tolerance = 1e-14)
    }
  }
  # corr > 0 (possible under assume_efficient = FALSE): the lookup uses
  # |corr| -- identical cvs to the mirrored negative-corr input. (MissAdapt's
  # signed-grid convention would silently extrapolate ~0.5-1.9 off-grid here.)
  # The pair shares (YO, VO, VU) and flips the sign of VUO: corr = -/+ 0.3.
  core_neg <- did:::.edid_aks_core(YR = 0.5, VR = 1.4, YU = 0.2, VU = 1, VUR = 0.7,
                                   tables = tab)
  core_pos <- did:::.edid_aks_core(YR = 0.5, VR = 2.6, YU = 0.2, VU = 1, VUR = 1.3,
                                   tables = tab)
  expect_lt(core_neg$corr, 0)
  expect_gt(core_pos$corr, 0)
  expect_equal(abs(core_pos$corr), abs(core_neg$corr), tolerance = 1e-12)
  ci_neg <- did:::.edid_aks_ci(core_neg, tab, c(0, 1), st_cv = "missadapt")
  ci_pos <- did:::.edid_aks_ci(core_pos, tab, c(0, 1), st_cv = "missadapt")
  expect_equal(ci_pos$cv, ci_neg$cv, tolerance = 1e-12)
})

test_that("|corr| outside the tabulated grid is clamped (with the core's warning), not extrapolated", {
  tab <- did:::.edid_aks_lookup()
  # |corr| ~ 0.0316 < grid minimum 0.04996: VUR = VR convention with VR/VU close to 1
  expect_warning(
    core <- did:::.edid_aks_core(YR = 0.1, VR = 0.999, YU = 0, VU = 1, VUR = 0.999,
                                 tables = tab),
    "outside the tabulated grid"
  )
  expect_lt(abs(core$corr), min(tab$corr_grid))
  expect_equal(core$acorr_eval, min(tab$corr_grid), tolerance = 1e-12)
  # the CI layer reuses the clamped value: cv equals the grid-edge lookup
  ci <- did:::.edid_aks_ci(core, tab, 0, st_cv = "missadapt")
  edge <- did:::.edid_aks_flci_cv(tab$flci_cv_adaptive, 1L, tab$corr_grid,
                                  min(tab$corr_grid))
  expect_equal(ci$cv[ci$variant == "adaptive"], edge, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 5. edid_adaptive() integration: fits, both covariance conventions, API
# ---------------------------------------------------------------------------

make_panel_aksci <- function(seed, n = 320L) {
  set.seed(seed)
  Tt  <- 6L
  coh <- sample(c(3, 5, Inf), n, replace = TRUE, prob = c(.3, .3, .4))
  df  <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  ufe  <- rnorm(n)
  df$y <- ufe[df$id] + 0.2 * df$time + rnorm(nrow(df), 0, 0.5) + 1 * (df$time >= df$g)
  df
}

test_that("edid_adaptive attaches B-FLCIs on simulated fits (AUTO TRUE and FALSE paths)", {
  skip_on_cran()
  df <- make_panel_aksci(20260611L)
  fit_R <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
                aggregate = "event_study", cband = FALSE)
  fit_U <- edid(df, "y", "id", "time", "g", pt_assumption = "post",
                aggregate = "event_study", cband = FALSE)

  # AUTO -> TRUE (efficient restricted fit, no covariates)
  ad <- edid_adaptive(fit_U, fit_R)
  expect_true(ad$assume_efficient)
  expect_true(ad$assume_efficient_auto)
  expect_s3_class(ad$ci, "data.frame")
  expect_identical(nrow(ad$ci), 4L)   # {adaptive, soft_threshold} x {0, Inf}
  expect_identical(ad$ci$variant, rep(c("adaptive", "soft_threshold"), each = 2L))
  expect_identical(ad$ci$B, rep(c(0, Inf), 2L))
  expect_equal(ad$ci$B_tilde, rep(c(0.01, 9), 2L), tolerance = 1e-12)
  expect_identical(ad$st_cv, "exact")
  expect_identical(ad$ci$cv_source, c("table", "table", "exact", "exact"))
  expect_identical(ad$ci_level, 0.95)
  expect_match(ad$ci_note, "uniformly over biases")
  # geometry: every interval is center +- cv * sigma_U with sigma_U = sqrt(VU)
  expect_equal(ad$sigma_U, sqrt(ad$VU), tolerance = 1e-12)
  expect_equal(ad$ci$sigma_U, rep(ad$sigma_U, 4L), tolerance = 1e-12)
  expect_equal(ad$ci$center[ad$ci$variant == "adaptive"], rep(ad$adaptive, 2L),
               tolerance = 1e-12)
  expect_equal(ad$ci$center[ad$ci$variant == "soft_threshold"],
               rep(ad$adaptive_st, 2L), tolerance = 1e-12)
  expect_equal(ad$ci$lower, ad$ci$center - ad$ci$cv * ad$ci$sigma_U, tolerance = 1e-12)
  expect_equal(ad$ci$upper, ad$ci$center + ad$ci$cv * ad$ci$sigma_U, tolerance = 1e-12)
  # the cv is nondecreasing in B for each variant (0-FLCI vs Inf-FLCI bracket)
  cv_ad <- ad$ci$cv[ad$ci$variant == "adaptive"]
  expect_gte(cv_ad[2L], cv_ad[1L])
  expect_true(all(ad$ci$cv > 0))

  # user-supplied B adds rows (B = 0.3 exercises the tolerance match end-to-end)
  ad3 <- edid_adaptive(fit_U, fit_R, B = c(0.3, 1))
  expect_identical(nrow(ad3$ci), 8L)
  expect_equal(sort(unique(ad3$ci$B)), c(0, 0.3, 1, Inf))
  # default rows are unchanged by adding B values
  sub <- ad3$ci[ad3$ci$B %in% c(0, Inf), ]
  expect_equal(sub$cv, ad$ci$cv, tolerance = 1e-12)

  # st_cv = "missadapt" replicates the shipped-table source
  adm <- edid_adaptive(fit_U, fit_R, st_cv = "missadapt")
  expect_identical(adm$ci$cv_source, c("table", "table", "missadapt_table", "missadapt_table"))
  expect_equal(adm$ci$cv[adm$ci$variant == "adaptive"],
               ad$ci$cv[ad$ci$variant == "adaptive"], tolerance = 1e-12)
  expect_identical(adm$st_cv, "missadapt")

  # ci = FALSE suppresses the layer
  ad0 <- edid_adaptive(fit_U, fit_R, ci = FALSE)
  expect_null(ad0$ci)
  expect_null(ad0$ci_note)

  # AUTO -> FALSE (uniform-weights restricted fit is not bound-attaining)
  fit_Ru <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
                 weight_scheme = "uniform", aggregate = "event_study",
                 cband = FALSE)
  adF <- edid_adaptive(fit_U, fit_Ru)
  expect_false(adF$assume_efficient)
  expect_true(adF$assume_efficient_auto)
  expect_s3_class(adF$ci, "data.frame")
  expect_identical(nrow(adF$ci), 4L)
  expect_equal(adF$ci$lower, adF$ci$center - adF$ci$cv * adF$ci$sigma_U,
               tolerance = 1e-12)
  # explicit FALSE on the efficient pair also runs the empirical-covariance path
  adE <- edid_adaptive(fit_U, fit_R, assume_efficient = FALSE)
  expect_false(adE$assume_efficient)
  expect_s3_class(adE$ci, "data.frame")

  # event-study parameter: per-e rows
  ade <- edid_adaptive(fit_U, fit_R, parameter = "event_study")
  expect_s3_class(ade$ci, "data.frame")
  expect_identical(nrow(ade$ci), length(ade$e_set) * 4L)
  expect_true(all(c("e", "variant", "B", "B_tilde", "center", "sigma_U",
                    "cv", "lower", "upper", "cv_source") %in% names(ade$ci)))
  for (e in ade$e_set) {
    rows_e <- ade$ci[ade$ci$e == e & ade$ci$variant == "adaptive", ]
    tab_e  <- ade$table[ade$table$e == e, ]
    expect_equal(unique(rows_e$center), tab_e$adaptive, tolerance = 1e-12)
    expect_equal(unique(rows_e$sigma_U), sqrt(tab_e$VU), tolerance = 1e-12)
  }

  # print methods surface the FLCIs and the validity note
  expect_output(print(ad), "95% adaptive FLCIs", fixed = TRUE)
  expect_output(print(ad), "uniformly over PT violations", fixed = TRUE)
  expect_output(print(ade), "95% adaptive FLCIs", fixed = TRUE)
  expect_output(print(ad0), "ci = TRUE", fixed = TRUE)
})

test_that("level != 0.95 errors loudly (the tables are 95%-only)", {
  skip_on_cran()
  df <- make_panel_aksci(20260612L, n = 240L)
  fit_R <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
                aggregate = "event_study", cband = FALSE)
  fit_U <- edid(df, "y", "id", "time", "g", pt_assumption = "post",
                aggregate = "event_study", cband = FALSE)
  expect_error(edid_adaptive(fit_U, fit_R, level = 0.90), "95% level only", fixed = TRUE)
  expect_error(edid_adaptive(fit_U, fit_R, level = 0.99), "level must be 0.95")
  expect_error(edid_adaptive(fit_U, fit_R, level = NA_real_), "level must be 0.95")
  # invalid B propagates the row-resolution error through the main entry point
  expect_error(edid_adaptive(fit_U, fit_R, B = 0.05), "tabulated")
})
