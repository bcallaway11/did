library(testthat)

# ===========================================================================
# Section-5 toolkit of Chen, Sant'Anna & Xie (2025):
#   edid_hausman()  -- Theorem 5.1 joint + eqn (5.5) scalar Hausman tests
#   edid_sargan()   -- Section 5.1 incremental Sargan procedure (Holm step-down)
#   edid_frontier() -- Theorem 5.2 / eqn (5.6) robustness frontier
#   edid_adaptive() -- Proposition 5.1 / eqn (5.4) AKS adaptive estimator
# plus the internal `moment_set` pair-restriction mechanism in edid().
# ===========================================================================

# Staggered panel generator for the toolkit tests. Under viol = 0 the DGP
# satisfies PT-All exactly (two-way structure + iid noise + post-treatment
# effects). viol != 0 is the declared robustness deviation: it shifts cohort
# g = 3's outcomes in period 1 ONLY, which violates PT-All (the t'' = 1
# baseline moments and the cross-cohort pre-period moments are contaminated)
# while leaving PT-Post intact (all post-treatment comparisons use base period
# g - 1 > 1, untouched).
make_panel_toolkit <- function(seed, n = 400L, viol = 0) {
  set.seed(seed)
  Tt  <- 6L
  coh <- sample(c(3, 5, Inf), n, replace = TRUE, prob = c(.3, .3, .4))
  df  <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  ufe  <- rnorm(n)
  df$y <- ufe[df$id] + 0.2 * df$time + rnorm(nrow(df), 0, 0.5) + 1 * (df$time >= df$g)
  if (viol != 0) df$y <- df$y + viol * (df$g == 3 & df$time == 1)
  df
}

fit_pair_toolkit <- function(df, ...) {
  list(
    R = edid(df, "y", "id", "time", "g", pt_assumption = "all",
             aggregate = "event_study", cband = FALSE, ...),
    U = edid(df, "y", "id", "time", "g", pt_assumption = "post",
             aggregate = "event_study", cband = FALSE, ...)
  )
}

# ===========================================================================
# edid_hausman: Monte Carlo size under the null (PT-All true)
# ===========================================================================
test_that("Hausman test has approximately correct size under PT-All and power under violation", {
  skip_on_cran()
  R_null <- 200L
  p_null <- vapply(seq_len(R_null), function(r) {
    df <- make_panel_toolkit(20260610L + r)
    ff <- fit_pair_toolkit(df)
    edid_hausman(ff$U, ff$R)$p_value
  }, numeric(1L))
  rej_null <- mean(p_null < 0.05)
  expect_gte(rej_null, 0.01)
  expect_lte(rej_null, 0.12)
  # p-values roughly uniform: median not far from 0.5
  expect_gt(median(p_null), 0.30)
  expect_lt(median(p_null), 0.70)

  # Power under a PT-All violation that preserves PT-Post
  R_alt <- 60L
  p_alt <- vapply(seq_len(R_alt), function(r) {
    df <- make_panel_toolkit(40000L + r, viol = 0.7)
    ff <- fit_pair_toolkit(df)
    edid_hausman(ff$U, ff$R)$p_value
  }, numeric(1L))
  rej_alt <- mean(p_alt < 0.05)
  expect_gt(rej_alt, rej_null)   # power > size
  expect_gt(rej_alt, 0.5)
})

# ===========================================================================
# edid_hausman: structure, scalar version by hand, unit-order invariance
# ===========================================================================
test_that("edid_hausman returns the documented structure and matches by-hand computation", {
  df <- make_panel_toolkit(1L, n = 300L)
  ff <- fit_pair_toolkit(df)
  h  <- edid_hausman(ff$U, ff$R)

  expect_s3_class(h, "edid_hausman")
  expect_true(all(c("statistic", "df", "p_value", "d", "D", "scalar") %in% names(h)))
  expect_identical(h$parameter, "event_study")
  expect_identical(h$e_set, c(0, 1, 2, 3))
  expect_equal(dim(h$D), c(4L, 4L))
  expect_true(is.finite(h$statistic) && h$statistic >= 0)
  expect_true(h$df >= 1L && h$df <= 4L)
  expect_output(print(h), "Hausman test of PT-All vs PT-Post")

  # By-hand reconstruction from the AGGTEobj influence functions (iid case):
  # d = ES_U - ES_R over post e's; D = E_n[xi xi']; H = n d' D^{-1} d.
  n   <- ff$R$n
  aR  <- ff$R$event_study; aU <- ff$U$event_study
  pos_R <- match(h$e_set, aR$egt); pos_U <- match(h$e_set, aU$egt)
  d_hand  <- aU$att.egt[pos_U] - aR$att.egt[pos_R]
  xi_hand <- aU$inf.function$dynamic.inf.func.e[, pos_U, drop = FALSE] -
             aR$inf.function$dynamic.inf.func.e[, pos_R, drop = FALSE]
  D_hand  <- crossprod(xi_hand) / n
  H_hand  <- as.numeric(n * t(d_hand) %*% solve(D_hand) %*% d_hand)
  expect_equal(unname(h$d), d_hand, tolerance = 1e-10)
  expect_equal(h$D, D_hand, tolerance = 1e-10)
  if (h$df == length(h$e_set)) expect_equal(h$statistic, H_hand, tolerance = 1e-8)

  # Scalar eqn (5.5) statistics by hand: H_e = n d_e^2 / mean(xi_e^2)
  for (j in seq_along(h$e_set)) {
    H_j <- n * d_hand[j]^2 / mean(xi_hand[, j]^2)
    expect_equal(h$scalar$H[j], H_j, tolerance = 1e-10)
    expect_equal(h$scalar$p_value[j],
                 pchisq(H_j, df = 1, lower.tail = FALSE), tolerance = 1e-12)
  }
  # ES_avg row present
  expect_identical(h$scalar$parameter[nrow(h$scalar)], "ES_avg")

  # overall parameter: scalar test with df <= 1
  h_ov <- edid_hausman(ff$U, ff$R, parameter = "overall")
  expect_true(h_ov$df %in% c(0L, 1L))
  expect_equal(h_ov$statistic,
               h$scalar$H[h$scalar$parameter == "ES_avg"], tolerance = 1e-10)
})

test_that("edid_hausman scalar version matches by-hand on a 2-parameter toy", {
  # One cohort g = 3 with T = 4 -> exactly two post e's {0, 1}
  set.seed(33)
  n <- 150L; Tt <- 4L
  coh <- sample(c(3, Inf), n, replace = TRUE)
  df <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  df$y <- rnorm(n)[df$id] + 0.1 * df$time + rnorm(nrow(df), 0, 0.4) + (df$time >= df$g)
  ff <- fit_pair_toolkit(df)
  h  <- edid_hausman(ff$U, ff$R)
  expect_identical(h$e_set, c(0, 1))
  expect_equal(dim(h$D), c(2L, 2L))
  aR <- ff$R$event_study; aU <- ff$U$event_study
  for (j in 1:2) {
    e  <- h$e_set[j]
    dU <- aU$att.egt[match(e, aU$egt)] - aR$att.egt[match(e, aR$egt)]
    xi <- aU$inf.function$dynamic.inf.func.e[, match(e, aU$egt)] -
          aR$inf.function$dynamic.inf.func.e[, match(e, aR$egt)]
    expect_equal(h$scalar$H[j], n * dU^2 / mean(xi^2), tolerance = 1e-10)
  }
})

test_that("edid_hausman is invariant to unit relabeling/order", {
  df <- make_panel_toolkit(2L, n = 250L)
  ff1 <- fit_pair_toolkit(df)
  h1  <- edid_hausman(ff1$U, ff1$R)

  # Relabel units with a permutation: same data, different internal unit order
  set.seed(99)
  perm <- sample(250L)
  df2  <- df
  df2$id <- perm[df$id]
  ff2 <- fit_pair_toolkit(df2)
  h2  <- edid_hausman(ff2$U, ff2$R)

  expect_equal(h1$statistic, h2$statistic, tolerance = 1e-8)
  expect_identical(h1$df, h2$df)
  expect_equal(h1$p_value, h2$p_value, tolerance = 1e-8)
  expect_equal(h1$scalar$H, h2$scalar$H, tolerance = 1e-8)
})

test_that("edid_hausman cluster-aware path runs and validation catches mismatches", {
  df <- make_panel_clustered(n_clusters_treat = 8L, n_clusters_never = 8L,
                             units_per_cluster = 5L, n_periods = 5L, seed = 5L)
  fR <- edid(df, "outcome", "unit", "time", "first_treat", pt_assumption = "all",
             aggregate = "event_study", cband = FALSE, clustervars = "cluster_id")
  fU <- edid(df, "outcome", "unit", "time", "first_treat", pt_assumption = "post",
             aggregate = "event_study", cband = FALSE, clustervars = "cluster_id")
  h <- edid_hausman(fU, fR)
  expect_true(isTRUE(h$clustered))
  expect_true(is.finite(h$statistic) && h$statistic >= 0)
  expect_true(is.finite(h$p_value))

  # Mismatched clustering: one clustered fit, one not
  fU_noclust <- edid(df, "outcome", "unit", "time", "first_treat", pt_assumption = "post",
                     aggregate = "event_study", cband = FALSE)
  expect_error(edid_hausman(fU_noclust, fR), "cluster")

  # Mismatched samples: different n
  df_small <- df[df$unit <= 60L, ]
  fU_small <- edid(df_small, "outcome", "unit", "time", "first_treat",
                   pt_assumption = "post", aggregate = "event_study", cband = FALSE)
  expect_error(suppressWarnings(edid_hausman(fU_small, fR)), "sample sizes")

  # Swapped pt assumptions warn
  df2 <- make_panel_toolkit(3L, n = 200L)
  ff  <- fit_pair_toolkit(df2)
  expect_warning(edid_hausman(ff$R, ff$U), "pt_assumption")
})

# ===========================================================================
# edid_hausman: degenerate-contrast guard on the joint statistic
# ===========================================================================

# Panel whose treated cohorts are ALL below the default min_pair_units = 5, so
# the thin-cohort guard pins every cell to its just-identified moment and the
# PT-All fit coincides with the PT-Post fit up to float dust.
make_panel_pinned <- function(seed, n = 60L) {
  set.seed(seed)
  coh <- c(rep(3, 3L), rep(4, 2L), rep(Inf, n - 5L))
  df  <- data.frame(id = rep(seq_len(n), each = 6L), time = rep(1:6, n))
  df$g <- coh[df$id]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + rnorm(nrow(df), 0, 0.5) + 1 * (df$time >= df$g)
  df
}

test_that(".edid_if_diff_quadform guards degenerate contrasts (no rank found in float dust)", {
  # The PDMP/Johnson gate-run failure mode: a contrast that is pure numerical
  # noise (d ~ 1e-19, D entries ~ 1e-36) must NOT be ranked by the relative
  # eigenvalue threshold into a spurious large H with p ~ 0.
  set.seed(99L)
  n  <- 24L
  d  <- rnorm(4L) * 1e-19
  xi <- matrix(rnorm(n * 4L), n, 4L) * 1e-18
  qf <- .edid_if_diff_quadform(d, xi, n, NULL, v_scale = 1)
  expect_identical(qf$statistic, 0)
  expect_identical(qf$df, 0L)
  expect_identical(qf$p_value, 1)
  expect_true(isTRUE(qf$degenerate))

  # A genuine contrast on the same scale as v_scale is untouched by the guard.
  xi2 <- matrix(rnorm(n * 4L), n, 4L)
  d2  <- colMeans(xi2)
  qf2 <- .edid_if_diff_quadform(d2, xi2, n, NULL, v_scale = 1)
  expect_gt(qf2$df, 0L)
  expect_false(isTRUE(qf2$degenerate))
  expect_true(is.finite(qf2$statistic) && qf2$statistic >= 0)
})

test_that("edid_hausman joint test is degenerate (p = 1), not spurious, when the guard pins both fits", {
  df <- make_panel_pinned(20260612L)
  fR <- suppressWarnings(edid(df, "y", "id", "time", "g", pt_assumption = "all",
                              aggregate = "event_study", cband = FALSE))
  fU <- suppressWarnings(edid(df, "y", "id", "time", "g", pt_assumption = "post",
                              aggregate = "event_study", cband = FALSE))
  expect_identical(nrow(fR$thin_cohorts), 2L)   # both cohorts pinned

  # Joint contrast: previously the rank-thresholded pseudoinverse "found" rank
  # in the 1e-18 noise (e.g. H = 18.5, df = 4, p = 0.001); now degenerate.
  expect_message(h <- edid_hausman(fU, fR), "degenerate")
  expect_identical(h$statistic, 0)
  expect_identical(h$df, 0L)
  expect_identical(h$p_value, 1)
  expect_true(isTRUE(h$degenerate))
  expect_output(print(h), "degenerate contrast")

  # The scalar path already guarded this case; joint and scalar now agree.
  expect_true(all(h$scalar$H == 0))
  expect_true(all(h$scalar$p_value == 1))

  # The overall (scalar) parameter goes through the same joint code path.
  h_ov <- suppressMessages(edid_hausman(fU, fR, parameter = "overall"))
  expect_identical(h_ov$df, 0L)
  expect_identical(h_ov$p_value, 1)
})

# ===========================================================================
# edid_frontier (Theorem 5.2 / eqn 5.6)
# ===========================================================================
test_that("edid_frontier: H >= 0, radii monotone in tau, frontier centered at theta_R", {
  df <- make_panel_toolkit(4L, n = 300L)
  ff <- fit_pair_toolkit(df)
  fr <- edid_frontier(ff$U, ff$R)
  expect_s3_class(fr, "edid_frontier")
  tab <- fr$table

  expect_true(all(tab$H >= 0))
  expect_true(all(is.finite(tab$radius)) && all(tab$radius >= 0))
  expect_equal(tab$radius, tab$tau * sqrt(tab$H) * tab$se_R, tolerance = 1e-12)
  expect_equal(tab$frontier_low,  tab$theta_R - tab$radius, tolerance = 1e-12)
  expect_equal(tab$frontier_high, tab$theta_R + tab$radius, tolerance = 1e-12)

  # Radii monotone (strictly increasing when H > 0) in tau within parameter
  for (p in unique(tab$parameter)) {
    sub <- tab[tab$parameter == p, ]
    sub <- sub[order(sub$tau), ]
    expect_true(all(diff(sub$radius) >= 0))
    if (sub$H[1] > 0) expect_true(all(diff(sub$radius) > 0))
  }

  # Default rows: ES(e) for each shared post e + ES_avg, each at 3 tau values
  expect_identical(sort(unique(tab$tau)), c(0.25, 0.5, 1))
  expect_true("ES_avg" %in% tab$parameter)
  expect_equal(nrow(tab), (length(fr$e_set) + 1L) * 3L)

  # Frontier H agrees with the Hausman scalar statistics
  h <- edid_hausman(ff$U, ff$R)
  for (j in seq_along(h$e_set)) {
    expect_equal(unique(tab$H[tab$e == h$e_set[j] & !is.na(tab$e)]),
                 h$scalar$H[j], tolerance = 1e-10)
  }
  expect_output(print(fr), "Robustness frontier")
})

test_that("edid_frontier degenerate-D path returns H = 0 (same fit twice)", {
  df <- make_panel_toolkit(5L, n = 200L)
  fR <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
             aggregate = "event_study", cband = FALSE)
  # Passing the same fit twice: xi = 0 exactly -> guard reports H = 0, p = 1,
  # zero radius (frontier collapses to the point estimate), never NaN.
  fr <- suppressWarnings(edid_frontier(fR, fR))   # warning: pt pairing not (post, all)
  tab <- fr$table
  expect_true(all(tab$H == 0))
  expect_true(all(tab$p_value == 1))
  expect_true(all(tab$radius == 0))
  expect_equal(tab$frontier_low,  tab$theta_R, tolerance = 1e-15)
  expect_equal(tab$frontier_high, tab$theta_R, tolerance = 1e-15)
  expect_false(any(is.nan(tab$H)))
})

# ===========================================================================
# edid_adaptive (Proposition 5.1 / eqn 5.4)
# ===========================================================================
test_that("adaptive core reproduces the Application implementation on a fixture", {
  # Fixture generated by running Application/utils/aks_adaptive.R
  # (aks_adaptive_estimate) on the SAME lookup tables with
  # YR = 0.8, VR = 0.04, YU = 1.0, VU = 0.09, VUR = 0.041.
  r1 <- .edid_aks_core(YR = 0.8, VR = 0.04, YU = 1.0, VU = 0.09, VUR = 0.041)
  expect_equal(r1$YO,  -0.2,   tolerance = 1e-12)
  expect_equal(r1$VO,   0.048, tolerance = 1e-12)
  expect_equal(r1$VUO, -0.049, tolerance = 1e-12)
  expect_equal(r1$tO,         -0.912870929175, tolerance = 1e-9)
  expect_equal(r1$corr,       -0.745511258826, tolerance = 1e-9)
  expect_equal(r1$rho_aks_sq,  0.555787037037, tolerance = 1e-9)
  expect_equal(r1$GMM,         0.795833333333, tolerance = 1e-9)
  expect_equal(r1$V_GMM,       0.039979166667, tolerance = 1e-9)
  expect_equal(r1$adaptive,    0.855174444719, tolerance = 1e-8)
  expect_equal(r1$adaptive_st, 0.860515005906, tolerance = 1e-8)
  expect_equal(r1$adaptive_ht, 0.795833333333, tolerance = 1e-8)
  expect_equal(r1$pretest,     0.795833333333, tolerance = 1e-8)
  expect_equal(r1$erm,         0.888636363636, tolerance = 1e-9)
  expect_equal(r1$adaptive_erm, 0.865226962671, tolerance = 1e-8)
  expect_equal(r1$soft_threshold, 0.623665940400, tolerance = 1e-8)
  expect_equal(r1$hard_threshold, 1.392253712519, tolerance = 1e-8)
  expect_equal(r1$erm_lambda,     1.618460736429, tolerance = 1e-8)

  # Efficient case (Hausman identity VUR = VR): GMM reduces to YR exactly
  r2 <- .edid_aks_core(YR = 0.8, VR = 0.04, YU = 1.0, VU = 0.09, VUR = 0.04)
  expect_equal(r2$VO,  0.05, tolerance = 1e-12)
  expect_equal(r2$VUO, -0.05, tolerance = 1e-12)
  expect_equal(r2$GMM, 0.8,  tolerance = 1e-12)              # = YR
  expect_equal(r2$rho_aks_sq, 1 - 0.04 / 0.09, tolerance = 1e-12)  # 1 - VR/VU
  expect_equal(r2$adaptive, 0.857090679392, tolerance = 1e-8)
})

test_that("adaptive core edge behavior: sigma_O ~ 0, clamps, and VO <= 0 assert", {
  # sigma_O ~ 0 with YO = 0: tO = 0, delta*(0) ~ 0 -> adaptive ~ GMM = YR
  r <- suppressWarnings(.edid_aks_core(YR = 1.0, VR = 0.04, YU = 1.0, VU = 0.0401, VUR = 0.04))
  expect_lt(abs(r$adaptive - 1.0), 1e-3)
  expect_equal(r$GMM, 1.0, tolerance = 1e-12)

  # |corr| below the tabulated grid (rho^2 -> 0): clamped with a warning
  expect_warning(
    .edid_aks_core(YR = 0.8, VR = 0.10, YU = 1.0, VU = 0.09, VUR = 0.0901),
    "outside the tabulated grid")

  # tO outside the tabulated y-grid: clamped with a warning
  expect_warning(
    .edid_aks_core(YR = 10, VR = 0.04, YU = 0, VU = 0.09, VUR = 0.01),
    "outside the tabulated y-grid")

  # VO <= 0: no valid over-identification direction -> error
  expect_error(.edid_aks_core(YR = 0.8, VR = 0.04, YU = 1.0, VU = 0.09, VUR = 0.07),
               "not positive")
})

test_that("edid_adaptive's VO <= 0 error is diagnostic (coinciding fits / thin-cohort guard)", {
  # The error must explain the generic cause -- the efficient and conservative
  # fits coincide, the common outcome when the thin-cohort guard pins every
  # cell -- and point at the guard diagnostics, not just state VO <= 0.
  expect_error(.edid_aks_core(YR = 0.8, VR = 0.04, YU = 1.0, VU = 0.09, VUR = 0.07),
               "no over-identification direction to adapt over")
  expect_error(.edid_aks_core(YR = 0.8, VR = 0.04, YU = 1.0, VU = 0.09, VUR = 0.07),
               "thin-cohort guard")
  expect_error(.edid_aks_core(YR = 0.8, VR = 0.04, YU = 1.0, VU = 0.09, VUR = 0.07),
               "edid_weights")

  # Guard-pinned scenario: identical estimator pair (VU = VR = VUR exactly).
  # AUTO resolves assume_efficient = TRUE for the no-covariate restricted fit,
  # so VO = VU - VR = 0 exactly -> the documented hard error, with the hint.
  df <- make_panel_pinned(20260613L)
  fR <- suppressWarnings(edid(df, "y", "id", "time", "g", pt_assumption = "all",
                              aggregate = "event_study", cband = FALSE))
  expect_error(suppressWarnings(edid_adaptive(fR, fR)), "thin-cohort guard")
})

test_that("edid_adaptive runs on fits, overall and per-e, with components consistent", {
  df <- make_panel_toolkit(6L, n = 300L)
  ff <- fit_pair_toolkit(df)
  # The component identities below are the EMPIRICAL-covariance definitions, so pin the convention
  # explicitly (the NULL default is AUTO and resolves to TRUE here: the restricted fit is a no-covariate
  # non-uniform fit, hence bound-attaining).
  ad <- edid_adaptive(ff$U, ff$R, assume_efficient = FALSE)
  expect_s3_class(ad, "edid_adaptive")
  expect_identical(ad$parameter, "overall")
  expect_false(ad$assume_efficient_auto)   # explicit override recorded as non-AUTO

  # Components reproduce the definition from the aggregation IFs (iid case)
  n  <- ff$R$n
  aR <- ff$R$event_study; aU <- ff$U$event_study
  YR <- aR$overall.att; YU <- aU$overall.att
  pR <- aR$inf.function$dynamic.inf.func
  pU <- aU$inf.function$dynamic.inf.func
  expect_equal(ad$YR, YR, tolerance = 1e-12)
  expect_equal(ad$YU, YU, tolerance = 1e-12)
  expect_equal(ad$VR, sum(pR^2) / n^2, tolerance = 1e-12)
  expect_equal(ad$VU, sum(pU^2) / n^2, tolerance = 1e-12)
  expect_equal(ad$VUR, sum(pU * pR) / n^2, tolerance = 1e-12)
  expect_equal(ad$YO, YR - YU, tolerance = 1e-12)
  expect_equal(ad$tO, ad$YO / sqrt(ad$VO), tolerance = 1e-12)
  expect_equal(ad$rho_aks_sq, ad$corr^2, tolerance = 1e-12)
  expect_true(is.finite(ad$adaptive))
  expect_output(print(ad), "Adaptive event-study estimator")
  expect_output(print(ad), "95% adaptive FLCIs", fixed = TRUE)   # the B-FLCI block is printed
  expect_output(print(ad), "uniformly over PT violations")      # with its validity note

  # Per-e variant
  ad_e <- edid_adaptive(ff$U, ff$R, parameter = "event_study", assume_efficient = FALSE)
  expect_identical(ad_e$parameter, "event_study")
  expect_equal(nrow(ad_e$table), length(ad_e$e_set))
  expect_true(all(is.finite(ad_e$table$adaptive)))
})

test_that("assume_efficient = TRUE imposes the Hausman identity VUR = VR (GMM = YR exactly)", {
  df <- make_panel_toolkit(6L, n = 300L)
  ff <- fit_pair_toolkit(df)

  ad_emp <- edid_adaptive(ff$U, ff$R, assume_efficient = FALSE) # empirical IF covariance (explicit)
  ad_eff <- edid_adaptive(ff$U, ff$R, assume_efficient = TRUE)  # imposed identity

  # The identity holds EXACTLY (identical, not all.equal): VUR = VR, GMM = YR,
  # V_GMM = VR, and corr = -sqrt(1 - VR/VU).
  expect_true(ad_eff$assume_efficient)
  expect_identical(ad_eff$VUR, ad_eff$VR)
  expect_identical(ad_eff$GMM, ad_eff$YR)
  expect_identical(ad_eff$V_GMM, ad_eff$VR)
  expect_equal(ad_eff$corr, -sqrt(1 - ad_eff$VR / ad_eff$VU), tolerance = 1e-12)
  expect_equal(ad_eff$rho_aks_sq, 1 - ad_eff$VR / ad_eff$VU,  tolerance = 1e-12)
  # YR / YU / VR / VU are convention-independent
  expect_identical(ad_eff$YR, ad_emp$YR)
  expect_identical(ad_eff$YU, ad_emp$YU)
  expect_identical(ad_eff$VR, ad_emp$VR)
  expect_identical(ad_eff$VU, ad_emp$VU)
  # explicit FALSE honored
  expect_false(ad_emp$assume_efficient)

  # AUTO default (assume_efficient = NULL): the restricted fit here is a no-covariate, non-uniform
  # (default "efficient") fit, hence bound-attaining -> AUTO resolves to TRUE and reproduces the
  # efficient-variant numbers EXACTLY (the example.R VUR = VR convention).
  ad_auto <- edid_adaptive(ff$U, ff$R)
  expect_true(ad_auto$assume_efficient)
  expect_true(ad_auto$assume_efficient_auto)
  expect_identical(ad_auto$VUR,      ad_eff$VUR)
  expect_identical(ad_auto$GMM,      ad_eff$GMM)
  expect_identical(ad_auto$adaptive, ad_eff$adaptive)
  expect_identical(ad_auto$tO,       ad_eff$tO)
  expect_output(print(ad_auto), "AUTO")
  # ... while a non-bound-attaining restricted fit (uniform weights) AUTO-resolves to FALSE
  fitRu <- edid(df, "y", "id", "time", "g", pt_assumption = "all", weight_scheme = "uniform",
                aggregate = "event_study", cband = FALSE)
  ad_auto_u <- edid_adaptive(ff$U, fitRu)
  expect_false(ad_auto_u$assume_efficient)
  expect_true(ad_auto_u$assume_efficient_auto)

  # The two conventions coincide asymptotically when the restricted fit is
  # efficient (the empirical IF covariance VUR consistently estimates VR);
  # in finite samples they differ only by the sampling noise in VUR. On this
  # DGP (seed 6, n = 300): VUR/VR = 0.99848, overall adaptive gap = 3.99e-5
  # (estimates ~1.0, se_U ~ 0.067, i.e. ~6e-4 of an SE). Tolerance generous.
  expect_lt(abs(ad_eff$adaptive - ad_emp$adaptive), 0.02)
  expect_lt(abs(ad_eff$tO - ad_emp$tO), 0.2)

  # Per-e branch: the identity holds row-wise; gaps stay small (max observed
  # 3.78e-4 on this DGP).
  ade_eff <- edid_adaptive(ff$U, ff$R, parameter = "event_study", assume_efficient = TRUE)
  ade_emp <- edid_adaptive(ff$U, ff$R, parameter = "event_study", assume_efficient = FALSE)
  expect_identical(ade_eff$table$GMM, ade_eff$table$YR)
  expect_identical(ade_eff$table$VUR, ade_eff$table$VR)
  expect_lt(max(abs(ade_eff$table$adaptive - ade_emp$table$adaptive)), 0.02)

  # The convention is announced by print()
  expect_output(print(ad_eff), "assume_efficient = TRUE")

  # Validation: assume_efficient must be a single non-NA logical
  expect_error(edid_adaptive(ff$U, ff$R, assume_efficient = NA), "is.na")
})

test_that("shipped rds lookup tables match the vendored MissAdapt .mat files", {
  skip_if_not_installed("R.matlab")
  dir <- system.file("extdata", "aks_lookup", package = "did")
  expect_true(nzchar(dir))
  tab <- .edid_aks_lookup()
  policy     <- R.matlab::readMat(file.path(dir, "policy.mat"))
  thresholds <- R.matlab::readMat(file.path(dir, "thresholds.mat"))
  mse_data   <- R.matlab::readMat(file.path(dir, "emse_corr.mat"))
  expect_equal(tab$y_grid,  as.numeric(policy$y.grid),  tolerance = 0)
  expect_equal(tab$psi_mat, unname(policy$psi.mat),     tolerance = 0)
  expect_equal(tab$st, as.numeric(thresholds$st.mat),   tolerance = 0)
  expect_equal(tab$ht, as.numeric(thresholds$ht.mat),   tolerance = 0)
  expect_equal(tab$mse_lambda, as.numeric(mse_data$MSE.lambda.mat), tolerance = 0)
  expect_equal(tab$corr_grid, abs(tanh(seq(-3, -0.05, 0.05))), tolerance = 1e-12)
  # The corr grid is also stored (signed) by the authors inside emse_corr.mat
  expect_equal(abs(as.numeric(mse_data$Sigma.UO.grid)), tab$corr_grid, tolerance = 1e-12)
  # Provenance attributes embedded by data-raw/aks_lookup.R
  expect_identical(attr(tab, "source_repo"), "https://github.com/lsun20/MissAdapt")
  expect_identical(attr(tab, "source_commit"), "98d823a0818eebbec37ce7d1acf9ca0b78aee46b")
  expect_match(attr(tab, "license"), "MIT License")
  expect_match(attr(tab, "license"), "Copyright \\(c\\) 2023 Sophie Sun")
  expect_match(attr(tab, "grid_conventions")[["psi_mat"]], "rows = y_grid")
})

# ===========================================================================
# moment_set mechanism + edid_sargan
# ===========================================================================
test_that("edid() with moment_set = NULL is byte-identical to the default call", {
  df <- make_panel_toolkit(7L, n = 200L)
  f0 <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
             aggregate = "event_study", cband = FALSE)
  f1 <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
             aggregate = "event_study", cband = FALSE, moment_set = NULL)
  expect_identical(f0$att_gt$att, f1$att_gt$att)
  expect_identical(f0$att_gt$se,  f1$att_gt$se)
  expect_identical(f0$eif,        f1$eif)
  expect_identical(f0$event_study$att.egt, f1$event_study$att.egt)
  expect_identical(f0$event_study$se.egt,  f1$event_study$se.egt)
})

test_that("moment_set restricted to the PT-Post base reproduces the pt_assumption='post' fit", {
  df <- make_panel_toolkit(8L, n = 250L)
  base_ms <- data.frame(g = c(3, 5), gp = c(3, 5), tpre = c(2, 4))
  f_base <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
                 aggregate = "event_study", cband = FALSE, moment_set = base_ms)
  f_post <- edid(df, "y", "id", "time", "g", pt_assumption = "post",
                 aggregate = "event_study", cband = FALSE)
  # Algebraically identical moments (Y_t - Y_{g-1} comparisons); only float
  # association order differs.
  expect_equal(f_base$att_gt$att, f_post$att_gt$att, tolerance = 1e-10)
  expect_equal(f_base$att_gt$se,  f_post$att_gt$se,  tolerance = 1e-10)
  expect_equal(f_base$event_study$att.egt, f_post$event_study$att.egt, tolerance = 1e-10)
  # Pair restriction is visible in the cells
  expect_true(all(f_base$att_gt$n_pairs == 1L))
})

test_that("moment_set validation rejects malformed input and rows outside the enumeration are ignored", {
  df <- make_panel_toolkit(9L, n = 150L)
  expect_error(edid(df, "y", "id", "time", "g", moment_set = data.frame(a = 1)),
               "columns `g`, `gp`, `tpre`")
  expect_error(edid(df, "y", "id", "time", "g",
                    moment_set = data.frame(g = numeric(0), gp = numeric(0), tpre = numeric(0))),
               "zero rows")
  expect_error(edid(df, "y", "id", "time", "g",
                    moment_set = data.frame(g = 3, gp = "x", tpre = 1)),
               "must be numeric")
  # A row that is not a valid pair is ignored (intersection semantics): cells
  # for cohort 5 then have no pairs -> NA, cohort 3 keeps its restricted pair.
  ms <- data.frame(g = c(3, 5), gp = c(3, 5), tpre = c(2, 99))
  f  <- suppressWarnings(edid(df, "y", "id", "time", "g", pt_assumption = "all",
                              aggregate = "none", cband = FALSE, moment_set = ms))
  expect_true(all(is.na(f$att_gt$att[f$att_gt$group == 5])))
  expect_true(all(is.finite(f$att_gt$att[f$att_gt$group == 3])))
})

test_that(".edid_holm implements the Holm-Bonferroni step-down (hand-check, L = 3)", {
  # All rejected: thresholds 0.05/3, 0.05/2, 0.05/1
  h1 <- .edid_holm(c(0.001, 0.02, 0.04), alpha = 0.05)
  expect_identical(h1$rejected, c(TRUE, TRUE, TRUE))
  expect_equal(sort(h1$threshold), c(0.05 / 3, 0.05 / 2, 0.05 / 1), tolerance = 1e-15)

  # Step-down stop: smallest rejected, second exceeds its threshold -> stop,
  # so the third is NOT rejected even though p3 < alpha/(L+1-3) would hold.
  h2 <- .edid_holm(c(0.0166, 0.026, 0.04), alpha = 0.05)
  expect_identical(h2$rejected, c(TRUE, FALSE, FALSE))
  # Thresholds aligned to the original order: p1 is smallest -> alpha/3, etc.
  expect_equal(h2$threshold, c(0.05 / 3, 0.05 / 2, 0.05 / 1), tolerance = 1e-15)

  # Unordered input: rejection set maps back to original positions
  h3 <- .edid_holm(c(0.018, 0.001, 0.2), alpha = 0.05)
  expect_identical(h3$rejected, c(TRUE, TRUE, FALSE))
  expect_equal(h3$threshold, c(0.05 / 2, 0.05 / 3, 0.05 / 1), tolerance = 1e-15)
})

test_that("edid_sargan produces the documented step-down table on a staggered design", {
  df  <- make_panel_toolkit(10L, n = 300L)
  fit <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
              aggregate = "event_study", cband = FALSE)
  sg <- edid_sargan(fit, data = df)
  expect_s3_class(sg, "edid_sargan")

  # Candidate enumeration for cohorts {3, 5}, T = 6: full PT-All pairs minus
  # each cohort's own base, unique over (gp, tpre), ordered by (gp, tpre):
  # (3,1), (3,2), (5,1), (5,2), (5,3), (5,4)  [(3,2)/(5,4) are bases for their
  # own cohort but candidates for the other].
  expect_identical(sg$L, 6L)
  expect_identical(sg$table$gp,   c(3, 3, 5, 5, 5, 5))
  expect_identical(sg$table$tpre, c(1, 2, 1, 2, 3, 4))
  expect_identical(names(sg$table),
                   c("gp", "tpre", "H_statistic", "df", "p_value", "holm_threshold", "rejected"))
  expect_identical(sg$base, data.frame(g = c(3, 5), gp = c(3, 5), tpre = c(2, 4)))

  # Statistics well-formed; df = rank of the IF-difference covariance >= 1
  expect_true(all(is.finite(sg$table$H_statistic)) && all(sg$table$H_statistic >= 0))
  expect_true(all(sg$table$df >= 1L))
  expect_true(all(sg$table$p_value >= 0 & sg$table$p_value <= 1))

  # Holm thresholds: sorted p-values get alpha/(L+1-l)
  ord <- order(sg$table$p_value)
  expect_equal(sg$table$holm_threshold[ord], 0.05 / (6:1), tolerance = 1e-15)

  # Under the null DGP nothing should be (overwhelmingly) rejected; the
  # admissible set is the complement of the rejected set.
  expect_identical(sg$admissible,
                   sg$table[!sg$table$rejected, c("gp", "tpre"), drop = FALSE])
  expect_output(print(sg), "Incremental Sargan")
})

test_that("edid_sargan detects the violated moment under a PT-All violation", {
  skip_on_cran()
  # Period-1 shift for cohort 3: the contaminated restrictions are those using
  # t'' = 1 with cohort 3's pre-periods or its baseline; the base set (g-1
  # comparisons) is clean. The candidate (gp = 3, tpre = 1)-type restrictions
  # should be rejected far more often than clean ones.
  df  <- make_panel_toolkit(11L, n = 600L, viol = 1.0)
  fit <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
              aggregate = "event_study", cband = FALSE)
  sg <- edid_sargan(fit, data = df)
  expect_true(any(sg$table$rejected))
  # the directly contaminated candidate (3, 1) must be among the rejected
  expect_true(sg$table$rejected[sg$table$gp == 3 & sg$table$tpre == 1])
})

test_that("edid_sargan returns NULL with a message when the model is just-identified", {
  # One cohort at g = 2 with T = 3: PT-All pairs are self (2, 1) only -> no
  # candidate beyond the base.
  set.seed(12)
  n <- 80L
  coh <- sample(c(2, Inf), n, replace = TRUE)
  df <- data.frame(id = rep(seq_len(n), each = 3L), time = rep(1:3, n))
  df$g <- coh[df$id]
  df$y <- rnorm(n)[df$id] + 0.1 * df$time + rnorm(nrow(df), 0, 0.3) + (df$time >= df$g)
  fit <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
              aggregate = "event_study", cband = FALSE)
  expect_message(out <- edid_sargan(fit, data = df), "just-identified")
  expect_null(out)
})

# ===========================================================================
# edid_sargan refits: the fit's stored argument snapshot, not call replay
# ===========================================================================
test_that("edid_sargan refits use the fit's stored arguments, not the caller's mutated variables", {
  # Panel with two candidate covariates; y depends on x1 only.
  set.seed(77L)
  n   <- 80L
  coh <- sample(c(3, Inf), n, replace = TRUE)
  x1u <- rnorm(n); x2u <- rnorm(n)
  df  <- data.frame(id = rep(seq_len(n), each = 4L), time = rep(1:4, n))
  df$g  <- coh[df$id]
  df$x1 <- x1u[df$id]; df$x2 <- x2u[df$id]
  df$y  <- rnorm(n)[df$id] + 0.3 * df$time + 0.5 * df$x1 +
    1 * (df$time >= df$g) + rnorm(nrow(df), 0, 0.4)

  xf  <- ~ x1
  fit <- edid(df, "y", "id", "time", "g", xformla = xf, pt_assumption = "all",
              aggregate = "event_study", cband = FALSE, misspec_robust = FALSE)
  ref <- edid(df, "y", "id", "time", "g", xformla = ~ x1, pt_assumption = "all",
              aggregate = "event_study", cband = FALSE, misspec_robust = FALSE)
  s_ref <- edid_sargan(ref, data = df, inference = "plugin_fast")

  # Mutate the caller's variable AFTER fitting: previously the refits
  # re-evaluated `xformla = xf` in this environment and silently used ~ x2.
  xf <- ~ x2
  s_fit <- edid_sargan(fit, data = df, inference = "plugin_fast")
  expect_equal(s_fit$table, s_ref$table, tolerance = 1e-10)

  # Even with the variable gone, the stored snapshot carries the formula.
  rm(xf)
  expect_no_error(edid_sargan(fit, data = df, inference = "plugin_fast"))

  # Wrong data is loud, not silent: the refit sample must match the fit.
  expect_error(edid_sargan(fit, data = df[df$id <= 70L, ], inference = "plugin_fast"),
               "does not match the fitted sample")
})

test_that("edid_sargan runs on programmatically built fits (wrapper / do.call / lapply)", {
  df <- make_panel_toolkit(13L, n = 200L)

  # `...`-forwarding wrapper: the stored call holds `..1`-style promises that
  # cannot be re-evaluated ("..3 used in an incorrect context" previously).
  wrap  <- function(...) edid(...)
  fit_w <- wrap(df, "y", "id", "time", "g", pt_assumption = "all",
                aggregate = "event_study", cband = FALSE)
  expect_s3_class(edid_sargan(fit_w, data = df), "edid_sargan")

  # do.call-built fit
  fit_d <- do.call(edid, list(data = df, yname = "y", idname = "id", tname = "time",
                              gname = "g", pt_assumption = "all",
                              aggregate = "event_study", cband = FALSE))
  expect_s3_class(edid_sargan(fit_d, data = df), "edid_sargan")

  # Fit built inside lapply: the call references a lambda-local variable.
  fit_l <- lapply(list("y"), function(yn) {
    edid(df, yn, "id", "time", "g", pt_assumption = "all",
         aggregate = "event_study", cband = FALSE)
  })[[1L]]
  expect_s3_class(edid_sargan(fit_l, data = df), "edid_sargan")
})

test_that(".edid_refit_args returns the stored snapshot and falls back to call replay for legacy fits", {
  df  <- make_panel_toolkit(15L, n = 150L)
  fit <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
              aggregate = "event_study", cband = FALSE)
  args <- .edid_refit_args(fit)
  expect_identical(args, fit$args)
  expect_false("data" %in% names(args))
  expect_identical(args$yname, "y")
  expect_identical(args$pt_assumption, "all")
  expect_identical(args$weight_scheme, "efficient")
  expect_identical(args$min_pair_units, 5L)

  # Legacy fit (no $args, e.g. loaded from an old .rds): call re-evaluation,
  # which for a literal call recovers the supplied arguments.
  legacy <- fit
  legacy$args <- NULL
  args_legacy <- .edid_refit_args(legacy, envir = environment())
  expect_identical(args_legacy$pt_assumption, "all")
  expect_identical(args_legacy$aggregate, "event_study")
  expect_false("data" %in% names(args_legacy))
})

# ===========================================================================
# Conservative path reconciliation (the package's pt_assumption = "post" fit
# vs an inline port of the Application's conservative estimator + ES
# aggregation with the weight-IF correction)
# ===========================================================================
test_that("edid(pt='post') event study matches the Application conservative estimator", {
  df <- make_panel_2cohort(n_g3 = 30L, n_g5 = 30L, n_never = 40L, n_periods = 7L, seed = 14L)
  fU <- edid(df, "outcome", "unit", "time", "first_treat", pt_assumption = "post",
             aggregate = "event_study", cband = FALSE)

  # ---- inline port of Application/utils/conservative_did.R (no covariates) ----
  Y <- df$outcome; Time <- df$time; G <- df$first_treat; id <- df$unit
  G[!is.finite(G)] <- 0
  n_units <- length(unique(id))
  tlist <- sort(unique(Time)); glist <- sort(unique(G)); g_treated <- glist[glist > 0]
  ord1 <- order(id[Time == tlist[1]])
  G_cs <- (G[Time == tlist[1]])[ord1]
  pi_g <- sapply(glist, function(g) mean(G_cs == g)); names(pi_g) <- as.character(glist)
  gt <- expand.grid(g = g_treated, tt = tlist); gt <- gt[gt$tt >= gt$g, ]
  gt <- gt[order(gt$g, gt$tt), ]; rownames(gt) <- NULL
  est <- numeric(nrow(gt)); IF <- matrix(0, n_units, nrow(gt))
  Yw <- sapply(tlist, function(s) (Y[Time == s])[order(id[Time == s])])  # n x T wide, unit-sorted
  for (k in seq_len(nrow(gt))) {
    g <- gt$g[k]; tt <- gt$tt[k]; tb <- g - 1
    dg <- Yw[G_cs == g, which(tlist == tt)] - Yw[G_cs == g, which(tlist == tb)]
    d0 <- Yw[G_cs == 0, which(tlist == tt)] - Yw[G_cs == 0, which(tlist == tb)]
    est[k] <- mean(dg) - mean(d0)
    IF[G_cs == g, k] <- (Yw[G_cs == g, which(tlist == tt)] - Yw[G_cs == g, which(tlist == tb)] - mean(dg)) / pi_g[as.character(g)]
    IF[G_cs == 0, k] <- -(Yw[G_cs == 0, which(tlist == tt)] - Yw[G_cs == 0, which(tlist == tb)] - mean(d0)) / pi_g["0"]
  }
  # ES aggregation with the estimated-weight (wif) correction
  eseq <- sort(unique(gt$tt - gt$g))
  pg_gt <- pi_g[as.character(gt$g)]
  es_att <- numeric(length(eseq)); es_IF <- matrix(0, n_units, length(eseq))
  for (j in seq_along(eseq)) {
    idx <- which(gt$tt - gt$g == eseq[j])
    pge <- pg_gt[idx] / sum(pg_gt[idx])
    es_att[j] <- sum(est[idx] * pge)
    sum_pg <- sum(pg_gt[idx])
    if1 <- sapply(idx, function(k) ((G_cs == gt$g[k]) - pi_g[as.character(gt$g[k])]) / sum_pg)
    if2s <- rowSums(sapply(idx, function(k) (G_cs == gt$g[k]) - pi_g[as.character(gt$g[k])]))
    if2 <- if2s %*% t(pg_gt[idx] / sum_pg^2)
    es_IF[, j] <- IF[, idx, drop = FALSE] %*% pge + (if1 - if2) %*% est[idx]
  }
  epos <- which(eseq >= 0)
  overall <- mean(es_att[epos])
  overall_IF <- as.numeric(es_IF[, epos, drop = FALSE] %*% rep(1 / length(epos), length(epos)))
  # ---- end inline port ----

  aU <- fU$event_study
  pos <- match(eseq[epos], aU$egt)
  expect_equal(aU$att.egt[pos], unname(es_att[epos]), tolerance = 1e-10)
  expect_equal(aU$overall.att, overall, tolerance = 1e-10)
  expect_equal(as.numeric(aU$inf.function$dynamic.inf.func), overall_IF, tolerance = 1e-8)
  expect_equal(unname(as.matrix(aU$inf.function$dynamic.inf.func.e[, pos])),
               unname(es_IF[, epos]), tolerance = 1e-8)
})
