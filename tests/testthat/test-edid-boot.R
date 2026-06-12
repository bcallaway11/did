library(testthat)

# ===========================================================================
# Finite-sample bootstrap tools for edid fits (R/edid-boot.R):
#   edid_refit_bootstrap()        -- nuisance-refitting nonparametric cluster bootstrap
#   edid_perturbation_bootstrap() -- sieve-coefficient perturbation bootstrap (no refit)
# Calibration provenance: the Chen-Sant'Anna-Xie efficient-DiD inference study
# (refit bootstrap = the small-n remedy on weak-overlap long-horizon cells;
# perturbation = ~88% of that gain at matmul cost).
# ===========================================================================

# Small staggered panel with a covariate; well-behaved overlap. Cohorts 2, 3,
# never-treated; T = 4. True ATT = 1 for every post cell.
make_panel_boot <- function(seed, n = 300L, Tt = 4L, clustered = FALSE) {
  set.seed(seed)
  coh <- sample(c(2, 3, Inf), n, replace = TRUE, prob = c(.3, .3, .4))
  df  <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  x    <- rnorm(n); df$x1 <- x[df$id]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + 0.3 * df$x1 * df$time +
    1 * (df$time >= df$g) + rnorm(n * Tt, 0, .5)
  if (clustered) df$cl <- ((df$id - 1L) %% 30L) + 1L
  df
}

# The study's weak-overlap DGP (edid_inference_tests, sim_panel): nonlinear
# propensity index in x1 => strong overlap deterioration in the tails; the
# long-horizon cell ATT(2,4) is the documented hard case. True ATT = 1.
make_panel_weak_overlap <- function(seed, n = 300L, Tn = 4L) {
  set.seed(seed)
  x1  <- runif(n, -2, 2)
  eta <- 1.1 * x1 + 0.7 * x1^2 - 0.5
  P   <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gcat  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alpha <- rnorm(n, 0.5 * x1, 1)
  do.call(rbind, lapply(seq_len(Tn), function(tt) {
    ht  <- (tt - 1) * (0.5 * x1 + 0.45 * x1^2)
    tau <- ifelse(is.finite(gcat) & tt >= gcat, 1.0, 0)
    data.frame(id = seq_len(n), time = tt, g = gcat, x1 = x1,
               y = alpha + 0.3 * tt + ht + tau + rnorm(n))
  }))
}

fit_boot_uniform <- function(df, ...) {
  edid(df, "y", "id", "time", "g", xformla = ~ x1, weight_scheme = "uniform",
       aggregate = "event_study", cband = FALSE, ...)
}

# ===========================================================================
# edid_refit_bootstrap: clean run, sane SEs, documented structure
# ===========================================================================
test_that("refit bootstrap runs clean on a small staggered sim with SEs near the analytic SEs", {
  df  <- make_panel_boot(1L)
  fit <- fit_boot_uniform(df)
  rb  <- expect_no_warning(
    edid_refit_bootstrap(fit, data = df, B = 29L, seed = 42L, agg = "event_study"))

  expect_s3_class(rb, "edid_refit_bootstrap")
  expect_true(all(c("att_gt", "aggregates", "B", "n_failed", "failed_messages",
                    "resample", "n_resample_units", "alpha", "seed", "call") %in% names(rb)))
  expect_identical(rb$B, 29L)
  expect_identical(rb$n_failed, 0L)
  expect_identical(rb$resample, "unit")
  expect_identical(rb$n_resample_units, 300L)
  expect_identical(rb$alpha, fit$alpha)

  # cell table aligned on the fit, full draw count, positive bootstrap SEs
  expect_identical(rb$att_gt$group, fit$att_gt$group)
  expect_identical(rb$att_gt$time,  fit$att_gt$time)
  expect_equal(rb$att_gt$att, fit$att_gt$att)
  expect_equal(rb$att_gt$se_analytic, fit$att_gt$se)
  expect_true(all(rb$att_gt$n_boot == 29L))
  expect_true(all(is.finite(rb$att_gt$se_boot) & rb$att_gt$se_boot > 0))

  # bootstrap SEs in a 0.4x-3x band of the analytic SEs (sanity, not calibration)
  ratio <- rb$att_gt$se_boot / rb$att_gt$se_analytic
  expect_true(all(ratio > 0.4 & ratio < 3))

  # symmetric normal-quantile CI around att (the validated convention) + percentile CI
  z <- qnorm(1 - rb$alpha / 2)
  expect_equal(rb$att_gt$ci_lower, rb$att_gt$att - z * rb$att_gt$se_boot)
  expect_equal(rb$att_gt$ci_upper, rb$att_gt$att + z * rb$att_gt$se_boot)
  expect_true(all(rb$att_gt$pct_lower <= rb$att_gt$pct_upper))

  # event-study aggregate: rows e<egt> plus overall, aligned with the fit's AGGTEobj
  es <- rb$aggregates$event_study
  expect_identical(es$parameter, c(paste0("e", fit$event_study$egt), "overall"))
  expect_equal(es$att, c(fit$event_study$att.egt, fit$event_study$overall.att))
  expect_true(all(is.finite(es$se_boot) & es$se_boot > 0))
  expect_true(all(es$n_boot == 29L))

  expect_output(print(rb), "Nuisance-refitting cluster bootstrap")
})

# ===========================================================================
# Reproducibility: same seed identical; cores = 1 vs cores = 2 identical
# ===========================================================================
test_that("refit bootstrap is seed-reproducible and cores-invariant", {
  df  <- make_panel_boot(2L, n = 150L)
  fit <- fit_boot_uniform(df)

  r1 <- edid_refit_bootstrap(fit, data = df, B = 9L, seed = 5L, agg = "event_study")
  r2 <- edid_refit_bootstrap(fit, data = df, B = 9L, seed = 5L, agg = "event_study")
  expect_identical(r1$att_gt, r2$att_gt)
  expect_identical(r1$aggregates, r2$aggregates)

  # cores > 1 is numerically identical (per-draw seeds seed + b; draws whose
  # forked worker dies -- e.g. fork-unsafe BLAS -- are recomputed serially with
  # the same seed, with a warning, so identity holds regardless)
  r3 <- suppressWarnings(
    edid_refit_bootstrap(fit, data = df, B = 9L, seed = 5L, cores = 2L, agg = "event_study"))
  expect_identical(r1$att_gt, r3$att_gt)
  expect_identical(r1$aggregates, r3$aggregates)

  # different seeds give different draws
  r4 <- edid_refit_bootstrap(fit, data = df, B = 9L, seed = 6L, agg = "event_study")
  expect_false(isTRUE(all.equal(r1$att_gt$se_boot, r4$att_gt$se_boot)))

  # a seeded call does not disturb the caller's RNG stream
  set.seed(99L); a <- runif(1)
  set.seed(99L)
  invisible(edid_refit_bootstrap(fit, data = df, B = 2L, seed = 5L, agg = "overall"))
  expect_identical(a, runif(1))
})

# ===========================================================================
# Clustered fits resample whole clusters
# ===========================================================================
test_that("refit bootstrap resamples clusters for a clustered fit", {
  df  <- make_panel_boot(3L, n = 150L, clustered = TRUE)
  fit <- edid(df, "y", "id", "time", "g", clustervars = "cl",
              aggregate = "event_study", cband = FALSE)
  rb  <- edid_refit_bootstrap(fit, data = df, B = 9L, seed = 7L, agg = "event_study")

  expect_identical(rb$resample, "cluster")
  expect_identical(rb$n_resample_units, 30L)   # 30 clusters, not 150 units
  post <- !rb$att_gt$is_pre
  expect_true(all(is.finite(rb$att_gt$se_boot[post])))
})

# ===========================================================================
# Failed-draw accounting: a constructed degenerate design must warn
# ===========================================================================
test_that("refit bootstrap accounts for failed draws and warns above 5 percent", {
  # exactly ONE never-treated unit: a resample that drops it has no comparison
  # group at all, so the inner edid() refit fails for that draw
  set.seed(4L)
  n <- 40L; Tt <- 3L
  coh <- c(rep(2, n - 1L), Inf)
  df  <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) + rnorm(n * Tt, 0, .5)
  # the base fit itself warns about the single never-treated unit -- by construction
  fit <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE))

  expect_warning(
    rb <- edid_refit_bootstrap(fit, data = df, B = 29L, seed = 1L, agg = "overall"),
    "failed entirely")
  expect_gt(rb$n_failed, 0.05 * rb$B)
  expect_true(length(rb$failed_messages) >= 1L)
  expect_match(rb$failed_messages[1L], "never-treated", ignore.case = TRUE)
  # surviving draws still feed the SEs; the per-coordinate draw count reflects the losses
  expect_true(all(rb$att_gt$n_boot <= rb$B - rb$n_failed))
})

# ===========================================================================
# Argument validation + data recovery
# ===========================================================================
test_that("refit bootstrap validates inputs and recovers data from the fit's call", {
  df  <- make_panel_boot(5L, n = 120L)
  fit <- fit_boot_uniform(df)

  expect_error(edid_refit_bootstrap(list()), "edid_fit")
  expect_error(edid_refit_bootstrap(fit, data = df, B = 1L), "`B` must be")
  expect_error(edid_refit_bootstrap(fit, data = df, B = 9L, cores = 0L), "`cores` must be")
  expect_error(edid_refit_bootstrap(fit, data = df, B = 9L, seed = 1.5), "seed")
  expect_error(edid_refit_bootstrap(fit, data = df, B = 9L, seed = 3e9), "seed")

  # update() idiom: data = NULL re-evaluates the fit call's data in the caller
  r0 <- edid_refit_bootstrap(fit, B = 5L, seed = 11L, agg = "overall")
  r1 <- edid_refit_bootstrap(fit, data = df, B = 5L, seed = 11L, agg = "overall")
  expect_identical(r0$att_gt, r1$att_gt)
})

test_that("bootstrap tools refit from the fit's stored args (wrapper-built fits, mutated variables)", {
  df <- make_panel_boot(9L, n = 120L)

  # `...`-forwarding wrapper: the stored call holds `..N` promises, which the
  # legacy call re-evaluation could not handle ("..3 used in an incorrect
  # context"); the per-draw refits now consume fit$args.
  wrap <- function(...) edid(...)
  fit_w <- wrap(df, "y", "id", "time", "g", xformla = ~ x1, weight_scheme = "uniform",
                aggregate = "event_study", cband = FALSE)
  rb <- edid_refit_bootstrap(fit_w, data = df, B = 5L, seed = 3L, agg = "overall")
  expect_s3_class(rb, "edid_refit_bootstrap")
  pb <- edid_perturbation_bootstrap(fit_w, data = df, B = 9L, seed = 3L, agg = "overall")
  expect_s3_class(pb, "edid_perturbation_bootstrap")
  expect_identical(pb$weight_scheme, "uniform")

  # Mutating a caller variable used in the original call must not change the
  # refit configuration (previously it was silently re-evaluated).
  xf  <- ~ x1
  fit <- edid(df, "y", "id", "time", "g", xformla = xf, weight_scheme = "uniform",
              aggregate = "event_study", cband = FALSE)
  ref <- edid_refit_bootstrap(fit, data = df, B = 5L, seed = 13L, agg = "overall")
  xf <- ~ I(x1^3)
  mut <- edid_refit_bootstrap(fit, data = df, B = 5L, seed = 13L, agg = "overall")
  expect_identical(ref$att_gt, mut$att_gt)
  expect_identical(ref$aggregates, mut$aggregates)
})

# ===========================================================================
# edid_perturbation_bootstrap: uniform and efficient covariate fits
# ===========================================================================
test_that("perturbation bootstrap runs clean on a uniform-weight covariate fit", {
  df  <- make_panel_boot(6L)
  fit <- fit_boot_uniform(df)
  pb  <- expect_no_warning(
    edid_perturbation_bootstrap(fit, data = df, B = 99L, seed = 8L, agg = "event_study"))

  expect_s3_class(pb, "edid_perturbation_bootstrap")
  expect_identical(pb$B, 99L)
  expect_identical(pb$n_failed, 0L)
  expect_identical(pb$weight_scheme, "uniform")
  expect_gt(pb$n_nuisances, 0L)

  # cells aligned on the fit; combined SE = sqrt(plug^2 + pert^2) >= plug
  expect_equal(pb$att_gt$att, fit$att_gt$att)
  expect_equal(pb$att_gt$se_analytic, fit$att_gt$se)
  expect_true(all(pb$att_gt$n_pert == 99L))
  expect_true(all(is.finite(pb$att_gt$se_plug) & pb$att_gt$se_plug > 0))
  expect_true(all(is.finite(pb$att_gt$se_pert) & pb$att_gt$se_pert > 0))
  expect_equal(pb$att_gt$se_combined,
               sqrt(pb$att_gt$se_plug^2 + pb$att_gt$se_pert^2))
  expect_true(all(pb$att_gt$se_combined >= pb$att_gt$se_plug))
  # the perturbation simulates a higher-order channel: smaller than the
  # first-order SE, but a real correction (sanity band)
  expect_true(all(pb$att_gt$se_pert < pb$att_gt$se_plug))

  # Wald CI on the combined SE
  z <- qnorm(1 - pb$alpha / 2)
  expect_equal(pb$att_gt$ci_lower, pb$att_gt$att - z * pb$att_gt$se_combined)
  expect_equal(pb$att_gt$ci_upper, pb$att_gt$att + z * pb$att_gt$se_combined)

  # aggregates: exact linear map of the cells, rows e<egt> plus overall
  es <- pb$aggregates$event_study
  expect_identical(es$parameter, c(paste0("e", fit$event_study$egt), "overall"))
  expect_equal(es$att, c(fit$event_study$att.egt, fit$event_study$overall.att),
               tolerance = 1e-6)
  expect_true(all(es$se_combined >= es$se_plug))

  expect_output(print(pb), "perturbation bootstrap")
})

test_that("perturbation bootstrap supports the efficient scheme and is reproducible", {
  df  <- make_panel_boot(7L, n = 200L)
  fit <- edid(df, "y", "id", "time", "g", xformla = ~ x1, weight_scheme = "efficient",
              aggregate = "event_study", cband = FALSE)

  p1 <- edid_perturbation_bootstrap(fit, data = df, B = 99L, seed = 9L, agg = "event_study")
  expect_identical(p1$weight_scheme, "efficient")
  expect_identical(p1$n_failed, 0L)
  expect_true(all(is.finite(p1$att_gt$se_pert) & p1$att_gt$se_pert > 0))
  expect_true(all(p1$att_gt$se_combined >= p1$att_gt$se_plug))

  # same seed identical; cores-invariant
  p2 <- edid_perturbation_bootstrap(fit, data = df, B = 99L, seed = 9L, agg = "event_study")
  expect_identical(p1$att_gt, p2$att_gt)
  expect_identical(p1$aggregates, p2$aggregates)
  p3 <- suppressWarnings(
    edid_perturbation_bootstrap(fit, data = df, B = 99L, seed = 9L, cores = 2L,
                                agg = "event_study"))
  expect_identical(p1$att_gt, p3$att_gt)
})

test_that("perturbation bootstrap errors informatively where its construction does not apply", {
  df <- make_panel_boot(8L, n = 150L)

  # unsupported weight scheme -> informative error pointing to the refit bootstrap
  fa <- edid(df, "y", "id", "time", "g", xformla = ~ x1, weight_scheme = "averaged",
             aggregate = "none", cband = FALSE)
  expect_error(edid_perturbation_bootstrap(fa, data = df, B = 9L, seed = 1L),
               "supports weight_scheme = 'efficient' or 'uniform'")
  expect_error(edid_perturbation_bootstrap(fa, data = df, B = 9L, seed = 1L),
               "edid_refit_bootstrap")

  # no covariates -> no sieve coefficients to perturb
  fn <- edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE)
  expect_error(edid_perturbation_bootstrap(fn, data = df, B = 9L, seed = 1L),
               "requires a covariate fit")

  # changed data -> the exactness guard refuses to mis-center the perturbation
  fit <- fit_boot_uniform(df)
  df2 <- df; df2$y <- df2$y + rnorm(nrow(df2), 0, 0.3)
  expect_error(edid_perturbation_bootstrap(fit, data = df2, B = 9L, seed = 1L),
               "could not reproduce")

  expect_error(edid_perturbation_bootstrap(list()), "edid_fit")
  expect_error(edid_perturbation_bootstrap(fit, data = df, B = 1L), "`B` must be")
})

# ===========================================================================
# Calibration smoke (skip_on_cran): on the study's weak-overlap design the
# refit bootstrap CI is wider than the plug-in analytic CI on the long-horizon
# cell ATT(2,4) -- the documented small-n under-coverage remedy.
# ===========================================================================
test_that("refit bootstrap widens the long-horizon CI on a weak-overlap design", {
  skip_on_cran()
  ratios <- vapply(1:5, function(s) {
    df  <- make_panel_weak_overlap(20260600L + s)
    fit <- suppressWarnings(
      edid(df, "y", "id", "time", "g", xformla = ~ x1, weight_scheme = "uniform",
           aggregate = "none", cband = FALSE, misspec_robust = FALSE))
    rb  <- suppressWarnings(
      edid_refit_bootstrap(fit, data = df, B = 39L, seed = 1000L + s, agg = "overall"))
    i24 <- which(rb$att_gt$group == 2 & rb$att_gt$time == 4)
    ciw_boot     <- rb$att_gt$ci_upper[i24] - rb$att_gt$ci_lower[i24]
    ciw_analytic <- 2 * qnorm(1 - rb$alpha / 2) * rb$att_gt$se_analytic[i24]
    ciw_boot / ciw_analytic
  }, numeric(1L))
  # per-design ratios fluctuate; the average widening is the documented effect
  # (study: bootstrap-SE/plug-in-SE ~ 1.2-1.4 on this cell at n in the hundreds)
  expect_gt(mean(ratios), 1.1)
  expect_true(all(ratios > 0.4 & ratios < 3))
})
