# test-edid-thin-cohort.R
# The thin-cohort guard (`min_pair_units`).
#
# Provenance (audit repo Efficient_DiD_Claude, quality_reports/drafts/):
#   verify_imputation_nesting.R Section 7 and verify_dominance_did2s.R design 10
# established that under pt_assumption = "all" a 3-unit cohort at n = 2000 produces
# overidentified cells whose analytic SE understates the true sampling SD by up to
# 9-25x (cell coverage 0.10-0.71; sup-t 0.08), that the "efficient" cell is NOISIER
# than the just-identified one, and that with a 1-unit cohort the contamination
# spills over into healthy cohorts' cells (ATT(3,3) = -2.07, se 0.04, truth 1).
# Restricting to the just-identified moment restores calibration; uniform weights
# do not. The guard implements exactly that restriction:
#   * thin TARGET cohort  -> cells pinned to the just-identified (pt = "post") moment;
#   * thin COMPARISON cohort -> its cross pairs excised from other cohorts' cells.

# ---------------------------------------------------------------------------
# DGP helpers (no-covariate (3,4,4) design of nesting Section 7; tau = 1)
# ---------------------------------------------------------------------------
gen_thin_outcomes <- function(G, T_per = 4L) {
  n <- length(G)
  eta <- rnorm(n) + 0.4 * (G == 3L) - 0.3 * (G == 4L)
  Y <- eta + matrix(0.3 * seq_len(T_per), n, T_per, byrow = TRUE) +
    matrix(rnorm(n * T_per), n, T_per)
  for (g in c(3L, 4L)) for (t in seq_len(T_per)) if (t >= g) Y[G == g, t] <- Y[G == g, t] + 1
  data.frame(id   = rep(seq_len(n), each = T_per),
             time = rep(seq_len(T_per), times = n),
             y    = as.vector(t(Y)),
             gvar = rep(G, each = T_per))
}
make_thin_panel <- function(n = 400L, n_thin = 3L, seed = 20260611L, T_per = 4L) {
  set.seed(seed)
  G <- sample(c(3L, 0L), n, TRUE, c(0.4, 0.6))
  G[seq_len(n_thin)] <- 4L
  gen_thin_outcomes(G, T_per)
}
fit_thin <- function(df, pt = "all", ...) {
  edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
       pt_assumption = pt, weight_scheme = "efficient",
       aggregate = "none", cband = FALSE, ...)
}
cell_of <- function(fit, g, t) fit$att_gt[fit$att_gt$group == g & fit$att_gt$time == t, ]
expect_identical_unless_covr <- function(object, expected) {
  if (identical(Sys.getenv("R_COVR"), "true")) {
    expect_equal(object, expected, tolerance = 0, ignore_attr = TRUE)
  } else {
    expect_identical(object, expected)
  }
}

# Legacy fingerprints: fit_thin() output of the PRE-GUARD package (commit 5133060,
# the tip of pr/260 before the thin-cohort guard) on the seeded designs below,
# captured with dput(control = c("all", "hexNumeric")) so the doubles round-trip
# exactly. Row order: (3,2) (3,3) (3,4) (4,2) (4,3) (4,4).
legacy_thin3_att <- c(0x1.bc0944d642ecep-3, 0x1.41a35a08e8cf4p+0, 0x1.1092891263p+0,
                      0x1.082cefd2a1fcp+0, -0x1.f2f3a4e9b852ap-3, 0x1.e331f106a6ad6p-1)
legacy_thin3_se  <- c(0x1.2cb1e3b314381p-3, 0x1.dd919addf68d4p-4, 0x1.f97d8e96ecd4bp-4,
                      0x1.f0e2a421f38adp-1, 0x1.a7aecb2240933p-3, 0x1.385b3537f9999p-4)
legacy_healthy_att <- c(0x1.80a8bb02c2241p-4, 0x1.30430f4584176p+0, 0x1.1df68f8021b34p+0,
                        0x1.83194cef211d4p-5, 0x1.5c0c5b152a07ep-3, 0x1.0bc69aa543607p+0)
legacy_healthy_se  <- c(0x1.3063c054803a6p-3, 0x1.e7f74ff8545eep-4, 0x1.060dabbbaad95p-3,
                        0x1.c501ee53d809p-3, 0x1.b39842b70f8b9p-3, 0x1.696af518bcd55p-3)

# ===========================================================================
# Guard helper: pure unit tests
# ===========================================================================
test_that("apply_thin_cohort_guard_edid pins a thin target to the just-identified self pair", {
  pairs <- data.frame(gp = c(4, 4, 4, 3), tpre = c(1, 2, 3, 2))   # target 4 in a (3,4,4) design
  sizes <- c("3" = 100, "4" = 3)
  out <- apply_thin_cohort_guard_edid(4, pairs, sizes, 5L, "all")
  expect_true(out$degraded)
  expect_identical(out$pairs, data.frame(gp = 4, tpre = 3))       # self pair, max tpre = g - 1
  expect_length(out$excised_gp, 0L)
})

test_that("apply_thin_cohort_guard_edid excises thin comparison cohorts from a healthy target", {
  pairs <- data.frame(gp = c(3, 3, 4, 4), tpre = c(1, 2, 2, 3))   # target 3 in a (3,4,4) design
  sizes <- c("3" = 100, "4" = 3)
  out <- apply_thin_cohort_guard_edid(3, pairs, sizes, 5L, "all")
  expect_false(out$degraded)
  expect_identical(out$excised_gp, 4)
  expect_identical(out$pairs, data.frame(gp = c(3, 3), tpre = c(1, 2)))
})

test_that("apply_thin_cohort_guard_edid is inert under pt = 'post' and above the threshold", {
  pairs_post <- data.frame(gp = Inf, tpre = 3)
  sizes <- c("3" = 100, "4" = 1)
  out <- apply_thin_cohort_guard_edid(4, pairs_post, sizes, 5L, "post")
  expect_identical(out$pairs, pairs_post)                          # untouched (same object)
  expect_false(out$degraded)
  pairs_all <- data.frame(gp = c(3, 3, 4, 4), tpre = c(1, 2, 2, 3))
  out2 <- apply_thin_cohort_guard_edid(3, pairs_all, c("3" = 5, "4" = 5), 5L, "all")
  expect_identical(out2$pairs, pairs_all)
  expect_false(out2$degraded)
  expect_length(out2$excised_gp, 0L)
})

# ===========================================================================
# 3-unit cohort: guard semantics, flags, warnings, pt = "post" equivalence
# ===========================================================================
test_that("3-unit cohort: thin target pinned to just-identified moment, pairs excised elsewhere", {
  df <- make_thin_panel(400L, 3L)
  w <- capture_warnings(fit <- fit_thin(df))

  # (a) loud, specific warnings: degraded target + excised comparison, both naming
  # cohort 4 with its unit count and recommending edid_refit_bootstrap()
  expect_true(any(grepl("Thin-cohort guard: treated cohort\\(s\\) 4 \\(3 units\\)", w) &
                  grepl("just-identified moment", w) &
                  grepl("analytic SEs for overidentified efficient weighting are unreliable", w) &
                  grepl("edid_refit_bootstrap", w)))
  expect_true(any(grepl("Thin-cohort guard: comparison cohort\\(s\\) 4 \\(3 units\\)", w) &
                  grepl("excised from the moment sets of target cohort\\(s\\) 3", w) &
                  grepl("edid_refit_bootstrap", w)))
  # the legacy "<2 units" warning does not apply here (3 units) and must not fire
  expect_false(any(grepl("fewer than 2 units", w)))

  # moment sets: cohort-4 cells just-identified (self pair (4,3)); cohort-3 cells
  # keep only their self pairs (cross pairs with the thin cohort excised)
  expect_identical(fit$att_gt$n_pairs, c(2L, 2L, 2L, 1L, 1L, 1L))
  k44 <- which(vapply(fit$cells, function(x) x$group == 4 && x$time == 4, logical(1L)))
  expect_identical(fit$cells[[k44]]$pairs, data.frame(gp = 4, tpre = 3))
  k33 <- which(vapply(fit$cells, function(x) x$group == 3 && x$time == 3, logical(1L)))
  expect_identical(fit$cells[[k33]]$pairs, data.frame(gp = c(3, 3), tpre = c(1, 2)))

  # flags: per-cell thin_cohort_degraded on the thin cohort only; fit-level record
  flags <- vapply(fit$cells, function(x) isTRUE(x$thin_cohort_degraded), logical(1L))
  expect_identical(flags, fit$att_gt$group == 4)
  expect_identical(
    fit$thin_cohorts,
    data.frame(cohort = 4, n_units = 3L, degraded_target = TRUE, excised_comparison = TRUE))

  # the degraded cells ARE the pt = "post" cells (never-treated comparison, base
  # period g-1) up to floating-point reassociation -- the calibrated moment
  fpost <- suppressWarnings(fit_thin(df, pt = "post"))
  i4 <- fit$att_gt$group == 4
  expect_equal(fit$att_gt$att[i4], fpost$att_gt$att[fpost$att_gt$group == 4], tolerance = 1e-10)
  expect_equal(fit$att_gt$se[i4],  fpost$att_gt$se[fpost$att_gt$group == 4],  tolerance = 1e-10)
})

test_that("the legacy '<2 units' warning is subsumed under pt='all' and kept under pt='post'", {
  df1 <- make_thin_panel(400L, 1L)
  w_all <- capture_warnings(fit_all <- fit_thin(df1))
  expect_false(any(grepl("fewer than 2 units", w_all)))            # subsumed by the guard
  expect_true(any(grepl("Thin-cohort guard: treated cohort\\(s\\) 4 \\(1 unit\\)", w_all)))
  w_post <- capture_warnings(fit_post <- fit_thin(df1, pt = "post"))
  expect_true(any(grepl("fewer than 2 units", w_post)))            # guard inert: legacy warning stays
  expect_false(any(grepl("Thin-cohort guard", w_post)))
  expect_null(fit_post$thin_cohorts)
})

# ===========================================================================
# (b) 1-unit cohort: healthy cells byte-identical to excluding the singleton's pairs
# ===========================================================================
test_that("1-unit cohort: healthy cells byte-identical to a fit excluding the singleton's pairs", {
  df <- make_thin_panel(400L, 1L)
  w <- capture_warnings(fg <- fit_thin(df))                        # guarded fit (default 5)
  expect_true(any(grepl("Thin-cohort guard: treated cohort\\(s\\) 4 \\(1 unit\\)", w)))
  expect_true(any(grepl("comparison cohort\\(s\\) 4 \\(1 unit\\)", w) &
                  grepl("target cohort\\(s\\) 3", w)))

  # comparison fit: the SAME exclusion expressed through the documented moment_set
  # mechanism -- cohort 3 keeps only its self pairs (no gp = 4 pairs anywhere),
  # cohort 4 the just-identified (4,3) moment
  ms <- rbind(data.frame(g = 3, gp = 3, tpre = c(1, 2)),
              data.frame(g = 4, gp = 4, tpre = 3))
  fms <- suppressWarnings(fit_thin(df, moment_set = ms, min_pair_units = 2L))
  expect_identical(fg$att_gt, fms$att_gt)                          # entire cell table, byte-identical
  expect_identical(fg$eif,    fms$eif)                             # influence functions too

  # spillover gone: the audited 1-unit failure had ATT(3,3) at -2.07 (truth 1, se
  # 0.04); post-guard the healthy cohort's cell must sit within 4 SEs of the truth
  c33 <- cell_of(fg, 3, 3)
  expect_lt(abs(c33$att - 1), 4 * c33$se)
  # and the healthy cells are free of the singleton: no gp = 4 pair anywhere
  for (k in seq_along(fg$cells)) {
    if (fg$att_gt$group[k] == 3) expect_false(any(fg$cells[[k]]$pairs$gp == 4))
  }
})

# ===========================================================================
# (c) healthy design: the guard is inert (byte-identical to the legacy package)
# ===========================================================================
test_that("healthy design (all cohorts >= 5): guard inert, byte-identical to the legacy fit", {
  df <- make_thin_panel(400L, 40L)                                 # cohort 4 has 40 units
  # nocov_shrink = FALSE: the legacy fingerprints were captured on the pre-guard
  # package (commit 5133060), which had no pole-target shrinkage either; the
  # legacy-reproduction contract is the unshrunk pipeline (nocov_shrink = FALSE
  # is documented to reproduce it bit-for-bit).
  expect_no_warning(f5 <- fit_thin(df, nocov_shrink = FALSE))      # no guard, no legacy warnings
  expect_identical_unless_covr(f5$att_gt$att, legacy_healthy_att)  # pre-guard package, bit-for-bit
  expect_identical_unless_covr(f5$att_gt$se,  legacy_healthy_se)
  expect_null(f5$thin_cohorts)
  expect_false(any(vapply(f5$cells, function(x) isTRUE(x$thin_cohort_degraded), logical(1L))))
  # min_pair_units = 2 is equally inert here: identical fit
  f2 <- fit_thin(df, min_pair_units = 2L, nocov_shrink = FALSE)
  expect_identical(f5$att_gt, f2$att_gt)
  expect_identical(f5$eif,    f2$eif)
  # the guard is equally inert under the (default) shrinkage: same pair sets and
  # thin-cohort bookkeeping, only the weight regularization differs
  expect_no_warning(f5s <- fit_thin(df))
  expect_null(f5s$thin_cohorts)
  expect_identical(f5s$att_gt$n_pairs, f5$att_gt$n_pairs)
})

# ===========================================================================
# (d) min_pair_units = 2 reproduces the legacy behavior bit-for-bit
# ===========================================================================
test_that("min_pair_units = 2 reproduces the pre-guard fit bit-for-bit on the 3-unit design", {
  df <- make_thin_panel(400L, 3L)
  # nocov_shrink = FALSE: legacy fingerprints predate the pole-target shrinkage
  # (see the healthy-design test above for the contract).
  w <- capture_warnings(f2 <- fit_thin(df, min_pair_units = 2L, nocov_shrink = FALSE))
  expect_false(any(grepl("Thin-cohort guard", w)))                 # guard never fires at 2 for 3 units
  expect_identical_unless_covr(f2$att_gt$att, legacy_thin3_att)    # pre-guard package, bit-for-bit
  expect_identical_unless_covr(f2$att_gt$se,  legacy_thin3_se)
  expect_identical(f2$att_gt$n_pairs, rep(4L, 6L))                 # full overidentified moment sets
  expect_null(f2$thin_cohorts)
  expect_identical(f2$min_pair_units, 2L)
})

test_that("min_pair_units is validated (integer scalar >= 2)", {
  df <- make_thin_panel(120L, 6L)
  expect_error(fit_thin(df, min_pair_units = 1L),  "min_pair_units")
  expect_error(fit_thin(df, min_pair_units = 2.5), "min_pair_units")
  expect_error(fit_thin(df, min_pair_units = NA),  "min_pair_units")
})

# ===========================================================================
# guarded fits keep working downstream: aggregations + refit bootstrap call
# ===========================================================================
test_that("aggregations run on a guarded fit and recommend-able tools accept it", {
  df <- make_thin_panel(400L, 3L)
  fit <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", aggregate = "all", cband = TRUE, seed = 1L))
  expect_true(is.finite(fit$overall$overall.att))
  expect_true(is.finite(fit$simple$overall.att))
  expect_true(all(is.finite(fit$event_study$att.egt[fit$event_study$egt >= 0])))
})

# ===========================================================================
# (a) 3-unit cohort, n = 1500: refit-bootstrap agreement + MC calibration
# ===========================================================================
test_that("3-unit cohort: degraded cells' analytic SEs within 25% of the refit-bootstrap SE", {
  skip_on_cran()
  df <- make_thin_panel(1500L, 3L)
  # min_pair_units = 10 ensures the guard binds in essentially every bootstrap
  # resample as well (cohort-4 resample counts are ~Poisson(3); at the default 5
  # about 18% of draws cross back above the threshold and deliberately re-enter
  # the audited overidentified failure mode inside those draws, which the boot SE
  # then correctly reports as extra noise). With the guard binding throughout,
  # analytic and refit-bootstrap SEs must agree.
  # (Literal edid() calls: the refit bootstrap re-evaluates the matched call, so a
  # wrapper's local variables would be out of scope there.)
  fit <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE, min_pair_units = 10L))
  rb  <- suppressWarnings(edid_refit_bootstrap(fit, data = df, B = 200L, seed = 20260611L))
  tab <- rb$att_gt
  ok  <- is.finite(tab$se_analytic) & tab$se_analytic > 0          # drops the degenerate (4,3) placebo
  ratio <- tab$se_boot[ok] / tab$se_analytic[ok]
  deg <- tab$group[ok] == 4
  expect_true(all(ratio[deg] >= 0.75 & ratio[deg] <= 1.25))        # degraded cells: within 25%
  expect_true(all(ratio[!deg] >= 0.75 & ratio[!deg] <= 1.25))      # healthy cells too
  # regression vs the audited failure: pre-guard the (4,4) boot/analytic ratio was
  # ~8x even at the DEFAULT threshold; post-guard it must stay far below that
  fit5 <- suppressWarnings(
    edid(df, yname = "y", idname = "id", tname = "time", gname = "gvar",
         pt_assumption = "all", weight_scheme = "efficient",
         aggregate = "none", cband = FALSE))
  rb5  <- suppressWarnings(edid_refit_bootstrap(fit5, data = df, B = 200L, seed = 20260611L))
  r44  <- rb5$att_gt[rb5$att_gt$group == 4 & rb5$att_gt$time == 4, ]
  expect_lt(r44$se_boot / r44$se_analytic, 2.5)
})

test_that("3-unit cohort MC (n = 1500): spillover gone, just-identified calibration restored", {
  skip_on_cran()
  # Design held FIXED across reps (outcomes redrawn), as in nesting Section 7.
  set.seed(20260611L)
  G0 <- sample(c(3L, 0L), 1500L, TRUE, c(0.4, 0.6)); G0[1:3] <- 4L
  z <- qnorm(0.975)
  reps <- 150L
  M <- matrix(NA_real_, reps, 8L,
              dimnames = list(NULL, c("a44", "s44", "p44", "sp44", "a33", "s33", "a34", "s34")))
  for (r in seq_len(reps)) {
    set.seed(54000L + r)
    d  <- gen_thin_outcomes(G0)
    fa <- suppressWarnings(fit_thin(d))
    fp <- suppressWarnings(fit_thin(d, pt = "post"))
    M[r, ] <- c(unlist(cell_of(fa, 4, 4)[, c("att", "se")]),
                unlist(cell_of(fp, 4, 4)[, c("att", "se")]),
                unlist(cell_of(fa, 3, 3)[, c("att", "se")]),
                unlist(cell_of(fa, 3, 4)[, c("att", "se")]))
  }
  cover <- function(a, s) mean(abs(M[, a] - 1) <= z * M[, s])

  # healthy cohort's cells: spillover gone, coverage back at/above 0.90
  # (pre-guard, at n = 2000, these covered 0.85 / 0.88 with analytic/MC 0.28 / 0.67)
  expect_gte(cover("a33", "s33"), 0.90)
  expect_gte(cover("a34", "s34"), 0.90)
  expect_gte(mean(M[, "s33"]) / sd(M[, "a33"]), 0.85)              # analytic SE honest again
  expect_gte(mean(M[, "s34"]) / sd(M[, "a34"]), 0.85)

  # thin cohort's cell: numerically the just-identified pt = "post" cell rep by rep,
  # hence exactly its (calibrated) coverage. NOTE: with 3 units the analytic z-CI of
  # ANY estimator covers ~0.83, not 0.95 (a 2-df t statistic; the audit's
  # just-identified benchmark showed the same 0.70-0.83 at n_thin = 3) -- the guard
  # restores that benchmark from the pre-guard 0.10, and edid_refit_bootstrap() is
  # recommended (and warned about) for inference this thin.
  expect_equal(M[, "a44"], M[, "p44"], tolerance = 1e-10)
  expect_equal(M[, "s44"], M[, "sp44"], tolerance = 1e-10)
  expect_equal(cover("a44", "s44"), cover("p44", "sp44"))
  expect_gte(cover("a44", "s44"), 0.75)                            # vs 0.10 pre-guard
  expect_gte(mean(M[, "s44"]) / sd(M[, "a44"]), 0.60)              # vs 0.05 pre-guard
})

# ===========================================================================
# guarded covariate fit: perturbation bootstrap reproduces the guarded cells
# ===========================================================================
test_that("perturbation bootstrap reproduces a guarded covariate fit (exactness check)", {
  skip_on_cran()
  set.seed(7L)
  n <- 300L; Tt <- 4L
  coh <- sample(c(2, Inf), n, replace = TRUE, prob = c(.45, .55))
  coh[1:3] <- 3                                                    # thin cohort 3 (3 units)
  df <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  x <- rnorm(n); df$x1 <- x[df$id]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + 0.3 * df$x1 * df$time +
    1 * (df$time >= df$g) + rnorm(n * Tt, 0, .5)
  fit <- suppressWarnings(
    edid(df, "y", "id", "time", "g", xformla = ~ x1, weight_scheme = "uniform",
         aggregate = "event_study", cband = FALSE, seed = 2L))
  expect_identical(fit$thin_cohorts$cohort, 3)
  # the rebuild applies the fit's min_pair_units, so the exactness guard passes
  pb <- suppressWarnings(
    edid_perturbation_bootstrap(fit, data = df, B = 59L, seed = 3L, agg = "event_study"))
  expect_s3_class(pb, "edid_perturbation_bootstrap")
  expect_identical(pb$n_failed, 0L)
})
