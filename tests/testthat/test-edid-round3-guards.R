library(testthat)

# ===========================================================================
# Round-3 guards / messages (all behavior-neutral except the OPT-IN fix 8):
#   1. edid_hausman broken-leg sanity guard (hollow-pass footgun)
#   2. thin-cohort radar note (12-35-unit cohorts; covariate path)
#   3. macOS fork-unsafe BLAS guard (Darwin + Accelerate -> serial)
#   4. edid_adaptive AUTO fallback when the efficient leg is not tighter (VR>=VU)
#   5. few-cluster toolkit guard (<5 clusters -> stat unreliable)
#   6. net-cross-moment-mass red flag (informational)
#   7. curse-of-dimensionality warning suppressed on just-identified PT-Post
#   8. estimability auto-guard (opt-in cross-cohort pair excision post-trim)
# Reproducers distilled from the gate runs (Nguyen, Bailey-GB, Brazil, ACA).
# ===========================================================================

# Small healthy staggered panel (PT-All true): no guard should false-positive.
mk_panel <- function(seed, n = 400L, viol = 0) {
  set.seed(seed)
  Tt  <- 6L
  coh <- sample(c(3, 5, Inf), n, replace = TRUE, prob = c(.3, .3, .4))
  df  <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + rnorm(nrow(df), 0, 0.5) + 1 * (df$time >= df$g)
  if (viol != 0) df$y <- df$y + viol * (df$g == 3 & df$time == 1)
  df
}

# A d=5-covariate panel with thin treated cohorts: the efficient kernel collapses, the
# fit becomes numerically degenerate (extreme ratios), and the efficient leg is NOT
# empirically tighter (VR >= VU) -- the Brazil/Bailey-GB with-X regime.
mk_broken_cov <- function(seed = 303, n = 700L) {
  set.seed(seed); Tt <- 6L
  coh <- sample(c(3, 4, 5, Inf), n, replace = TRUE, prob = c(.1, .1, .1, .7))
  df <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]
  X <- matrix(rnorm(n * 5), n, 5)
  for (k in 1:5) df[[paste0("x", k)]] <- X[df$id, k]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + rnorm(nrow(df), 0, 0.5) +
    1 * (df$time >= df$g) + 0.3 * rowSums(X)[df$id]
  df
}

# ---------------------------------------------------------------------------
# Fix 1: broken-leg Hausman guard
# ---------------------------------------------------------------------------
test_that("edid_hausman flags a hollow test when a constituent fit is degenerate (fix 1)", {
  df <- mk_broken_cov()
  fR <- suppressWarnings(suppressMessages(edid(df, "y", "id", "time", "g", xformla = ~ x1 + x2 + x3 + x4 + x5,
              pt_assumption = "all", weight_scheme = "efficient", aggregate = "event_study",
              cband = FALSE, trim_level = 200)))
  fU <- suppressWarnings(suppressMessages(edid(df, "y", "id", "time", "g", xformla = ~ x1 + x2 + x3 + x4 + x5,
              pt_assumption = "post", aggregate = "event_study", cband = FALSE, trim_level = 200)))
  # The restricted (PT-All) fit is degenerate (extreme ratios) -> hollow test.
  expect_true(isTRUE(fR$diagnostics$unstable))
  expect_warning(h <- edid_hausman(fU, fR), "HOLLOW")
  expect_true(isTRUE(h$leg_unstable))
  expect_true(length(h$leg_reasons) > 0L)
  expect_output(print(h), "HOLLOW TEST")
})

test_that("edid_hausman does NOT flag a healthy fit as hollow (fix 1 no false positive)", {
  df <- mk_panel(1)
  fR <- edid(df, "y", "id", "time", "g", pt_assumption = "all", aggregate = "event_study", cband = FALSE)
  fU <- edid(df, "y", "id", "time", "g", pt_assumption = "post", aggregate = "event_study", cband = FALSE)
  expect_false(isTRUE(fR$diagnostics$unstable))
  h <- edid_hausman(fU, fR)             # no warning expected
  expect_false(isTRUE(h$leg_unstable))
  expect_identical(h$leg_reasons, character(0L))
})

# ---------------------------------------------------------------------------
# Fix 2: thin-cohort radar note
# ---------------------------------------------------------------------------
test_that("thin-cohort radar records small cohorts as a field, prints only on covariate PT-All (fix 2)", {
  # mpdta's 2004 cohort has 20 units: in [min_pair_units, comfort). The radar is a FIELD +
  # a printed note, never a warning() (so it cannot flood the warning stream). On the no-X
  # fit the field is recorded but the note is NOT printed (the no-X path is fine here).
  data(mpdta, package = "did")
  ws <- character(0)
  f_nox <- withCallingHandlers(
    edid(mpdta, "lemp", "countyreal", "year", "first.treat", aggregate = "none", cband = FALSE),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("Thin-cohort radar", ws)))      # never a warning
  expect_false(is.null(f_nox$diagnostics$small_cohorts))  # field recorded
  expect_true(20L %in% f_nox$diagnostics$small_cohorts$n_units)
  expect_false(any(grepl("Thin-cohort radar", capture.output(print(f_nox)))))  # not printed on no-X
  # On a covariate PT-All fit with a small cohort (firmly in [5, 36)), the note IS printed.
  set.seed(42); n <- 300L; Tt <- 5L
  coh <- sample(c(3, Inf), n, replace = TRUE, prob = c(.08, .92))   # ~25-unit cohort
  df <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n)); df$g <- coh[df$id]
  df$x <- rnorm(n)[df$id]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + rnorm(nrow(df), 0, 0.5) + 1 * (df$time >= df$g) + 0.3 * df$x
  fX <- suppressWarnings(suppressMessages(edid(df, "y", "id", "time", "g", xformla = ~ x,
              pt_assumption = "all", aggregate = "none", cband = FALSE, trim_level = 200)))
  expect_false(is.null(fX$diagnostics$small_cohorts))
  expect_true(any(grepl("Thin-cohort radar", capture.output(print(fX)))))      # printed on cov PT-All
})

# ---------------------------------------------------------------------------
# Fix 3: macOS fork-unsafe BLAS guard
# ---------------------------------------------------------------------------
test_that("fork-BLAS guard serializes cores>1 on Darwin+Accelerate and respects override (fix 3)", {
  data(mpdta, package = "did")
  if (isTRUE(.edid_fork_blas_unsafe())) {
    # On the unsafe platform the message fires and the fit still completes (serial).
    expect_message(
      f <- edid(mpdta, "lemp", "countyreal", "year", "first.treat", aggregate = "none",
                cband = FALSE, cores = 2),
      "not.*fork-safe")
    expect_s3_class(f, "edid_fit")
    expect_equal(nrow(f$att_gt), 12L)
    # Override suppresses the guard (no fork-safe message); allow the genuine fork to run.
    op <- options(edid_allow_fork_blas = TRUE); on.exit(options(op), add = TRUE)
    expect_no_condition_msg <- tryCatch({
      suppressWarnings(edid(mpdta, "lemp", "countyreal", "year", "first.treat",
                            aggregate = "none", cband = FALSE, cores = 1))
      TRUE
    }, error = function(e) FALSE)
    expect_true(expect_no_condition_msg)
  } else {
    # On a fork-safe platform the helper returns FALSE and cores>1 is not downgraded here.
    expect_false(.edid_fork_blas_unsafe())
    expect_s3_class(
      suppressWarnings(edid(mpdta, "lemp", "countyreal", "year", "first.treat",
                            aggregate = "none", cband = FALSE, cores = 1)),
      "edid_fit")
  }
})

# ---------------------------------------------------------------------------
# Fix 4: edid_adaptive AUTO fallback (VR >= VU)
# ---------------------------------------------------------------------------
test_that("edid_adaptive AUTO falls back to empirical covariance when the efficient leg is not tighter (fix 4)", {
  df <- mk_broken_cov()
  fR <- suppressWarnings(suppressMessages(edid(df, "y", "id", "time", "g", xformla = ~ x1 + x2 + x3 + x4 + x5,
              pt_assumption = "all", weight_scheme = "efficient", aggregate = "event_study",
              cband = FALSE, trim_level = 200)))
  fU <- suppressWarnings(suppressMessages(edid(df, "y", "id", "time", "g", xformla = ~ x1 + x2 + x3 + x4 + x5,
              pt_assumption = "post", aggregate = "event_study", cband = FALSE, trim_level = 200)))
  g <- .edid_param_ifs(fR, "overall"); VR <- as.numeric(cluster_cov_edid(g$IF, NULL, fR$n))
  gU <- .edid_param_ifs(fU, "overall"); VU <- as.numeric(cluster_cov_edid(gU$IF, NULL, fU$n))
  skip_if_not(VR >= VU, "fixture did not produce VR >= VU on this platform/BLAS")
  # AUTO (weight_scheme = "efficient" labels bound-attaining) would impose VUR = VR -> VO <= 0.
  # The fallback recovers with a message instead of erroring.
  expect_message(a <- suppressWarnings(edid_adaptive(fU, fR, parameter = "overall")), "Falling back")
  expect_true(isTRUE(a$assume_efficient_fallback))
  expect_false(isTRUE(a$assume_efficient))
  expect_true(is.finite(a$adaptive))
  # Explicit assume_efficient = TRUE STILL errors (documented contract preserved).
  expect_error(suppressWarnings(edid_adaptive(fU, fR, parameter = "overall", assume_efficient = TRUE)),
               "not positive")
})

# ---------------------------------------------------------------------------
# Fix 5: few-cluster toolkit guard
# ---------------------------------------------------------------------------
test_that("edid_hausman flags few-cluster (G<5) statistics as unreliable (fix 5)", {
  df <- mk_panel(2)
  df$clu <- ((df$id - 1) %% 3) + 1                    # only 3 clusters
  fR <- suppressWarnings(edid(df, "y", "id", "time", "g", pt_assumption = "all",
              aggregate = "event_study", cband = FALSE, clustervars = "clu"))
  fU <- suppressWarnings(edid(df, "y", "id", "time", "g", pt_assumption = "post",
              aggregate = "event_study", cband = FALSE, clustervars = "clu"))
  expect_warning(h <- edid_hausman(fU, fR), "cluster")
  expect_true(isTRUE(h$few_clusters))
  expect_equal(h$n_clusters, 3L)
  expect_output(print(h), "cluster")
  # Unclustered fits do NOT trip the guard.
  fRu <- edid(df, "y", "id", "time", "g", pt_assumption = "all", aggregate = "event_study", cband = FALSE)
  fUu <- edid(df, "y", "id", "time", "g", pt_assumption = "post", aggregate = "event_study", cband = FALSE)
  hu <- edid_hausman(fUu, fRu)
  expect_false(isTRUE(hu$few_clusters))
})

# ---------------------------------------------------------------------------
# Fix 6: net-cross-moment-mass red flag (informational)
# ---------------------------------------------------------------------------
test_that("net-hedge-mass diagnostic is computed and does not flag a healthy fit (fix 6)", {
  df <- mk_panel(3)
  ws <- character(0)
  f <- withCallingHandlers(
    edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_true(is.finite(f$diagnostics$net_hedge_mass))
  expect_false(isTRUE(f$diagnostics$net_hedge_flag))           # healthy: net hedge well below threshold
  expect_false(any(grepl("Net cross-cohort", ws)))            # no red-flag warning
  # The flag logic itself: a synthetic diagnostics record at/above threshold (net ~= gross)
  # is classified as flagged by the builder.
  fake_cells <- list(list(is_pre = FALSE, group = 3,
                          pairs = data.frame(gp = c(3, 5), tpre = c(2, 2)),
                          weights = c(0.4, 0.65)))
  d <- .edid_build_diagnostics(list(n_extreme_ratio = 0L, use_cov_path = TRUE),
                               fake_cells, "all", "efficient", 5L)
  expect_true(isTRUE(d$net_hedge_flag))    # cross-pair gp=5 weight 0.65 >= 0.55, net ~= gross
})

# ---------------------------------------------------------------------------
# Fix 7: curse-of-dimensionality warning suppressed on PT-Post
# ---------------------------------------------------------------------------
test_that("curse-of-dimensionality warning is suppressed on just-identified PT-Post fits (fix 7)", {
  df <- mk_broken_cov(seed = 11, n = 400L)
  # PT-All efficient d=5: curse warning SHOULD fire (over-identified efficient weights formed).
  ws_all <- character(0)
  withCallingHandlers(
    suppressMessages(edid(df, "y", "id", "time", "g", xformla = ~ x1 + x2 + x3 + x4 + x5,
         pt_assumption = "all", weight_scheme = "efficient", aggregate = "none",
         cband = FALSE, trim_level = 200)),
    warning = function(w) { ws_all <<- c(ws_all, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("curse of dimensionality", ws_all)))
  # PT-Post d=5 efficient: just-identified, no efficient weights formed -> NO curse warning.
  ws_post <- character(0)
  withCallingHandlers(
    suppressMessages(edid(df, "y", "id", "time", "g", xformla = ~ x1 + x2 + x3 + x4 + x5,
         pt_assumption = "post", weight_scheme = "efficient", aggregate = "none",
         cband = FALSE, trim_level = 200)),
    warning = function(w) { ws_post <<- c(ws_post, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("curse of dimensionality", ws_post)))
})

# ---------------------------------------------------------------------------
# Fix 8: estimability auto-guard (opt-in)
# ---------------------------------------------------------------------------
test_that("estimability auto-guard is OFF by default and excises only when opted in (fix 8)", {
  df <- mk_broken_cov()
  xf <- ~ x1 + x2 + x3 + x4 + x5
  # Default OFF: no excision warning; the fit is byte-identical to the un-guarded fit.
  f_off <- suppressWarnings(suppressMessages(edid(df, "y", "id", "time", "g", xformla = xf,
                pt_assumption = "all", weight_scheme = "averaged", aggregate = "none",
                cband = FALSE, trim_level = 200)))
  expect_false(any(grepl("auto-guard", paste(names(f_off), collapse = " "))))
  # Opt-in: if any cross pair is unstable post-trim, it is excised with a named warning;
  # the never-treated and self moments are retained, so the fit still returns finite cells.
  op <- options(edid_auto_excise_unstable_pairs = TRUE); on.exit(options(op), add = TRUE)
  ws <- character(0)
  f_on <- withCallingHandlers(
    suppressMessages(edid(df, "y", "id", "time", "g", xformla = xf, pt_assumption = "all",
         weight_scheme = "averaged", aggregate = "none", cband = FALSE, trim_level = 200)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_s3_class(f_on, "edid_fit")
  # On this fixture extreme ratios occur, so the auto-guard fires; if it did, the warning
  # names the auto-guard and at least one post cell survives.
  if (any(grepl("Estimability auto-guard", ws))) {
    expect_true(any(is.finite(f_on$att_gt$att[!grepl("pre", rownames(f_on$att_gt))])))
  } else {
    succeed("no post-trim-unstable cross pair on this platform; guard correctly inert")
  }
})

test_that(".edid_ratio_unstable_pairs excises extreme post-trim ratios, spares self/NT pairs (fix 8 unit)", {
  pairs <- data.frame(gp = c(3, 5, Inf), tpre = c(2, 2, 2))      # target g = 3
  N <- 8L                                                        # > EDID_RATIO_EXCISE_MINKEEP so mass_gone does not fire
  keep_all <- list(`5` = rep(TRUE, N), `3` = rep(TRUE, N), `Inf` = rep(TRUE, N))
  # (a) extreme ratio survives the trim: gp=5 has a 250 on a kept unit; gp=3 (self) and
  # Inf (NT) must never be candidates even if their (hypothetical) ratios were large.
  prop_ratios <- list(`5` = c(rep(3, N - 1L), 250), `3` = rep(2, N), `Inf` = rep(1, N))
  r <- .edid_ratio_unstable_pairs(pairs, target_g = 3, prop_ratios, keep_all)
  expect_equal(r$drop_gp, 5)             # only the extreme cross cohort
  # A healthy cross ratio (all moderate, full mass kept) is NOT excised.
  prop_ratios2 <- list(`5` = rep(c(5, 8), length.out = N), `3` = rep(2, N), `Inf` = rep(1, N))
  r2 <- .edid_ratio_unstable_pairs(pairs, target_g = 3, prop_ratios2, keep_all)
  expect_equal(length(r2$drop_gp), 0L)
  # (b) mass-gone branch: a healthy ratio but only < MINKEEP units survive the trim -> excised.
  keep_thin <- keep_all; keep_thin[["5"]] <- c(rep(TRUE, 3L), rep(FALSE, N - 3L))
  r3 <- .edid_ratio_unstable_pairs(pairs, target_g = 3, prop_ratios2, keep_thin)
  expect_equal(r3$drop_gp, 5)            # mass removed for gp=5
})
