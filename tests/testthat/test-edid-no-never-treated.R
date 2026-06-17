library(testthat)

# ============================================================
# No never-treated group: edid() coerces the last-treated cohort into the
# comparison group (mirrors att_gt's control_group = "nevertreated"
# pre-processing). When every cohort is eventually treated (no g == Inf and
# no g == 0), edid() drops all periods at/after the last cohort's effective
# onset g_max - anticipation and recasts g_max as never-treated.
#
# CORRECTNESS ORACLE: the internal transform must produce results
# BYTE-IDENTICAL to running edid() on the SAME panel transformed by hand
# (drop t >= g_max - anticipation; relabel g == g_max to Inf). When a
# never-treated group is already present (Inf, or the att_gt 0 convention),
# the transform is a NO-OP.
# ============================================================

# Balanced panel from an EXPLICIT cohort vector. (Do NOT use sample(c(k), ...):
# a length-1 first arg triggers R's sample(1:k) convenience gotcha.)
mk_nnt <- function(gvec, Tn = 6L, seed = 7L) {
  set.seed(seed)
  n   <- length(gvec)
  x1  <- runif(n, -2, 2)
  alph <- rnorm(n, 0.4 * x1, 1)
  wv  <- exp(0.5 * rnorm(n)); wv <- wv / mean(wv)            # dispersed mean-1 weights
  cl  <- rep_len(1:12, n)
  do.call(rbind, lapply(seq_len(Tn), function(tt) {
    tau <- ifelse(is.finite(gvec) & tt >= gvec, 1, 0)
    data.frame(id = seq_len(n), t = tt, g = gvec, x1 = x1, w = wv, cl = cl,
               y = alph + 0.3 * tt + (tt - 1) * 0.4 * x1 + tau + rnorm(n))
  }))
}

# Hand transform = exactly what edid() should do internally.
hand_transform <- function(d, anticipation = 0L) {
  g_max  <- max(d$g[is.finite(d$g)])
  cutoff <- g_max - anticipation
  d2 <- d[d$t < cutoff, , drop = FALSE]
  d2$g[d2$g == g_max] <- Inf
  d2
}

fit <- function(d, ...) suppressWarnings(
  edid(d, "y", "id", "t", "g", aggregate = "none", bstrap = FALSE, seed = 1L, ...))

# A typical no-never-treated panel: cohorts {2,3,4}, periods 1..6, no Inf/0.
gv  <- rep_len(c(2, 3, 4), 360L)
D   <- mk_nnt(gv)

# ---- (1) the transform fires (warning) and produces a valid fit -----------------
test_that("no never-treated: edid() warns and coerces the last cohort as comparison", {
  expect_warning(
    edid(D, "y", "id", "t", "g", aggregate = "none", bstrap = FALSE, seed = 1L),
    "No never-treated group is available")
  f <- fit(D)
  expect_false(is.null(f$att_gt))
  expect_true(all(is.finite(f$att_gt$att)))
  expect_true(all(is.finite(f$att_gt$se)))
})

# ---- (2) CORRECTNESS ORACLE: internal transform == hand transform ---------------
test_that("internal no-never-treated transform == edid() on the hand-transformed panel", {
  configs <- list(
    nocov          = list(xformla = ~1),
    nocov_post     = list(xformla = ~1, pt_assumption = "post"),
    cov            = list(xformla = ~x1),
    cov_weighted   = list(xformla = ~x1, weightsname = "w"),
    cov_averaged   = list(xformla = ~x1, weight_scheme = "averaged"),
    cov_misspec    = list(xformla = ~x1, misspec_robust = TRUE)
  )
  Dman <- hand_transform(D)
  for (nm in names(configs)) {
    fa <- do.call(fit, c(list(D),    configs[[nm]]))   # internal transform
    fm <- do.call(fit, c(list(Dman), configs[[nm]]))   # hand-transformed input (already has Inf)
    expect_equal(fa$att_gt$att, fm$att_gt$att, tolerance = 1e-10, info = nm)
    expect_equal(fa$att_gt$se,  fm$att_gt$se,  tolerance = 1e-10, info = nm)
    expect_equal(as.numeric(fa$eif), as.numeric(fm$eif), tolerance = 1e-9, info = nm)
  }
})

# ---- (3) anticipation: cutoff is g_max - anticipation ---------------------------
test_that("anticipation shifts the cutoff to g_max - anticipation (internal == hand)", {
  gv2 <- rep_len(c(3, 4, 5), 360L)
  D2  <- mk_nnt(gv2, Tn = 8L)
  # a in {0,1}: both leave every retained cohort estimable; the two cutoffs differ
  # (keep t<5 vs t<4), so matching the hand transform for each proves the shift.
  for (a in c(0L, 1L)) {
    fa <- fit(D2, xformla = ~x1, anticipation = a)
    fm <- fit(hand_transform(D2, anticipation = a), xformla = ~x1, anticipation = a)
    expect_equal(fa$att_gt$att, fm$att_gt$att, tolerance = 1e-10, info = paste0("anticipation=", a))
    expect_equal(fa$att_gt$se,  fm$att_gt$se,  tolerance = 1e-10, info = paste0("anticipation=", a))
  }
})

# ---- (4) aggregation + clustering propagate through the transform ----------------
test_that("no-never-treated transform is consistent across aggregations and clustering", {
  Dman <- hand_transform(D)
  for (ag in c("group", "event_study", "calendar", "overall")) {
    fa <- suppressWarnings(edid(D,    "y", "id", "t", "g", xformla = ~x1, aggregate = ag,
                                bstrap = FALSE, seed = 1L))
    fm <- suppressWarnings(edid(Dman, "y", "id", "t", "g", xformla = ~x1, aggregate = ag,
                                bstrap = FALSE, seed = 1L))
    la <- fa[[ag]]; lm <- fm[[ag]]
    expect_equal(la$att.egt,     lm$att.egt,     tolerance = 1e-9, info = ag)
    expect_equal(la$overall.att, lm$overall.att, tolerance = 1e-9, info = ag)
  }
  fa <- fit(D,             xformla = ~x1, clustervars = "cl")
  fm <- fit(hand_transform(D), xformla = ~x1, clustervars = "cl")
  expect_equal(fa$att_gt$se, fm$att_gt$se, tolerance = 1e-9)
})

# ---- (5) NO-OP when a never-treated group already exists -------------------------
test_that("never-treated already present (Inf): transform is a no-op (no warning)", {
  gv3 <- rep_len(c(2, 3, Inf), 360L)
  D3  <- mk_nnt(gv3)
  expect_no_warning(
    suppressMessages(edid(D3, "y", "id", "t", "g", xformla = ~x1, aggregate = "none",
                          bstrap = FALSE, seed = 1L)),
    message = "No never-treated group is available")
})

test_that("att_gt 0-convention never-treated is recoded to Inf -> transform is a no-op", {
  gv4 <- rep_len(c(2, 3, 0), 360L)               # 0 = att_gt never-treated convention
  D4  <- mk_nnt(gv4)
  fired <- FALSE
  f <- withCallingHandlers(
    suppressMessages(edid(D4, "y", "id", "t", "g", xformla = ~x1, aggregate = "none",
                          bstrap = FALSE, seed = 1L)),
    warning = function(w) {
      if (grepl("No never-treated group is available", conditionMessage(w))) fired <<- TRUE
      invokeRestart("muffleWarning")
    })
  expect_false(fired)
  expect_true(all(is.finite(f$att_gt$se)))
})

# ---- (6) edge-case guards -------------------------------------------------------
test_that("single treated cohort + no never-treated errors clearly (nothing to compare)", {
  D5 <- mk_nnt(rep(3, 300L))                     # ONE cohort, no Inf/0
  expect_error(
    edid(D5, "y", "id", "t", "g", xformla = ~x1, bstrap = FALSE),
    "only one treated cohort")
})

test_that("anticipation that leaves < 2 retained periods errors clearly", {
  # cohorts {2,3}, periods 1..4 -> g_max=3; anticipation=2 -> cutoff=1 -> 0 periods kept
  D6 <- mk_nnt(rep_len(c(2, 3), 200L), Tn = 4L)
  expect_error(
    edid(D6, "y", "id", "t", "g", xformla = ~x1, anticipation = 2L, bstrap = FALSE),
    "only .* period")
})

# ---- (7) the coerced comparison cohort actually anchors estimable cells ----------
test_that("the recast last cohort yields the same estimable (g,t) cells as the hand panel", {
  fa <- fit(D, xformla = ~x1)
  fm <- fit(hand_transform(D), xformla = ~x1)
  expect_equal(fa$att_gt[, c("group", "time")], fm$att_gt[, c("group", "time")])
  expect_true(all(fa$att_gt$group < max(gv)))   # the recast cohort (g_max) is no longer a target group
})

# ---- (8) BOTH bootstraps flow through (refit re-runs edid; perturbation reuses the
#          shared coercion helper so its inline panel rebuild matches the fit) -------
test_that("edid_perturbation_bootstrap() works on a no-never-treated fit", {
  f <- fit(D, xformla = ~x1)
  # (data passed explicitly: the test's fit() wrapper hides `data` from the captured
  # call, so data=NULL recovery is unavailable here -- that is a call-capture detail
  # unrelated to the no-never-treated coercion, which is what this test exercises.)
  pb0 <- suppressWarnings(edid_perturbation_bootstrap(f, data = D, B = 30L, seed = 1L))
  expect_false(is.null(pb0$att_gt))
  expect_true(all(is.finite(pb0$att_gt$se[is.finite(pb0$att_gt$se)])))
  # AUTO == MANUAL: same seeded perturbation SEs as on the hand-transformed panel
  fm  <- fit(hand_transform(D), xformla = ~x1)
  pba <- suppressWarnings(edid_perturbation_bootstrap(f,  data = D,                 B = 60L, seed = 7L))
  pbm <- suppressWarnings(edid_perturbation_bootstrap(fm, data = hand_transform(D), B = 60L, seed = 7L))
  expect_equal(pba$att_gt$se, pbm$att_gt$se, tolerance = 1e-8)
})

# ---- (9) NA in gname defers to validation (no spurious coercion / max() warning) --
test_that("NA in gname does not trigger the coercion path (deferred to validation)", {
  Dna <- D; Dna$g[Dna$id == 1L] <- NA
  w <- character(0)
  e <- withCallingHandlers(
    tryCatch(edid(Dna, "y", "id", "t", "g", bstrap = FALSE), error = conditionMessage),
    warning = function(wn) { w <<- c(w, conditionMessage(wn)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("coerced as the never-treated", w)))   # no spurious coercion notice
  expect_match(e, "NA|missing")                                 # validation reports the real cause
  Dall <- D; Dall$g <- NA_real_                                 # all-NA: no base max() warning either
  w2 <- character(0)
  e2 <- withCallingHandlers(
    tryCatch(edid(Dall, "y", "id", "t", "g", bstrap = FALSE), error = conditionMessage),
    warning = function(wn) { w2 <<- c(w2, conditionMessage(wn)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("no non-missing", w2)))
})

# ---- (10) anticipation-aware is_pre: the effective-onset cell is POST, not pre -----
test_that("is_pre respects anticipation (the g - anticipation cell is labeled post)", {
  D2 <- mk_nnt(rep_len(c(3, 4, 5), 360L), Tn = 8L)
  f  <- suppressWarnings(edid(D2, "y", "id", "t", "g", pt_assumption = "post",
                              aggregate = "none", bstrap = FALSE, anticipation = 1L))
  skip_if(is.null(f$att_gt$is_pre))
  row <- f$att_gt[f$att_gt$group == 4 & f$att_gt$time == 3, ]   # effective onset 4 - 1 = 3 -> post
  skip_if(nrow(row) == 0L)
  expect_false(isTRUE(row$is_pre))
})
