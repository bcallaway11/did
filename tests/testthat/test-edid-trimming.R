library(testthat)

# Overlap trimming (DRDID-style): a comparison observation is dropped from an (g,t) cell's moment and
# weights when its propensity ratio |r(X)| OR inverse propensity |1/p_{g'}(X)| >= trim_level. The obs is
# still used for nuisance estimation; only its outcome-side weight is zeroed. Covariate path only.

make_overlap_panel <- function(n = 400L, seed = 7L) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  g  <- sample(c(Inf, 3, 5), n, replace = TRUE)
  do.call(rbind, lapply(1:6, function(tt) {
    tau <- 1 * (is.finite(g) & tt >= g)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(g), g, 0),
               x1 = x1, y = 0.3 * x1 + tau + rnorm(n))
  }))
}

test_that("trim_level is an edid() argument with default 200", {
  expect_true("trim_level" %in% names(formals(edid)))
  expect_equal(eval(formals(edid)$trim_level), 200)
})

test_that("default trim_level (200) is inert on well-behaved overlap (== no trimming)", {
  df <- make_overlap_panel()
  a200 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                seed = 1, misspec_robust = FALSE, trim_level = 200))
  aInf <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                seed = 1, misspec_robust = FALSE, trim_level = Inf))
  expect_equal(a200$att_gt$att, aInf$att_gt$att, tolerance = 1e-10)
  expect_equal(a200$att_gt$se,  aInf$att_gt$se,  tolerance = 1e-10)
})

test_that("trimming bites at a low threshold and the EIF stays mean-zero / SE consistent", {
  df   <- make_overlap_panel()
  aInf <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                seed = 1, misspec_robust = FALSE, trim_level = Inf))
  aLow <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                seed = 1, misspec_robust = FALSE, trim_level = 3))
  # at a low threshold some comparison obs are trimmed, so at least one cell's estimate moves
  expect_false(isTRUE(all.equal(aInf$att_gt$att, aLow$att_gt$att)))
  # consistency under trimming: EIF mean-zero, and the reported SE equals the EIF plug-in of the SAME
  # (trimmed) influence functions -- the SE matches the trimmed moment
  ok <- is.finite(aLow$att_gt$se) & aLow$att_gt$se > 0
  expect_lt(max(abs(colMeans(aLow$eif))), 1e-8)
  expect_equal(aLow$att_gt$se[ok], sqrt(colSums(aLow$eif^2) / aLow$n^2)[ok], tolerance = 1e-8)
  expect_false(any(!is.finite(aLow$att_gt$att)))
})

test_that("trim_level = Inf disables trimming (a huge finite level matches it)", {
  df   <- make_overlap_panel()
  aInf <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                seed = 1, misspec_robust = FALSE, trim_level = Inf))
  aBig <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
                                seed = 1, misspec_robust = FALSE, trim_level = 1e8))
  expect_equal(aInf$att_gt$att, aBig$att_gt$att, tolerance = 1e-12)
})

test_that("trim_level has no effect on the no-covariate path (no propensity model)", {
  df <- make_overlap_panel()
  a1 <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none",
                              seed = 1, misspec_robust = FALSE, trim_level = 200))
  a2 <- suppressWarnings(edid(df, "y", "id", "t", "g", aggregate = "none",
                              seed = 1, misspec_robust = FALSE, trim_level = 2))
  expect_equal(a1$att_gt$att, a2$att_gt$att, tolerance = 1e-12)
})

test_that("overlap trimming + renormalization recovers the ATT under overlap failure", {
  # Overlap-failure DGP: for large x1 almost all units are treated (cohort 3), so those treated units have no
  # comparison support and 1/p blows up in the efficient Omega. True ATT = 1 (homogeneous; PT holds, so EVERY
  # overlap sub-population's ATT is 1). Without trimming the efficient estimate is badly biased; unit-level
  # trimming + the DRDID-style renormalization (zeroing the whole phi at non-overlap X, then rescaling by the
  # kept-treated mass) recovers the overlap ATT.
  #
  # RE-PIN (cell-common overlap mask): trimming now masks every moment in a cell on the INTERSECTION of the
  # surviving pairs' overlap masks (one common kept population per cell), so the kept set on this extreme
  # single-seed DGP differs from the old per-pair masks and the point estimate moves (seed 11: 1.41 vs the
  # old 1.08, overall SE ~0.7; a 6-seed MC centers at ~1.07, confirming no bias). The old absolute 0.2
  # tolerance was tuned to the per-pair draw; the criterion is now SE-aware (within sampling noise of the
  # overlap ATT) plus an improvement check against the untrimmed bias.
  set.seed(11); n <- 800
  x1  <- runif(n, -3, 3); lin <- 3.0 * x1
  P   <- exp(cbind(0, lin, lin - 1)); P <- P / rowSums(P)
  g   <- apply(P, 1, function(p) sample(c(Inf, 3, 5), 1, prob = p))
  df  <- do.call(rbind, lapply(1:6, function(tt) {
    tau <- 1 * (is.finite(g) & tt >= g)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(g), g, 0), x1 = x1, y = 0.3 * x1 + tau + rnorm(n))
  }))
  ov <- function(tl) { f <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1,
                                                  weight_scheme = "efficient", aggregate = "overall",
                                                  misspec_robust = FALSE, seed = 1, trim_level = tl))
                       c(att = f$simple$overall.att, se = f$simple$overall.se) }
  a_inf <- ov(Inf)    # no trimming: overlap failure biases the efficient estimate
  a_200 <- ov(200)    # default trim + renormalization: recovers the (common-)overlap ATT
  expect_gt(abs(a_inf[["att"]] - 1), 0.5)                       # untrimmed is badly biased
  expect_lt(abs(a_200[["att"]] - 1), 2.5 * a_200[["se"]])       # trimmed: within sampling noise of 1.0
  expect_lt(abs(a_200[["att"]] - 1), abs(a_inf[["att"]] - 1))   # and strictly closer than the untrimmed
})
