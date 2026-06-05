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
