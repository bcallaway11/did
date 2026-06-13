library(testthat)

# Regression test for the PT-Post covariate generated-outcome base-period bug.
# DGP: cohorts {Inf, 3}, periods 1-4. A cohort-3 level shift `delta` starting at period 2 (= g-1) breaks
# PT-All (base period 1) but NOT PT-Post (base period g-1 = 2). True ATT(3,t) = 1 for t >= 3.
# Before the fix the covariate path used base period 1 even under pt="post" (the gp=Inf cross-pair formula
# collapses to the PT-All moment for g >= 3), biasing ATT(3,t) by ~delta. After the fix pt="post" is unbiased.

make_ptpost_panel <- function(seed = 1, n = 4000, delta = 0.6) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  g  <- ifelse(runif(n) < plogis(0.8 * x1), 3, Inf)
  alpha <- rnorm(n, 0.3 * x1, 1)
  do.call(rbind, lapply(1:4, function(t) {
    trend <- 0.3 * t + 0.5 * x1 * (t - 1)        # parallel trend (both cohorts)
    viol  <- delta * (g == 3) * (t >= 2)         # cohort-3 shift from period 2: breaks PT-All, not PT-Post
    tau   <- 1.0 * (g == 3) * (t >= 3)           # true ATT = 1 from period 3
    data.frame(id = 1:n, tt = t, g = ifelse(is.finite(g), g, 0), x1 = x1,
               y = alpha + trend + viol + tau + rnorm(n))
  }))
}

test_that("PT-Post covariate path uses base period g-1 (unbiased when PT-Post holds, PT-All fails)", {
  df <- make_ptpost_panel(seed = 1, n = 4000, delta = 0.6)
  fp <- edid(df, "y", "id", "tt", "g", xformla = ~ x1, pt_assumption = "post", aggregate = "none")
  a  <- fp$att_gt[fp$att_gt$group == 3 & fp$att_gt$time >= 3, ]
  # PT-Post holds => unbiased; the bug would give ~1 + delta = 1.6
  expect_lt(abs(a$att[a$time == 3] - 1), 0.2)
  expect_lt(abs(a$att[a$time == 4] - 1), 0.2)

  # sanity: PT-All is genuinely violated here, so pt="all" should be materially biased away from 1
  fa <- edid(df, "y", "id", "tt", "g", xformla = ~ x1, pt_assumption = "all", aggregate = "none")
  aa <- fa$att_gt[fa$att_gt$group == 3 & fa$att_gt$time == 3, ]
  expect_gt(aa$att - 1, 0.15)
})
