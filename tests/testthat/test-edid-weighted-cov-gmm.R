library(testthat)

# ============================================================
# Phase 3c (full weight propagation): the non-default weight schemes
# (gmm: invert the unconditional sample covariance; averaged: invert the
# pooled Omega-bar) must also propagate observation weights -- the gmm
# weight inverts the WEIGHTED covariance, the gmm/averaged weight-estimation
# corrections carry the per-unit weight.
#
#   (A) constant weight column => weight_scheme in {gmm, averaged} fit
#       (att + misspec_robust SE) == the unweighted fit.
#   (B) dispersed weights run + finite SE under misspec_robust.
# ============================================================

make_gmm_panel <- function(n = 320, seed = 9) {
  set.seed(n + seed); Tn <- 4L
  x1u <- runif(n, -2, 2); eta <- 0.5 * x1u
  P   <- exp(cbind(0, eta, 0.55 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.4 * x1u, 1); wu <- exp(0.6 * rnorm(n)); wu <- wu / mean(wu)
  rows <- lapply(1:Tn, function(tt) {
    ht <- (tt - 1) * (0.4 * x1u); tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

fit_ws <- function(d, wt, ws) suppressWarnings(edid(d, "y","id","t","g", xformla = ~ x1,
  weightsname = wt, weight_scheme = ws, aggregate = "none", bstrap = FALSE, seed = 1L,
  misspec_robust = TRUE))

test_that("constant weight column => gmm / averaged weighted fit == unweighted fit (att + misspec SE)", {
  df <- make_gmm_panel(); df_c <- df; df_c$w <- 6
  for (ws in c("gmm", "averaged")) {
    fu <- fit_ws(df,   NULL, ws)
    fc <- fit_ws(df_c, "w",  ws)
    ok <- is.finite(fu$att_gt$att) & is.finite(fc$att_gt$att)
    expect_equal(fc$att_gt$att[ok], fu$att_gt$att[ok], tolerance = 1e-7, info = paste(ws, "att"))
    okse <- ok & is.finite(fu$att_gt$se) & is.finite(fc$att_gt$se)
    expect_equal(fc$att_gt$se[okse], fu$att_gt$se[okse], tolerance = 1e-7, info = paste(ws, "se"))
  }
})

test_that("dispersed weights: gmm / averaged run with finite misspec SE", {
  df <- make_gmm_panel(seed = 13)
  for (ws in c("gmm", "averaged")) {
    f <- fit_ws(df, "w", ws)
    expect_true(any(is.finite(f$att_gt$se) & f$att_gt$se > 0), info = ws)
  }
})
