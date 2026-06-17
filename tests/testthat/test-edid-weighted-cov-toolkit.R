library(testthat)

# ============================================================
# Phase 3c toolkit gate: the post-estimation toolkit
# (edid_weights / edid_sargan / edid_hausman / edid_frontier /
# edid_adaptive) must consume weighted-covariate fits correctly.
#
# Invariant: a CONSTANT weight column normalizes to unit_weights == 1,
# so every weighted-covariate fit is byte-identical to the unweighted
# covariate fit (validated in test-edid-weighted-cov-e2e). Therefore
# every toolkit function applied to the constant-weight fits must return
# output byte-identical to the unweighted-fit output. Plus: under a
# DISPERSED weight, every tool runs and returns finite output.
# ============================================================

make_tk_panel <- function(n = 360, seed = 4) {
  set.seed(seed)
  Tn  <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 0.5 * x1u
  P   <- exp(cbind(0, eta, 0.55 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.4 * x1u, 1)
  wu   <- exp(0.6 * rnorm(n)); wu <- wu / mean(wu)
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.4 * x1u)
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

# numeric leaves of a (possibly nested) result, for a structure-agnostic comparison
numleaves <- function(x) {
  out <- c()
  rec <- function(z) {
    if (is.numeric(z)) out <<- c(out, as.numeric(z))
    else if (is.list(z)) for (el in z) rec(el)
  }
  rec(x); out[is.finite(out)]
}

mk <- function(df, wt, pt) suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1,
  weightsname = wt, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE,
  seed = 1L, misspec_robust = FALSE, pt_assumption = pt))

df   <- make_tk_panel()
df_c <- df; df_c$w <- 4                                   # constant weight column -> unit_weights == 1

# unrestricted = PT-Post, restricted = PT-All (the Section-5 conservative-vs-efficient pairing)
Uu_post <- mk(df,   NULL, "post"); Uu_all <- mk(df,   NULL, "all")
Cc_post <- mk(df_c, "w",  "post"); Cc_all <- mk(df_c, "w",  "all")

tk_pairs <- list(
  weights  = list(u = quote(edid_weights(Uu_all)),                       c = quote(edid_weights(Cc_all))),
  sargan   = list(u = quote(edid_sargan(Uu_all, data = df)),             c = quote(edid_sargan(Cc_all, data = df_c))),
  # The toolkit refits the legs in plug-in mode; pass `data` explicitly (the fits were built inside mk(),
  # whose call captures the symbol `df`, so automatic recovery would pick the wrong panel for df_c).
  hausman  = list(u = quote(edid_hausman(Uu_post, Uu_all, data = df)),   c = quote(edid_hausman(Cc_post, Cc_all, data = df_c))),
  frontier = list(u = quote(edid_frontier(Uu_post, Uu_all, data = df)),  c = quote(edid_frontier(Cc_post, Cc_all, data = df_c))),
  adaptive = list(u = quote(edid_adaptive(Uu_post, Uu_all, data = df)),  c = quote(edid_adaptive(Cc_post, Cc_all, data = df_c)))
)

test_that("toolkit on a constant-weight covariate fit == toolkit on the unweighted covariate fit", {
  for (nm in names(tk_pairs)) {
    ru <- tryCatch(suppressWarnings(eval(tk_pairs[[nm]]$u)), error = function(e) e)
    rc <- tryCatch(suppressWarnings(eval(tk_pairs[[nm]]$c)), error = function(e) e)
    expect_false(inherits(ru, "error"), info = paste(nm, "unweighted errored"))
    expect_false(inherits(rc, "error"), info = paste(nm, "const-weight errored"))
    nu <- numleaves(ru); nc <- numleaves(rc)
    expect_equal(length(nu), length(nc), info = paste(nm, "shape"))
    expect_equal(nc, nu, tolerance = 1e-8, info = paste(nm, "numeric leaves"))
  }
})

test_that("toolkit runs + returns finite output on a DISPERSED-weight covariate fit", {
  Wp <- mk(df, "w", "post"); Wa <- mk(df, "w", "all")
  expect_true(length(numleaves(suppressWarnings(edid_weights(Wa))))  > 0L)
  expect_true(length(numleaves(suppressWarnings(edid_sargan(Wa))))   > 0L)
  expect_true(length(numleaves(suppressWarnings(edid_hausman(Wp, Wa))))  > 0L)
  expect_true(length(numleaves(suppressWarnings(edid_frontier(Wp, Wa)))) > 0L)
  expect_true(length(numleaves(suppressWarnings(edid_adaptive(Wp, Wa)))) > 0L)
})
