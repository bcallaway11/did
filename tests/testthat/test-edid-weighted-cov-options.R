library(testthat)

# ============================================================
# FULL option-matrix weight-compatibility gate (full-adaptation rule):
# observation weights must be compatible with EVERY option. A CONSTANT
# weight column normalizes to unit_weights == 1, so for EVERY option
# combination the weighted-covariate fit MUST collapse to the unweighted
# fit (att + SE + aggregate). This sweeps the estimation surface and
# asserts weighted-const == unweighted for each; a single failing
# combination flags an option through which weights do not propagate.
#
# (Correctness under DISPERSED weights -- FD oracles + MC coverage -- is in
# test-edid-weighted-cov-{e2e,higher-order} and quality_reports/wcov; this
# file is the COMPATIBILITY matrix.)
# ============================================================

mk_opt_panel <- function(n = 320, seed = 5) {
  set.seed(seed); Tn <- 4L
  x1u <- runif(n, -2, 2); x2u <- rnorm(n)
  eta <- 0.5 * x1u + 0.2 * x2u
  P   <- exp(cbind(0, eta, 0.55 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.4 * x1u, 1); wu <- exp(0.6 * rnorm(n)); wu <- wu / mean(wu)
  cl  <- sample(1:24, n, replace = TRUE)
  rows <- lapply(1:Tn, function(tt) {
    ht <- (tt - 1) * (0.4 * x1u); tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, x2 = x2u, w = wu, cl = cl, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

DF   <- mk_opt_panel()
DF_C <- DF; DF_C$w <- 5                                       # constant weight column -> unit_weights == 1

# Run edid with explicit args + options(); return finite numeric leaves of att_gt (+ aggregate).
fit_num <- function(d, wt, args, .opts = list(), agg = "none") {
  op <- options(); on.exit(options(op), add = TRUE)
  if (length(.opts)) do.call(options, .opts)
  call_args <- c(list(data = d, yname = "y", idname = "id", tname = "t", gname = "g",
                      xformla = ~ x1 + x2, weightsname = wt, aggregate = agg, seed = 1L), args)
  if (is.null(call_args$bstrap)) call_args$bstrap <- FALSE
  f <- tryCatch(suppressWarnings(do.call(edid, call_args)), error = function(e) e)
  if (inherits(f, "error")) return(list(err = conditionMessage(f)))
  out <- c(f$att_gt$att, f$att_gt$se)
  if (!is.null(f[[agg]]$att.egt))     out <- c(out, f[[agg]]$att.egt, f[[agg]]$se.egt)
  if (!is.null(f[[agg]]$overall.att)) out <- c(out, f[[agg]]$overall.att, f[[agg]]$overall.se)
  out[is.finite(out)]
}

cmp <- function(lab, args, .opts = list(), agg = "none") {
  u  <- fit_num(DF,   NULL, args, .opts, agg)
  cc <- fit_num(DF_C, "w",  args, .opts, agg)
  uerr <- if (is.list(u))  u$err  else NULL          # fit_num returns list(err=) on error, else a numeric vector
  cerr <- if (is.list(cc)) cc$err else NULL
  expect_null(uerr, info = paste(lab, "unweighted error:",  uerr))
  expect_null(cerr, info = paste(lab, "const-weight error:", cerr))
  if (is.null(uerr) && is.null(cerr)) {
    expect_equal(length(cc), length(u), info = paste(lab, "shape"))
    expect_equal(cc, u, tolerance = 1e-7, info = lab)
  }
}

K  <- list(edid_omega_method = "kernel")
KO <- list(edid_omega_method = "kernel_orig")
S  <- list(edid_omega_method = "sieve")

test_that("constant weight column == unweighted across the full option matrix", {
  cmp("kernel|exp|efficient|plugin",  list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=FALSE, estimation_effect=FALSE), K)
  cmp("kernel|exp|efficient|ee",      list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=FALSE, estimation_effect=TRUE),  K)
  cmp("kernel|exp|efficient|misspec", list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=TRUE),  K)
  cmp("kernel|exp|efficient|higher",  list(weight_scheme="efficient", ratio_method="exp",    higher_order=TRUE),    K)
  cmp("kernel|direct|eff|misspec",    list(weight_scheme="efficient", ratio_method="direct", misspec_robust=TRUE),  K)
  cmp("korig|exp|eff|misspec",        list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=TRUE),  KO)
  cmp("sieve|exp|eff|misspec",        list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=TRUE),  S)
  cmp("sieve|exp|eff|ee",             list(weight_scheme="efficient", ratio_method="exp",    estimation_effect=TRUE, misspec_robust=FALSE), S)
  cmp("kernel|exp|averaged|misspec",  list(weight_scheme="averaged",  ratio_method="exp",    misspec_robust=TRUE),  K)
  cmp("kernel|exp|gmm|misspec",       list(weight_scheme="gmm",       ratio_method="exp",    misspec_robust=TRUE),  K)
  cmp("kernel|exp|uniform",           list(weight_scheme="uniform",   ratio_method="exp"),   K)
  cmp("shrink-none|misspec",          list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=TRUE, omega_cov_shrink="none"),       K)
  cmp("shrink-LW|misspec",            list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=TRUE, omega_cov_shrink="ledoit_wolf"),K)
  cmp("bsdf-ic|ee",                   list(weight_scheme="efficient", ratio_method="exp",    estimation_effect=TRUE, misspec_robust=FALSE, bs_df="ic"), S)
  cmp("ptpost|eff",                   list(weight_scheme="efficient", ratio_method="exp",    pt_assumption="post"), K)
  cmp("trim50|misspec",               list(weight_scheme="efficient", ratio_method="exp",    misspec_robust=TRUE, trim_level=50), K)
})

test_that("constant weight column == unweighted across aggregations + clustering", {
  base <- list(weight_scheme = "efficient", ratio_method = "exp", misspec_robust = TRUE)
  for (ag in c("group", "event_study", "calendar", "overall"))
    cmp(paste("agg", ag), base, K, agg = ag)
  cmp("clustered", c(base, list(clustervars = "cl")), K)
})

test_that("constant weight column == unweighted under the multiplier bootstrap (seeded band)", {
  cmp("mult-bootstrap",
      list(weight_scheme = "efficient", ratio_method = "exp", misspec_robust = FALSE,
           bstrap = TRUE, biters = 199L, cband = TRUE), K)
})
