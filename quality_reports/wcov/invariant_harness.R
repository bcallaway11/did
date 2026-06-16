## =====================================================================
## quality_reports/wcov/invariant_harness.R
##
## Byte-identical INVARIANT guard for the weighted-covariate-path rollout.
## Captures the CURRENT (pre-change) numbers on the two paths that MUST NOT
## move, then re-runs and asserts byte-identical after every phase:
##   (1) UNWEIGHTED covariate path  (xformla = ~x1, no weightsname)
##   (2) NO-COV WEIGHTED path        (xformla = NULL, weightsname = "w")
## across 3 smoothers x ratio_methods x EE-configs x cov-shrink.
##
## Usage (from the package repo root):
##   Rscript quality_reports/wcov/invariant_harness.R baseline   # capture ref
##   Rscript quality_reports/wcov/invariant_harness.R check      # compare to ref
## =====================================================================
suppressPackageStartupMessages({ library(pkgload) })
pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE)
options(edid_mc_cores = 1L)

REF <- "/Users/pcostag/Documents/GitHub/did/quality_reports/wcov/baseline_ref.rds"
mode <- (function() { a <- commandArgs(TRUE); if (length(a)) a[1] else "check" })()

## ---- covariate DGP (from tests/testthat/test-edid-ach-correction.R) + a weight col ----
make_panel <- function(n = 300, seed = 1) {
  set.seed(seed)
  Tn <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gcat <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  ## dispersed, time-INVARIANT per-unit weight (Kish n_eff well below n)
  wu <- exp(0.8 * rnorm(n)); wu <- wu / mean(wu)
  rows <- lapply(1:Tn, function(tt) {
    ht <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gcat) & tt >= gcat, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gcat), gcat, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

## ---- one fit -> compact numeric fingerprint ----
fingerprint <- function(args, opts = list()) {
  old <- options(); on.exit(options(old), add = TRUE)
  if (length(opts)) do.call(options, opts)
  fit <- tryCatch(suppressWarnings(do.call(edid, args)), error = function(e) e)
  if (inherits(fit, "error")) return(list(err = conditionMessage(fit)))
  eif <- fit$eif
  list(att = fit$att_gt$att, se = fit$att_gt$se,
       eif_ss = if (is.null(eif)) NA_real_ else sum(eif * eif, na.rm = TRUE),
       eif_dim = if (is.null(eif)) NA else paste(dim(as.matrix(eif)), collapse = "x"))
}

## ---- config grid (the must-stay-identical set) ----
base_cov <- function(df, extra = list(), wt = NULL)
  c(list(data = df, yname = "y", idname = "id", tname = "t", gname = "g",
         xformla = ~ x1, weight_scheme = "efficient", aggregate = "none",
         bstrap = FALSE, seed = 1L, weightsname = wt), extra)
base_nocov <- function(df, extra = list(), wt = NULL)
  c(list(data = df, yname = "y", idname = "id", tname = "t", gname = "g",
         xformla = NULL, weight_scheme = "efficient", aggregate = "none",
         bstrap = FALSE, seed = 1L, weightsname = wt), extra)

build_grid <- function(df) {
  g <- list()
  smoothers <- c("kernel", "kernel_orig", "sieve")
  for (sm in smoothers) for (rm in c("exp", "direct")) {
    op <- list(edid_omega_method = sm)
    ## plug-in
    g[[sprintf("cov|%s|%s|plugin", sm, rm)]] <- list(
      args = base_cov(df, list(ratio_method = rm, misspec_robust = FALSE, estimation_effect = FALSE)), opts = op)
    ## estimation_effect (ACH) only
    g[[sprintf("cov|%s|%s|ee", sm, rm)]] <- list(
      args = base_cov(df, list(ratio_method = rm, misspec_robust = FALSE, estimation_effect = TRUE)), opts = op)
    ## misspec_robust default (the headline path)
    g[[sprintf("cov|%s|%s|misspec", sm, rm)]] <- list(
      args = base_cov(df, list(ratio_method = rm, misspec_robust = TRUE)), opts = op)
  }
  ## higher_order + cov-shrink variants on the default smoother
  g[["cov|kernel|exp|higher"]] <- list(args = base_cov(df, list(ratio_method = "exp", higher_order = TRUE)), opts = list(edid_omega_method = "kernel"))
  g[["cov|kernel|exp|shrink_none"]] <- list(args = base_cov(df, list(ratio_method = "exp", omega_cov_shrink = "none")), opts = list(edid_omega_method = "kernel"))
  ## NO-COV WEIGHTED (must stay identical: the existing weighted path)
  for (rm in c("exp", "direct")) {
    g[[sprintf("nocovW|%s|plugin", rm)]] <- list(args = base_nocov(df, list(ratio_method = rm), wt = "w"), opts = list())
    g[[sprintf("nocovW|%s|ee", rm)]] <- list(args = base_nocov(df, list(ratio_method = rm, estimation_effect = TRUE), wt = "w"), opts = list())
  }
  g
}

df <- make_panel(n = 300, seed = 21)
grid <- build_grid(df)
res <- lapply(grid, function(cfg) fingerprint(cfg$args, cfg$opts))

if (identical(mode, "baseline")) {
  dir.create(dirname(REF), showWarnings = FALSE, recursive = TRUE)
  saveRDS(res, REF)
  ne <- sum(vapply(res, function(r) is.null(r$err), TRUE))
  cat(sprintf("[BASELINE] captured %d configs (%d ok, %d errored) -> %s\n",
              length(res), ne, length(res) - ne, REF))
  for (k in names(res)) if (!is.null(res[[k]]$err)) cat("  ERR ", k, ": ", substr(res[[k]]$err, 1, 70), "\n", sep = "")
} else {
  stopifnot(file.exists(REF)); ref <- readRDS(REF)
  worst <- 0; nfail <- 0
  for (k in names(ref)) {
    r0 <- ref[[k]]; r1 <- res[[k]]
    if (!is.null(r0$err) || !is.null(r1$err)) {
      if (!identical(is.null(r0$err), is.null(r1$err))) { cat("FAIL(err-state) ", k, "\n"); nfail <- nfail + 1 }
      next
    }
    d <- max(max(abs(r1$att - r0$att), na.rm = TRUE),
             max(abs(r1$se  - r0$se ), na.rm = TRUE),
             abs(r1$eif_ss - r0$eif_ss))
    worst <- max(worst, d)
    if (!is.finite(d) || d > 1e-12) { cat(sprintf("FAIL %-28s max|diff|=%.3e\n", k, d)); nfail <- nfail + 1 }
  }
  cat(sprintf("\n[CHECK] %d configs | worst max|diff| = %.3e | %s\n",
              length(ref), worst, if (nfail == 0) "ALL BYTE-IDENTICAL (PASS)" else sprintf("%d FAILED", nfail)))
  quit(status = if (nfail == 0) 0 else 1)
}
