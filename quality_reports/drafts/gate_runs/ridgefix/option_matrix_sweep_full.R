#!/usr/bin/env Rscript
# ===========================================================================
# FULL option-matrix smoke sweep (full-adaptation audit item 4e) for the
# effective-n ridge/LW fix. Independent re-run + extension of
# option_matrix_sweep.R. Axes exercised (task spec):
#   ratio_method {exp,direct}
#   x weight_scheme {efficient,averaged,gmm,uniform}
#   x omega smoother {kernel, kernel_orig, sieve}  (cov path, via edid_omega_method)
#   x misspec_robust / estimation_effect / higher_order
#   x nocov_shrink {none, ledoit_wolf, ridge}
#   x weightsname {NULL, weighted col "w"}
#   x clustervars x both bootstraps (analytic + multiplier cband; bstrap)
#   x toolkit {edid_hausman, edid_sargan, edid_frontier, edid_adaptive, edid_weights}
#   x print / summary x aggregate variants x moment_set x bs_df="ic"
# EVERY combination must run without error and (for fits) return finite att/se>0.
# FAIL LOUDLY: nbad>0 => non-zero summary line; per-combo [ERR]/[BAD] printed.
# ===========================================================================
suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE)))

`%||%` <- function(a, b) if (is.null(a)) b else a
BASE <- list(yname = "outcome", tname = "time", idname = "unit", gname = "first_treat")
nbad <- 0L
ncombo <- 0L

# --- data factories --------------------------------------------------------
mk_w <- function(n_g3 = 18L, n_g5 = 16L, n_never = 24L, np = 7L, seed = 5L, disp = 1.0,
                 with_cluster = FALSE) {
  set.seed(seed); n <- n_g3 + n_g5 + n_never
  uid <- rep(seq_len(n), each = np); tid <- rep(seq_len(np), times = n)
  ft  <- c(rep(3L, n_g3 * np), rep(5L, n_g5 * np), rep(Inf, n_never * np))
  ufe <- rep(rnorm(n, 0, 1), each = np); tfe <- rep(seq(0, 1, length.out = np), times = n)
  e   <- rnorm(n * np, 0, 0.6)
  t3  <- (uid <= n_g3) & (tid >= 3L); t5 <- (uid > n_g3 & uid <= n_g3 + n_g5) & (tid >= 5L)
  df  <- data.frame(unit = uid, time = tid, outcome = ufe + tfe + e + 1.5 * t3 + 2.5 * t5,
                    first_treat = ft, w = rep(exp(rnorm(n, 0, disp)), each = np))
  if (with_cluster) df$cl <- rep(sample.int(8L, n, TRUE), each = np)
  df
}
mk_cov <- function(seed = 7L, ...) {
  df <- mk_w(seed = seed, ...)
  set.seed(seed + 100L); uids <- sort(unique(df$unit)); xx <- rnorm(length(uids))
  df$X <- xx[match(df$unit, uids)]
  # make outcome mildly depend on X so the cov path is non-degenerate
  df$outcome <- df$outcome + 0.5 * df$X
  df
}

dfw  <- mk_w()
dfwc <- mk_w(with_cluster = TRUE)
dfx  <- mk_cov()
dfxc <- mk_cov(seed = 7L); dfxc$cl <- rep(sample.int(8L, length(unique(dfxc$unit)), TRUE),
                                          each = 7L)[seq_len(nrow(dfxc))]

# --- runners ---------------------------------------------------------------
run <- function(desc, ...) {
  ncombo <<- ncombo + 1L
  args <- modifyList(BASE, list(...))
  agg_none <- isTRUE(args$aggregate == "none")   # aggregate="none" returns only att_gt by design
  out <- tryCatch({
    fit <- suppressWarnings(suppressMessages(do.call(edid, args)))
    if (agg_none) {
      a  <- fit$att_gt$att; s <- fit$att_gt$se
      ok <- length(a) > 0L && all(is.finite(a)) && all(is.finite(s)) && all(s > 0)
      sprintf("  [%s] %-62s att_gt cells=%d finite", if (ok) "OK " else "BAD", desc, length(a))
    } else {
      agg <- fit$overall %||% fit$event_study %||% fit$group %||% fit$calendar %||% fit$simple
      att <- agg$overall.att; se <- agg$overall.se
      ok  <- length(att) == 1L && is.finite(att) && length(se) == 1L && is.finite(se) && se > 0
      sprintf("  [%s] %-62s att=%+.4f se=%.4f", if (ok) "OK " else "BAD", desc,
              if (length(att)) att else NA_real_, if (length(se)) se else NA_real_)
    }
  }, error = function(e) sprintf("  [ERR] %-62s %s", desc, conditionMessage(e)))
  cat(out, "\n")
  if (!startsWith(trimws(out), "[OK")) nbad <<- nbad + 1L
  invisible(NULL)
}

# ---------------------------------------------------------------------------
cat("=== [A] WEIGHTED no-cov: omega_cov_shrink x weight_scheme (the fix's live cells) ===\n")
for (sh in c("ridge", "ledoit_wolf", "none"))
  for (ws in c("efficient", "averaged", "gmm", "uniform"))
    run(sprintf("w | shrink=%s scheme=%s", sh, ws),
        data = dfw, weightsname = "w", omega_cov_shrink = sh, weight_scheme = ws)

cat("\n=== [B] UNWEIGHTED no-cov: same matrix (byte-identity / no-op leg) ===\n")
for (sh in c("ridge", "ledoit_wolf", "none"))
  for (ws in c("efficient", "averaged", "gmm", "uniform"))
    run(sprintf("nw | shrink=%s scheme=%s", sh, ws),
        data = dfw, omega_cov_shrink = sh, weight_scheme = ws)

cat("\n=== [C] WEIGHTED no-cov: ratio_method x shrink ===\n")
for (rm in c("exp", "direct"))
  for (sh in c("ridge", "ledoit_wolf"))
    run(sprintf("w | ratio_method=%s shrink=%s", rm, sh),
        data = dfw, weightsname = "w", omega_cov_shrink = sh, ratio_method = rm)

cat("\n=== [D] WEIGHTED no-cov: misspec_robust / estimation_effect / higher_order ===\n")
run("w | ridge estimation_effect=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", estimation_effect = TRUE)
run("w | ledoit_wolf estimation_effect=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ledoit_wolf", estimation_effect = TRUE)
run("w | ridge misspec_robust=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", misspec_robust = TRUE)
run("w | ridge misspec_robust=FALSE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", misspec_robust = FALSE)
# higher_order is meaningless without xformla (documented no-cov guard: zero higher-order
# variance by construction) -- the valid weighted no-cov axis is ee only; ho is exercised
# on the cov path in section [I]. Keep an explicit higher_order=FALSE no-cov to pin the guard.
run("w | ridge ee=TRUE ho=FALSE (explicit)",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge",
    estimation_effect = TRUE, higher_order = FALSE)

cat("\n=== [E] WEIGHTED no-cov: clustervars + both bootstraps + cband_method ===\n")
run("w | ridge clustervars",
    data = dfwc, weightsname = "w", omega_cov_shrink = "ridge", clustervars = "cl")
run("w | ridge bstrap=TRUE (multiplier)",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", bstrap = TRUE, biters = 199L)
run("w | ledoit_wolf bstrap=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ledoit_wolf", bstrap = TRUE, biters = 199L)
run("w | ridge cband_method=multiplier",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", cband_method = "multiplier")
run("w | ridge cband_method=analytic",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", cband_method = "analytic")
run("w | ridge clustervars + bstrap",
    data = dfwc, weightsname = "w", omega_cov_shrink = "ridge",
    clustervars = "cl", bstrap = TRUE, biters = 199L)

cat("\n=== [F] WEIGHTED no-cov: aggregate variants + moment_set + pt_assumption ===\n")
for (ag in c("event_study", "group", "calendar", "overall", "none"))
  run(sprintf("w | ridge aggregate=%s", ag),
      data = dfw, weightsname = "w", omega_cov_shrink = "ridge", aggregate = ag)
ms_valid <- rbind(data.frame(g = 3, gp = 3, tpre = c(1, 2)),
                  data.frame(g = 3, gp = 0, tpre = 2),
                  data.frame(g = 5, gp = 5, tpre = c(1, 2, 3, 4)),
                  data.frame(g = 5, gp = 0, tpre = 4))
run("w | ridge moment_set=<valid df>",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", moment_set = ms_valid)
run("w | ridge pt_assumption=post",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", pt_assumption = "post")
run("w | ridge pt_assumption=post scheme=uniform",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge",
    pt_assumption = "post", weight_scheme = "uniform")

cat("\n=== [G] COV path (unweighted): omega smoother {kernel,kernel_orig,sieve} x shrink ===\n")
for (sm in c("kernel", "kernel_orig", "sieve")) {
  options(edid_omega_method = sm)
  for (sh in c("ridge", "ledoit_wolf", "none"))
    run(sprintf("cov | smoother=%-11s shrink=%s", sm, sh),
        data = dfx, xformla = ~X, omega_cov_shrink = sh)
}
options(edid_omega_method = NULL)

cat("\n=== [H] COV path (unweighted): smoother x weight_scheme x ratio_method ===\n")
for (sm in c("kernel", "sieve")) {
  options(edid_omega_method = sm)
  for (ws in c("efficient", "averaged", "gmm"))
    for (rm in c("exp", "direct"))
      run(sprintf("cov | sm=%-6s scheme=%-9s rm=%s", sm, ws, rm),
          data = dfx, xformla = ~X, weight_scheme = ws, ratio_method = rm)
}
options(edid_omega_method = NULL)

cat("\n=== [I] COV path: misspec_robust / estimation_effect / higher_order / clustervars / bstrap / bs_df ===\n")
run("cov | kernel estimation_effect=TRUE",
    data = dfx, xformla = ~X, estimation_effect = TRUE)
run("cov | kernel higher_order=TRUE",
    data = dfx, xformla = ~X, higher_order = TRUE)
run("cov | kernel misspec_robust=TRUE",
    data = dfx, xformla = ~X, misspec_robust = TRUE)
run("cov | kernel clustervars",
    data = dfxc, xformla = ~X, clustervars = "cl")
run("cov | kernel bstrap=TRUE",
    data = dfx, xformla = ~X, bstrap = TRUE, biters = 199L)
run("cov | kernel bs_df=ic",
    data = dfx, xformla = ~X, bs_df = "ic")
options(edid_omega_method = "sieve")
run("cov | sieve estimation_effect=TRUE",
    data = dfx, xformla = ~X, estimation_effect = TRUE)
options(edid_omega_method = NULL)

cat("\n=== [J] Weighted-cov contract: weightsname + xformla must error by design (scoped out) ===\n")
ncombo <- ncombo + 1L
out_j <- tryCatch({
  suppressWarnings(suppressMessages(do.call(edid, modifyList(BASE,
    list(data = dfx, xformla = ~X, weightsname = "w")))))
  nbad <<- nbad + 1L
  "  [BAD] weighted-cov did NOT error (contract says it must)"
}, error = function(e) sprintf("  [OK ] weighted-cov errors by design: %s",
                               substr(conditionMessage(e), 1, 60)))
cat(out_j, "\n")

cat("\n=== [K] TOOLKIT on a WEIGHTED ridge fit (unrestricted + PT-Post restricted) ===\n")
fitw   <- suppressWarnings(suppressMessages(
  do.call(edid, modifyList(BASE, list(data = dfw, weightsname = "w", omega_cov_shrink = "ridge")))))
fitw_r <- suppressWarnings(suppressMessages(
  do.call(edid, modifyList(BASE, list(data = dfw, weightsname = "w",
          omega_cov_shrink = "ridge", pt_assumption = "post")))))
tk <- function(nm, expr) {
  ncombo <<- ncombo + 1L
  out <- tryCatch({ v <- force(expr); sprintf("  [OK ] %s", nm) },
                  error = function(e) { nbad <<- nbad + 1L; sprintf("  [ERR] %s : %s", nm, conditionMessage(e)) })
  cat(out, "\n")
}
tk("edid_weights",   edid_weights(fitw))
tk("edid_sargan",    suppressWarnings(edid_sargan(fitw_r, data = dfw)))
tk("edid_hausman",   suppressWarnings(edid_hausman(fitw, fitw_r)))
tk("edid_frontier",  suppressWarnings(edid_frontier(fitw, fitw_r)))
tk("edid_adaptive",  suppressWarnings(edid_adaptive(fitw, fitw_r)))
tk("print",          { z <- capture.output(print(fitw)); invisible(z) })
tk("summary",        { z <- capture.output(summary(fitw)); invisible(z) })
tk("$args refit snapshot present", stopifnot(!is.null(fitw$args)))

cat("\n=== [L] TOOLKIT on an UNWEIGHTED cov ridge fit (invariant leg) ===\n")
fitx   <- suppressWarnings(suppressMessages(
  do.call(edid, modifyList(BASE, list(data = dfx, xformla = ~X, omega_cov_shrink = "ridge")))))
fitx_r <- suppressWarnings(suppressMessages(
  do.call(edid, modifyList(BASE, list(data = dfx, xformla = ~X,
          omega_cov_shrink = "ridge", pt_assumption = "post")))))
tk("cov edid_weights",  edid_weights(fitx))
tk("cov edid_sargan",   suppressWarnings(edid_sargan(fitx_r, data = dfx)))
tk("cov edid_hausman",  suppressWarnings(edid_hausman(fitx, fitx_r)))
tk("cov edid_frontier", suppressWarnings(edid_frontier(fitx, fitx_r)))
tk("cov print",         { z <- capture.output(print(fitx)); invisible(z) })
tk("cov summary",       { z <- capture.output(summary(fitx)); invisible(z) })

cat(sprintf("\n>>> FULL OPTION-MATRIX SWEEP: %s  (%d combinations, %d bad)\n",
            if (nbad == 0L) "ALL PASS" else "FAILURES", ncombo, nbad))
if (nbad > 0L) quit(status = 1L)
