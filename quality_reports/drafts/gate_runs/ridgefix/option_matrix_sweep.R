#!/usr/bin/env Rscript
# Option-matrix smoke sweep (full-adaptation audit item 4e). Representative combinations
# across the edid option matrix, on WEIGHTED no-cov data (where the n_eff fix is active),
# plus UNWEIGHTED cov combinations (where the cov-path lift is a structural no-op). Each
# combination must run without error and return finite att/se.
suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE)))

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

`%||%` <- function(a, b) if (is.null(a)) b else a
BASE <- list(yname = "outcome", tname = "time", idname = "unit", gname = "first_treat")
run <- function(desc, ...) {
  args <- modifyList(BASE, list(...))
  out <- tryCatch({
    fit <- suppressWarnings(suppressMessages(do.call(edid, args)))
    # pick whichever aggregate object this call produced (overall may be NULL for group/calendar)
    agg <- fit$overall %||% fit$event_study %||% fit$group %||% fit$calendar %||% fit$simple
    att <- agg$overall.att; se <- agg$overall.se
    ok  <- length(att) == 1L && is.finite(att) && length(se) == 1L && is.finite(se) && se > 0
    sprintf("  [%s] %-58s att=%+.4f se=%.4f", if (ok) "OK " else "BAD", desc,
            if (length(att)) att else NA_real_, if (length(se)) se else NA_real_)
  }, error = function(e) sprintf("  [ERR] %-58s %s", desc, conditionMessage(e)))
  cat(out, "\n")
  startsWith(trimws(out), "[OK")
}

dfw  <- mk_w()
dfwc <- mk_w(with_cluster = TRUE)
base <- list(yname = "outcome", tname = "time", idname = "unit", gname = "first_treat")
nbad <- 0L
chk  <- function(...) { if (!run(...)) nbad <<- nbad + 1L }

cat("=== WEIGHTED no-cov: omega_cov_shrink x weight_scheme ===\n")
for (sh in c("ridge", "ledoit_wolf", "none"))
  for (ws in c("efficient", "averaged", "gmm", "uniform"))
    chk(sprintf("w | shrink=%s scheme=%s", sh, ws),
        data = dfw, weightsname = "w", omega_cov_shrink = sh, weight_scheme = ws)

cat("\n=== WEIGHTED no-cov: ratio_method ===\n")
for (rm in c("exp", "direct"))
  chk(sprintf("w | ratio_method=%s ridge", rm),
      data = dfw, weightsname = "w", omega_cov_shrink = "ridge", ratio_method = rm)

cat("\n=== WEIGHTED no-cov: estimation_effect / misspec_robust ===\n")
chk("w | ridge estimation_effect=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", estimation_effect = TRUE)
chk("w | ledoit_wolf estimation_effect=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ledoit_wolf", estimation_effect = TRUE)
chk("w | ridge misspec_robust=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", misspec_robust = TRUE)

cat("\n=== WEIGHTED no-cov: clustervars + bootstraps ===\n")
chk("w | ridge clustervars",
    data = dfwc, weightsname = "w", omega_cov_shrink = "ridge", clustervars = "cl")
chk("w | ridge bstrap=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ridge", bstrap = TRUE, biters = 199L)
chk("w | ledoit_wolf bstrap=TRUE",
    data = dfw, weightsname = "w", omega_cov_shrink = "ledoit_wolf", bstrap = TRUE, biters = 199L)

cat("\n=== WEIGHTED no-cov: aggregate variants ===\n")
for (ag in c("event_study", "group", "calendar", "overall"))
  chk(sprintf("w | ridge aggregate=%s", ag),
      data = dfw, weightsname = "w", omega_cov_shrink = "ridge", aggregate = ag)

cat("\n=== UNWEIGHTED cov (cov-path ridge lift = structural no-op): smoothers ===\n")
dfx <- mk_w(); dfx$X <- rep(rnorm(length(unique(dfx$unit))), each = 7L)[seq_len(nrow(dfx))]
# rebuild a clean per-unit X
set.seed(9); uids <- sort(unique(dfx$unit)); xx <- rnorm(length(uids))
dfx$X <- xx[match(dfx$unit, uids)]
for (sh in c("ridge", "ledoit_wolf", "none"))
  chk(sprintf("cov(unweighted) | shrink=%s", sh),
      data = dfx, xformla = ~X, omega_cov_shrink = sh)

cat("\n=== TOOLKIT on a weighted ridge fit ===\n")
fitw <- suppressWarnings(suppressMessages(
  do.call(edid, modifyList(BASE, list(data = dfw, weightsname = "w", omega_cov_shrink = "ridge")))))
tk <- function(nm, expr) {
  out <- tryCatch({ force(expr); sprintf("  [OK ] %s", nm) },
                  error = function(e) { nbad <<- nbad + 1L; sprintf("  [ERR] %s : %s", nm, conditionMessage(e)) })
  cat(out, "\n")
}
# restricted (PT-Post) weighted ridge fit for the two-fit toolkit members
fitw_r <- suppressWarnings(suppressMessages(
  do.call(edid, modifyList(BASE, list(data = dfw, weightsname = "w",
          omega_cov_shrink = "ridge", pt_assumption = "post")))))
tk("edid_weights",  edid_weights(fitw))
tk("edid_sargan",   suppressWarnings(edid_sargan(fitw, data = dfw)))
tk("edid_hausman",  suppressWarnings(edid_hausman(fitw, fitw_r)))
tk("edid_frontier", suppressWarnings(edid_frontier(fitw, fitw_r)))
tk("print",         print(fitw))
tk("summary",       summary(fitw))

cat(sprintf("\n>>> OPTION-MATRIX SWEEP: %s  (%d bad combinations)\n",
            if (nbad == 0L) "ALL PASS" else "FAILURES", nbad))
