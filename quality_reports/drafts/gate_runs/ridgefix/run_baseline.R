#!/usr/bin/env Rscript
# Effective-n ridge fix: baseline + self-check harness.
# Usage: Rscript run_baseline.R <out.rds> [label]
# Runs a fixed battery of unweighted fits (no-cov + cov; none/ledoit_wolf/ridge;
# small staggered designs) plus weighted ridge fits, and records per-cell
# fingerprints (att, se, lambda, cond#). The UNWEIGHTED battery must be
# byte-identical before/after the n_eff patch (the fix is a no-op at w-equal).

suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE, export_all = FALSE)))

args   <- commandArgs(trailingOnly = TRUE)
outrds <- if (length(args) >= 1L) args[[1]] else "/tmp/gate_runs/ridgefix/baseline.rds"
label  <- if (length(args) >= 2L) args[[2]] else "baseline"

# ---- data factories (mirror tests/testthat/helper-edid.R; self-contained) ----
make_panel_1cohort <- function(n_treat = 20L, n_never = 20L, n_periods = 5L, seed = 42L) {
  set.seed(seed); n <- n_treat + n_never
  unit_ids <- rep(seq_len(n), each = n_periods)
  time_ids <- rep(seq_len(n_periods), times = n)
  ft <- c(rep(3L, n_treat * n_periods), rep(Inf, n_never * n_periods))
  unit_fe <- rep(rnorm(n, 0, 1), each = n_periods)
  time_fe <- rep(seq(0, 0.5, length.out = n_periods), times = n)
  noise   <- rnorm(n * n_periods, 0, 0.5)
  tp <- (unit_ids <= n_treat) & (time_ids >= 3L)
  data.frame(unit = unit_ids, time = time_ids,
             outcome = unit_fe + time_fe + noise + 2 * tp,
             first_treat = ft)
}
make_panel_2cohort <- function(n_g3 = 15L, n_g5 = 15L, n_never = 20L, n_periods = 7L, seed = 123L) {
  set.seed(seed); n <- n_g3 + n_g5 + n_never
  unit_ids <- rep(seq_len(n), each = n_periods)
  time_ids <- rep(seq_len(n_periods), times = n)
  ft <- c(rep(3L, n_g3 * n_periods), rep(5L, n_g5 * n_periods), rep(Inf, n_never * n_periods))
  unit_fe <- rep(rnorm(n, 0, 1), each = n_periods)
  time_fe <- rep(seq(0, 1, length.out = n_periods), times = n)
  noise   <- rnorm(n * n_periods, 0, 0.5)
  t3 <- (unit_ids <= n_g3) & (time_ids >= 3L)
  t5 <- (unit_ids > n_g3 & unit_ids <= n_g3 + n_g5) & (time_ids >= 5L)
  data.frame(unit = unit_ids, time = time_ids,
             outcome = unit_fe + time_fe + noise + 1.5 * t3 + 2.5 * t5,
             first_treat = ft)
}
# covariate panel: time-invariant X driving the outcome (so cov path is non-trivial)
make_panel_cov <- function(n_g3 = 18L, n_g5 = 0L, n_never = 22L, n_periods = 6L, seed = 321L) {
  set.seed(seed); n <- n_g3 + n_g5 + n_never
  unit_ids <- rep(seq_len(n), each = n_periods)
  time_ids <- rep(seq_len(n_periods), times = n)
  ft <- c(rep(3L, n_g3 * n_periods),
          if (n_g5 > 0) rep(5L, n_g5 * n_periods) else integer(0),
          rep(Inf, n_never * n_periods))
  x_unit <- rnorm(n, 0, 1)
  x       <- rep(x_unit, each = n_periods)
  unit_fe <- rep(rnorm(n, 0, 0.7), each = n_periods)
  time_fe <- rep(seq(0, 0.6, length.out = n_periods), times = n)
  noise   <- rnorm(n * n_periods, 0, 0.5)
  tp <- (unit_ids <= n_g3) & (time_ids >= 3L)
  data.frame(unit = unit_ids, time = time_ids,
             outcome = unit_fe + time_fe + 0.8 * x * time_ids / n_periods + noise + 2 * tp,
             first_treat = ft, X = x)
}

# ---- weighted variant: attach a dispersed per-unit weight column -------------
add_weights <- function(df, seed = 999L, dispersion = 2.0) {
  set.seed(seed)
  uids <- sort(unique(df$unit))
  w_unit <- exp(rnorm(length(uids), 0, dispersion))     # heavy-tailed => low Kish ESS
  df$w <- w_unit[match(df$unit, uids)]
  df
}

# ---- fingerprint extractor ----------------------------------------------------
fp_cells <- function(fit) {
  cells <- fit$cells
  if (is.null(cells)) return(NULL)
  do.call(rbind, lapply(seq_along(cells), function(k) {
    cl <- cells[[k]]
    g  <- if (!is.null(cl$group)) cl$group else NA_real_
    t  <- if (!is.null(cl$time))  cl$time  else NA_real_
    lam <- if (!is.null(cl$nocov_shrink_lambda))
             suppressWarnings(as.numeric(cl$nocov_shrink_lambda)[1]) else NA_real_
    cn  <- if (!is.null(cl$condition_num))
             suppressWarnings(as.numeric(cl$condition_num)[1]) else NA_real_
    data.frame(k = k, g = g, t = t,
               att = suppressWarnings(as.numeric(cl$att)[1]),
               se  = suppressWarnings(as.numeric(cl$se)[1]),
               lambda = lam, cond = cn)
  }))
}

fit_one <- function(df, weightsname = NULL, xformla = NULL, shrink = "ridge", tag = "") {
  out <- tryCatch({
    fit <- suppressWarnings(suppressMessages(
      edid(yname = "outcome", tname = "time", idname = "unit",
           gname = "first_treat", data = df,
           xformla = xformla, weightsname = weightsname,
           omega_cov_shrink = shrink)))
    oa <- if (!is.null(fit$overall)) fit$overall$overall.att else NA_real_
    os <- if (!is.null(fit$overall)) fit$overall$overall.se  else NA_real_
    list(ok = TRUE,
         att = suppressWarnings(as.numeric(oa)[1]),
         se  = suppressWarnings(as.numeric(os)[1]),
         cells = fp_cells(fit))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  out$tag <- tag
  out
}

# ---- battery ------------------------------------------------------------------
designs <- list(
  d1 = make_panel_1cohort(n_treat = 20L, n_never = 20L, n_periods = 5L, seed = 42L),
  d2 = make_panel_2cohort(n_g3 = 15L, n_g5 = 15L, n_never = 20L, n_periods = 7L, seed = 123L),
  d3 = make_panel_2cohort(n_g3 = 10L, n_g5 = 12L, n_never = 14L, n_periods = 6L, seed = 7L)
)
cov_designs <- list(
  c1 = make_panel_cov(n_g3 = 18L, n_never = 22L, n_periods = 6L, seed = 321L),
  c2 = make_panel_cov(n_g3 = 14L, n_g5 = 12L, n_never = 18L, n_periods = 7L, seed = 88L)
)

results <- list()

## UNWEIGHTED no-cov: none / ledoit_wolf / ridge
for (dn in names(designs)) for (sh in c("none", "ledoit_wolf", "ridge")) {
  key <- paste0("nocov_", dn, "_", sh)
  results[[key]] <- fit_one(designs[[dn]], shrink = sh, tag = key)
}
## UNWEIGHTED cov: none / ledoit_wolf / ridge
for (dn in names(cov_designs)) for (sh in c("none", "ledoit_wolf", "ridge")) {
  key <- paste0("cov_", dn, "_", sh)
  results[[key]] <- fit_one(cov_designs[[dn]], xformla = ~X, shrink = sh, tag = key)
}
## WEIGHTED ridge (no-cov + cov): these MOVE under the fix (recorded, not invariant)
for (dn in names(designs)) {
  dfw <- add_weights(designs[[dn]], seed = 1000L + as.integer(factor(dn)))
  key <- paste0("wnocov_", dn, "_ridge")
  results[[key]] <- fit_one(dfw, weightsname = "w", shrink = "ridge", tag = key)
}
# NOTE: weighted covariate path is UNCONDITIONALLY scoped out in edid() (errors by
# design: weighted kernel/sieve EE not yet derived). So the cov-path ridge lift only
# ever runs UNWEIGHTED -> n_eff == n there always -> site (3) is a structural no-op
# today (still wired to n_eff so it is correct-by-construction if the guard is lifted).
# We therefore do NOT attempt weighted-cov fits here (they would error by design).
## WEIGHTED ledoit_wolf (no-cov) -- LW intensity also moves under the fix
for (dn in names(designs)) {
  dfw <- add_weights(designs[[dn]], seed = 3000L + as.integer(factor(dn)))
  key <- paste0("wnocov_", dn, "_ledoit_wolf")
  results[[key]] <- fit_one(dfw, weightsname = "w", shrink = "ledoit_wolf", tag = key)
}

saveRDS(list(label = label, results = results, ts = Sys.time()), outrds)
cat(sprintf("[%s] wrote %s  (%d configs)\n", label, outrds, length(results)))
# quick console summary
for (k in names(results)) {
  r <- results[[k]]
  if (isTRUE(r$ok)) cat(sprintf("  %-28s att=%+.8f se=%.8f\n", k, r$att, r$se))
  else              cat(sprintf("  %-28s ERROR: %s\n", k, r$msg))
}
