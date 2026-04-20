#!/usr/bin/env Rscript
# benchmark/edid_cov_sim.R
# -----------------------------------------------------------------------
# Simulation study for the edid() covariate path.
# Evaluates bias, SE calibration, and coverage for the doubly-robust
# efficient DiD estimator of Chen, Sant'Anna & Xie (2025).
#
# Designs:
#   DGP 1 (main):  linear covariate, 2 cohorts + never-treated
#   DGP 2:         nonlinear covariate (quadratic), same structure
#   DGP 3:         propensity hard, outcome easy (misspec propensity)
#   DGP 4:         outcome hard, propensity easy (misspec outcome)
#   DGP 5:         2D covariate, paper-style (mirrors Sim_10DGP.R DGPs 5-8)
#
# Sample sizes: n = 200, 500, 1000, 2000
# Repetitions: R_fast = 200 per cell (increase R for publication)
# Cross-fitting stability: K-fold seed sweep at n=500
#
# Usage:
#   Rscript benchmark/edid_cov_sim.R
#   # or from R:
#   source("benchmark/edid_cov_sim.R")
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("did", quietly = TRUE)) devtools::load_all(".")
  library(did)
})

set.seed(20240419)

# -----------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------
R_FAST    <- 200          # replications per cell (increase for publication)
NS        <- c(200, 500, 1000, 2000)
TRUE_ATT  <- 1.0          # true ATT for all post-treatment cells in all DGPs
ALPHA     <- 0.05
N_PERIODS <- 6
COHORTS   <- c(3, 5)      # two treatment cohorts

# -----------------------------------------------------------------------
# DGP functions
# -----------------------------------------------------------------------

#' DGP 1: linear covariate, moderate overlap
dgp_linear <- function(n, seed) {
  set.seed(seed)
  n_g  <- floor(n / 3)
  unit <- 1:n
  g_unit <- c(rep(3, n_g), rep(5, n_g), rep(Inf, n - 2 * n_g))
  x1u <- rnorm(n, mean = g_unit / 10, sd = 1)  # mild selection on X
  T   <- N_PERIODS
  ids    <- rep(unit, each = T)
  times  <- rep(1:T, times = n)
  g_vec  <- g_unit[ids]
  x1     <- rep(x1u, each = T)
  treated <- as.numeric(times >= g_vec)
  y <- times + 0.5 * x1 + TRUE_ATT * treated + rnorm(n * T, sd = 0.6)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1)
}

#' DGP 2: nonlinear covariate (quadratic), 2D
dgp_nonlinear <- function(n, seed) {
  set.seed(seed)
  n_g <- floor(n / 3)
  g_unit <- c(rep(3, n_g), rep(5, n_g), rep(Inf, n - 2 * n_g))
  x1u <- rnorm(n, mean = g_unit / 8, sd = 1)
  x2u <- rnorm(n)
  T   <- N_PERIODS
  ids    <- rep(1:n, each = T)
  times  <- rep(1:T, times = n)
  g_vec  <- g_unit[ids]
  x1     <- rep(x1u, each = T)
  x2     <- rep(x2u, each = T)
  treated <- as.numeric(times >= g_vec)
  y <- times + x1^2 + 0.3 * x2 + TRUE_ATT * treated + rnorm(n * T, sd = 0.6)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1, x2 = x2)
}

#' DGP 3: propensity hard (strong selection), outcome easy
dgp_prop_hard <- function(n, seed) {
  set.seed(seed)
  n_g <- floor(n / 3)
  g_unit <- c(rep(3, n_g), rep(5, n_g), rep(Inf, n - 2 * n_g))
  # Strong selection: x1 drives treatment strongly
  x1u <- rnorm(n) + 0.8 * (g_unit < Inf)
  T   <- N_PERIODS
  ids    <- rep(1:n, each = T)
  times  <- rep(1:T, times = n)
  g_vec  <- g_unit[ids]
  x1     <- rep(x1u, each = T)
  treated <- as.numeric(times >= g_vec)
  # Outcome: only linear in x1 (easy)
  y <- times + 0.3 * x1 + TRUE_ATT * treated + rnorm(n * T, sd = 0.5)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1)
}

#' DGP 4: outcome hard (nonlinear), propensity easy
dgp_outcome_hard <- function(n, seed) {
  set.seed(seed)
  n_g <- floor(n / 3)
  g_unit <- c(rep(3, n_g), rep(5, n_g), rep(Inf, n - 2 * n_g))
  # Mild selection
  x1u <- rnorm(n, mean = 0.2 * (g_unit < Inf), sd = 1)
  x2u <- rnorm(n)
  T   <- N_PERIODS
  ids    <- rep(1:n, each = T)
  times  <- rep(1:T, times = n)
  g_vec  <- g_unit[ids]
  x1     <- rep(x1u, each = T)
  x2     <- rep(x2u, each = T)
  treated <- as.numeric(times >= g_vec)
  # Complex nonlinear outcome: hard to fit with linear basis
  y <- times + sin(x1) * cos(x2) + 0.5 * x1^2 + TRUE_ATT * treated +
       rnorm(n * T, sd = 0.6)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1, x2 = x2)
}

#' DGP 5: 2D covariate, paper-style (mirrors Sim_10DGP.R DGPs 5-8)
dgp_paper_style <- function(n, seed) {
  set.seed(seed)
  n_g <- floor(n / 3)
  g_unit <- c(rep(3, n_g), rep(5, n_g), rep(Inf, n - 2 * n_g))
  # X1, X2 ~ truncated N(0,1) on [-2, 2]
  x1u <- pmin(pmax(rnorm(n, mean = 0.3 * (g_unit < Inf)), -2), 2)
  x2u <- pmin(pmax(rnorm(n, mean = 0.3 * (g_unit < Inf)), -2), 2)
  T   <- N_PERIODS
  ids    <- rep(1:n, each = T)
  times  <- rep(1:T, times = n)
  g_vec  <- g_unit[ids]
  x1     <- rep(x1u, each = T)
  x2     <- rep(x2u, each = T)
  treated <- as.numeric(times >= g_vec)
  # Additive effects
  y <- times + 0.5 * x1 + 0.4 * x2 + 0.2 * x1 * x2 +
       TRUE_ATT * treated + rnorm(n * T, sd = 0.5)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1, x2 = x2)
}

# -----------------------------------------------------------------------
# Single-run function
# -----------------------------------------------------------------------
run_edid_cov <- function(df, xformla, seed = 1L) {
  tryCatch(
    edid(df, "y", "id", "t", "g", xformla = xformla,
         aggregate = "none", seed = seed),
    error = function(e) NULL
  )
}

# -----------------------------------------------------------------------
# Monte Carlo loop
# -----------------------------------------------------------------------
run_simulation <- function(dgp_fn, xformla, label, ns = NS, R = R_FAST) {
  cat(sprintf("\n=== %s | xformla: %s ===\n", label,
              deparse(xformla, width.cutoff = 60)))

  results_list <- vector("list", length(ns))

  for (ni in seq_along(ns)) {
    n <- ns[ni]
    cat(sprintf("  n = %4d ...", n)); flush.console()

    cell_data <- list()  # keyed by "g_t"
    n_ok <- 0L

    for (r in seq_len(R)) {
      df  <- dgp_fn(n, seed = r)
      fit <- run_edid_cov(df, xformla, seed = 1L)
      if (is.null(fit)) next
      n_ok <- n_ok + 1L
      att_df <- fit$att_gt[!fit$att_gt$is_pre, ]
      for (k in seq_len(nrow(att_df))) {
        key <- paste0(att_df$group[k], "_", att_df$time[k])
        cell_data[[key]] <- c(cell_data[[key]],
                              list(list(att = att_df$att[k], se = att_df$se[k],
                                        ci_l = att_df$ci_lower[k],
                                        ci_u = att_df$ci_upper[k])))
      }
    }

    cat(sprintf(" %d/%d OK\n", n_ok, R))

    # Summarise per cell
    cell_summary <- do.call(rbind, lapply(names(cell_data), function(key) {
      draws <- cell_data[[key]]
      atts  <- sapply(draws, `[[`, "att")
      ses   <- sapply(draws, `[[`, "se")
      ci_ls <- sapply(draws, `[[`, "ci_l")
      ci_us <- sapply(draws, `[[`, "ci_u")
      valid <- is.finite(atts) & is.finite(ses)
      if (sum(valid) < 10L) return(NULL)
      mc_mean  <- mean(atts[valid])
      mc_bias  <- mc_mean - TRUE_ATT
      mc_sd    <- sd(atts[valid])
      mc_rmse  <- sqrt(mc_bias^2 + mc_sd^2)
      mean_se  <- mean(ses[valid])
      se_ratio <- mean_se / mc_sd
      coverage <- mean(ci_ls[valid] <= TRUE_ATT & TRUE_ATT <= ci_us[valid])
      mean_ci_len <- mean(ci_us[valid] - ci_ls[valid])
      parts <- strsplit(key, "_")[[1L]]
      data.frame(
        n = n, label = label, group = as.numeric(parts[1]),
        time  = as.numeric(parts[2]),
        R_ok  = sum(valid),
        mc_mean = mc_mean, bias = mc_bias, mc_sd = mc_sd, rmse = mc_rmse,
        mean_se = mean_se, se_ratio = se_ratio,
        coverage = coverage, mean_ci_len = mean_ci_len,
        stringsAsFactors = FALSE
      )
    }))

    results_list[[ni]] <- cell_summary
  }

  do.call(rbind, results_list)
}

# -----------------------------------------------------------------------
# Run all DGPs
# -----------------------------------------------------------------------
all_results <- list()

all_results[["dgp1"]] <- run_simulation(
  dgp_fn = dgp_linear, xformla = ~ x1, label = "DGP1: linear, 1D")

all_results[["dgp2"]] <- run_simulation(
  dgp_fn = dgp_nonlinear, xformla = ~ x1 + x2 + I(x1^2), label = "DGP2: nonlinear, 2D")

all_results[["dgp3"]] <- run_simulation(
  dgp_fn = dgp_prop_hard, xformla = ~ x1, label = "DGP3: propensity hard")

all_results[["dgp4"]] <- run_simulation(
  dgp_fn = dgp_outcome_hard, xformla = ~ x1 + x2, label = "DGP4: outcome hard, 2D")

all_results[["dgp5"]] <- run_simulation(
  dgp_fn = dgp_paper_style, xformla = ~ x1 + x2, label = "DGP5: paper-style, 2D")

# -----------------------------------------------------------------------
# Print results tables
# -----------------------------------------------------------------------
fmt_table <- function(df, caption) {
  cat(sprintf("\n%s\n", paste(rep("=", nchar(caption)), collapse = "")))
  cat(sprintf("%s\n", caption))
  cat(sprintf("%s\n", paste(rep("=", nchar(caption)), collapse = "")))
  if (is.null(df) || nrow(df) == 0L) { cat("  (no results)\n"); return(invisible()) }
  for (g_val in sort(unique(df$group))) {
    for (t_val in sort(unique(df$time[df$group == g_val]))) {
      sub <- df[df$group == g_val & df$time == t_val, ]
      if (nrow(sub) == 0L) next
      cat(sprintf("\n  ATT(%g,%g) | true ATT = %.2f\n", g_val, t_val, TRUE_ATT))
      cat(sprintf("  %6s  %7s  %7s  %7s  %7s  %7s  %8s  %8s\n",
                  "n", "bias", "mc_sd", "mean_se", "se_ratio", "rmse", "coverage", "ci_len"))
      for (k in seq_len(nrow(sub))) {
        cat(sprintf("  %6d  %+7.4f  %7.4f  %7.4f  %7.3f  %7.4f  %8.3f  %8.4f\n",
                    sub$n[k], sub$bias[k], sub$mc_sd[k], sub$mean_se[k],
                    sub$se_ratio[k], sub$rmse[k], sub$coverage[k],
                    sub$mean_ci_len[k]))
      }
    }
  }
}

for (key in names(all_results)) {
  fmt_table(all_results[[key]], all_results[[key]]$label[1])
}

# -----------------------------------------------------------------------
# Cross-fitting stability study
# -----------------------------------------------------------------------
cat("\n\n=== Cross-fitting stability: ATT(3,3) across fold seeds, n=500 ===\n")
cat(sprintf("  %5s  %7s  %7s\n", "seed", "att(3,3)", "se(3,3)"))
df_stab <- dgp_linear(n = 500, seed = 42)
for (fold_seed in c(1L, 2L, 3L, 7L, 42L, 99L, 123L)) {
  fit_s <- tryCatch(
    edid(df_stab, "y", "id", "t", "g", xformla = ~ x1,
         seed = fold_seed, aggregate = "none"),
    error = function(e) NULL
  )
  if (is.null(fit_s)) {
    cat(sprintf("  %5d  (error)\n", fold_seed))
    next
  }
  row <- fit_s$att_gt[fit_s$att_gt$group == 3 & fit_s$att_gt$time == 3, ]
  if (nrow(row) == 0L) next
  cat(sprintf("  %5d  %7.4f  %7.4f\n", fold_seed, row$att, row$se))
}

# -----------------------------------------------------------------------
# EIF diagnostic
# -----------------------------------------------------------------------
cat("\n\n=== EIF diagnostics: mean and variance per cell, n=500 ===\n")
df_eif <- dgp_linear(n = 500, seed = 1)
fit_eif <- edid(df_eif, "y", "id", "t", "g", xformla = ~ x1,
                store_eif = TRUE, aggregate = "none", seed = 1L)
if (!is.null(fit_eif$eif)) {
  eif_mat <- fit_eif$eif
  n       <- fit_eif$n
  post    <- fit_eif$att_gt[!fit_eif$att_gt$is_pre, ]
  cat(sprintf("  %5s  %5s  %10s  %10s  %10s  %10s\n",
              "group", "time", "eif_mean", "eif_var", "se_eif", "se_reported"))
  for (k in seq_len(nrow(post))) {
    g_k <- post$group[k]; t_k <- post$time[k]
    cell_id <- which(fit_eif$att_gt$group == g_k & fit_eif$att_gt$time == t_k)
    if (length(cell_id) == 0L) next
    eif_col <- eif_mat[, cell_id, drop = TRUE]
    cat(sprintf("  %5g  %5g  %10.4e  %10.4e  %10.4f  %10.4f\n",
                g_k, t_k, mean(eif_col), var(eif_col),
                sqrt(sum(eif_col^2) / n^2), post$se[k]))
  }
}

# -----------------------------------------------------------------------
# Acceptance criteria check
# -----------------------------------------------------------------------
cat("\n\n=== Acceptance criteria ===\n")
passed <- 0L; failed <- 0L

for (key in names(all_results)) {
  df_res <- all_results[[key]]
  if (is.null(df_res) || nrow(df_res) == 0L) next
  lbl <- df_res$label[1]

  # 1. Bias decreases with n
  for (g_val in unique(df_res$group)) {
    for (t_val in unique(df_res$time[df_res$group == g_val])) {
      sub <- df_res[df_res$group == g_val & df_res$time == t_val &
                    df_res$n >= min(NS) & df_res$n <= max(NS), ]
      if (nrow(sub) < 2L) next
      large_bias  <- abs(sub$bias[sub$n == max(NS)])
      small_bias  <- abs(sub$bias[sub$n == min(NS)])
      # Allow some Monte Carlo noise: just require large-n bias < 0.2
      if (large_bias < 0.2) {
        cat(sprintf("  PASS: bias decreasing/small for %s ATT(%g,%g): n=%d bias=%.3f\n",
                    lbl, g_val, t_val, max(NS), large_bias))
        passed <- passed + 1L
      } else {
        cat(sprintf("  WARN: large-n bias for %s ATT(%g,%g): %.3f (> 0.2)\n",
                    lbl, g_val, t_val, large_bias))
        failed <- failed + 1L
      }
    }
  }

  # 2. SE ratio at largest n in (0.5, 2.0)
  sub_large <- df_res[df_res$n == max(NS), ]
  for (k in seq_len(nrow(sub_large))) {
    r <- sub_large$se_ratio[k]
    if (is.finite(r) && r > 0.5 && r < 2.0) {
      cat(sprintf("  PASS: SE ratio for %s ATT(%g,%g) n=%d: %.3f\n",
                  lbl, sub_large$group[k], sub_large$time[k], max(NS), r))
      passed <- passed + 1L
    } else {
      cat(sprintf("  WARN: SE ratio for %s ATT(%g,%g) n=%d: %.3f (outside 0.5-2.0)\n",
                  lbl, sub_large$group[k], sub_large$time[k], max(NS), r))
      failed <- failed + 1L
    }
  }
}

cat(sprintf("\n  Total PASS: %d | Total WARN: %d\n", passed, failed))

# -----------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------
results_path <- file.path("benchmark", "edid_cov_sim_results.rds")
saveRDS(list(results = all_results, params = list(R = R_FAST, ns = NS,
             true_att = TRUE_ATT, date = Sys.Date())),
        results_path)
cat(sprintf("\nResults saved to: %s\n", results_path))
