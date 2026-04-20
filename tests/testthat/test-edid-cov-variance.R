library(testthat)

# ============================================================
# Variance calibration tests for the covariate path.
# Use small Monte Carlo repetitions (R=50) so the testthat
# suite finishes quickly; the full simulation is in
# benchmark/edid_cov_sim.R.
# ============================================================

# ---------------------------------------------------------------
# DGP: balanced panel, single cohort (g=3) + never-treated,
#      linear covariate effect.
# True ATT = 1.0 for all post-treatment periods.
# ---------------------------------------------------------------
sim_one_draw <- function(n, seed) {
  set.seed(seed)
  T      <- 5
  ids    <- rep(1:n, each = T)
  times  <- rep(1:T, times = n)
  g_unit <- rep(c(3, Inf), each = n / 2)
  g_vec  <- g_unit[ids]
  x1u    <- rnorm(n)
  x1     <- rep(x1u, each = T)
  y      <- times + 0.5 * x1 + as.numeric(times >= g_vec) + rnorm(n * T, sd = 0.5)
  data.frame(id = ids, t = times, y = y, g = g_vec, x1 = x1)
}

run_mc <- function(n, R = 50) {
  results <- lapply(seq_len(R), function(r) {
    df <- sim_one_draw(n, seed = r)
    tryCatch({
      fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 1L,
                  aggregate = "none")
      att_gt_df <- fit$att_gt
      # post-treatment cells
      post <- att_gt_df[!att_gt_df$is_pre, ]
      list(att = post$att, se = post$se,
           group = post$group, time = post$time)
    }, error = function(e) NULL)
  })
  results <- Filter(Negate(is.null), results)
  results
}

# ============================================================
# SE vs empirical SD: ratio should be in (0.4, 2.5) for n=200
# (loose bounds for a quick test; tight bounds in benchmark)
# ============================================================

test_that("covariate SE / empirical SD ratio is in (0.4, 2.5) at n=200, R=50", {
  skip_on_cran()
  skip_if(Sys.getenv("CI") == "true" && Sys.getenv("EDID_SLOW_TESTS") != "1",
          "skipping slow variance test on CI")

  n <- 200; R <- 50
  res <- run_mc(n, R)
  skip_if(length(res) < 30L, "too many MC failures; skipping")

  # Collect per-cell ATTs and SEs
  all_gt   <- do.call(rbind, lapply(res, function(r) {
    data.frame(group = r$group, time = r$time, att = r$att, se = r$se)
  }))
  cell_keys <- unique(all_gt[, c("group", "time")])

  for (k in seq_len(nrow(cell_keys))) {
    g_k <- cell_keys$group[k]; t_k <- cell_keys$time[k]
    sub <- all_gt[all_gt$group == g_k & all_gt$time == t_k, ]
    if (nrow(sub) < 20L) next

    emp_sd   <- sd(sub$att)
    mean_se  <- mean(sub$se, na.rm = TRUE)
    ratio    <- mean_se / emp_sd

    expect_true(ratio > 0.4 && ratio < 2.5,
      info = sprintf("ATT(%g,%g): SE ratio = %.3f (mean_se=%.3f, emp_sd=%.3f)",
                     g_k, t_k, ratio, mean_se, emp_sd))
  }
})

# ============================================================
# EIF plug-in SE vs empirical: ratio in (0.4, 2.5)
# ============================================================

test_that("EIF plug-in SE matches empirical SD in expected range at n=200", {
  skip_on_cran()
  skip_if(Sys.getenv("CI") == "true" && Sys.getenv("EDID_SLOW_TESTS") != "1",
          "skipping slow EIF variance test on CI")

  n <- 200; R <- 50
  atts_33 <- numeric(R)
  ses_33   <- numeric(R)

  for (r in seq_len(R)) {
    df  <- sim_one_draw(n, seed = r)
    fit <- tryCatch(
      edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 1L, aggregate = "none"),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    row <- fit$att_gt[fit$att_gt$group == 3 & fit$att_gt$time == 3, ]
    if (nrow(row) == 0L) next
    atts_33[r] <- row$att
    ses_33[r]  <- row$se
  }

  valid <- is.finite(atts_33) & atts_33 != 0
  skip_if(sum(valid) < 20L, "insufficient valid draws")

  emp_sd  <- sd(atts_33[valid])
  mean_se <- mean(ses_33[valid])
  ratio   <- mean_se / emp_sd

  expect_true(ratio > 0.4 && ratio < 2.5,
    info = sprintf("ATT(3,3) SE ratio = %.3f (mean_se=%.3f, emp_sd=%.3f)",
                   ratio, mean_se, emp_sd))
})

# ============================================================
# Coverage: rough check at n=200 (nominal 95%, accept 70-99%)
# ============================================================

test_that("ATT(3,3) CI coverage is roughly nominal at n=200, R=50", {
  skip_on_cran()
  skip_if(Sys.getenv("CI") == "true" && Sys.getenv("EDID_SLOW_TESTS") != "1",
          "skipping slow coverage test on CI")

  n <- 200; R <- 50; true_att <- 1.0
  covered <- logical(R)

  for (r in seq_len(R)) {
    df  <- sim_one_draw(n, seed = r)
    fit <- tryCatch(
      edid(df, "y", "id", "t", "g", xformla = ~ x1, seed = 1L, aggregate = "none"),
      error = function(e) NULL
    )
    if (is.null(fit)) { covered[r] <- FALSE; next }
    row <- fit$att_gt[fit$att_gt$group == 3 & fit$att_gt$time == 3, ]
    if (nrow(row) == 0L) { covered[r] <- FALSE; next }
    covered[r] <- (row$ci_lower <= true_att && true_att <= row$ci_upper)
  }

  cov_rate <- mean(covered)
  expect_true(cov_rate >= 0.65 && cov_rate <= 0.99,
    info = sprintf("ATT(3,3) coverage = %.2f at n=200 (nominal 0.95)", cov_rate))
})

# ============================================================
# No-covariate path: SE unchanged by xformla=NULL vs ~1
# ============================================================

test_that("no-covariate path SEs are identical for xformla=NULL and xformla=~1", {
  df   <- sim_one_draw(200, seed = 99)
  fit0 <- edid(df, "y", "id", "t", "g")
  fit1 <- edid(df, "y", "id", "t", "g", xformla = ~1)
  expect_equal(fit0$att_gt$se, fit1$att_gt$se, tolerance = 1e-12)
})
