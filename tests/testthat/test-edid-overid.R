# Tests for edid_overid (omnibus over-identification J).
# Structural + regression + honesty-safeguard checks. The size/power calibration lives in the MC
# study (quality_reports/overid_J_mc/); these are fast, deterministic structural guards.

make_stg <- function(n_per = 120, cohorts = c(3, 4, 5, Inf), Tt = 5, with_x = FALSE, seed = 20260619) {
  set.seed(seed)
  g <- sample(cohorts, n_per * length(cohorts), replace = TRUE); N <- length(g)
  do.call(rbind, lapply(seq_len(N), function(i) {
    gi <- g[i]; u <- rnorm(1); x <- rnorm(1)
    data.frame(id = i, time = 1:Tt, g = gi, x = x,
               y = u + (if (with_x) 0.3 * x else 0) + 0.2 * (1:Tt) + 1 * ((1:Tt) >= gi) + rnorm(Tt, 0, 0.5))
  }))
}

test_that("edid_overid returns a well-formed object with sensible statistics", {
  df  <- make_stg()
  fit <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
              aggregate = "event_study", cband = FALSE)
  ov  <- suppressWarnings(suppressMessages(edid_overid(fit, data = df)))
  expect_s3_class(ov, "edid_overid")
  expect_true(is.data.frame(ov$table))
  expect_true(all(c("parameter", "J_statistic", "df", "p_value", "n_moments",
                    "n_params", "nominal_df") %in% names(ov$table)))
  ov_row <- ov$table[ov$table$parameter == "overall", , drop = FALSE]
  expect_equal(nrow(ov_row), 1L)
  expect_true(is.finite(ov_row$J_statistic) && ov_row$J_statistic >= 0)
  expect_true(ov_row$df >= 1L)
  expect_true(ov_row$p_value >= 0 && ov_row$p_value <= 1)
  # engine df is the detected rank, never above the naive moment count Q - p
  expect_true(ov_row$df <= ov_row$nominal_df)
  # per-cell breakdown present
  expect_true(is.data.frame(ov$cells) && nrow(ov$cells) >= 1L)
})

test_that("edid_overid is invariant to the fit's reporting options (over-id refits in bare plug-in)", {
  df <- make_stg()
  f_eff <- edid(df, "y", "id", "time", "g", weight_scheme = "efficient",
                pt_assumption = "all", aggregate = "event_study", cband = FALSE)
  f_avg <- edid(df, "y", "id", "time", "g", weight_scheme = "averaged",
                pt_assumption = "all", aggregate = "event_study", cband = FALSE)
  j_eff <- suppressWarnings(suppressMessages(edid_overid(f_eff, data = df)))$table
  j_avg <- suppressWarnings(suppressMessages(edid_overid(f_avg, data = df)))$table
  je <- j_eff$J_statistic[j_eff$parameter == "overall"]
  ja <- j_avg$J_statistic[j_avg$parameter == "overall"]
  expect_equal(je, ja, tolerance = 1e-6)
})

test_that("rel_tol is a no-op for the no-covariate spectrum and reduces the covariate rank", {
  # No-covariate: the adaptive rel_tol = "auto" = r_bare/n_eff must NOT drop any genuine direction
  # (exact zeros after the true rank), so it reproduces the bare numerical-rank cut (rel_tol = 0).
  d0 <- make_stg(with_x = FALSE)
  f0 <- edid(d0, "y", "id", "time", "g", pt_assumption = "all",
             aggregate = "event_study", cband = FALSE)
  o_adapt <- suppressWarnings(suppressMessages(edid_overid(f0, data = d0, parameter = "overall")))
  o_bare  <- suppressWarnings(suppressMessages(edid_overid(f0, data = d0, parameter = "overall", rel_tol = 0)))
  expect_equal(o_adapt$table$df, o_bare$table$df)
  expect_equal(o_adapt$table$J_statistic, o_bare$table$J_statistic, tolerance = 1e-8)

  # Covariate: the decaying spectrum makes the bare cut over-count; the adaptive floor must reduce the
  # detected rank below the bare-cut rank (recovering the effective over-id dimension).
  dX <- make_stg(with_x = TRUE)
  fX <- edid(dX, "y", "id", "time", "g", xformla = ~ x, pt_assumption = "all",
             aggregate = "event_study", cband = FALSE)
  oX_adapt <- suppressWarnings(suppressMessages(edid_overid(fX, data = dX, parameter = "overall")))
  oX_bare  <- suppressWarnings(suppressMessages(edid_overid(fX, data = dX, parameter = "overall", rel_tol = 0)))
  expect_true(oX_adapt$table$df <= oX_bare$table$df)
  expect_true(is.finite(oX_adapt$table$J_statistic) && oX_adapt$table$df >= 1L)
  # the legacy "not yet calibrated" warning is gone (covariate path is now validated)
  ws <- character(0)
  withCallingHandlers(
    suppressMessages(edid_overid(fX, data = dX, parameter = "overall")),
    warning = function(w) { ws[[length(ws) + 1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("not yet calibrated", ws, ignore.case = TRUE)))
})

test_that("edid_overid errors informatively on a non-edid object", {
  expect_error(edid_overid(list(1, 2, 3)), "edid_fit")
})

test_that("extreme over-id rank-deficiency (q >= n_eff) returns NA, not a misleading p = 1", {
  # The engine: rel_tol floor that drops every direction => NA (uncomputable), rank_deficient = TRUE.
  set.seed(11); n <- 200; xi <- matrix(rnorm(n * 8), n, 8); d <- rnorm(8) * 0.1
  qf <- .edid_if_diff_quadform(d, xi, n, NULL, rel_tol = 5)   # rel_tol > 1 drops all directions
  expect_true(is.na(qf$statistic)); expect_true(is.na(qf$p_value)); expect_true(isTRUE(qf$rank_deficient))
  # Genuine coincidence (D ~ 0) still reports p = 1, NOT NA.
  qf0 <- .edid_if_diff_quadform(rep(0, 4), matrix(1e-12 * rnorm(n * 4), n, 4), n, NULL, rel_tol = 0.01)
  expect_equal(qf0$p_value, 1)
  # rel_tol = 0 (edid_hausman/edid_sargan convention) never hits the NA branch.
  qfb <- .edid_if_diff_quadform(rep(0, 4), matrix(1e-12 * rnorm(n * 4), n, 4), n, NULL, rel_tol = 0)
  expect_equal(qfb$p_value, 1); expect_false(isTRUE(qfb$rank_deficient))
})

test_that("cluster-rank SATURATION (r_bare = G_eff-1) returns NA on the auto path, never under rel_tol=0", {
  # The Bailey-GB mechanism: a cluster-robust D-hat whose numerical rank reaches the cluster ceiling
  # (G_eff - 1) is full-rank-for-its-budget -> no null space -> the joint over-id is uncomputable.
  set.seed(21); G <- 6L; per <- 12L; n <- G * per
  ci <- rep(seq_len(G), each = per)                 # G_eff = 6 clusters
  xi <- matrix(rnorm(n * 9), n, 9)                  # cluster-robust D-hat has rank <= G_eff-1 = 5
  d  <- rnorm(9) * 0.1
  # auto path: r_bare reaches G_eff-1 -> saturation -> NA (rank_deficient), NOT a misleading number.
  qa <- .edid_if_diff_quadform(d, xi, n, ci, rel_tol = "auto")
  expect_true(is.na(qa$statistic)); expect_true(is.na(qa$p_value)); expect_true(isTRUE(qa$rank_deficient))
  # rel_tol = 0 (edid_hausman/edid_sargan) is UNAFFECTED by the saturation guard (gated on auto): it
  # computes its usual rank-aware statistic, never NA-by-saturation.
  qb <- .edid_if_diff_quadform(d, xi, n, ci, rel_tol = 0)
  expect_false(isTRUE(qb$rank_deficient)); expect_true(is.finite(qb$statistic))
})
