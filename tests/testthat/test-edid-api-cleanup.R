library(testthat)

# Tests for the Tier 2 (API) and Tier 3 (correctness) edid changes:
#   - weight_scheme rename; control_group removed; balance_e honored; multi-value aggregate
#   - aggte_edid() bootstrap reproducibility (seed); analytic t_stat/p_value consistency
#   - degenerate-data guards (single cluster, single-unit cohort)
#   - covariate-path SE == EIF plug-in identity across ALL weight schemes (coverage gap)

# Covariate panel with cohorts {never, 3, 5}, a single covariate x1, plenty of units per cohort.
make_cov_panel <- function(n = 400L, seed = 11L) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  P  <- exp(cbind(0, 0.4 * x1, 0.2 * x1)); P <- P / rowSums(P)
  g  <- apply(P, 1, function(pr) sample(c(Inf, 3, 5), 1, prob = pr))
  alpha <- rnorm(n, 0.3 * x1, 1)
  do.call(rbind, lapply(1:6, function(t) {
    tau <- 1 * (is.finite(g) & t >= g)
    data.frame(id = 1:n, t = t, g = ifelse(is.finite(g), g, 0),
               x1 = x1, y = alpha + 0.3 * t + 0.4 * x1 * (t - 1) + tau + rnorm(n))
  }))
}

# ---------------------------------------------------------------------------
# Tier 2: weight_scheme rename + control_group removal
# ---------------------------------------------------------------------------
test_that("weight_scheme is the weighting argument; the old `weights` arg is gone", {
  df <- make_panel_1cohort()
  expect_true("weight_scheme" %in% names(formals(edid)))
  expect_false("weights" %in% names(formals(edid)))
  fit <- edid(df, "outcome", "unit", "time", "first_treat", weight_scheme = "averaged",
              aggregate = "none")
  expect_s3_class(fit, "edid_fit")
  expect_error(
    edid(df, "outcome", "unit", "time", "first_treat", weights = "efficient"),
    "unused argument"
  )
})

test_that("control_group is removed from edid() and the family always uses never-treated", {
  df <- make_panel_1cohort()
  expect_false("control_group" %in% names(formals(edid)))
  expect_false("control_group" %in% names(formals(prepare_edid_panel)))
  expect_error(
    edid(df, "outcome", "unit", "time", "first_treat", control_group = "nevertreated"),
    "unused argument"
  )
  fit <- edid(df, "outcome", "unit", "time", "first_treat", aggregate = "none")
  expect_false("control_group" %in% names(fit))
})

# ---------------------------------------------------------------------------
# Tier 2: multi-value aggregate no longer crashes; only requested slots populated
# ---------------------------------------------------------------------------
test_that("aggregate accepts a vector of more than one type without error", {
  df <- make_panel_2cohort()
  fit <- expect_no_error(
    edid(df, "outcome", "unit", "time", "first_treat", aggregate = c("group", "calendar"))
  )
  expect_false(is.null(fit$group))
  expect_false(is.null(fit$calendar))
  expect_null(fit$event_study)
  expect_null(fit$simple)
})

# ---------------------------------------------------------------------------
# Tier 2: balance_e is honored (forwarded to the dynamic aggregation), not a no-op
# ---------------------------------------------------------------------------
test_that("balance_e / max_e restrict the dynamic aggregation", {
  df  <- make_panel_2cohort()
  fit <- edid(df, "outcome", "unit", "time", "first_treat", aggregate = "event_study")
  full <- fit$event_study
  expect_gt(max(full$egt), 1)                               # unrestricted spans beyond e = 1

  # via aggte_edid(max_e=): restriction takes effect
  restricted <- aggte_edid(fit, type = "dynamic", max_e = 1, na.rm = TRUE)
  expect_lte(max(restricted$egt), 1)
  expect_lt(length(restricted$egt), length(full$egt))

  # via edid(balance_e=): forwarded into .agg (was previously silently ignored)
  fit_b <- edid(df, "outcome", "unit", "time", "first_treat",
                aggregate = "event_study", balance_e = 1)
  expect_lte(max(fit_b$event_study$egt), 1)
  expect_false(identical(sort(unique(fit_b$event_study$egt)),
                         sort(unique(full$egt))))
})

# ---------------------------------------------------------------------------
# Tier 3: aggte_edid() multiplier bootstrap is reproducible (seed)
# ---------------------------------------------------------------------------
test_that("aggte_edid() bootstrap is reproducible via the fit's seed", {
  df  <- make_panel_2cohort()
  fit <- edid(df, "outcome", "unit", "time", "first_treat",
              aggregate = "none", bstrap = TRUE, biters = 199L,
              cband_method = "multiplier", seed = 5L)
  a1 <- aggte_edid(fit, type = "dynamic", na.rm = TRUE)
  a2 <- aggte_edid(fit, type = "dynamic", na.rm = TRUE)
  expect_equal(a1$se.egt,        a2$se.egt)
  expect_equal(a1$crit.val.egt,  a2$crit.val.egt)
  # explicit seed override is also reproducible
  b1 <- aggte_edid(fit, type = "dynamic", na.rm = TRUE, seed = 999L)
  b2 <- aggte_edid(fit, type = "dynamic", na.rm = TRUE, seed = 999L)
  expect_equal(b1$se.egt, b2$se.egt)
})

# ---------------------------------------------------------------------------
# Tier 3: analytic t_stat / p_value are consistent with the reported SE
# ---------------------------------------------------------------------------
test_that("analytic-path t_stat and p_value match the reported SE (default and higher_order)", {
  df  <- make_cov_panel()
  fit <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none", seed = 1L)
  ok  <- is.finite(fit$att_gt$se) & fit$att_gt$se > 0
  expect_equal(fit$att_gt$t_stat[ok], (fit$att_gt$att / fit$att_gt$se)[ok])
  expect_equal(fit$att_gt$p_value[ok],
               2 * stats::pnorm(-abs(fit$att_gt$t_stat[ok])))

  fitH <- edid(df, "y", "id", "t", "g", xformla = ~ x1, aggregate = "none",
               higher_order = TRUE, seed = 1L)
  okH  <- is.finite(fitH$att_gt$se) & fitH$att_gt$se > 0
  # higher_order inflates SE; t/p must reflect the inflated SE, not the plug-in SE
  expect_equal(fitH$att_gt$t_stat[okH], (fitH$att_gt$att / fitH$att_gt$se)[okH])
  expect_true(all(fitH$att_gt$se[okH] >= fit$att_gt$se[okH] - 1e-10))
})

# ---------------------------------------------------------------------------
# Tier 3: degenerate-data guards warn loudly
# ---------------------------------------------------------------------------
test_that("a single distinct cluster warns that cluster-robust SEs are undefined", {
  df <- make_panel_1cohort()
  df$cl <- 1L                                               # one cluster for everyone
  expect_warning(
    edid(df, "outcome", "unit", "time", "first_treat", clustervars = "cl", aggregate = "none"),
    "only 1 distinct cluster"
  )
})

test_that("cohorts with fewer than 2 units warn that SEs are degenerate", {
  df <- make_degenerate_panel()                            # 1 treated + 1 never-treated unit
  expect_warning(
    edid(df, "outcome", "unit", "time", "first_treat", aggregate = "none"),
    "fewer than 2 units|never-treated unit"
  )
})

# ---------------------------------------------------------------------------
# Tier 4 coverage: covariate-path reported SE == EIF plug-in identity, ALL schemes
# ---------------------------------------------------------------------------
test_that("covariate-path SE equals sqrt(colSums(eif^2)/n^2) for every weight scheme", {
  df <- make_cov_panel()
  for (ws in c("efficient", "averaged", "gmm", "uniform")) {
    fit <- suppressWarnings(
      edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = ws,
           aggregate = "none", seed = 1L, misspec_robust = FALSE)   # plug-in EIF-SE identity
    )
    ok  <- is.finite(fit$att_gt$se) & fit$att_gt$se > 0
    eif_se <- sqrt(colSums(fit$eif^2) / fit$n^2)
    expect_equal(fit$att_gt$se[ok], eif_se[ok], tolerance = 1e-8,
                 info = paste("weight_scheme =", ws))
  }
})
