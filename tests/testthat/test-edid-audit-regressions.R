library(testthat)

# ===========================================================================
# Audited-fix regression batch (Phase 2 audit program). Each test pins a fix
# from the audit sessions; the fixed file / behavior is cited in a comment so
# a future failure points straight at the regressed change.
# ===========================================================================

# Collect every warning a call emits (so tests can assert exact warning sets
# without testthat swallowing or re-signalling extras).
.collect_warnings <- function(expr) {
  ws <- character(0L)
  val <- withCallingHandlers(expr,
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = val, warnings = ws)
}

# Covariate panel used by several tests below (mild covariate effect: the
# default fit and the bs_df = "ic" fit are warning-free on it).
make_panel_cov_audit <- function(seed = 11L, n = 200L, Tt = 4L) {
  set.seed(seed)
  coh <- sample(c(3, Inf), n, replace = TRUE, prob = c(.45, .55))
  x1  <- rnorm(n)
  df  <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  df$g <- coh[df$id]; df$x1 <- x1[df$id]
  ufe <- rnorm(n)
  df$y <- ufe[df$id] + 0.3 * df$time + 0.3 * df$x1 +
    rnorm(nrow(df), 0, 0.5) + 1 * (df$time >= df$g)
  df
}

# ---------------------------------------------------------------------------
# 1. clustervars coinciding with a design column, and reserved column names
#    Fix: edid-validate.R -- clustervars %in% c(yname, idname, tname, gname)
#    now errors (clustering on a design column corrupted as_MP_edid's per-unit
#    frame), and ".w"/".edid_cluster" are reserved internal names.
# ---------------------------------------------------------------------------
test_that("clustervars == gname/yname/idname rejected; reserved column names rejected", {
  df <- make_panel_1cohort()

  for (cv in c("first_treat", "outcome", "unit")) {
    expect_error(
      edid(df, "outcome", "unit", "time", "first_treat", clustervars = cv),
      "coincides with yname/idname/tname/gname")
  }

  # reserved internal names: as_MP_edid() writes .w (sampling weight) and
  # .edid_cluster (cluster codes) into its per-unit frame
  df_w <- df; df_w$.w <- df_w$outcome
  expect_error(edid(df_w, ".w", "unit", "time", "first_treat"),
               "reserved for internal use")
  df_c <- df; df_c$.edid_cluster <- df_c$unit
  expect_error(edid(df_c, "outcome", ".edid_cluster", "time", "first_treat"),
               "reserved for internal use")
})

# ---------------------------------------------------------------------------
# 2. The C1 corruption scenario: cohort-level clustering via a duplicate of
#    gname under a different name.
#    Fix: as_MP_edid() (edid-mp.R) stores the EIF-aligned cluster codes under
#    the reserved ".edid_cluster" column instead of overwriting the caller's
#    column in the per-unit frame -- previously a cluster column duplicating
#    gname overwrote the cohort column, silently corrupting compute.aggte()'s
#    group shares (wrong POINT estimates downstream); it also now sets
#    DIDparams$nG/nT/est_method so broom::glance() works on the MP.
# ---------------------------------------------------------------------------
test_that("cohort-level clustering changes only SEs; MP aggregation and glance() are sound", {
  # UNEQUAL cohort sizes (30 / 90 / 60) so a corrupted group share cannot cancel
  set.seed(202606)
  n_g3 <- 30L; n_g5 <- 90L; n_nv <- 60L; Tt <- 6L
  n  <- n_g3 + n_g5 + n_nv
  df <- data.frame(id = rep(seq_len(n), each = Tt), time = rep(seq_len(Tt), n))
  g_u <- c(rep(3, n_g3), rep(5, n_g5), rep(Inf, n_nv))
  df$g <- g_u[df$id]
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
    rnorm(nrow(df), 0, 0.5)
  # legit duplicate-of-gname cluster column under a DIFFERENT name (finite codes)
  df$cohort_cl <- ifelse(is.finite(df$g), df$g, 999)

  f0 <- edid(df, "y", "id", "time", "g", aggregate = "all", cband = FALSE)
  fc <- edid(df, "y", "id", "time", "g", aggregate = "all", cband = FALSE,
             clustervars = "cohort_cl")

  # (i) clustering must change ONLY the standard errors
  expect_identical(fc$att_gt$att, f0$att_gt$att)
  expect_identical(fc$event_study$att.egt, f0$event_study$att.egt)
  expect_identical(fc$group$att.egt,       f0$group$att.egt)
  expect_identical(fc$overall$overall.att, f0$overall$overall.att)
  expect_false(identical(fc$att_gt$se, f0$att_gt$se))   # the cluster path did engage

  # (ii) glance() works on the edid-backed MP (needs DIDparams$nG/nT/est_method)
  gl <- generics::glance(as_MP_edid(fc))
  expect_s3_class(gl, "data.frame")
  expect_identical(gl$nobs, n)
  expect_identical(gl$ngroup, 2L)
  expect_identical(gl$ntime, Tt)
  expect_identical(gl$est.method, "edid")

  # (iii) aggte through the clustered MP reproduces the unclustered att.egt
  a_mp <- aggte(as_MP_edid(fc, bstrap = FALSE, cband = FALSE),
                type = "dynamic", na.rm = TRUE, bstrap = FALSE)
  expect_equal(a_mp$att.egt, f0$event_study$att.egt, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. NA unit id
#    Fix: edid-validate.R check 5b -- an NA id used to slip through the balance
#    arithmetic and silently drop the unit in prepare_edid_panel (or trip a
#    misleading duplicate-rows error); it is now rejected explicitly.
# ---------------------------------------------------------------------------
test_that("NA unit id gives an informative error", {
  df <- make_panel_1cohort()
  df$unit[1] <- NA
  expect_error(edid(df, "outcome", "unit", "time", "first_treat"),
               "idname.*contains NA values")
})

# ---------------------------------------------------------------------------
# 4. Full overlap trim: every treated unit removed
#    Fix: fit_edid_cells (edid-fit.R) -- a trim_level at/below the smallest
#    inverse propensity used to return a confident-looking exact att = 0 with
#    an NA SE; the cell is now NA with a single explanatory warning.
# ---------------------------------------------------------------------------
test_that("trim_level = 1.0001 yields NA cells plus the full-trim warning (not a silent att = 0)", {
  df <- make_panel_cov_audit()
  res <- .collect_warnings(
    edid(df, "y", "id", "time", "g", xformla = ~ x1,
         trim_level = 1.0001, aggregate = "none", cband = FALSE))
  fit <- res$value
  post <- !fit$att_gt$is_pre
  expect_true(any(post))
  expect_true(all(is.na(fit$att_gt$att[post])))         # NA, not 0
  expect_true(any(grepl("Overlap trimming removed every treated unit", res$warnings)))
})

# ---------------------------------------------------------------------------
# 5. balance_e beyond the feasible event-time span
#    Fix: aggte_edid (edid-aggte.R) -- previously an opaque subscript error out
#    of compute.aggte; now a feasibility error naming balance_e and the bound.
# ---------------------------------------------------------------------------
test_that("balance_e beyond the feasible span errors informatively", {
  df <- make_panel_2cohort()   # cohorts 3 and 5, periods 1..7 -> max feasible e = 4
  expect_error(
    edid(df, "outcome", "unit", "time", "first_treat",
         aggregate = "event_study", balance_e = 10, cband = FALSE),
    "balance_e")
})

# ---------------------------------------------------------------------------
# 6. Band-label truth
#    Fixes: summary.AGGTEobj (AGGTEobj.R) now labels from the effective
#    DIDparams$cband (edid's analytic sup-t path sets it TRUE when it installs
#    a simultaneous crit; the old bstrap && cband rule printed "Pointwise" for
#    those bands), and print.edid_fit (edid-methods.R) labels "Pointwise" for
#    bstrap = TRUE with cband = FALSE.
# ---------------------------------------------------------------------------
test_that("analytic simultaneous bands print 'Simult.'; bootstrap pointwise prints 'Pointwise'", {
  df <- make_panel_1cohort()

  f_an <- edid(df, "outcome", "unit", "time", "first_treat",
               aggregate = "event_study", cband = TRUE)   # default analytic sup-t
  out_an <- capture.output(summary(f_an$event_study))
  expect_true(any(grepl("Simult.", out_an, fixed = TRUE)))
  expect_false(any(grepl("Pointwise", out_an, fixed = TRUE)))

  f_pw <- edid(df, "outcome", "unit", "time", "first_treat",
               aggregate = "none", bstrap = TRUE, biters = 60L,
               cband = FALSE, seed = 5L)
  out_pw <- capture.output(print(f_pw))
  expect_true(any(grepl("Pointwise", out_pw, fixed = TRUE)))
  expect_false(any(grepl("Simult.", out_pw, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 7. Clustered consistency of the aggregate SEs
#    Fix: .edid_analytic_cband_agg (edid-aggte.R) -- the aggregate SEs now carry
#    the same G/(G-1) finite-cluster factor as the cell-level SEs and vcov()
#    (did's getSE() does not apply it), so the two conventions agree exactly.
# ---------------------------------------------------------------------------
test_that("clustered event-study SEs equal sqrt(diag(vcov())) (G/(G-1) alignment)", {
  df  <- make_panel_clustered()
  fit <- edid(df, "outcome", "unit", "time", "first_treat",
              clustervars = "cluster_id", aggregate = "event_study",
              cband = FALSE)
  v  <- vcov(fit, "event_study")
  expect_equal(unname(sqrt(diag(v))), unname(fit$event_study$se.egt),
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 8. PT-Post H = 1 invariance across weight schemes
#    Fix: compute_generated_outcomes_cov_edid (edid-cov-eif.R) routes the
#    PT-Post single moment to the self/two-period branch (base period g-1, not
#    period_1), and the no-covariate path honors weight_method. Under PT-Post
#    each cell has exactly one pair, so all four schemes must coincide.
# ---------------------------------------------------------------------------
test_that("pt_assumption = 'post': all four weight schemes give identical att (H = 1)", {
  df <- make_panel_cov_audit(seed = 42L)
  atts <- lapply(c("efficient", "averaged", "gmm", "uniform"), function(sch) {
    # gmm legitimately warns about its two-step bias on the covariate path
    suppressWarnings(
      edid(df, "y", "id", "time", "g", xformla = ~ x1, pt_assumption = "post",
           weight_scheme = sch, aggregate = "none", cband = FALSE)$att_gt$att)
  })
  for (a in atts[-1]) expect_equal(a, atts[[1]], tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 9. anticipation = 1 with covariates
#    Fix: the covariate nuisance/m-cache build under anticipation (edid-fit.R /
#    edid-pairs.R effective-onset handling): the run must complete cleanly and
#    the e = -1 (anticipation-period) cells must be estimated.
# ---------------------------------------------------------------------------
test_that("anticipation = 1 with xformla runs clean and e = -1 cells exist", {
  df  <- make_panel_cov_audit(seed = 77L, Tt = 5L)
  res <- .collect_warnings(
    edid(df, "y", "id", "time", "g", xformla = ~ x1, anticipation = 1L,
         aggregate = "event_study", cband = FALSE))
  expect_identical(res$warnings, character(0L))
  fit <- res$value
  expect_true(-1 %in% fit$event_study$egt)
  expect_true(is.finite(fit$event_study$att.egt[fit$event_study$egt == -1]))
  # the (g, g-1) anticipation cell itself is estimated
  expect_true(any(fit$att_gt$group == 3 & fit$att_gt$time == 2 &
                    is.finite(fit$att_gt$att)))
})

# ---------------------------------------------------------------------------
# 10. alp threading into the pointwise CIs
#     Fix: alp is threaded through fit_edid_cells / analytic_bands_edid, so the
#     pointwise CI width scales exactly by the normal quantile ratio.
# ---------------------------------------------------------------------------
test_that("alp = 0.10 pointwise CI width ratio equals qnorm(.95)/qnorm(.975)", {
  df  <- make_panel_2cohort()
  f05 <- edid(df, "outcome", "unit", "time", "first_treat",
              alp = 0.05, cband = FALSE, aggregate = "none")
  f10 <- edid(df, "outcome", "unit", "time", "first_treat",
              alp = 0.10, cband = FALSE, aggregate = "none")
  ok  <- is.finite(f05$att_gt$se) & f05$att_gt$se > 0
  expect_true(any(ok))
  w05 <- (f05$att_gt$ci_upper - f05$att_gt$ci_lower)[ok]
  w10 <- (f10$att_gt$ci_upper - f10$att_gt$ci_lower)[ok]
  expect_equal(w10 / w05,
               rep(qnorm(0.95) / qnorm(0.975), sum(ok)), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 11. aggregate = "calendar" standalone leaves $overall NULL
#     Fix: coef.edid_fit / as.data.frame.edid_fit (edid-methods.R) are NULL-safe
#     for $overall (shape-stable empty returns instead of an error).
# ---------------------------------------------------------------------------
test_that("calendar-only aggregation: NULL $overall is safe in coef() and as.data.frame()", {
  df  <- make_panel_2cohort()
  fit <- edid(df, "outcome", "unit", "time", "first_treat",
              aggregate = "calendar", cband = FALSE)
  expect_null(fit$overall)
  expect_s3_class(fit$calendar, "AGGTEobj")

  co <- coef(fit, "overall")
  expect_identical(co, numeric(0L))

  dd <- as.data.frame(fit, which = "overall")
  expect_s3_class(dd, "data.frame")
  expect_identical(nrow(dd), 0L)
  expect_true(all(c("att", "se", "ci_lower", "ci_upper") %in% names(dd)))
})

# ---------------------------------------------------------------------------
# 12. PT-Post baseline on a non-integer time grid
#     Fix: enumerate_valid_pairs_edid (edid-pairs.R) uses the strict
#     "last observed period < g - anticipation" rule, so periods {1, 1.5, 2, 3}
#     with g = 2 take baseline 1.5 (the "<= g - 1 - anticipation" arithmetic
#     skipped back to 1).
# ---------------------------------------------------------------------------
test_that("PT-Post picks tpre = 1.5 on the grid {1, 1.5, 2, 3} with g = 2", {
  # direct enumeration
  pr <- enumerate_valid_pairs_edid(
    target_g = 2, treatment_groups = 2, time_periods = c(1, 1.5, 2, 3),
    period_1 = 1, pt_assumption = "post", anticipation = 0L)
  expect_identical(nrow(pr), 1L)
  expect_identical(pr$gp, Inf)
  expect_identical(pr$tpre, 1.5)

  # and through a full fit (the stored per-cell pairs carry the same key)
  set.seed(8)
  n  <- 80L; per <- c(1, 1.5, 2, 3)
  df <- data.frame(id = rep(seq_len(n), each = length(per)),
                   time = rep(per, n))
  df$g <- rep(sample(c(2, Inf), n, replace = TRUE), each = length(per))
  df$y <- rnorm(n)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
    rnorm(nrow(df), 0, 0.5)
  fit <- edid(df, "y", "id", "time", "g", pt_assumption = "post",
              aggregate = "none", cband = FALSE)
  cell <- Filter(function(cc) cc$group == 2 && cc$time == 2, fit$cells)[[1]]
  expect_identical(cell$pairs$tpre, 1.5)
  expect_identical(cell$pairs$gp, Inf)
  expect_true(is.finite(cell$att))
})

# ---------------------------------------------------------------------------
# 13. bs_df argument (this phase): default byte-identity, "ic" selection,
#     and rejection of df < 3.
#     Feature: edid(bs_df =) threads through fit_edid_cells to the sieve
#     nuisance estimators (edid-cov.R); "ic" selects per fit by the paper's
#     2*E_n[loss] + log(n)*K/n criterion over df 3:8 and stores the selected
#     dimensions on the fit.
# ---------------------------------------------------------------------------
test_that("bs_df: default is byte-identical to 4L; 'ic' runs clean and records selections; 2 errors", {
  df <- make_panel_cov_audit()

  f_def <- edid(df, "y", "id", "time", "g", xformla = ~ x1, seed = 1L)
  f_4   <- edid(df, "y", "id", "time", "g", xformla = ~ x1, seed = 1L, bs_df = 4L)
  expect_identical(f_def$att_gt, f_4$att_gt)            # att, se, ci, t, p -- all of it
  expect_identical(unname(f_def$eif), unname(f_4$eif))
  expect_null(f_def$bs_df_selected)

  res <- .collect_warnings(
    edid(df, "y", "id", "time", "g", xformla = ~ x1, seed = 1L, bs_df = "ic"))
  expect_identical(res$warnings, character(0L))
  f_ic <- res$value
  sel  <- f_ic$bs_df_selected
  expect_s3_class(sel, "data.frame")
  expect_true(all(c("g", "nuisance", "key", "bs_df") %in% names(sel)))
  expect_true(all(sel$nuisance %in% c("r", "s", "m")))
  expect_true(all(sel$bs_df %in% 3:8))
  expect_true(all(c("r", "s", "m") %in% sel$nuisance))
  # the downstream variance channels (misspec_robust default) produced finite SEs
  expect_true(all(is.finite(f_ic$att_gt$se[!f_ic$att_gt$is_pre])))

  expect_error(edid(df, "y", "id", "time", "g", bs_df = 2),
               "bs_df.*integer >= 3")
  expect_error(edid(df, "y", "id", "time", "g", bs_df = "aic"),
               "bs_df.*integer >= 3")
})

# ---------------------------------------------------------------------------
# 14. edid_weights() / edid_weight_plot() (this phase): the paper's
#     weight-decomposition diagnostic exposed as a tidy accessor + heatmap.
#     The stored cell weights are labeled by their (g', t_pre) pair and the
#     pair keys ride on the cell ($pairs).
# ---------------------------------------------------------------------------
test_that("edid_weights returns labeled tidy weights; uniform weights are equal; plot is a ggplot", {
  df  <- make_panel_2cohort()
  fit <- edid(df, "outcome", "unit", "time", "first_treat",
              aggregate = "none", cband = FALSE)

  w <- edid_weights(fit)
  expect_s3_class(w, "data.frame")
  expect_identical(names(w),
                   c("group", "time", "gp", "tpre", "weight", "n_pairs", "condition_num"))
  expect_s3_class(attr(w, "na_cells"), "data.frame")

  # keys match the enumeration for every cell, in enumeration order
  for (key in unique(paste(w$group, w$time))) {
    wk <- w[paste(w$group, w$time) == key, ]
    pr <- enumerate_valid_pairs_edid(
      target_g = wk$group[1], treatment_groups = fit$treatment_groups,
      time_periods = fit$time_periods, period_1 = min(fit$time_periods),
      pt_assumption = "all", anticipation = 0L)
    expect_identical(wk$gp, pr$gp)
    expect_identical(wk$tpre, pr$tpre)
    expect_identical(wk$n_pairs[1], nrow(pr))
    # weights sum to ~1 within the cell
    expect_equal(sum(wk$weight), 1, tolerance = 1e-8)
  }

  # the stored vector itself carries the "gp=..,tpre=.." names
  cell1 <- Filter(function(cc) !is.null(cc$weights), fit$cells)[[1]]
  expect_identical(names(cell1$weights),
                   paste0("gp=", cell1$pairs$gp, ",tpre=", cell1$pairs$tpre))

  # uniform scheme: all weights equal 1/n_pairs
  f_u <- edid(df, "outcome", "unit", "time", "first_treat",
              weight_scheme = "uniform", aggregate = "none", cband = FALSE)
  wu <- edid_weights(f_u)
  expect_equal(wu$weight, 1 / wu$n_pairs, tolerance = 1e-12)

  # heatmap
  p <- edid_weight_plot(fit)
  expect_s3_class(p, "ggplot")
  p1 <- edid_weight_plot(fit, cells = data.frame(group = 3, time = 4))
  expect_s3_class(p1, "ggplot")
  expect_error(edid_weight_plot(fit, cells = data.frame(group = 99, time = 99)),
               "None of the requested")

  # a fit stripped of its cells errors informatively
  fit_nc <- fit; fit_nc$cells <- NULL
  expect_error(edid_weights(fit_nc), "no stored cells")
})
