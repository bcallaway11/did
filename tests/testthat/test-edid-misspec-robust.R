# Tests for the opt-in misspecification-robust SE in edid (misspec_robust = TRUE): the weight-estimation
# influence-function channel psi_Omega folded into the EIF.
#
# Coverage:
#   (a) misspec_robust = FALSE is byte-identical to the default (no behavior change).
#   (b) misspec_robust = TRUE leaves point estimates unchanged and the augmented EIF mean-zero (efficient + averaged).
#   (c) the production fold equals the jackknife-validated diagnostic: eif(mr=T) - eif(mr=F) == psi_Omega (data - corr).
#   (d) composition with estimation_effect is additive: eif(ee=T, mr=T) - eif(ee=T, mr=F) == the same psi_Omega; mean-zero.
#   (e) aggregations + clustering inherit the channel (finite, ordered; att unchanged), no separate assembler.
#   (f) guards: gmm / uniform / no-covariate warn and fall back to the plug-in SE; cband_method is NOT coerced.
#   (g) the existing edid suite stays green (run separately).

make_mr_panel <- function(n = 320, seed = 7, bump = 0.6) {
  set.seed(seed)
  Tn   <- 4L
  x1u  <- runif(n, -2, 2)
  eta  <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P    <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gcat <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  is_inf <- is.infinite(gcat)
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)            # x^2 trend the ~x1 cond-mean cannot capture => misspec
    tau <- ifelse(is.finite(gcat) & tt >= gcat, 1, 0)
    bt  <- if (tt == 1L) bump * is_inf * (1 + 0.5 * x1u) else 0
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gcat), gcat, 0),
               x1 = x1u, y = alph + 0.3 * tt + ht + tau + bt + rnorm(n))
  })
  do.call(rbind, rows)
}

fit_mr <- function(df, weights = "efficient", ...) {
  suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1,
                        control_group = "nevertreated", aggregate = "none", weights = weights, ...))
}

# (a) -----------------------------------------------------------------------------------------------------
test_that("misspec_robust = FALSE is byte-identical to the default", {
  df <- make_mr_panel(n = 240, seed = 5)
  for (wm in c("efficient", "averaged", "gmm")) {
    f0 <- fit_mr(df, weights = wm)
    fF <- fit_mr(df, weights = wm, misspec_robust = FALSE)
    expect_identical(fF$att_gt$se, f0$att_gt$se)
    expect_identical(fF$att_gt$ci_lower, f0$att_gt$ci_lower)
    expect_identical(fF$att_gt$ci_upper, f0$att_gt$ci_upper)
    expect_identical(fF$eif, f0$eif)
  }
})

# (b) -----------------------------------------------------------------------------------------------------
test_that("misspec_robust = TRUE: point estimates unchanged, augmented EIF mean-zero, SEs finite", {
  df <- make_mr_panel(n = 320, seed = 7)
  for (wm in c("efficient", "averaged", "gmm")) {
    f0 <- fit_mr(df, weights = wm)
    f1 <- fit_mr(df, weights = wm, misspec_robust = TRUE)
    k  <- !f1$att_gt$is_pre
    expect_equal(f1$att_gt$att, f0$att_gt$att, tolerance = 1e-12)          # att UNCHANGED
    expect_true(all(is.finite(f1$att_gt$se[k])) && all(f1$att_gt$se[k] > 0))
    expect_lt(max(abs(colMeans(as.matrix(f1$eif)[, which(k), drop = FALSE]))), 1e-8)  # augmented EIF mean-zero
  }
})

# (c) the production fold IS the jackknife-validated diagnostic psi_Omega -----------------------------------
test_that("eif(misspec_robust=TRUE) - eif(FALSE) equals the diagnostic psi_Omega (data - corr)", {
  df <- make_mr_panel(n = 320, seed = 7)
  for (wm in c("efficient", "averaged", "gmm")) {
    f0 <- fit_mr(df, weights = wm)                                          # plug-in eif (no fold)
    f1 <- fit_mr(df, weights = wm, misspec_robust = TRUE)                   # folded eif
    old <- options(edid_store_psiomega = TRUE, edid_psiomega_acc = list())  # diagnostic accumulator
    on.exit(options(old), add = TRUE)
    fd <- fit_mr(df, weights = wm)                                          # computes psi but does NOT fold
    pa <- getOption("edid_psiomega_acc"); options(edid_store_psiomega = NULL, edid_psiomega_acc = NULL)
    cn <- paste0(fd$att_gt$group, "_", fd$att_gt$time)
    e0 <- as.matrix(f0$eif); e1 <- as.matrix(f1$eif)
    checked <- 0L
    for (ci in which(!fd$att_gt$is_pre)) {
      p <- pa[[cn[ci]]]
      if (is.null(p)) {                                                    # fallback cell: psi = NULL => fold adds +0
        expect_equal(e1[, ci], e0[, ci], tolerance = 1e-12)               # ...so the EIF is unchanged (no stale psi)
        next
      }
      expect_equal(e1[, ci] - e0[, ci], p$data - p$corr, tolerance = 1e-9) # fold == validated psi exactly
      checked <- checked + 1L
    }
    expect_gt(checked, 0L)
  }
})

# (d) composition with estimation_effect -------------------------------------------------------------------
test_that("misspec_robust composes additively with estimation_effect (same psi added; combined EIF mean-zero)", {
  df <- make_mr_panel(n = 320, seed = 11)
  for (wm in c("efficient", "averaged", "gmm")) {                                  # assert BOTH fold sites' ACH+psi ordering
    f_e  <- fit_mr(df, weights = wm, estimation_effect = TRUE)
    f_em <- fit_mr(df, weights = wm, estimation_effect = TRUE, misspec_robust = TRUE)
    old <- options(edid_store_psiomega = TRUE, edid_psiomega_acc = list())
    fd <- fit_mr(df, weights = wm)
    pa <- getOption("edid_psiomega_acc"); options(old)
    cn <- paste0(fd$att_gt$group, "_", fd$att_gt$time)
    ee <- as.matrix(f_e$eif); eem <- as.matrix(f_em$eif); k <- !fd$att_gt$is_pre
    expect_equal(f_em$att_gt$att, f_e$att_gt$att, tolerance = 1e-12)        # att unchanged by misspec_robust
    expect_lt(max(abs(colMeans(eem[, which(k), drop = FALSE]))), 1e-8)      # combined (-ach + psi) EIF mean-zero
    for (ci in which(k)) {
      p <- pa[[cn[ci]]]; if (is.null(p)) next
      expect_equal(eem[, ci] - ee[, ci], p$data - p$corr, tolerance = 1e-9) # adds the SAME psi on top of the ACH eif
    }
  }
})

# (e) aggregations + clustering ride-through ---------------------------------------------------------------
test_that("misspec_robust: the channel rides through aggregations + clustering (aggregate SE MOVES, att fixed)", {
  df <- make_mr_panel(n = 360, seed = 13)
  df$cl <- (df$id %% 30L)                                                   # 30 clusters
  for (agg in c("group", "event_study")) {
    a0 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, control_group = "nevertreated",
                                weights = "efficient", aggregate = agg))
    a1 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, control_group = "nevertreated",
                                weights = "efficient", aggregate = agg, misspec_robust = TRUE))
    r0 <- a0[[agg]]; r1 <- a1[[agg]]                                        # the aggregate result (att.egt/se.egt/overall.*)
    expect_equal(r1$att.egt, r0$att.egt, tolerance = 1e-12)                # aggregate point estimates UNCHANGED
    expect_true(all(is.finite(r1$se.egt)) && is.finite(r1$overall.se))    # finite
    expect_gt(max(abs(r1$se.egt - r0$se.egt)), 1e-8)                       # the channel PROPAGATED to the aggregate SE
  }
  cl0 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, control_group = "nevertreated",
                               weights = "averaged", aggregate = "none", clustervars = "cl"))
  cl1 <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, control_group = "nevertreated",
                               weights = "averaged", aggregate = "none", clustervars = "cl", misspec_robust = TRUE))
  k <- !cl1$att_gt$is_pre
  expect_true(all(is.finite(cl1$att_gt$se[k])))                            # cluster-robust fold rides rowsum
  expect_gt(max(abs(cl1$att_gt$se[k] - cl0$att_gt$se[k])), 1e-8)           # ...and actually moves the clustered SE
})

# (f) guards ----------------------------------------------------------------------------------------------
test_that("misspec_robust warns + falls back to plug-in for uniform / no covariates; stays active for gmm", {
  df <- make_mr_panel(n = 240, seed = 5)
  catch_warnings <- function(expr) {                                       # collect ALL warnings (edid emits several)
    ws <- character(0)
    val <- withCallingHandlers(expr,
      warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(value = val, warnings = ws)
  }
  # uniform: fixed weights have no estimation channel -> warn-disable, SE unchanged
  f0u <- fit_mr(df, weights = "uniform")
  ru  <- catch_warnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, control_group = "nevertreated",
                             aggregate = "none", weights = "uniform", misspec_robust = TRUE))
  expect_true(any(grepl("misspec_robust", ru$warnings)))                   # the warn-disable warning fired
  expect_identical(ru$value$att_gt$se, f0u$att_gt$se)                      # SE unchanged (fell back to plug-in)
  expect_false(ru$value$misspec_robust)                                    # stored S3 flag reflects the downgrade
  # gmm: now ACTIVE (sample-cov channel + the ACH correction for u'Cw) -> SE changes, flag stays TRUE
  f0g <- fit_mr(df, weights = "gmm")
  fMg <- fit_mr(df, weights = "gmm", misspec_robust = TRUE)
  expect_true(fMg$misspec_robust)                                          # gmm is NOT warn-disabled
  expect_false(isTRUE(all.equal(fMg$att_gt$se[!fMg$att_gt$is_pre], f0g$att_gt$se[!f0g$att_gt$is_pre])))  # SE moves
  expect_equal(fMg$att_gt$att, f0g$att_gt$att, tolerance = 1e-12)          # point estimates unchanged
  # efficient stays active; no covariates warn-disables
  expect_true(fit_mr(df, weights = "efficient", misspec_robust = TRUE)$misspec_robust)
  r <- catch_warnings(edid(df, "y", "id", "t", "g", control_group = "nevertreated",   # no xformla
                           weights = "averaged", aggregate = "none", misspec_robust = TRUE))
  expect_true(any(grepl("without covariates", r$warnings)))
})

test_that("misspec_robust does NOT coerce cband_method (a real IF is carried by the multiplier bootstrap)", {
  df <- make_mr_panel(n = 220, seed = 9)
  fM <- suppressWarnings(edid(df, "y", "id", "t", "g", xformla = ~ x1, control_group = "nevertreated",
                              weights = "efficient", aggregate = "none", misspec_robust = TRUE,
                              cband_method = "multiplier", bstrap = TRUE, biters = 199L))
  expect_identical(fM$cband_method, "multiplier")
})
