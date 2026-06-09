library(testthat)

# ============================================================
# Tests for edid(estimation_effect = ...): the ACH (Ackerberg,
# Chen & Hahn 2012) first-step nuisance-estimation correction.
# ============================================================

make_cfs_panel <- function(n = 200, seed = 1) {
  set.seed(seed)
  Tn   <- 4L
  x1u  <- runif(n, -2, 2)
  eta  <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P    <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gcat <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gcat) & tt >= gcat, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gcat), gcat, 0),
               x1 = x1u, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}

fit_cfs <- function(df, cfs) {
  # misspec_robust = FALSE isolates the ACH (estimation_effect) channel. The misspec_robust master switch
  # now defaults TRUE and would otherwise fold the weight-estimation channel into every fit; these tests
  # target the ACH nuisance correction only.
  edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
       estimation_effect = cfs, misspec_robust = FALSE)
}

test_that("under misspec_robust = FALSE, estimation_effect defaults to FALSE (byte-identical EIF)", {
  df  <- make_cfs_panel(n = 200, seed = 11)
  # The misspec_robust master switch defaults TRUE; with it OFF the fine-grained estimation_effect defaults
  # FALSE, so the plug-in fit is byte-identical to an explicit estimation_effect = FALSE.
  fd  <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none",
              bstrap = FALSE, seed = 1L, misspec_robust = FALSE)
  fF  <- fit_cfs(df, FALSE)
  expect_false(isTRUE(fd$estimation_effect))
  expect_identical(fd$eif, fF$eif)
  expect_identical(fd$att_gt$se, fF$att_gt$se)
})

test_that("correction changes the EIF/SE but NOT the point estimates", {
  df <- make_cfs_panel(n = 200, seed = 12)
  fF <- fit_cfs(df, FALSE)
  fT <- fit_cfs(df, TRUE)
  expect_true(isTRUE(fT$estimation_effect))
  # point estimates are identical (the correction touches only the influence function)
  expect_equal(fF$att_gt$att, fT$att_gt$att, tolerance = 1e-12)
  # the EIF actually changed
  expect_false(isTRUE(all.equal(fF$eif, fT$eif)))
  # all corrected SEs are finite and positive
  ok <- is.finite(fT$att_gt$se)
  expect_true(all(fT$att_gt$se[ok] > 0))
})

test_that("corrected EIF stays mean-zero to machine precision", {
  df <- make_cfs_panel(n = 200, seed = 13)
  fT <- fit_cfs(df, TRUE)
  expect_true(all(abs(colMeans(fT$eif, na.rm = TRUE)) < 1e-10),
              info = paste("max col mean:", max(abs(colMeans(fT$eif, na.rm = TRUE)))))
})

test_that("correction propagates to the event-study aggregation (points equal, SE may differ)", {
  df   <- make_cfs_panel(n = 200, seed = 14)
  esF  <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "event_study", bstrap = FALSE, seed = 1L,
               estimation_effect = FALSE)$event_study
  esT  <- edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "event_study", bstrap = FALSE, seed = 1L,
               estimation_effect = TRUE)$event_study
  skip_if(is.null(esF) || is.null(esT))
  expect_equal(esF$att.egt, esT$att.egt, tolerance = 1e-10)  # ES point estimates unchanged
})

test_that("estimation_effect has no effect without covariates (warns, disabled)", {
  df <- make_cfs_panel(n = 200, seed = 15)
  expect_warning(
    edid(df, "y", "id", "t", "g", xformla = NULL, aggregate = "none", bstrap = FALSE, seed = 1L, estimation_effect = TRUE),
    "no effect without covariates"
  )
})

test_that("conditional-mean ACH correction has the CORRECT SIGN (matches the numerical two-step IF)", {
  # Guards against the OLS-Jacobian sign trap: the moment B*resid has Jacobian -E[BB'], so its
  # first-step IF is +H^{-1} s and the correction must be ADDED. With the shared subtract convention
  # the m-score must carry a leading minus. We verify the package's m-channel EIF change is POSITIVELY
  # correlated with an INDEPENDENT numerical two-step influence function (Gamma x phi_beta, phi_beta by
  # OLS reweighting). A regression of this kind catches a flipped sign (which would give cor = -1).
  df  <- make_cfs_panel(n = 300, seed = 21)
  pn  <- prepare_edid_panel(df, "y", "id", "t", "g", xformla = ~ x1,
                            anticipation = 0L)
  g <- 2L; t <- 4L
  prs <- enumerate_valid_pairs_edid(g, pn$treatment_groups, pn$time_periods, pn$period_1, "all", 0L)
  pfn <- prs; sc <- is.finite(pfn$gp) & (pfn$gp == g); if (any(sc)) pfn$gp[sc] <- Inf
  cr  <- prs[is.finite(prs$gp) & prs$gp != g, , drop = FALSE]
  if (nrow(cr) > 0L) pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cr$tpre))))
  fold <- rep(1L, pn$n)
  cm <- estimate_all_conditional_means(pn, pfn, t_val = t, bs_df = 4L, K_folds = 1L, fold_id = fold, return_aux = TRUE)
  pr <- estimate_all_propensity_ratios(pn, g, pfn, bs_df = 4L, K_folds = 1L, fold_id = fold, return_aux = TRUE)
  prop_ratios <- pr$predictions; cond_means <- cm$predictions; m_aux <- cm$aux
  H <- nrow(prs); w <- rep(1 / H, H)
  pkg_change <- -compute_ach_correction_cov_edid(pn, g, t, prs, prop_ratios, cond_means, w, m_aux, list())
  go0 <- compute_generated_outcomes_cov_edid(pn, g, t, prs, prop_ratios, cond_means, "all")
  m0  <- mean(as.vector(go0 %*% w)); true_corr <- numeric(pn$n)
  for (key in names(m_aux)) {
    a <- m_aux[[key]]; if (isTRUE(a$is_fallback) || is.null(a$B_test)) next
    B <- a$B_test; p <- ncol(B); base <- cond_means[[key]]; eps <- 1e-6 * (1 + max(abs(base)))
    Gamma <- vapply(seq_len(p), function(j) { cm2 <- cond_means; cm2[[key]] <- base + eps * B[, j]
      (mean(as.vector(compute_generated_outcomes_cov_edid(pn, g, t, prs, prop_ratios, cm2, "all") %*% w)) - m0) / eps }, numeric(1))
    parts <- strsplit(key, "_", fixed = TRUE)[[1]]; cohort <- parts[1]; t1 <- as.numeric(parts[2])
    mask <- if (cohort == "Inf") is.infinite(pn$unit_cohorts) else pn$unit_cohorts == as.numeric(cohort)
    col1 <- pn$period_to_col[[as.character(pn$period_1)]]; colt1 <- pn$period_to_col[[as.character(t1)]]
    yd <- pn$outcome_wide[, colt1] - pn$outcome_wide[, col1]; Bc <- B[mask, , drop = FALSE]; yc <- yd[mask]
    bhat <- solve(crossprod(Bc), crossprod(Bc, yc)); idxc <- which(mask); phib <- matrix(0, pn$n, p)
    for (ii in idxc) { wt <- rep(1, length(yc)); wt[match(ii, idxc)] <- 1 + 1e-4
      phib[ii, ] <- (solve(t(Bc) %*% (wt * Bc), t(Bc) %*% (wt * yc)) - bhat) / 1e-4 }
    true_corr <- true_corr + as.vector(phib %*% Gamma)
  }
  cc <- suppressWarnings(stats::cor(pkg_change, true_corr))
  # Under covr instrumentation the deep numerical-derivative internals can degenerate one input to a
  # constant vector, making cor() NA; a genuine sign flip would give a finite cor ~= -1 (not NA), so
  # skipping only on NA keeps the sign guard fully intact everywhere it can actually run.
  #
  # The m-channel ACH is Neyman-orthogonal under the UNIFORM weights used here, so on this fixture the
  # correction is ~0: the exact analytic Gamma (the default) returns ~1e-14, while the old finite difference
  # returned ~1e-8 of numerical residual. cor() of a ~0 vector has no defined sign, so skip the sign check when
  # the correction is negligible. The non-negligible (efficient-weight) ACH -- where a sign flip WOULD matter --
  # is validated directly against the FD oracle in the next test.
  skip_if(is.na(cc) || max(abs(pkg_change)) < 1e-6,
          "m-channel ACH is orthogonal (~0) under uniform weights; sign undefined (see FD-oracle test)")
  expect_gt(cc, 0.95)
})

test_that("analytic ACH reproduces the finite-difference oracle where the correction is non-negligible", {
  # The default analytic ACH Gamma is exact (closed form, no finite differences). On the efficient-weight path
  # the ACH is a LARGE, real correction (it moves the SEs substantially), so it is the right place to verify the
  # analytic against the forced finite-difference oracle -- any sign or coefficient error would show as a gross
  # mismatch. (Under uniform weights the m-channel is orthogonal and ~0, which is why the sign test above is
  # skipped there; here the correction is non-zero so the comparison is meaningful.)
  df  <- make_cfs_panel(n = 300, seed = 21)
  run <- function(ach, ee) {
    op <- options(edid_ach = ach); on.exit(options(op), add = TRUE)
    edid(df, "y", "id", "t", "g", xformla = ~ x1, weight_scheme = "efficient", aggregate = "none",
         bstrap = FALSE, seed = 1L, misspec_robust = FALSE, estimation_effect = ee)$att_gt$se
  }
  se_noee <- run("analytic", FALSE)   # plug-in SE (no ACH)
  se_an   <- run("analytic", TRUE)    # ACH via the exact analytic Gamma (default)
  se_fd   <- run("fd",       TRUE)    # ACH via the finite-difference oracle
  ok <- is.finite(se_noee) & is.finite(se_an) & is.finite(se_fd)
  skip_if(sum(ok) < 2L, "too few non-degenerate cells")
  expect_gt(max(abs(se_an[ok] - se_noee[ok])), 1e-3)   # the ACH is a real, non-negligible correction here
  expect_equal(se_an[ok], se_fd[ok], tolerance = 1e-5) # analytic reproduces the FD oracle (sign + magnitude)
})
