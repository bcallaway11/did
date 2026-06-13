# Tests for the NO-COVARIATE estimation_effect: the closed-form second-order
# weight-estimation variance correction (compute_nocov_ee_correction_edid).
#
#   (a) FD oracle: the analytic per-unit Jacobian directions D (plain map AND through the
#       nocov_shrink chain) reproduce brute-force finite differences of w(Omega_sh(Omega));
#       the assembled var_add matches the FD-rebuilt assembly.
#   (b) structural identities: Omega-hat == crossprod(psi)/n^2; sum_i d_i = 0; q_opt == -cov_lead;
#       eif == psi %*% w; lambda = 1 clamp => var_add == delta_df (pure Bessel).
#   (c) API: default no-covariate fits are bit-for-bit unchanged; estimation_effect = TRUE engages
#       the correction (att unchanged, SE^2 = plug-in SE^2 + var_add, sigma_nocov_ee stored);
#       explicit misspec_robust = TRUE on a no-covariate fit auto-enables it; uniform weights
#       warn-disable; PT-Post has no correction (H = 1).
#   (d) propagation: the event-study / overall aggregate SEs inherit the increment.
#   (e) asymptotic no-op: the correction is O(1/n) relative -- negligible at large n.

mk_ee_panel <- function(n, seed, rho = 0, Tn = 6L) {
  set.seed(seed)
  G <- sample(c(3L, 5L, 0L), n, TRUE, c(.35, .30, .35))
  e <- matrix(rnorm(n * Tn), n, Tn)
  if (rho > 0) for (s in 2:Tn) e[, s] <- rho * e[, s - 1L] + sqrt(1 - rho^2) * e[, s]
  Y <- rnorm(n) + matrix(0.3 * seq_len(Tn), n, Tn, byrow = TRUE) + e
  for (g in c(3L, 5L)) for (t in g:Tn) Y[G == g, t] <- Y[G == g, t] + 1 + 0.3 * (t - g)
  data.frame(id = rep(seq_len(n), each = Tn), time = rep(seq_len(Tn), n),
             y = as.vector(t(Y)), g = rep(G, each = Tn))
}

# internal-pieces fixture for one cell (uses Inf never-treated convention directly)
ee_cell_pieces <- function(n, seed, rho, gg, tt, use_shrink) {
  df <- mk_ee_panel(n, seed, rho)
  df$g <- ifelse(df$g == 0L, Inf, df$g)
  pn  <- prepare_edid_panel(df, "y", "id", "time", "g", anticipation = 0L)
  prs <- enumerate_valid_pairs_edid(gg, pn$treatment_groups, pn$time_periods,
                                    pn$period_1, "all", 0L)
  Om  <- compute_omega_star_nocov_edid(gg, tt, prs, pn, "all")
  psi <- compute_psi_moments_nocov_edid(gg, tt, prs, pn)
  lam <- NA_real_; Omu <- Om
  if (use_shrink) {
    sh <- shrink_omega_nocov_edid(Om, gg, tt, prs, pn)
    lam <- sh$lambda; Omu <- sh$omega
  }
  w  <- compute_efficient_weights_edid(Omu)
  ee <- compute_nocov_ee_correction_edid(gg, tt, prs, pn, omega_raw = Om, omega_used = Omu,
                                         weights = w, shrink_lambda = lam, return_D = TRUE)
  list(df = df, pn = pn, prs = prs, Om = Om, Omu = Omu, psi = psi, w = w, lam = lam, ee = ee,
       q4 = sum(rowSums(psi * psi)^2))
}

# the full weight map Omega -> (optional shrink, with S, n, q4 fixed) -> w, for FD
w_of_omega_ref <- function(Om, shrink = NULL) {
  Os <- Om
  if (!is.null(shrink)) {
    S <- shrink$S; ss <- sum(S * S)
    sigma2 <- sum(Om * S) / ss
    target <- sigma2 * S
    d2 <- sum((Om - target)^2)
    b2 <- (shrink$q4 / shrink$n^2 - shrink$n * sum(Om * Om)) / shrink$n^2
    lam <- min(1, max(0, b2) / d2)
    Os <- (1 - lam) * Om + lam * target
  }
  A <- solve(Os); u <- drop(A %*% rep(1, nrow(Os))); u / sum(u)
}

test_that("FD oracle: analytic Jacobian directions and var_add match finite differences (plain + shrink chain)", {
  skip_on_cran()
  for (cfg in list(list(seed = 101, n = 60,  rho = 0.0, shrink = FALSE, g = 3L, t = 4L),
                   list(seed = 202, n = 80,  rho = 0.7, shrink = FALSE, g = 5L, t = 6L),
                   list(seed = 505, n = 100, rho = 0.5, shrink = TRUE,  g = 5L, t = 5L),
                   list(seed = 606, n = 150, rho = 0.9, shrink = TRUE,  g = 3L, t = 4L))) {
    p <- ee_cell_pieces(cfg$n, cfg$seed, cfg$rho, cfg$g, cfg$t, cfg$shrink)
    skip_if(!isTRUE(p$ee$applied), "correction not applied on this draw")
    # the lambda = 1 clamp is a kink (one-sided derivative; correction identically the Bessel
    # piece there) -- FD comparison is only meaningful at interior lambda / no shrink
    skip_if(cfg$shrink && (!is.finite(p$lam) || p$lam <= 0 || p$lam >= 1), "lambda clamped")
    shr <- if (cfg$shrink) list(S = compute_pole_structure_nocov_edid(cfg$g, cfg$t, p$prs, p$pn),
                                n = p$pn$n, q4 = p$q4) else NULL
    n_u <- p$pn$n; H <- nrow(p$prs)
    h <- 1e-6 * max(abs(p$Om))
    Dfd <- matrix(NA_real_, n_u, H)
    for (i in seq_len(n_u)) {
      Vi <- tcrossprod(p$psi[i, ]) / n_u - p$Om
      Dfd[i, ] <- (w_of_omega_ref(p$Om + h * Vi, shr) - w_of_omega_ref(p$Om - h * Vi, shr)) / (2 * h)
    }
    expect_lt(max(abs(p$ee$D - Dfd)) / max(abs(Dfd)), 1e-5)
    # rebuild var_add from the FD directions with the same assembly
    a   <- drop(p$psi %*% p$w)
    qfd <- -sum(a * rowSums(Dfd * p$psi)) / n_u^3
    expect_lt(abs(p$ee$var_add - (p$ee$delta_df + 2 * qfd)) / max(abs(p$ee$var_add), 1e-300), 1e-5)
  }
})

test_that("structural identities: psi/Omega, sum d_i = 0, q_opt = -cov_lead, eif = psi w, B-orthogonality", {
  p <- ee_cell_pieces(80, 7, 0.5, 3L, 5L, use_shrink = FALSE)
  expect_lt(max(abs(p$Om - crossprod(p$psi) / p$pn$n^2)), 1e-12 * max(abs(p$Om)))  # exact identity
  expect_true(isTRUE(p$ee$applied))
  expect_lt(max(abs(colSums(p$ee$D))), 1e-10 * max(abs(p$ee$D)))                   # sum_i d_i = 0 exactly
  expect_equal(p$ee$q_opt, -p$ee$cov_lead, tolerance = 1e-12)
  expect_gte(p$ee$q_opt, 0)    # exact path: B PSD => optimism is non-negative
  expect_gt(p$ee$delta_df, 0)
  # the cell EIF is psi %*% w
  m   <- compute_generated_outcomes_nocov_edid(3L, 5L, p$prs, p$pn, "all")
  eif <- compute_eif_nocov_edid(3L, 5L, p$prs, p$w, p$pn, sum(p$w * m), "all")
  expect_equal(drop(p$psi %*% p$w), eif, tolerance = 1e-12)
})

test_that("lambda = 1 clamp: the Jacobian vanishes and var_add is the pure Bessel piece", {
  p <- ee_cell_pieces(80, 7, 0.5, 3L, 5L, use_shrink = FALSE)
  # force the clamp: omega_used = sigma2 * S (what shrinkage at lambda = 1 inverts)
  S <- compute_pole_structure_nocov_edid(3L, 5L, p$prs, p$pn)
  sigma2 <- sum(p$Om * S) / sum(S * S)
  w1 <- compute_efficient_weights_edid(sigma2 * S)
  ee1 <- compute_nocov_ee_correction_edid(3L, 5L, p$prs, p$pn, omega_raw = p$Om,
                                          omega_used = sigma2 * S, weights = w1,
                                          shrink_lambda = 1, return_D = TRUE)
  skip_if(!isTRUE(ee1$applied), "pole matrix not invertible on this draw")
  expect_lt(max(abs(ee1$D)), 1e-10)                       # B S w = B Omega_sh w / sigma2 = 0
  expect_equal(ee1$var_add, ee1$delta_df, tolerance = 1e-12)
})

test_that("default no-covariate fits are bit-for-bit unchanged; estimation_effect = TRUE engages the correction", {
  df <- mk_ee_panel(120, 11)
  f_def <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE))
  f_off <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                                 estimation_effect = FALSE, misspec_robust = FALSE))
  expect_identical(f_def$att_gt$se,  f_off$att_gt$se)     # default == flags off (bit-for-bit)
  expect_identical(f_def$att_gt$att, f_off$att_gt$att)
  expect_null(f_def$sigma_nocov_ee)
  expect_false(isTRUE(f_def$estimation_effect))

  f_ee <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                                estimation_effect = TRUE))
  expect_true(isTRUE(f_ee$estimation_effect))
  expect_equal(f_ee$att_gt$att, f_def$att_gt$att, tolerance = 1e-12)   # point estimates unchanged
  va <- vapply(f_ee$cells, function(cc) {
    if (is.null(cc$nocov_ee) || !isTRUE(cc$nocov_ee$applied)) NA_real_ else cc$nocov_ee$var_add
  }, numeric(1L))
  expect_true(any(is.finite(va)))                          # the correction applied somewhere
  k <- which(is.finite(va))
  # fold identity: corrected SE^2 = plug-in SE^2 + var_add, cell by cell
  expect_equal(f_ee$att_gt$se[k]^2, f_def$att_gt$se[k]^2 + va[k], tolerance = 1e-10)
  expect_true(all(f_ee$att_gt$se[k] > f_def$att_gt$se[k]))  # delta_df > 0 and exact-path q >= 0
  expect_false(is.null(f_ee$sigma_nocov_ee))
  expect_equal(diag(f_ee$sigma_nocov_ee)[k], unname(va[k]), tolerance = 1e-12)
  # no obsolete downgrade warning
  w <- character(0)
  withCallingHandlers(
    edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE, estimation_effect = TRUE),
    warning = function(ww) { w <<- c(w, conditionMessage(ww)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("no effect without covariates", w)))
})

test_that("explicit misspec_robust = TRUE on a no-covariate fit engages the weight-estimation correction", {
  df <- mk_ee_panel(100, 13)
  w <- character(0)
  f_mr <- withCallingHandlers(
    edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE, misspec_robust = TRUE),
    warning = function(ww) { w <<- c(w, conditionMessage(ww)); invokeRestart("muffleWarning") })
  expect_true(isTRUE(f_mr$estimation_effect))              # auto-enabled by the master switch
  expect_false(isTRUE(f_mr$misspec_robust))                # the psi_Omega fold itself stays covariate-only
  expect_true(any(grepl("without covariates", w)))         # the reroute is announced
  f_ee <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                                estimation_effect = TRUE))
  expect_identical(f_mr$att_gt$se, f_ee$att_gt$se)         # same correction either way
})

test_that("uniform weights warn-disable; PT-Post has no correction", {
  df <- mk_ee_panel(100, 17)
  expect_warning(
    f_u <- edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                weight_scheme = "uniform", estimation_effect = TRUE),
    "has no effect for weights = 'uniform'")
  f_u0 <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                                weight_scheme = "uniform"))
  expect_identical(f_u$att_gt$se, f_u0$att_gt$se)
  f_p <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                               pt_assumption = "post", estimation_effect = TRUE))
  expect_true(all(vapply(f_p$cells, function(cc) is.null(cc$nocov_ee), logical(1L))))
  expect_null(f_p$sigma_nocov_ee)
})

test_that("the correction propagates to the event-study and overall aggregate SEs", {
  skip_on_cran()
  df <- mk_ee_panel(120, 19, rho = 0.5)
  f_def <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE, seed = 5L))
  f_ee  <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE, seed = 5L,
                                 estimation_effect = TRUE))
  skip_if(is.null(f_ee$sigma_nocov_ee), "no applied correction on this draw")
  es0 <- aggte_edid(f_def, type = "dynamic"); es1 <- aggte_edid(f_ee, type = "dynamic")
  ov0 <- aggte_edid(f_def, type = "simple");  ov1 <- aggte_edid(f_ee, type = "simple")
  expect_equal(es0$att.egt, es1$att.egt, tolerance = 1e-12)          # point estimates unchanged
  expect_gt(max(es1$se.egt - es0$se.egt), 0)                          # some horizon inherits the increment
  expect_true(all(es1$se.egt >= es0$se.egt - 1e-12))                  # exact-path increments are >= 0
  expect_gt(ov1$overall.se, ov0$overall.se - 1e-12)
})

test_that("multiplier-bootstrap path warns that it cannot carry the correction", {
  skip_on_cran()
  df <- mk_ee_panel(100, 23)
  w <- character(0)
  withCallingHandlers(
    edid(df, "y", "id", "time", "g", aggregate = "none", estimation_effect = TRUE,
         bstrap = TRUE, biters = 20L, cband = FALSE, seed = 2L),
    warning = function(ww) { w <<- c(w, conditionMessage(ww)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("multiplier bootstrap", w)))
})

test_that("asymptotic no-op: the correction is negligible at large n", {
  skip_on_cran()
  df <- mk_ee_panel(20000, 29, rho = 0.7)
  f_def <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE))
  f_ee  <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                                 estimation_effect = TRUE))
  ok <- is.finite(f_def$att_gt$se) & is.finite(f_ee$att_gt$se)
  expect_lt(max(abs(f_ee$att_gt$se[ok] / f_def$att_gt$se[ok] - 1)), 0.005)  # < 0.5% at n = 20k
})

test_that("the K x K increment: diagonal == per-cell var_add; cross entries match a direct recompute", {
  df <- mk_ee_panel(100, 31, rho = 0.4)
  f_ee <- suppressWarnings(edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE,
                                estimation_effect = TRUE))
  skip_if(is.null(f_ee$sigma_nocov_ee), "no applied correction on this draw")
  va <- vapply(f_ee$cells, function(cc) {
    if (is.null(cc$nocov_ee) || !isTRUE(cc$nocov_ee$applied)) 0 else cc$nocov_ee$var_add
  }, numeric(1L))
  expect_equal(unname(diag(f_ee$sigma_nocov_ee)), unname(va), tolerance = 1e-12)
  # direct recompute of one cross entry from internals
  ks <- which(va != 0)
  skip_if(length(ks) < 2L, "need two corrected cells")
  k1 <- ks[1]; k2 <- ks[2]
  dfI <- df; dfI$g <- ifelse(dfI$g == 0L, Inf, dfI$g)
  pn <- prepare_edid_panel(dfI, "y", "id", "time", "g", anticipation = 0L)
  piece <- function(k) {
    cc <- f_ee$cells[[k]]
    prs <- cc$pairs
    Om  <- compute_omega_star_nocov_edid(cc$group, cc$time, prs, pn, "all")
    sh  <- shrink_omega_nocov_edid(Om, cc$group, cc$time, prs, pn)
    w   <- compute_efficient_weights_edid(sh$omega)
    ee  <- compute_nocov_ee_correction_edid(cc$group, cc$time, prs, pn, Om, sh$omega, w, sh$lambda)
    psi <- compute_psi_moments_nocov_edid(cc$group, cc$time, prs, pn)
    list(a = drop(psi %*% w), s = ee$s_vec)
  }
  p1 <- piece(k1); p2 <- piece(k2)
  m_i <- as.numeric(table(pn$unit_cohorts)[as.character(pn$unit_cohorts)])
  n <- pn$n
  cross <- sum(p1$a * p2$a / pmax(m_i - 1, 1)) / n^2 -
           (sum(p1$s * p2$a) + sum(p1$a * p2$s)) / n^3
  expect_equal(f_ee$sigma_nocov_ee[k1, k2], cross, tolerance = 1e-10)
  expect_equal(f_ee$sigma_nocov_ee[k1, k2], f_ee$sigma_nocov_ee[k2, k1], tolerance = 1e-12)
})
