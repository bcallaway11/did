# Tests for the analytic sup-t (MOPM) uniform confidence bands in edid (cband_method = "analytic").

test_that("supt_crit_edid matches the analytic Sidak crit for independent coordinates", {
  for (p in c(3L, 7L, 20L)) {
    sidak <- qnorm(1 - (1 - 0.95^(1 / p)) / 2)
    crit  <- supt_crit_edid(diag(p), alp = 0.05, B = 2e5L, seed = 7L)
    expect_equal(crit, sidak, tolerance = 0.02)   # Monte Carlo, ~B-noise
  }
})

test_that("build_kernel_weights_edid handles huge finite covariates without non-finite kernel distances", {
  X <- matrix(rep(1e308, 40L * 3L), nrow = 40L, ncol = 3L)
  K <- build_kernel_weights_edid(X)$K
  expect_true(all(is.finite(K)))
  expect_true(all(K >= 0) && all(!is.nan(K)))
  expect_true(all(abs(K - K[1L, 1L]) < .Machine$double.eps^0.5))
})

test_that("compute_pointwise_weights_edid falls back to uniform when local Omega has non-finite entries", {
  omega_array <- array(0, dim = c(200L, 2L, 2L))
  omega_array[1L, , ] <- matrix(c(NaN, 0.1, 0.1, 1), ncol = 2L)
  for (i in 2:200) omega_array[i, , ] <- matrix(c(4, 1, 1, 3), ncol = 2L)

  gen_out <- matrix(rep(1, 400L), nrow = 200L, ncol = 2L)
  out <- compute_pointwise_weights_edid(omega_array, d = 1L, gen_out_mat = gen_out, need_coup = TRUE)

  expect_equal(out$W[1L, ], c(0.5, 0.5))
  # for the finite block, weights match Ω^{-1}1/(1'Ω^{-1}1) for Ω = [[4,1],[1,3]]
  expect_equal(out$W[-1L, ], matrix(c(0.4, 0.6), nrow = 199L, ncol = 2L, byrow = TRUE))
  expect_equal(out$Q[1L, ], c(0, 0))
  expect_equal(out$C[1L, , ], matrix(0, 2L, 2L))
})

test_that("supt_crit_edid and analytic_bands_edid skip coordinates with non-finite covariance rows", {
  Sigma <- matrix(c(2, NA, 0.2,
                    NA, 1, 0.3,
                    0.2, 0.3, 1), nrow = 3L, byrow = TRUE)
  att <- c(0.1, 0.2, 0.3)
  expect_equal(supt_crit_edid(Sigma, seed = 1L), qnorm(0.975))
  expect_warning(
    bands <- analytic_bands_edid(att, Sigma, seed = 1L),
    "degenerate variance"
  )

  # only coordinate 3 has finite row/column, so the simultaneous crit is pointwise.
  expect_equal(bands$crit, qnorm(0.975))
  expect_true(all(is.na(bands$ci_lower[1L:2])) && all(is.na(bands$ci_upper[1L:2]))
              && all(is.finite(bands$ci_lower[3L])) && all(is.finite(bands$ci_upper[3L])))
})

test_that("supt_crit_edid is between the pointwise z and the Bonferroni bound, and rises with correlation", {
  p <- 6L
  zc  <- qnorm(0.975)
  bon <- qnorm(1 - 0.05 / (2 * p))
  c_ind <- supt_crit_edid(diag(p), seed = 1L)
  expect_gt(c_ind, zc); expect_lt(c_ind, bon)
  # higher positive correlation -> coordinates move together -> smaller sup-t crit
  R_hi <- matrix(0.7, p, p); diag(R_hi) <- 1
  expect_lt(supt_crit_edid(R_hi, seed = 1L), c_ind)
  # degenerate (single coordinate) -> pointwise
  expect_equal(supt_crit_edid(matrix(1, 1, 1)), zc)
})

test_that("cluster_cov_edid: iid form and sqrt(diag) == safe_inference_edid SE", {
  set.seed(1); n <- 200L; M <- matrix(rnorm(n * 3L), n, 3L)
  v <- cluster_cov_edid(M, NULL, n)
  expect_equal(v, crossprod(M) / n^2)
  expect_equal(sqrt(diag(v))[1], safe_inference_edid(M[, 1L], NULL, att = 0)$se)
  # cluster-robust: G/(G-1) sandwich and matches safe_inference_edid clustered SE
  ci <- rep(1:40, length.out = n)
  vc <- cluster_cov_edid(M, ci, n)
  expect_equal(sqrt(diag(vc))[1], safe_inference_edid(M[, 1L], ci, att = 0)$se)
})

test_that("edid default produces analytic SIMULTANEOUS bands (crit > z); SE unchanged vs pointwise", {
  data(mpdta, package = "did")
  fa <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year", gname = "first.treat",
             xformla = ~lpop, seed = 1L)
  expect_identical(fa$cband_method, "analytic")
  ok <- is.finite(fa$att_gt$se) & fa$att_gt$se > 0
  crit_impl <- (fa$att_gt$ci_upper[ok] - fa$att_gt$att[ok]) / fa$att_gt$se[ok]
  expect_true(all(crit_impl > qnorm(0.975) - 1e-8))      # simultaneous, wider than pointwise
  expect_true(all(abs(diff(crit_impl)) < 1e-6))          # one common crit across cells
  expect_true(all(is.finite(fa$att_gt$ci_lower[ok])) && all(fa$att_gt$ci_lower[ok] < fa$att_gt$ci_upper[ok]))
  # SE is the analytic first-order SE -> identical to cband = FALSE (only the band differs)
  fp <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year", gname = "first.treat",
             xformla = ~lpop, cband = FALSE)
  expect_equal(fa$att_gt$se, fp$att_gt$se)
  pw <- (fp$att_gt$ci_upper[ok] - fp$att_gt$att[ok]) / fp$att_gt$se[ok]
  expect_true(all(abs(pw - qnorm(0.975)) < 1e-8))        # cband = FALSE -> pointwise
})

test_that("cband_method = 'multiplier' preserves the bootstrap path and honors cband", {
  data(mpdta, package = "did")
  # multiplier without bootstrap is coerced to analytic instead of silently producing pointwise bands
  expect_warning(
    fm0 <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year", gname = "first.treat",
                xformla = ~lpop, cband_method = "multiplier"),
    "requires bstrap = TRUE"
  )
  expect_identical(fm0$cband_method, "analytic")
  ok <- is.finite(fm0$att_gt$se) & fm0$att_gt$se > 0
  pw <- (fm0$att_gt$ci_upper[ok] - fm0$att_gt$att[ok]) / fm0$att_gt$se[ok]
  expect_true(all(pw > qnorm(0.975) - 1e-8))
  # multiplier + bstrap + cband = FALSE -> bootstrap SEs with pointwise intervals
  fm_pw <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year", gname = "first.treat",
                xformla = ~lpop, bstrap = TRUE, biters = 199L, cband = FALSE,
                cband_method = "multiplier", seed = 1L)
  ok_pw <- is.finite(fm_pw$att_gt$se) & fm_pw$att_gt$se > 0
  crit_pw <- (fm_pw$att_gt$ci_upper[ok_pw] - fm_pw$att_gt$att[ok_pw]) / fm_pw$att_gt$se[ok_pw]
  expect_true(all(abs(crit_pw - qnorm(0.975)) < 1e-8))
  # multiplier + bstrap = TRUE -> the did multiplier bootstrap (simultaneous), finite & ordered
  set.seed(1)
  fm1 <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year", gname = "first.treat",
              xformla = ~lpop, bstrap = TRUE, biters = 199L, cband_method = "multiplier", seed = 1L)
  ok1 <- is.finite(fm1$att_gt$se)
  expect_true(all(fm1$att_gt$ci_lower[ok1] < fm1$att_gt$ci_upper[ok1]))
  crit_sim <- (fm1$att_gt$ci_upper[ok1] - fm1$att_gt$att[ok1]) / fm1$att_gt$se[ok1]
  expect_true(all(crit_sim > qnorm(0.975) - 1e-8))
})

test_that("bstrap = TRUE selects the multiplier bootstrap under the default cband_method, and inference_type is honest", {
  data(mpdta, package = "did")
  base <- list(data = mpdta, yname = "lemp", idname = "countyreal", tname = "year",
               gname = "first.treat", xformla = ~lpop)
  f_def <- do.call(edid, base)                                          # default: analytic, no bootstrap
  expect_identical(f_def$cband_method, "analytic")
  expect_identical(f_def$inference_type, "analytical")
  f_bs <- do.call(edid, c(base, list(bstrap = TRUE, biters = 199L, seed = 1L)))   # legacy contract restored
  expect_identical(f_bs$cband_method, "multiplier")
  expect_identical(f_bs$inference_type, "bootstrap")
  expect_false(isTRUE(all.equal(f_bs$att_gt$se, f_def$att_gt$se)))      # the bootstrap actually ran
  f_an <- do.call(edid, c(base, list(bstrap = TRUE, cband_method = "analytic")))  # explicit analytic wins
  expect_identical(f_an$cband_method, "analytic")
  expect_identical(f_an$inference_type, "analytical")                  # not misreported as bootstrap
  f_ho <- do.call(edid, c(base, list(bstrap = TRUE, higher_order = TRUE)))        # higher_order forces analytic
  expect_identical(f_ho$cband_method, "analytic")
  expect_identical(f_ho$inference_type, "analytical")
  expect_false(as_MP_edid(f_ho)$DIDparams$bstrap)
  expect_false(aggte_edid(f_ho, type = "dynamic", na.rm = TRUE)$DIDparams$bstrap)
  f_mb <- do.call(edid, c(base, list(bstrap = TRUE, biters = 99L, cband_method = "multiplier",
                                    cband = TRUE, seed = 1L)))
  expect_true(as_MP_edid(f_mb)$DIDparams$bstrap)
  expect_true(as_MP_edid(f_mb)$DIDparams$cband)
})

test_that("default analytic cband (seed = NULL) does not perturb the caller's RNG stream", {
  data(mpdta, package = "did")
  set.seed(123); a <- runif(1L)
  set.seed(123); invisible(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                                gname = "first.treat", xformla = ~lpop))
  b <- runif(1L)
  expect_equal(a, b)
})

test_that("analytic bands are cluster-robust and the aggregations carry the analytic sup-t crit", {
  data(mpdta, package = "did")
  set.seed(2); mpdta$clu <- as.integer(factor(mpdta$countyreal)) %% 30L
  f_iid <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year", gname = "first.treat",
                xformla = ~lpop, seed = 1L)
  f_cl  <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year", gname = "first.treat",
                xformla = ~lpop, clustervars = "clu", seed = 1L)
  expect_false(isTRUE(all.equal(f_iid$att_gt$se, f_cl$att_gt$se)))   # cluster-robust SE differs
  es <- f_cl$event_study
  expect_true(is.finite(es$crit.val.egt) && es$crit.val.egt > qnorm(0.975) - 1e-8)
})
