# =============================================================================
# mboot() post-processing: the rewritten critical-value computation (bSigma via
# a single two-probability quantile call, bT via an abs/NA->-Inf/pmax column
# fold) must reproduce the previous row-wise apply() reference, and the new
# all-degenerate branch (ncol(bres) == 0L after dropping degenerate dimensions)
# must return NA crit.val / all-NA se without error -- matching the old
# behavior where every row's max was -Inf and got dropped by is.finite().
# mboot() runs directly on a DIDparams-like list: with clustervars = NULL it
# never touches the data slot (see test-pretest-vectorization.R).
# =============================================================================

dp_mboot <- list(idname = "id", clustervars = NULL, biters = 100,
                 tname = "period", alp = 0.05, panel = TRUE,
                 true_repeated_cross_sections = FALSE,
                 allow_unbalanced_panel = FALSE, faster_mode = FALSE)

test_that("mboot on an all-degenerate influence function returns NA crit.val and all-NA se", {
  n <- 80
  ifz <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(n, 3))   # all-zero sparse n x 3
  set.seed(101)
  out <- mboot(ifz, dp_mboot)
  expect_identical(ncol(out$bres), 0L)          # every dimension dropped as degenerate
  expect_true(is.na(out$crit.val))
  expect_length(out$se, 3)
  expect_true(all(is.na(out$se)))
})

test_that("mboot pmax-fold critical value equals the old row-wise apply reference", {
  n <- 80
  set.seed(202)
  infm <- cbind(matrix(rnorm(n * 3), n, 3), 0)  # 3 finite columns + one all-zero
  set.seed(303)
  out <- mboot(infm, dp_mboot)

  # Degenerate column dropped before post-processing; its se slot stays NA.
  expect_identical(ncol(out$bres), 3L)
  expect_true(all(is.finite(out$se[1:3])))
  expect_true(is.na(out$se[4]))

  # Old reference (pre-rewrite mboot) recomputed on the same returned bres:
  # per-column IQR scale, then per-row max(abs(b / bSigma), na.rm = TRUE).
  bres <- out$bres
  iqr_norm <- qnorm(.75) - qnorm(.25)
  bSigma_old <- apply(bres, 2, function(b) {
    (quantile(b, .75, type = 1, na.rm = TRUE) -
       quantile(b, .25, type = 1, na.rm = TRUE)) / iqr_norm
  })
  bT_old <- suppressWarnings(apply(bres, 1, function(b) max(abs(b / bSigma_old), na.rm = TRUE)))
  bT_old <- bT_old[is.finite(bT_old)]
  crit_old <- unname(quantile(bT_old, 1 - dp_mboot$alp, type = 1, na.rm = TRUE))

  expect_equal(unname(out$crit.val), crit_old, tolerance = 1e-12)
  expect_equal(unname(out$se[1:3]), unname(bSigma_old) / sqrt(n), tolerance = 1e-12)
})
