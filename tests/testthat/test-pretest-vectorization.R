# =============================================================================
# Conditional pre-test internals: the vectorized indicator() and the vectorized
# multiplier bootstrap test.mboot() must reproduce their explicit (loop) reference
# definitions. test.mboot() is exercised both unclustered and clustered (the
# clustered path was previously untested).
# =============================================================================

test_that("indicator() vectorization matches the row-wise definition exactly", {
  old_indicator <- function(X, u) {
    cond <- t(apply(X, 1, function(x) x <= u))
    1 * apply(cond, 1, all)
  }
  # X always carries an intercept column in the pre-test (model.matrix output), so
  # use >= 2 columns -- exactly the shape indicator() is called with in practice.
  set.seed(3)
  X <- matrix(rnorm(500), 100, 5)
  for (i in 1:20) expect_identical(indicator(X, X[i, ]), old_indicator(X, X[i, ]))
  X2 <- cbind(1, matrix(rnorm(200), 100, 2))   # intercept + 2 covariates
  for (i in 1:10) expect_identical(indicator(X2, X2[i, ]), old_indicator(X2, X2[i, ]))
})

test_that("test.mboot equals the explicit per-draw multiplier bootstrap (unclustered)", {
  set.seed(1)
  n <- 120; k <- 2; nX <- 120
  arr <- array(rnorm(n * k * nX), c(n, k, nX))
  dp <- list(data = data.frame(id = 1:n, period = 1L), idname = "id",
             clustervars = NULL, biters = 200, tname = "period", alp = 0.05, panel = TRUE)
  ref_boot <- function(inf.func, n, biters) {
    sapply(1:biters, function(b) {
      Ub <- sample(c(-1, 1), n, replace = TRUE)
      Jb <- t(apply(Ub * inf.func, c(2, 3), mean))
      n * sum(apply(Jb^2, 2, mean))
    })
  }
  set.seed(42); ref <- ref_boot(arr, n, 200)
  set.seed(42); got <- test.mboot(arr, dp)$bres
  expect_equal(unname(got), ref, tolerance = 1e-9)
})

test_that("test.mboot equals the explicit clustered multiplier bootstrap", {
  set.seed(2)
  n <- 100; k <- 2; nX <- 100; ncl <- 20
  arr <- array(rnorm(n * k * nX), c(n, k, nX))
  cl <- rep(1:ncl, length.out = n)
  dp <- list(data = data.frame(id = 1:n, period = 1L, cl = cl), idname = "id",
             clustervars = "cl", biters = 200, tname = "period", alp = 0.05, panel = TRUE)
  ref_boot_cl <- function(inf.func, cl, n, biters) {
    uc <- unique(cl)
    sapply(1:biters, function(b) {
      Vb <- sample(c(-1, 1), length(uc), replace = TRUE)
      Ub <- Vb[match(cl, uc)]
      Jb <- t(apply(Ub * inf.func, c(2, 3), mean))
      n * sum(apply(Jb^2, 2, mean))
    })
  }
  set.seed(7); ref <- ref_boot_cl(arr, cl, n, 200)
  set.seed(7); got <- test.mboot(arr, dp)$bres
  expect_equal(unname(got), ref, tolerance = 1e-9)
})

test_that("test.mboot multi-chunk tiling accumulates across chunks correctly", {
  # test.mboot tiles over the X dimension with chunk = max(1L, 5e6 %/% (n*k))
  # and accumulates bout across iterations. Every other test in the suite has
  # n*k small enough that the loop runs exactly once; this is the only test
  # where the multi-chunk accumulation (chunk = 1, three iterations) is hit --
  # the regime real large-n pre-test runs land in.
  skip_on_cran()   # ~60 MB transient influence array
  set.seed(123)
  n <- 50; k <- 50010; nX <- 3
  expect_identical(max(1L, as.integer(5e6 %/% (as.numeric(n) * k))), 1L)  # multi-chunk path is exercised
  arr <- array(rnorm(n * k * nX), c(n, k, nX))
  dp <- list(data = data.frame(id = 1:n, period = 1L), idname = "id",
             clustervars = NULL, biters = 50, tname = "period", alp = 0.05, panel = TRUE)
  # Per-draw reference: Jb[x,j] = mean_i(Ub[i] * inf[i,j,x]); CvMb = (n/nX) * sum(Jb^2).
  # colMeans on the n x (k*nX) reshape is the same statistic as the apply-style
  # reference above, just fast enough for k = 50010.
  M <- matrix(arr, nrow = n)
  ref_boot <- function(M, n, nX, biters) {
    sapply(1:biters, function(b) {
      Ub <- sample(c(-1, 1), n, replace = TRUE)
      (n / nX) * sum(colMeans(Ub * M)^2)
    })
  }
  set.seed(42); ref <- ref_boot(M, n, nX, 50)
  set.seed(42); got <- test.mboot(arr, dp)$bres
  expect_equal(unname(got), ref, tolerance = 1e-9)
})
