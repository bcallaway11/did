# Cluster-robust multiplier bootstrap (Callaway & Sant'Anna 2021, Remark 10): the clustered standard
# error aggregates the influence function to cluster *sums*. For equal-sized clusters this coincides with
# the previous cluster-mean aggregation; for unbalanced clusters and repeated cross-sections it follows
# the cluster-sum form that aligns with the (1/n) empirical average defining the estimator.
library(testthat)

# small staggered panel with within-cluster correlation; bal controls equal vs unequal cluster sizes
.make_clustered <- function(seed, bal, G = 40L) {
  set.seed(seed)
  sz <- if (bal) rep(4L, G) else rep(c(1L, 2L, 3L, 8L), length.out = G)
  N <- sum(sz); cl <- rep(seq_len(G), times = sz)
  a <- rnorm(G, 0, 1)[cl]; nu <- rnorm(N, 0, 1); g <- ifelse(runif(N) < 0.5, 2L, 0L)
  d <- data.frame(id = rep(seq_len(N), each = 3L), t = rep(1:3, N), cl = rep(cl, each = 3L),
                  g = rep(g, each = 3L), a = rep(a, each = 3L), nu = rep(nu, each = 3L))
  d$y <- d$a + d$nu + 0.5 * d$t + (d$g == 2L & d$t >= 2L) + rnorm(nrow(d))
  d[order(d$id, d$t), ]
}

.cluster_targets <- function(res) {  # analytical cluster-SUM and cluster-MEAN for att(2,2)
  k <- which(res$group == 2L & res$t == 2L)
  cv <- res$DIDparams$cluster_vector
  inf <- as.matrix(res$inffunc)[, k]; n <- nrow(res$inffunc)
  S <- as.numeric(rowsum(matrix(inf, ncol = 1L), cv))
  nc <- as.numeric(table(cv)); Gc <- length(S)
  list(k = k, sum = sqrt(sum(S^2)) / n, mean = sqrt(sum((S / nc)^2)) / Gc, ok = length(cv) == n)
}

test_that("clustered mboot SE matches the cluster-sum (Remark 10) for UNBALANCED clusters", {
  d <- .make_clustered(101L, bal = FALSE)
  res <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                control_group = "nevertreated", bstrap = TRUE, biters = 5000L,
                clustervars = "cl", pl = FALSE, cband = FALSE, base_period = "varying")
  tg <- .cluster_targets(res)
  skip_if_not(tg$ok, "cluster_vector not aligned with influence function")
  # cluster-sum and cluster-mean genuinely differ when clusters are unbalanced
  expect_gt(abs(tg$sum - tg$mean) / tg$mean, 0.05)
  # reported SE tracks the cluster-SUM target (within bootstrap Monte Carlo tolerance) ...
  expect_equal(unname(res$se[tg$k]), tg$sum, tolerance = 0.08)
  # ... and is clearly NOT the cluster-mean
  expect_gt(abs(res$se[tg$k] - tg$mean) / tg$mean, 0.05)
})

test_that("clustered mboot SE is unchanged for BALANCED clusters (cluster-sum == cluster-mean)", {
  d <- .make_clustered(202L, bal = TRUE)
  res <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                control_group = "nevertreated", bstrap = TRUE, biters = 5000L,
                clustervars = "cl", pl = FALSE, cband = FALSE, base_period = "varying")
  tg <- .cluster_targets(res)
  skip_if_not(tg$ok, "cluster_vector not aligned with influence function")
  expect_equal(tg$sum, tg$mean, tolerance = 1e-8)            # identical targets when balanced
  expect_equal(unname(res$se[tg$k]), tg$sum, tolerance = 0.08)
})

test_that("clustering validation is preserved (at most one cluster variable beyond idname)", {
  d <- .make_clustered(303L, bal = FALSE); d$cl2 <- d$cl
  # two clustering variables (neither is idname) must still be rejected after the cluster-sum change
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
           control_group = "nevertreated", bstrap = TRUE, biters = 99L,
           clustervars = c("cl", "cl2"), pl = FALSE, cband = FALSE, base_period = "varying"),
    regexp = "length 1|character scalar|clustervars|[Aa]t most one cluster"
  )
})
