# Analytical (no-bootstrap) cluster-robust standard errors. With a cluster variable and bstrap = FALSE,
# att_gt() and aggte() report the cluster-robust variance (cluster sums of the influence function), the
# same cluster-sum aggregation as the multiplier bootstrap (Callaway & Sant'Anna 2021, Remark 10).
library(testthat)

# panel with cluster-by-time shocks -> genuine within-cluster correlation in the (first-differenced) IF
.make_clustered_shocks <- function(seed, G = 50L) {
  set.seed(seed)
  sz <- rep(c(1L, 2L, 4L, 10L), length.out = G); N <- sum(sz); cl <- rep(seq_len(G), times = sz)
  alpha <- rnorm(G, 0, 1)[cl]; nu <- rnorm(N, 0, 1); per <- 1:4
  eta <- matrix(rnorm(G * 4L, 0, 1.5), G, 4L)                     # cluster-by-time shocks
  g <- sample(c(2L, 3L, 0L), G, replace = TRUE, prob = c(.34, .33, .33))[cl]  # cluster-level treatment
  d <- data.frame(id = rep(seq_len(N), each = 4L), t = rep(per, N), cl = rep(cl, each = 4L),
                  g = rep(g, each = 4L), a = rep(alpha, each = 4L), nu = rep(nu, each = 4L))
  d$y <- d$a + d$nu + 0.3 * d$t + eta[cbind(d$cl, d$t)] + (d$g != 0L & d$t >= d$g) + rnorm(nrow(d))
  d[order(d$id, d$t), ]
}

test_that("analytical (no-bootstrap) clustered SE for att_gt equals the cluster-sum CRVE", {
  d <- .make_clustered_shocks(404L)
  res <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                control_group = "nevertreated", bstrap = FALSE, clustervars = "cl", base_period = "varying")
  cv <- res$DIDparams$cluster_vector
  skip_if_not(!is.null(cv) && length(cv) == nrow(res$inffunc), "cluster_vector not available/aligned")
  n <- nrow(res$inffunc)
  S <- rowsum(as.matrix(res$inffunc), cv)            # n_clusters x k cluster sums
  se_target <- sqrt(colSums(S^2)) / n                 # cluster-sum CRVE per ATT(g,t)
  ok <- is.finite(res$se) & is.finite(se_target)
  # analytical => deterministic => exact match (no Monte Carlo tolerance needed)
  expect_equal(unname(res$se[ok]), unname(se_target[ok]), tolerance = 1e-8)

  # and it differs from the i.i.d. SE under within-cluster correlation
  res_iid <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                    control_group = "nevertreated", bstrap = FALSE, base_period = "varying")
  k <- which(res$group == 2L & res$t == 2L)
  expect_gt(abs(res$se[k] - res_iid$se[k]) / res_iid$se[k], 0.05)
})

test_that("analytical cluster SE agrees with the bootstrap cluster SE (multiple DGPs)", {
  skip_on_cran()  # bootstrap-heavy
  # cluster-by-time shocks; treatment at the cluster level (within = FALSE) or within clusters (TRUE).
  mk <- function(seed, within, G = 60L) {
    set.seed(seed); sz <- rep(c(1L, 2L, 4L, 10L), length.out = G); N <- sum(sz); cl <- rep(seq_len(G), times = sz)
    alpha <- rnorm(G, 0, 1)[cl]; nu <- rnorm(N, 0, 1); eta <- matrix(rnorm(G * 4L, 0, 1.5), G, 4L)
    g <- if (within) sample(c(2L, 3L, 0L), N, replace = TRUE, prob = c(.34, .33, .33))
         else sample(c(2L, 3L, 0L), G, replace = TRUE, prob = c(.34, .33, .33))[cl]
    d <- data.frame(id = rep(seq_len(N), each = 4L), t = rep(1:4, N), cl = rep(cl, each = 4L),
                    g = rep(g, each = 4L), a = rep(alpha, each = 4L), nu = rep(nu, each = 4L))
    d$y <- d$a + d$nu + 0.3 * d$t + eta[cbind(d$cl, d$t)] + (d$g != 0L & d$t >= d$g) + rnorm(nrow(d))
    d[order(d$id, d$t), ]
  }
  for (within in c(FALSE, TRUE)) {
    d <- mk(11L, within)
    ra <- att_gt(yname="y", tname="t", idname="id", gname="g", data=d, control_group="nevertreated",
                 bstrap=FALSE, clustervars="cl", base_period="varying")
    rb <- att_gt(yname="y", tname="t", idname="id", gname="g", data=d, control_group="nevertreated",
                 bstrap=TRUE, biters=3000L, clustervars="cl", pl=FALSE, cband=FALSE, base_period="varying", seed=11L)
    k <- which(ra$group == 2L & ra$t == 2L)
    # the reported bootstrap SE (IQR scale) is within a few percent of the analytical at this G
    expect_lt(abs(rb$se[k] - ra$se[k]) / ra$se[k], 0.15)
    # the bootstrap's own SD scale matches the analytical closely (same cluster variance)
    dp <- ra$DIDparams; dp$biters <- 3000L; dp$pl <- FALSE
    set.seed(11L); bo <- mboot(as.matrix(ra$inffunc)[, k, drop = FALSE], dp)
    sd_se <- sd(as.numeric(bo$bres)) * sqrt(length(unique(dp$cluster_vector))) / nrow(ra$inffunc)
    expect_lt(abs(sd_se - ra$se[k]) / ra$se[k], 0.06)
  }
})

test_that("analytical clustered SE flows through aggte at all levels (simple/group/dynamic)", {
  d <- .make_clustered_shocks(505L)
  res_cl  <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                    control_group = "nevertreated", bstrap = FALSE, clustervars = "cl", base_period = "varying")
  res_iid <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                    control_group = "nevertreated", bstrap = FALSE, base_period = "varying")
  skip_if_not(!is.null(res_cl$DIDparams$cluster_vector), "cluster_vector not available")
  for (ty in c("simple", "group", "dynamic")) {
    a_cl  <- suppressWarnings(aggte(res_cl,  type = ty, bstrap = FALSE))
    a_iid <- suppressWarnings(aggte(res_iid, type = ty, bstrap = FALSE))
    expect_true(is.finite(a_cl$overall.se) && a_cl$overall.se > 0)
    # clustered overall SE differs from the i.i.d. one under within-cluster correlation
    expect_gt(abs(a_cl$overall.se - a_iid$overall.se) / a_iid$overall.se, 0.02)
  }
})
