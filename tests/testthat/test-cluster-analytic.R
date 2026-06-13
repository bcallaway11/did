# Analytical (no-bootstrap) cluster-robust standard errors. With a cluster variable and bstrap = FALSE,
# att_gt() and aggte() report the cluster-robust variance (cluster sums of the influence function), the
# same cluster-sum aggregation as the multiplier bootstrap (Callaway & Sant'Anna 2021, Remark 10).
# panel in which units within a cluster are dependent: they share a common shock each period (e.g. a
# state-by-year shock), so the cluster -- not the individual unit -- is the independent sampling unit
.make_clustered_shocks <- function(seed, G = 50L) {
  set.seed(seed)
  sz <- rep(c(1L, 2L, 4L, 10L), length.out = G); N <- sum(sz); cl <- rep(seq_len(G), times = sz)
  alpha <- rnorm(G, 0, 1)[cl]; nu <- rnorm(N, 0, 1); per <- 1:4
  eta <- matrix(rnorm(G * 4L, 0, 1.5), G, 4L)                     # common shock per cluster, per period
  g <- sample(c(2L, 3L, 0L), G, replace = TRUE, prob = c(.34, .33, .33))[cl]  # cluster-level treatment
  d <- data.frame(id = rep(seq_len(N), each = 4L), t = rep(per, N), cl = rep(cl, each = 4L),
                  g = rep(g, each = 4L), a = rep(alpha, each = 4L), nu = rep(nu, each = 4L))
  d$y <- d$a + d$nu + 0.3 * d$t + eta[cbind(d$cl, d$t)] + (d$g != 0L & d$t >= d$g) + rnorm(nrow(d))
  d[order(d$id, d$t), ]
}

test_that("analytical (no-bootstrap) clustered SE for att_gt equals the cluster-sum CRVE (faster_mode TRUE and FALSE)", {
  d <- .make_clustered_shocks(404L)
  for (fm in c(TRUE, FALSE)) {
    res <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                  control_group = "nevertreated", bstrap = FALSE, clustervars = "cl",
                  base_period = "varying", faster_mode = fm)
    cv <- res$DIDparams$cluster_vector
    # the per-unit cluster vector must be available in BOTH modes (derived from data in slower mode)
    expect_true(!is.null(cv) && length(cv) == nrow(res$inffunc))
    n <- nrow(res$inffunc)
    S <- rowsum(as.matrix(res$inffunc), cv)            # n_clusters x k cluster sums
    se_target <- sqrt(colSums(S^2)) / n                 # cluster-sum CRVE per ATT(g,t)
    ok <- is.finite(res$se) & is.finite(se_target)
    # analytical => deterministic => exact match (no Monte Carlo tolerance needed)
    expect_equal(unname(res$se[ok]), unname(se_target[ok]), tolerance = 1e-8)

    # and it differs from the i.i.d. SE when the units within a cluster are dependent
    res_iid <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                      control_group = "nevertreated", bstrap = FALSE, base_period = "varying", faster_mode = fm)
    k <- which(res$group == 2L & res$t == 2L)
    expect_gt(abs(res$se[k] - res_iid$se[k]) / res_iid$se[k], 0.05)
  }
})

test_that("analytical cluster SE agrees with the bootstrap cluster SE (multiple DGPs)", {
  skip_on_cran()  # bootstrap-heavy
  # units within a cluster share a common shock; treatment at the cluster level (within = FALSE) or within clusters (TRUE).
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
    set.seed(11L)
    rb <- att_gt(yname="y", tname="t", idname="id", gname="g", data=d, control_group="nevertreated",
                 bstrap=TRUE, biters=3000L, clustervars="cl", pl=FALSE, cband=FALSE, base_period="varying")
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

# true repeated cross-sections: a fresh cross-section is drawn each period within each cluster, the cluster
# shares a common per-period shock, and treatment is assigned at the cluster level. The cluster -- not the
# (non-existent) panel unit -- is the independent sampling unit.
.make_rcs <- function(seed, G = 60L, per_cell = 6L) {
  set.seed(seed); per <- 1:4
  gC <- sample(c(2L, 3L, 0L), G, replace = TRUE, prob = c(.34, .33, .33)); aC <- rnorm(G, 0, 1)
  eta <- matrix(rnorm(G * 4L, 0, 1.5), G, 4L); rows <- list(); k <- 0L
  for (cl in seq_len(G)) for (tt in per) {
    y <- aC[cl] + 0.3 * tt + eta[cl, tt] + (gC[cl] != 0L & tt >= gC[cl]) + rnorm(per_cell)
    k <- k + 1L; rows[[k]] <- data.frame(t = tt, cl = cl, g = gC[cl], y = y)
  }
  do.call(rbind, rows)[sample(G * 4L * per_cell), ]   # shuffle so row order is non-trivial
}

test_that("analytical clustered SE for repeated cross-sections (idname omitted/provided, faster_mode TRUE/FALSE)", {
  d    <- .make_rcs(909L)
  d_id <- d; d_id$uid <- seq_len(nrow(d_id))
  se_om <- list()
  for (fm in c(TRUE, FALSE)) {
    # idname OMITTED: pre_process_did() creates an internal .rowid; clustering must still aggregate to
    # cluster sums via that internal id (the function argument idname is NULL here).
    r_om <- att_gt(yname = "y", tname = "t", gname = "g", data = d, control_group = "nevertreated",
                   bstrap = FALSE, clustervars = "cl", panel = FALSE, base_period = "varying", faster_mode = fm)
    # idname PROVIDED (a per-row observation id)
    r_id <- att_gt(yname = "y", tname = "t", idname = "uid", gname = "g", data = d_id,
                   control_group = "nevertreated", bstrap = FALSE, clustervars = "cl",
                   panel = FALSE, base_period = "varying", faster_mode = fm)
    for (r in list(r_om, r_id)) {
      cv <- r$DIDparams$cluster_vector
      expect_true(!is.null(cv) && length(cv) == nrow(r$inffunc))   # must be present and aligned
      n <- nrow(r$inffunc); S <- rowsum(as.matrix(r$inffunc), cv)
      se_target <- sqrt(colSums(S^2)) / n                          # cluster-sum CRVE per ATT(g,t)
      ok <- is.finite(r$se) & is.finite(se_target)
      expect_equal(unname(r$se[ok]), unname(se_target[ok]), tolerance = 1e-8)
    }
    se_om[[as.character(fm)]] <- r_om$se[which(r_om$group == 2L & r_om$t == 2L)]
  }
  # the two code paths derive the cluster vector independently, yet must give the same analytical SE
  expect_equal(unname(se_om[["TRUE"]]), unname(se_om[["FALSE"]]), tolerance = 1e-6)
})

test_that("clustered bootstrap and analytical SE agree for repeated cross-sections (idname omitted)", {
  skip_on_cran()  # bootstrap-heavy
  d <- .make_rcs(910L)
  for (fm in c(TRUE, FALSE)) {
    ra <- att_gt(yname = "y", tname = "t", gname = "g", data = d, control_group = "nevertreated",
                 bstrap = FALSE, clustervars = "cl", panel = FALSE, base_period = "varying", faster_mode = fm)
    rb <- att_gt(yname = "y", tname = "t", gname = "g", data = d, control_group = "nevertreated",
                 bstrap = TRUE, biters = 5000L, clustervars = "cl", pl = FALSE, cband = FALSE,
                 panel = FALSE, base_period = "varying", faster_mode = fm)
    k <- which(ra$group == 2L & ra$t == 2L)
    expect_true(is.finite(rb$se[k]))
    expect_lt(abs(rb$se[k] - ra$se[k]) / ra$se[k], 0.12)
  }
})

test_that("analytical clustered SE flows through aggte at all levels (simple/group/dynamic)", {
  d <- .make_clustered_shocks(505L)
  res_cl  <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                    control_group = "nevertreated", bstrap = FALSE, clustervars = "cl", base_period = "varying")
  res_iid <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d,
                    control_group = "nevertreated", bstrap = FALSE, base_period = "varying")
  expect_true(!is.null(res_cl$DIDparams$cluster_vector))  # required for clustered SEs to propagate through aggte()
  for (ty in c("simple", "group", "dynamic")) {
    a_cl  <- suppressWarnings(aggte(res_cl,  type = ty, bstrap = FALSE))
    a_iid <- suppressWarnings(aggte(res_iid, type = ty, bstrap = FALSE))
    expect_true(is.finite(a_cl$overall.se) && a_cl$overall.se > 0)
    # clustered overall SE differs from the i.i.d. one when units within a cluster are dependent
    expect_gt(abs(a_cl$overall.se - a_iid$overall.se) / a_iid$overall.se, 0.02)
  }
})
