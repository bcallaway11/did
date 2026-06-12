# Regression tests for two faster_mode standard-error bugs on UNBALANCED panels with clustering.
#
# (1) att_gt(): the analytical clustered SE (bstrap = FALSE, clustervars coarser than idname) silently
#     fell back to the i.i.d. SE under faster_mode = TRUE, because the stored DIDparams$cluster_vector was
#     observation-length (the unbalanced time-invariant data keeps every row) rather than unit-length, so
#     it no longer aligned with the influence function. It is now rebuilt from the data to align with the
#     influence-function rows whenever it is missing or mis-sized.
#
# (2) aggte(): the estimated-weight influence term (wif) was built from aggregate(data, by id, mean),
#     which is sorted by id, but the influence function is in the data's first-appearance unit order.
#     get_agg_inf_func() adds the two row-wise without realigning, so under faster_mode = TRUE on an
#     unbalanced panel the wif was added to the WRONG units and the aggregated SE was wrong (the point
#     estimate is unaffected because wif is mean-zero). The aggregation data is now ordered to match the
#     influence function.
#
# In both cases faster_mode = TRUE must reproduce faster_mode = FALSE (the correct reference) to
# numerical precision, for balanced and unbalanced panels alike.
library(testthat)

# Units nested in unbalanced clusters that share a common shock each period (so the cluster, not the unit,
# is the independent sampling unit). Optionally drop rows to create an unbalanced panel and shuffle the
# row order, so that the internal first-appearance unit order differs from the sorted-by-id order -- the
# condition that triggered both bugs under faster_mode.
.make_ub_clustered <- function(seed, G = 40L, unbalanced = TRUE) {
  set.seed(seed)
  sz <- 1L + rpois(G, 4); N <- sum(sz); cl <- rep(seq_len(G), times = sz)
  alpha <- rnorm(G, 0, 1)[cl]; nu <- rnorm(N, 0, 1); Tt <- 6L
  eta <- matrix(rnorm(G * Tt, 0, 1), G, Tt)              # common shock per cluster, per period
  g <- sample(c(3L, 4L, 5L, 0L), N, replace = TRUE, prob = c(.2, .2, .2, .4))
  d <- data.frame(id = rep(seq_len(N), each = Tt), period = rep(seq_len(Tt), N),
                  cluster = rep(cl, each = Tt), g = rep(g, each = Tt),
                  a = rep(alpha, each = Tt), nu = rep(nu, each = Tt))
  d$y <- d$a + d$nu + 0.5 * d$period + eta[cbind(d$cluster, d$period)] +
    2 * (d$g != 0L & d$period >= d$g) + rnorm(nrow(d))
  d <- d[order(d$id, d$period), ]
  if (unbalanced) d <- d[sample(nrow(d), floor(0.85 * nrow(d))), ]   # unbalanced + shuffled
  d
}

.fit <- function(d, fm, clus, unb = TRUE) {
  suppressWarnings(suppressMessages(att_gt(
    yname = "y", tname = "period", idname = "id", gname = "g", data = d,
    control_group = "nevertreated", base_period = "universal", bstrap = FALSE, cband = FALSE,
    allow_unbalanced_panel = unb, faster_mode = fm,
    clustervars = if (clus) "cluster" else NULL)))
}

test_that("att_gt analytical clustered SE: faster_mode works (not i.i.d. fallback) on an unbalanced panel", {
  d <- .make_ub_clustered(202L)
  cT <- .fit(d, TRUE,  TRUE)
  cF <- .fit(d, FALSE, TRUE)
  iT <- .fit(d, TRUE,  FALSE)   # i.i.d. (no clustervars)

  # cluster vector now aligns with the influence function under faster_mode (the bug: length mismatch)
  cv <- cT$DIDparams$cluster_vector
  expect_true(!is.null(cv) && length(cv) == nrow(cT$inffunc))

  # analytical clustered SE equals the cluster-sum CRVE closed form (1/n) sqrt(sum_c S_c^2)
  n <- nrow(cT$inffunc)
  S <- rowsum(as.matrix(cT$inffunc), cv)
  se_cf <- sqrt(colSums(S^2)) / n
  ok <- is.finite(cT$se) & is.finite(se_cf)
  expect_equal(unname(cT$se[ok]), unname(se_cf[ok]), tolerance = 1e-8)

  # clustering is actually applied: differs from the i.i.d. SE (was a silent fallback to i.i.d.)
  expect_gt(max(abs(cT$se - iT$se), na.rm = TRUE), 1e-6)

  # faster_mode TRUE == FALSE
  expect_equal(unname(cT$se), unname(cF$se), tolerance = 1e-8)
})

test_that("aggte: faster_mode TRUE == FALSE on an unbalanced clustered panel (all aggregation types)", {
  d <- .make_ub_clustered(202L)
  for (clus in c(TRUE, FALSE)) {
    mT <- .fit(d, TRUE,  clus)
    mF <- .fit(d, FALSE, clus)
    for (ty in c("simple", "group", "dynamic", "calendar")) {
      aT <- suppressWarnings(suppressMessages(aggte(mT, type = ty, bstrap = FALSE, cband = FALSE)))
      aF <- suppressWarnings(suppressMessages(aggte(mF, type = ty, bstrap = FALSE, cband = FALSE)))
      lab <- paste(ty, if (clus) "clustered" else "iid")
      expect_equal(aT$overall.att, aF$overall.att, tolerance = 1e-10, label = paste(lab, "overall.att"))
      expect_equal(aT$overall.se,  aF$overall.se,  tolerance = 1e-8,  label = paste(lab, "overall.se"))
      expect_equal(unname(aT$se.egt), unname(aF$se.egt), tolerance = 1e-8, label = paste(lab, "se.egt"))
    }
  }
})

test_that("no regression: balanced panel still has faster_mode TRUE == FALSE (att_gt + aggte, clustered)", {
  d <- .make_ub_clustered(55L, unbalanced = FALSE)
  cT <- .fit(d, TRUE,  TRUE, unb = FALSE)
  cF <- .fit(d, FALSE, TRUE, unb = FALSE)
  expect_equal(unname(cT$se), unname(cF$se), tolerance = 1e-8)
  for (ty in c("simple", "group", "dynamic", "calendar")) {
    aT <- suppressWarnings(suppressMessages(aggte(cT, type = ty, bstrap = FALSE, cband = FALSE)))
    aF <- suppressWarnings(suppressMessages(aggte(cF, type = ty, bstrap = FALSE, cband = FALSE)))
    expect_equal(aT$overall.se, aF$overall.se, tolerance = 1e-8, label = paste(ty, "overall.se"))
    expect_equal(unname(aT$se.egt), unname(aF$se.egt), tolerance = 1e-8, label = paste(ty, "se.egt"))
  }
})
