# aggte()'s clustervars argument can only honor clustering that att_gt() actually set up (att_gt retains
# the cluster column and builds the per-unit cluster_vector only for the clustervars passed to it). When
# aggte() is asked to cluster on a variable that att_gt() did not use, the cluster information is
# unavailable, so aggte() must NOT silently return the i.i.d. SE (analytic path) or error in mboot()
# (bootstrap path) -- it warns and falls back to non-clustered standard errors. Inheriting clustering from
# the att_gt object, or overriding to the same variable att_gt used, work normally and emit no warning.
library(testthat)

.mk_two_cluster <- function(seed) {
  set.seed(seed)
  G <- 30L; sz <- rep(c(2L, 4L, 6L), length.out = G); N <- sum(sz); cl <- rep(seq_len(G), sz)
  region <- rep(sample(1:5, G, TRUE), times = sz)            # a second time-invariant grouping
  alpha <- rnorm(G, 0, 1)[cl]; nu <- rnorm(N); Tt <- 5L
  eta <- matrix(rnorm(G * Tt, 0, 1.5), G, Tt)                # shared within-cluster shock
  g <- sample(c(0L, 2L, 3L, 4L), N, TRUE, c(.4, .2, .2, .2))
  d <- data.frame(id = rep(seq_len(N), each = Tt), period = rep(seq_len(Tt), N),
                  cluster = rep(cl, each = Tt), region = rep(region, each = Tt),
                  g = rep(g, each = Tt), a = rep(alpha, each = Tt), nu = rep(nu, each = Tt))
  d$y <- d$a + d$nu + 0.4 * d$period + eta[cbind(d$cluster, d$period)] +
    (d$g != 0L & d$period >= d$g) + rnorm(nrow(d))
  d[order(d$id, d$period), ]
}

# did aggte() emit the "clustering unavailable" warning?
.warned <- function(expr) {
  flag <- FALSE
  withCallingHandlers(force(expr), warning = function(c) {
    if (grepl("cluster information needed is not available", conditionMessage(c))) flag <<- TRUE
    invokeRestart("muffleWarning")
  })
  flag
}
.agg <- function(m, ty, ...) suppressMessages(aggte(m, type = ty, bstrap = FALSE, cband = FALSE, ...))

test_that("inherited clustering and same-variable override are honored without warning", {
  d <- .mk_two_cluster(11)
  mc <- suppressWarnings(suppressMessages(att_gt(yname = "y", tname = "period", idname = "id", gname = "g",
        data = d, clustervars = "cluster", control_group = "nevertreated", bstrap = FALSE)))
  for (ty in c("simple", "dynamic")) {
    w_inh <- .warned(ai <- .agg(mc, ty))
    w_exp <- .warned(ae <- .agg(mc, ty, clustervars = "cluster"))
    expect_false(w_inh)
    expect_false(w_exp)
    expect_equal(ai$overall.se, ae$overall.se, tolerance = 1e-10)
  }
})

test_that("aggte warns and falls back to i.i.d. when att_gt was not clustered (analytic)", {
  d <- .mk_two_cluster(11)
  mi <- suppressWarnings(suppressMessages(att_gt(yname = "y", tname = "period", idname = "id", gname = "g",
        data = d, control_group = "nevertreated", bstrap = FALSE)))
  for (ty in c("simple", "dynamic", "group", "calendar")) {
    a_iid  <- .agg(mi, ty)                                    # no clustervars -> i.i.d.
    w <- .warned(a_over <- .agg(mi, ty, clustervars = "cluster"))
    expect_true(w)                                            # warning fired
    expect_equal(a_over$overall.se, a_iid$overall.se, tolerance = 1e-10)   # degraded to i.i.d.
  }
})

test_that("aggte bootstrap override warns and falls back instead of erroring", {
  d <- .mk_two_cluster(11)
  mi <- suppressWarnings(suppressMessages(att_gt(yname = "y", tname = "period", idname = "id", gname = "g",
        data = d, control_group = "nevertreated", bstrap = TRUE, biters = 200)))
  ab <- NULL
  w <- .warned(ab <- suppressMessages(aggte(mi, type = "dynamic", bstrap = TRUE, biters = 200,
                                            cband = FALSE, clustervars = "cluster")))
  expect_true(w)
  expect_false(is.null(ab))
  expect_true(all(is.finite(ab$se.egt)))                      # no "undefined columns selected" crash
})

test_that("aggte refuses to switch to a different cluster variable than att_gt used", {
  d <- .mk_two_cluster(11)
  mc <- suppressWarnings(suppressMessages(att_gt(yname = "y", tname = "period", idname = "id", gname = "g",
        data = d, clustervars = "cluster", control_group = "nevertreated", bstrap = FALSE)))
  a_inh <- .agg(mc, "dynamic")                                # correctly clustered on 'cluster'
  w <- .warned(a_sw <- .agg(mc, "dynamic", clustervars = "region"))
  expect_true(w)                                              # cannot switch to 'region'
  expect_false(isTRUE(all.equal(a_sw$overall.se, a_inh$overall.se)))   # degraded -> differs from clustered
})
