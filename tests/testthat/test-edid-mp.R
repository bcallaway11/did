# did-compatible MP construction: as_MP_edid() lets aggte() aggregate edid output unchanged.
library(testthat)

.mk_edid_panel <- function(seed, n = 400L, Tn = 4L) {
  set.seed(seed)
  g <- sample(c(0L, 2L, 3L), n, replace = TRUE, prob = c(.4, .3, .3))
  a <- rnorm(n); rows <- vector("list", Tn)
  for (tt in 1:Tn) {
    tau <- ifelse(g != 0L & tt >= g, 1 + 0.3 * (tt - g), 0)
    y <- a + 0.3 * tt + tau + rnorm(n)
    rows[[tt]] <- data.frame(id = seq_len(n), tt = tt, g = g, y = y)
  }
  do.call(rbind, rows)
}

test_that("as_MP_edid() builds a did MP and aggte aggregates edid cells (all types)", {
  d <- .mk_edid_panel(101L)
  fit <- edid(yname = "y", tname = "tt", idname = "id", gname = "g", data = d,
              aggregate = "none")
  mp <- as_MP_edid(fit)
  expect_s3_class(mp, "MP")
  expect_equal(length(mp$group), nrow(fit$att_gt))
  expect_equal(dim(mp$inffunc), c(fit$n, nrow(fit$att_gt)))
  for (ty in c("simple", "group", "dynamic", "calendar")) {
    a <- suppressWarnings(aggte(mp, type = ty, bstrap = FALSE, na.rm = TRUE))
    expect_true(is.finite(a$overall.att) && is.finite(a$overall.se) && a$overall.se > 0)
  }
})

test_that("aggte(dynamic) on the edid MP reproduces edid's per-cell event-time ATT", {
  d <- .mk_edid_panel(202L)
  fit <- edid(yname = "y", tname = "tt", idname = "id", gname = "g", data = d,
              aggregate = "none")
  ad  <- suppressWarnings(aggte(as_MP_edid(fit), type = "dynamic", bstrap = FALSE, na.rm = TRUE))
  # e = 2 has a single contributing cell (cohort 2 at time 4), so the event-study ATT equals that cell's
  agt <- fit$att_gt; k <- which(agt$group == 2 & agt$time == 4)
  e2  <- which(ad$egt == 2)
  expect_equal(unname(ad$att.egt[e2]), unname(agt$att[k]), tolerance = 1e-8)
})
