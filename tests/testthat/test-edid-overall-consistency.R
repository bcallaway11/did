# Regression: `$overall` is ALWAYS the dynamic event-study average (the average of the
# post-treatment event study), identical across aggregate = "all" / "event_study" / "overall".
# Previously aggregate = "overall" alone returned the cohort-share "simple" aggregate -- a
# different number. (overnight-aggfix-robust)

make_stag_panel <- function(n_per = 60L, Tn = 8L, seed = 1L) {
  set.seed(seed)
  cohorts <- c(0L, 4L, 6L)                      # never, g=4, g=6
  G <- rep(cohorts, each = n_per); I <- length(G)
  lam <- rnorm(I); del <- rnorm(Tn)
  id <- rep(seq_len(I), each = Tn); tt <- rep(seq_len(Tn), times = I)
  g  <- rep(G, each = Tn)
  tau <- ifelse(g != 0 & tt >= g, 0.4 * (tt - g + 1), 0)      # dynamic (grows with e): overall != simple
  y <- lam[id] + del[tt] + tau + rnorm(I * Tn)
  data.frame(id = id, t = tt, g = g, y = y)
}

test_that("aggregate='overall' returns the dynamic event-study average (= 'all' = 'event_study')", {
  df <- make_stag_panel()
  f_all <- edid(df, "y", "id", "t", "g", aggregate = "all",         bstrap = FALSE)
  f_es  <- edid(df, "y", "id", "t", "g", aggregate = "event_study", bstrap = FALSE)
  f_ov  <- edid(df, "y", "id", "t", "g", aggregate = "overall",     bstrap = FALSE)

  # headline overall identical across the three request modes
  expect_equal(f_ov$overall$overall.att, f_all$overall$overall.att, tolerance = 1e-10)
  expect_equal(f_ov$overall$overall.att, f_es$overall$overall.att,  tolerance = 1e-10)
  expect_equal(f_ov$overall$overall.se,  f_all$overall$overall.se,  tolerance = 1e-10)

  # overall carries the event-study breakdown (it IS the dynamic aggregation)
  expect_false(is.null(f_ov$overall$att.egt))
  expect_true(all(is.finite(f_ov$overall$att.egt)))

  # $simple still computed under "overall", and is a DIFFERENT estimand here (dynamic effects)
  expect_false(is.null(f_ov$simple))
  expect_false(isTRUE(all.equal(f_ov$overall$overall.att, f_ov$simple$overall.att)))

  # overall == the equal-weight mean of the post-treatment event-study coefficients
  egt <- f_ov$overall$egt; att <- f_ov$overall$att.egt
  expect_equal(f_ov$overall$overall.att, mean(att[egt >= 0]), tolerance = 1e-8)
})

test_that("calendar-/group-only requests still leave $overall NULL (shape-stable)", {
  df <- make_stag_panel()
  expect_null(edid(df, "y", "id", "t", "g", aggregate = "calendar", bstrap = FALSE)$overall)
  expect_null(edid(df, "y", "id", "t", "g", aggregate = "group",    bstrap = FALSE)$overall)
})
