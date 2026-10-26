test_that("bootstrap degeneracy does not depend on the units of the outcome", {
  set.seed(501)
  inf <- cbind(matrix(rnorm(120 * 3), nrow = 120), 0, NA_real_)
  dp <- list(idname = "id", clustervars = NULL, biters = 199,
             tname = "period", alp = 0.05, panel = TRUE,
             allow_unbalanced_panel = FALSE, faster_mode = FALSE)
  for (clustered in c(FALSE, TRUE)) {
    dp$clustervars <- if (clustered) "cluster" else NULL
    dp$cluster_vector <- if (clustered) rep(seq_len(30), each = 4) else NULL
    set.seed(502)
    reference <- mboot(inf, dp)
    expect_true(all(is.finite(reference$se[1:3])))
    expect_true(all(is.na(reference$se[4:5])))
    for (scale in c(1e-6, 1e-100, 1e6)) {
      set.seed(502)
      scaled <- mboot(inf * scale, dp)
      expect_identical(dim(scaled$bres), dim(reference$bres))
      expect_equal(scaled$se / scale, reference$se, tolerance = 1e-12)
      expect_equal(scaled$bres / scale, reference$bres, tolerance = 1e-12)
      expect_equal(scaled$crit.val, reference$crit.val, tolerance = 1e-12)
    }
  }
})

test_that("rescaling outcomes preserves att_gt and aggte bootstrap inference", {
  set.seed(813)
  d <- expand.grid(id = seq_len(120), t = seq_len(4))
  d$g <- rep(c(0, 2, 3, 4), each = 30)[d$id]
  d$y <- 0.1 * d$t + (d$g > 0 & d$t >= d$g) * d$g + rnorm(nrow(d))
  scale <- 3e-6
  for (faster in c(TRUE, FALSE)) {
    run <- function(multiplier) {
      dd <- d
      dd$y <- dd$y * multiplier
      set.seed(123)
      att_gt("y", "t", "id", "g", data = dd, bstrap = TRUE,
             cband = TRUE, biters = 1000, est_method = "reg",
             faster_mode = faster)
    }
    reference <- run(1)
    scaled <- run(scale)
    # These SEs exceed the separate near-zero-SE cutoff in att_gt/aggte.
    expect_true(all(reference$se * scale > sqrt(.Machine$double.eps) * 10))
    expect_equal(scaled$att / scale, reference$att, tolerance = 1e-10)
    expect_equal(scaled$se / scale, reference$se, tolerance = 1e-10)
    expect_equal(scaled$c, reference$c, tolerance = 1e-10)
    set.seed(503)
    a <- aggte(reference, "simple")
    set.seed(503)
    b <- aggte(scaled, "simple")
    expect_gt(a$overall.se * scale, sqrt(.Machine$double.eps) * 10)
    expect_equal(b$overall.se / scale, a$overall.se, tolerance = 1e-10)
  }
})
