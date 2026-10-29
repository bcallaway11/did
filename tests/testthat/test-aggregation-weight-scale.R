# Sampling weights specify relative masses: their common scale cannot change
# either the estimated aggregation shares or their empirical-law derivatives.
test_that("weight influence functions match ratio derivatives at extreme scales", {
  G <- rep(c(0, 2, 3, 4), each = 4)
  weights <- c(0, seq_len(15) / 10)
  group <- c(2, 2, 3, 3, 4)
  keepers <- c(1, 3, 4, 5) # repeated cohort masses occur in simple aggregation
  membership <- outer(G, group[keepers], `==`) * 1
  counts <- rowSums(membership)
  total <- sum(weights * counts)
  shares <- colSums(weights * membership) / total

  # Direct ratio derivative, without centering cohort means or squaring them.
  expected <- length(G) * weights *
    (membership - counts %o% shares) / total
  h <- 1e-5
  for (i in seq_along(G)) {
    perturbed.shares <- function(delta) {
      mass <- weights
      mass[i] <- mass[i] * (1 + delta)
      colSums(mass * membership) / sum(mass * counts)
    }
    derivative <- (perturbed.shares(h) - perturbed.shares(-h)) /
      (h / (length(G) + h) + h / (length(G) - h))
    expect_equal(derivative, expected[i, ], tolerance = 1e-8)
  }

  for (scale in c(1e-200, 1e-160, 1, 1e160, 1e200, 1e307)) {
    w <- weights * scale
    pg <- vapply(group, function(g) mean(w * (G == g)), numeric(1))
    actual <- did:::wif(keepers, pg, w, G, group)
    expect_true(all(is.finite(actual)))
    expect_equal(actual, expected, tolerance = 1e-12,
                 info = format(scale, scientific = TRUE))
  }
})

test_that("weight influence functions preserve ordinary arithmetic and invalid inputs", {
  old.wif <- function(keepers, pg, weights.ind, G, group) {
    Spg <- sum(pg[keepers])
    centered <- sapply(keepers, function(k) {
      weights.ind * 1 * BMisc::TorF(G == group[k]) - pg[k]
    })
    centered / Spg - rowSums(centered) %*% t(pg[keepers] / Spg^2)
  }
  G <- rep(c(0, 2, 3), 4)
  group <- c(2, 2, 3)
  for (scale in c(1e-100, 1, 1e100)) {
    weights <- seq_along(G) * scale
    pg <- vapply(group, function(g) mean(weights * (G == g)), numeric(1))
    expect_identical(did:::wif(1:3, pg, weights, G, group),
                     old.wif(1:3, pg, weights, G, group))
  }
  for (pg in list(c(0, 0, 0), c(Inf, 1, 1), c(-1, -1, 2))) {
    expect_identical(did:::wif(1:3, pg, rep(1, length(G)), G, group),
                     old.wif(1:3, pg, rep(1, length(G)), G, group))
  }
})

test_that("weight influence functions handle overflow of repeated selected masses", {
  G <- rep(c(0, 2, 3), each = 8)
  weights <- rep(1, length(G))
  group <- rep(c(2, 3), 40)
  pg <- vapply(group, function(g) mean(weights * (G == g)), numeric(1))
  expected <- did:::wif(seq_along(group), pg, weights, G, group)
  # sum(pg * scale) exceeds the largest finite double, while all supplied
  # weights and individual cohort masses remain finite.
  scale <- 1e307
  expect_true(is.infinite(sum(pg * scale)))
  actual <- did:::wif(seq_along(group), pg * scale, weights * scale, G, group)
  expect_true(all(is.finite(actual)))
  expect_equal(actual, expected, tolerance = 1e-12)
})

.aggregation_weight_scale_data <- function(panel) {
  set.seed(20261027)
  n <- 160L
  units <- data.frame(id = seq_len(n), G = rep(c(0, 2, 3, 4), each = n / 4),
                      X = rnorm(n), a = rnorm(n), weight = 0.5 + (seq_len(n) %% 7) / 10)
  d <- merge(units, data.frame(period = 1:4))
  d$Y <- d$a + 0.2 * d$period + 0.1 * d$X * d$period +
    (d$G > 0 & d$period >= d$G) * (0.5 + 0.2 * d$G) * (d$period - d$G + 1) +
    0.1 * rnorm(nrow(d))
  if (!panel) d$id <- seq_len(nrow(d))
  d
}

test_that("all public aggregations preserve full influence functions across weight scales", {
  scales <- c(1e-200, 1e-160, 1e160, 1e200, 1e307)
  types <- c("simple", "group", "dynamic", "calendar")
  for (panel in c(TRUE, FALSE)) {
    d <- .aggregation_weight_scale_data(panel)
    for (method in c("dr", "ipw", "reg")) {
      reference <- NULL
      for (faster in c(FALSE, TRUE)) {
        fit <- function(scale) {
          dd <- d
          dd$weight <- dd$weight * scale
          suppressMessages(att_gt(yname = "Y", tname = "period", idname = "id",
            gname = "G", data = dd, panel = panel, xformla = ~X,
            weightsname = "weight", est_method = method, faster_mode = faster,
            bstrap = FALSE, cband = FALSE))
        }
        ref <- fit(1)
        ref.agg <- lapply(types, function(type) aggte(ref, type = type))
        if (is.null(reference)) reference <- ref.agg
        for (scale in c(1, scales)) {
          out <- if (scale == 1) ref else fit(scale)
          info <- paste(panel, method, faster, format(scale, scientific = TRUE))
          expect_equal(out$att, ref$att, tolerance = 1e-9, info = info)
          expect_equal(out$inffunc, ref$inffunc, tolerance = 1e-9, info = info)
          expect_equal(out$se, ref$se, tolerance = 1e-9, info = info)
          for (j in seq_along(types)) {
            a <- aggte(out, type = types[j])
            expected <- reference[[j]]
            expect_true(is.finite(a$overall.se) && a$overall.se > 0, info = info)
            expect_equal(a$overall.att, expected$overall.att, tolerance = 1e-9, info = info)
            expect_equal(a$overall.se, expected$overall.se, tolerance = 1e-9, info = info)
            # Fast and slow fits store individuals in different row orders;
            # compare every IF entry to the same mode's ordinary-scale fit.
            expect_equal(a$inf.function, ref.agg[[j]]$inf.function,
                         tolerance = 1e-9, info = info)
            if (types[j] != "simple") {
              expect_equal(a$att.egt, expected$att.egt, tolerance = 1e-9, info = info)
              expect_equal(a$se.egt, expected$se.egt, tolerance = 1e-9, info = info)
            }
          }
        }
      }
    }
  }
})
