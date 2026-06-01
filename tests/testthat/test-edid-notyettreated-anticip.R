library(testthat)

# Regression test: notyettreated + anticipation must trim periods >= last_g - anticipation, so the
# last cohort (relabeled as the never-treated control) is never used as a control while it is
# anticipating. Before the fix the trim was `< last_g`, keeping the anticipation window [last_g-a, last_g)
# and biasing any cell whose control period falls in it (verified: ATT(3,4) = 0.62 vs true 1.0).

make_nyt_panel <- function(seed = 1, n = 600) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  P  <- exp(cbind(0, 0.5 * x1, 0.3 * x1)); P <- P / rowSums(P)
  g  <- apply(P, 1, function(p) sample(c(Inf, 3, 5), 1, prob = p))
  alpha <- rnorm(n, 0.3 * x1, 1)
  do.call(rbind, lapply(1:6, function(t) {
    tau <- 1 * (is.finite(g) & t >= g)
    data.frame(id = 1:n, tt = t, g = ifelse(is.finite(g), g, 0), x1 = x1,
               y = alpha + 0.3 * t + 0.4 * x1 * (t - 1) + tau + rnorm(n))
  }))
}

test_that("notyettreated trims periods >= last_g - anticipation", {
  df <- make_nyt_panel()
  # last finite cohort = 5; anticipation = 1 -> kept periods must all be < 5 - 1 = 4
  pn1 <- prepare_edid_panel(df, "y", "id", "tt", "g", xformla = ~ x1,
                            control_group = "notyettreated", anticipation = 1L)
  expect_true(all(pn1$time_periods < 4L),
              info = paste("kept periods:", paste(pn1$time_periods, collapse = ",")))
  expect_false(4L %in% pn1$time_periods)   # the anticipation-contaminated period is excluded

  # anticipation = 0 reduces to the original behavior: keep periods < last_g = 5
  pn0 <- prepare_edid_panel(df, "y", "id", "tt", "g", xformla = ~ x1,
                            control_group = "notyettreated", anticipation = 0L)
  expect_true(all(pn0$time_periods < 5L))
  expect_true(4L %in% pn0$time_periods)     # period 4 kept when anticipation = 0
})
