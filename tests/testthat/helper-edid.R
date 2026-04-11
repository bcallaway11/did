# helper-edid.R
# Shared test data factories for edid() tests.
# Auto-loaded by testthat before any test file runs.

#' Construct a minimal one-cohort balanced panel for testing
#'
#' @param n_treat  Number of treated units (cohort g=3)
#' @param n_never  Number of never-treated units
#' @param n_periods Number of time periods (periods 1..n_periods)
#' @param seed     RNG seed for outcome generation
#' @return data.frame with columns: unit, time, outcome, first_treat
make_panel_1cohort <- function(n_treat = 20L, n_never = 20L,
                                n_periods = 5L, seed = 42L) {
  set.seed(seed)
  n <- n_treat + n_never
  units <- seq_len(n)
  times <- seq_len(n_periods)

  unit_ids  <- rep(units, each = n_periods)
  time_ids  <- rep(times, times = n)
  first_treat_vals <- c(rep(3L, n_treat * n_periods),  # cohort g=3
                        rep(Inf, n_never * n_periods))  # never treated

  # unit fixed effects + time trend + noise
  unit_fe <- rep(rnorm(n, 0, 1), each = n_periods)
  time_fe <- rep(seq(0, 0.5, length.out = n_periods), times = n)
  noise   <- rnorm(n * n_periods, 0, 0.5)
  # Treatment effect of 2 for treated units in post periods
  treated_post <- (unit_ids <= n_treat) & (time_ids >= 3L)
  outcome <- unit_fe + time_fe + noise + 2 * treated_post

  data.frame(
    unit        = unit_ids,
    time        = time_ids,
    outcome     = outcome,
    first_treat = first_treat_vals,
    stringsAsFactors = FALSE
  )
}

#' Construct a two-cohort staggered balanced panel for testing
#'
#' @param n_g3  Units in cohort g=3
#' @param n_g5  Units in cohort g=5
#' @param n_never Never-treated units
#' @param n_periods Time periods (1..n_periods)
#' @param seed RNG seed
#' @return data.frame with columns: unit, time, outcome, first_treat
make_panel_2cohort <- function(n_g3 = 15L, n_g5 = 15L, n_never = 20L,
                                n_periods = 7L, seed = 123L) {
  set.seed(seed)
  n <- n_g3 + n_g5 + n_never
  units <- seq_len(n)
  times <- seq_len(n_periods)

  unit_ids <- rep(units, each = n_periods)
  time_ids <- rep(times, times = n)
  first_treat_vals <- c(
    rep(3L,  n_g3   * n_periods),
    rep(5L,  n_g5   * n_periods),
    rep(Inf, n_never * n_periods)
  )

  unit_fe <- rep(rnorm(n, 0, 1), each = n_periods)
  time_fe <- rep(seq(0, 1, length.out = n_periods), times = n)
  noise   <- rnorm(n * n_periods, 0, 0.5)
  treated_g3 <- (unit_ids <= n_g3) & (time_ids >= 3L)
  treated_g5 <- (unit_ids > n_g3 & unit_ids <= n_g3 + n_g5) & (time_ids >= 5L)
  outcome <- unit_fe + time_fe + noise + 1.5 * treated_g3 + 2.5 * treated_g5

  data.frame(
    unit        = unit_ids,
    time        = time_ids,
    outcome     = outcome,
    first_treat = first_treat_vals,
    stringsAsFactors = FALSE
  )
}

#' Construct a one-cohort panel with cluster variable
#'
#' @param n_clusters_treat Clusters among treated units
#' @param n_clusters_never Clusters among never-treated units
#' @param units_per_cluster Units per cluster
#' @param n_periods Time periods
#' @param seed RNG seed
#' @return data.frame with columns: unit, time, outcome, first_treat, cluster_id
make_panel_clustered <- function(n_clusters_treat = 5L,
                                  n_clusters_never  = 5L,
                                  units_per_cluster = 4L,
                                  n_periods = 5L,
                                  seed = 77L) {
  set.seed(seed)
  n_treat <- n_clusters_treat * units_per_cluster
  n_never <- n_clusters_never * units_per_cluster
  n       <- n_treat + n_never

  units   <- seq_len(n)
  times   <- seq_len(n_periods)
  unit_ids <- rep(units, each = n_periods)
  time_ids <- rep(times, times = n)

  cluster_ids <- c(
    rep(seq_len(n_clusters_treat), each = units_per_cluster * n_periods),
    rep(seq_len(n_clusters_never) + n_clusters_treat, each = units_per_cluster * n_periods)
  )

  first_treat_vals <- c(rep(3L, n_treat * n_periods),
                        rep(Inf, n_never * n_periods))

  cluster_fe <- rep(rnorm(n_clusters_treat + n_clusters_never, 0, 0.8),
                    each = units_per_cluster * n_periods)
  unit_fe    <- rep(rnorm(n, 0, 0.3), each = n_periods)
  noise      <- rnorm(n * n_periods, 0, 0.2)
  treated_post <- (unit_ids <= n_treat) & (time_ids >= 3L)
  outcome <- cluster_fe + unit_fe + noise + 1.8 * treated_post

  data.frame(
    unit        = unit_ids,
    time        = time_ids,
    outcome     = outcome,
    first_treat = first_treat_vals,
    cluster_id  = cluster_ids,
    stringsAsFactors = FALSE
  )
}

#' Minimal two-period, two-group panel (1 treated, 1 never-treated unit)
make_degenerate_panel <- function() {
  data.frame(
    unit        = c(1, 1, 2, 2),
    time        = c(1, 2, 1, 2),
    outcome     = c(0.5, 1.2, 0.3, 0.4),
    first_treat = c(2, 2, Inf, Inf),
    stringsAsFactors = FALSE
  )
}
