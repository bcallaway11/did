# edid-data.R
# Panel preparation and cluster indexing for the EDiD estimator.

#' Prepare the panel object used throughout edid estimation
#'
#' Reshapes the long-format input \code{data} into a wide outcome matrix and
#' builds all masks and maps needed by downstream functions.
#'
#' @param data data.frame (or data.table / tibble) already validated
#' @param outcome character scalar: outcome column name
#' @param unit character scalar: unit id column name
#' @param time character scalar: time column name
#' @param first_treat character scalar: first-treatment-period column name
#' @param covariates NULL (stub)
#' @param cluster character scalar or NULL
#' @param control_group \code{"never_treated"} or \code{"last_cohort"}
#' @param anticipation non-negative integer
#'
#' @return a \code{panel_obj} list; see spec Section 5.1
#' @keywords internal
prepare_edid_panel <- function(
  data, outcome, unit, time, first_treat,
  covariates = NULL, cluster = NULL,
  control_group = "never_treated",
  anticipation = 0L
) {

  # -----------------------------------------------------------------------
  # 1. Coerce to data.table and sort
  # -----------------------------------------------------------------------
  dt <- data.table::as.data.table(data)
  data.table::setkeyv(dt, c(unit, time))

  # -----------------------------------------------------------------------
  # 2-3. Extract sorted unique ids and time periods
  # -----------------------------------------------------------------------
  all_units   <- sort(unique(dt[[unit]]))
  time_periods <- sort(unique(dt[[time]]))
  n           <- length(all_units)
  T_periods   <- length(time_periods)

  # -----------------------------------------------------------------------
  # 4. period_1 and period_to_col map
  # -----------------------------------------------------------------------
  period_1     <- time_periods[1L]
  period_to_col <- stats::setNames(
    seq_along(time_periods),
    as.character(time_periods)
  )

  # -----------------------------------------------------------------------
  # 5. Pivot to wide outcome matrix (n x T_periods)
  #    Rows correspond to all_units (sorted), columns to time_periods (sorted)
  # -----------------------------------------------------------------------
  wide_dt <- data.table::dcast(
    dt,
    formula = stats::as.formula(paste(unit, "~ ", time)),
    value.var = outcome
  )
  # Ensure rows in same order as all_units
  setattr <- function(x, nm, val) { attr(x, nm) <- val; x }
  unit_order <- match(all_units, wide_dt[[unit]])
  wide_dt    <- wide_dt[unit_order, ]

  # Drop the unit id column; keep only the T_periods outcome columns
  # Column names after dcast are as.character(time_periods)
  col_order   <- as.character(time_periods)
  outcome_wide <- as.matrix(wide_dt[, col_order, with = FALSE])
  rownames(outcome_wide) <- NULL
  colnames(outcome_wide) <- col_order

  # -----------------------------------------------------------------------
  # 6. unit_cohorts: first_treat value per unit (Inf for never-treated)
  # -----------------------------------------------------------------------
  # Extract one first_treat per unit using base R tapply (avoids data.table NSE)
  ft_vals      <- dt[[first_treat]]
  unit_id_vals <- dt[[unit]]
  # Get first value of first_treat per unit (treatment is constant within unit)
  unit_ft_map  <- tapply(ft_vals, unit_id_vals, function(x) x[1L])
  # Map to all_units order
  unit_cohorts <- as.numeric(unit_ft_map[match(all_units, names(unit_ft_map))])

  # -----------------------------------------------------------------------
  # 7. Handle last_cohort control group
  # -----------------------------------------------------------------------
  if (control_group == "last_cohort") {
    finite_cohorts <- unit_cohorts[is.finite(unit_cohorts)]
    last_g         <- max(finite_cohorts)
    # Relabel last cohort as Inf (never-treated for estimation purposes)
    unit_cohorts[unit_cohorts == last_g] <- Inf
    # Trim time periods >= last_g
    keep_times  <- time_periods[time_periods < last_g]
    keep_cols   <- as.character(keep_times)
    outcome_wide <- outcome_wide[, keep_cols, drop = FALSE]
    time_periods <- keep_times
    T_periods    <- length(time_periods)
    period_1     <- time_periods[1L]
    period_to_col <- stats::setNames(
      seq_along(time_periods),
      as.character(time_periods)
    )
  }

  # -----------------------------------------------------------------------
  # 8. treatment_groups: sorted unique finite cohort values
  # -----------------------------------------------------------------------
  treatment_groups <- sort(unique(unit_cohorts[is.finite(unit_cohorts)]))

  # -----------------------------------------------------------------------
  # 9. cohort_masks: named list, one logical vector per cohort
  # -----------------------------------------------------------------------
  cohort_masks <- vector("list", length(treatment_groups))
  names(cohort_masks) <- as.character(treatment_groups)
  for (g_val in treatment_groups) {
    cohort_masks[[as.character(g_val)]] <- (unit_cohorts == g_val)
  }

  # -----------------------------------------------------------------------
  # 10. never_treated_mask
  # -----------------------------------------------------------------------
  never_treated_mask <- is.infinite(unit_cohorts)

  # -----------------------------------------------------------------------
  # 11. cohort_fractions: pi_g = n_g / n
  # -----------------------------------------------------------------------
  cohort_fractions <- stats::setNames(
    vapply(treatment_groups, function(g_val) sum(unit_cohorts == g_val) / n,
           numeric(1L)),
    as.character(treatment_groups)
  )

  # -----------------------------------------------------------------------
  # 12. Clustering
  # -----------------------------------------------------------------------
  cluster_indices <- NULL
  n_clusters      <- NULL
  if (!is.null(cluster)) {
    cluster_indices <- build_cluster_index(dt, unit, cluster, all_units)
    n_clusters      <- length(unique(cluster_indices))
  }

  # -----------------------------------------------------------------------
  # Assemble panel_obj
  # -----------------------------------------------------------------------
  panel_obj <- list(
    n                  = n,
    T_periods          = T_periods,
    outcome_wide       = outcome_wide,
    time_periods       = time_periods,
    period_1           = period_1,
    period_to_col      = period_to_col,
    all_units          = all_units,
    unit_cohorts       = unit_cohorts,
    treatment_groups   = treatment_groups,
    cohort_masks       = cohort_masks,
    never_treated_mask = never_treated_mask,
    cohort_fractions   = cohort_fractions,
    cluster_indices    = cluster_indices,
    n_clusters         = n_clusters,
    covariate_matrix   = NULL,   # deferred
    control_group      = control_group,
    anticipation       = as.integer(anticipation)
  )

  panel_obj
}

#' Build cluster integer index from cluster id column
#'
#' @param dt data.table (long format), sorted by unit then time
#' @param unit character scalar: unit id column name
#' @param cluster character scalar: cluster id column name
#' @param all_units sorted vector of unique unit ids
#'
#' @return integer vector length n (values 1..G)
#' @keywords internal
build_cluster_index <- function(dt, unit, cluster, all_units) {
  # Extract time-invariant cluster id per unit using base R tapply
  cl_vals      <- dt[[cluster]]
  unit_id_vals <- dt[[unit]]
  cl_map       <- tapply(cl_vals, unit_id_vals, function(x) x[1L])
  cl_ids       <- cl_map[match(all_units, names(cl_map))]

  # Map cluster id -> integer index
  sorted_cl <- sort(unique(cl_ids))
  cl_index  <- match(cl_ids, sorted_cl)
  cl_index
}
