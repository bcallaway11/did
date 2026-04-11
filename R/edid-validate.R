# edid-validate.R
# Input validation for edid(). All checks are performed before any computation.

#' Validate inputs to \code{edid()}
#'
#' Performs all structural and type checks on user-supplied arguments.
#' Returns invisibly \code{TRUE} on success; stops with an informative message
#' on any failure.
#'
#' @param data data.frame or coercible
#' @param outcome character scalar: outcome column name
#' @param unit character scalar: unit id column name
#' @param time character scalar: time column name
#' @param first_treat character scalar: first-treatment-period column name
#' @param covariates character vector or NULL
#' @param pt_assumption character scalar, already matched via \code{match.arg}
#' @param alpha numeric scalar in (0, 1)
#' @param cluster character scalar or NULL
#' @param control_group character scalar, already matched via \code{match.arg}
#' @param n_bootstrap non-negative integer
#' @param anticipation non-negative integer
#' @param survey_design always NULL (survey not yet implemented)
#'
#' @return invisibly TRUE
#' @keywords internal
validate_edid_inputs <- function(
  data, outcome, unit, time, first_treat, covariates,
  pt_assumption, alpha, cluster, control_group,
  n_bootstrap, anticipation, survey_design
) {

  # ------------------------------------------------------------------
  # 1. data is data.frame-like and has rows
  # ------------------------------------------------------------------
  if (!is.data.frame(data) && !inherits(data, "data.table") &&
      !inherits(data, "tbl_df")) {
    # try coercing
    tryCatch(
      data <- as.data.frame(data),
      error = function(e) stop("`data` must be a data.frame or coercible to one.")
    )
  }
  if (nrow(data) == 0L) {
    stop("`data` has no rows.")
  }

  # ------------------------------------------------------------------
  # 2. outcome / unit / time / first_treat are character scalars naming
  #    existing columns
  # ------------------------------------------------------------------
  .check_col <- function(arg, argname) {
    if (!is.character(arg) || length(arg) != 1L) {
      stop(sprintf("`%s` must be a character scalar (column name).", argname))
    }
    if (!arg %in% names(data)) {
      stop(sprintf("`%s` = \"%s\" is not a column in `data`.", argname, arg))
    }
  }
  .check_col(outcome,     "outcome")
  .check_col(unit,        "unit")
  .check_col(time,        "time")
  .check_col(first_treat, "first_treat")

  # Columns must be distinct
  col_names <- c(outcome, unit, time, first_treat)
  if (anyDuplicated(col_names)) {
    stop("`outcome`, `unit`, `time`, and `first_treat` must name distinct columns.")
  }

  # ------------------------------------------------------------------
  # 3. outcome column is numeric; no NA; all finite
  # ------------------------------------------------------------------
  y_col <- data[[outcome]]
  if (!is.numeric(y_col)) {
    stop(sprintf("Column `%s` (outcome) must be numeric.", outcome))
  }
  if (anyNA(y_col)) {
    stop(sprintf("Column `%s` (outcome) contains NA values. ",
                 outcome),
         "edid() requires a complete, balanced panel with no missing outcomes.")
  }
  if (!all(is.finite(y_col))) {
    stop(sprintf("Column `%s` (outcome) contains non-finite values (Inf/-Inf/NaN). ",
                 outcome),
         "edid() requires all outcomes to be finite.")
  }

  # ------------------------------------------------------------------
  # 4. time column is numeric; no NA
  # ------------------------------------------------------------------
  t_col <- data[[time]]
  if (!is.numeric(t_col)) {
    stop(sprintf("Column `%s` (time) must be numeric.", time))
  }
  if (anyNA(t_col)) {
    stop(sprintf("Column `%s` (time) contains NA values.", time))
  }

  # ------------------------------------------------------------------
  # 5. first_treat column is numeric; no NA
  # ------------------------------------------------------------------
  ft_col <- data[[first_treat]]
  if (!is.numeric(ft_col)) {
    stop(sprintf("Column `%s` (first_treat) must be numeric.", first_treat))
  }
  if (anyNA(ft_col)) {
    stop(sprintf("Column `%s` (first_treat) contains NA values. ",
                 first_treat),
         "Use Inf to denote never-treated units.")
  }

  # ------------------------------------------------------------------
  # 6. No duplicate (unit, time) rows
  # ------------------------------------------------------------------
  ut_key <- paste(data[[unit]], data[[time]], sep = "___")
  if (anyDuplicated(ut_key)) {
    stop("Duplicate (unit, time) pairs found in `data`. ",
         "edid() requires a balanced panel with exactly one observation per unit-period.")
  }

  # ------------------------------------------------------------------
  # 7. Panel is balanced: every unit appears in every time period
  # ------------------------------------------------------------------
  all_units_v  <- unique(data[[unit]])
  all_times_v  <- unique(data[[time]])
  n_units      <- length(all_units_v)
  n_times      <- length(all_times_v)
  expected_obs <- n_units * n_times
  if (nrow(data) != expected_obs) {
    stop(sprintf(
      "Panel is unbalanced: expected %d rows (%d units x %d periods) but found %d. ",
      expected_obs, n_units, n_times, nrow(data)),
      "edid() requires a balanced panel.")
  }

  # ------------------------------------------------------------------
  # 8. Treatment is absorbing: first_treat is constant within unit
  # ------------------------------------------------------------------
  ft_by_unit <- tapply(data[[first_treat]], data[[unit]], function(x) length(unique(x)))
  if (any(ft_by_unit > 1L)) {
    bad <- names(ft_by_unit)[ft_by_unit > 1L]
    stop(sprintf(
      "`%s` (first_treat) is not constant within unit for %d unit(s) (e.g., %s). ",
      first_treat, length(bad), bad[1]),
      "Treatment must be absorbing.")
  }

  # ------------------------------------------------------------------
  # 9-10. Control group availability
  # ------------------------------------------------------------------
  # Get time-invariant first_treat per unit (one row per unit)
  unit_ft <- tapply(data[[first_treat]], data[[unit]], `[`, 1L)

  if (control_group == "never_treated") {
    n_never <- sum(is.infinite(unit_ft))
    if (n_never == 0L) {
      stop("No never-treated units found (`first_treat == Inf`). ",
           "edid() requires never-treated units when `control_group = \"never_treated\"`.")
    }
  } else {
    # last_cohort
    finite_ft <- unit_ft[is.finite(unit_ft)]
    if (length(finite_ft) == 0L) {
      stop("No finite first-treatment values found; cannot determine last cohort.")
    }
    last_g <- max(finite_ft)
    n_last <- sum(finite_ft == last_g)
    if (n_last == 0L) {
      stop("No units in the last treated cohort found. ",
           "Cannot use `control_group = \"last_cohort\"`.")
    }
  }

  # ------------------------------------------------------------------
  # 11. covariates != NULL -> stop (stub)
  # ------------------------------------------------------------------
  if (!is.null(covariates)) {
    stop("covariate path not yet implemented in edid(). ",
         "Pass covariates = NULL or omit the covariates argument.")
  }

  # ------------------------------------------------------------------
  # 12. survey_design != NULL -> stop (stub)
  # ------------------------------------------------------------------
  if (!is.null(survey_design)) {
    stop("survey_design not yet implemented in edid(). ",
         "Pass survey_design = NULL or omit the argument.")
  }

  # ------------------------------------------------------------------
  # 13. alpha in (0, 1)
  # ------------------------------------------------------------------
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a numeric scalar strictly between 0 and 1.")
  }

  # ------------------------------------------------------------------
  # 14. n_bootstrap >= 0 integer
  # ------------------------------------------------------------------
  if (!is.numeric(n_bootstrap) || length(n_bootstrap) != 1L ||
      n_bootstrap < 0 || n_bootstrap != floor(n_bootstrap)) {
    stop("`n_bootstrap` must be a non-negative integer.")
  }

  # ------------------------------------------------------------------
  # 15. anticipation >= 0 integer; effective pre-treatment check
  # ------------------------------------------------------------------
  if (!is.numeric(anticipation) || length(anticipation) != 1L ||
      anticipation < 0 || anticipation != floor(anticipation)) {
    stop("`anticipation` must be a non-negative integer.")
  }
  min_ft <- min(unit_ft[is.finite(unit_ft)])
  min_t  <- min(all_times_v)
  if (min_ft - anticipation <= min_t) {
    stop(sprintf(
      "With `anticipation = %d`, the earliest treatment cohort (%g) would be treated at or before ",
      anticipation, min_ft),
      sprintf("the first observed period (%g). ", min_t),
      "There must be at least one pre-treatment period available.")
  }

  # ------------------------------------------------------------------
  # 16. cluster column checks
  # ------------------------------------------------------------------
  if (!is.null(cluster)) {
    if (!is.character(cluster) || length(cluster) != 1L) {
      stop("`cluster` must be a character scalar naming a column in `data`, or NULL.")
    }
    if (!cluster %in% names(data)) {
      stop(sprintf("`cluster` = \"%s\" is not a column in `data`.", cluster))
    }
    cl_col <- data[[cluster]]
    if (anyNA(cl_col)) {
      stop(sprintf("Cluster column `%s` contains NA values.", cluster))
    }
    # Time-invariant within unit
    cl_by_unit <- tapply(cl_col, data[[unit]], function(x) length(unique(x)))
    if (any(cl_by_unit > 1L)) {
      bad <- names(cl_by_unit)[cl_by_unit > 1L]
      stop(sprintf(
        "Cluster variable `%s` is not time-invariant for %d unit(s) (e.g., %s). ",
        cluster, length(bad), bad[1]),
        "Cluster variable must be constant within unit.")
    }
  }

  invisible(TRUE)
}
