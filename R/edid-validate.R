# edid-validate.R
# Input validation for edid(). All checks are performed before any computation.

#' Validate inputs to \code{edid()}
#'
#' Performs all structural and type checks on user-supplied arguments.
#' Returns invisibly \code{TRUE} on success; stops with an informative message
#' on any failure.
#'
#' @param data data.frame or coercible
#' @param yname character scalar: outcome column name
#' @param idname character scalar: unit id column name
#' @param tname character scalar: time column name
#' @param gname character scalar: first-treatment-period column name
#' @param covariates character vector or NULL
#' @param pt_assumption character scalar, already matched via \code{match.arg}
#' @param alp numeric scalar in (0, 1)
#' @param clustervars character scalar or NULL
#' @param control_group character scalar, already matched via \code{match.arg}
#' @param biters non-negative integer (internal bootstrap iterations)
#' @param anticipation non-negative integer
#' @param survey_design always NULL (survey not yet implemented)
#'
#' @return invisibly TRUE
#' @keywords internal
validate_edid_inputs <- function(
  data, yname, idname, tname, gname, xformla = NULL, covariates,
  pt_assumption, alp, clustervars, control_group,
  biters, anticipation, survey_design
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
  # 2. yname / idname / tname / gname are character scalars naming
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
  .check_col(yname,  "yname")
  .check_col(idname, "idname")
  .check_col(tname,  "tname")
  .check_col(gname,  "gname")

  # Columns must be distinct
  col_names <- c(yname, idname, tname, gname)
  if (anyDuplicated(col_names)) {
    stop("`yname`, `idname`, `tname`, and `gname` must name distinct columns.")
  }

  # ------------------------------------------------------------------
  # 3. yname column is numeric; no NA; all finite
  # ------------------------------------------------------------------
  y_col <- data[[yname]]
  if (!is.numeric(y_col)) {
    stop(sprintf("Column `%s` (yname) must be numeric.", yname))
  }
  if (anyNA(y_col)) {
    stop(sprintf("Column `%s` (yname) contains NA values. ", yname),
         "edid() requires a complete, balanced panel with no missing outcomes.")
  }
  if (!all(is.finite(y_col))) {
    stop(sprintf("Column `%s` (yname) contains non-finite values (Inf/-Inf/NaN). ", yname),
         "edid() requires all outcomes to be finite.")
  }

  # ------------------------------------------------------------------
  # 4. tname column is numeric; no NA
  # ------------------------------------------------------------------
  t_col <- data[[tname]]
  if (!is.numeric(t_col)) {
    stop(sprintf("Column `%s` (tname) must be numeric.", tname))
  }
  if (anyNA(t_col)) {
    stop(sprintf("Column `%s` (tname) contains NA values.", tname))
  }

  # ------------------------------------------------------------------
  # 5. gname column is numeric; no NA
  # ------------------------------------------------------------------
  ft_col <- data[[gname]]
  if (!is.numeric(ft_col)) {
    stop(sprintf("Column `%s` (gname) must be numeric.", gname))
  }
  if (anyNA(ft_col)) {
    stop(sprintf("Column `%s` (gname) contains NA values. ", gname),
         "Use Inf to denote never-treated units.")
  }

  # ------------------------------------------------------------------
  # 6. No duplicate (idname, tname) rows
  # ------------------------------------------------------------------
  ut_key <- paste(data[[idname]], data[[tname]], sep = "___")
  if (anyDuplicated(ut_key)) {
    stop("Duplicate (idname, tname) pairs found in `data`. ",
         "edid() requires a balanced panel with exactly one observation per unit-period.")
  }

  # ------------------------------------------------------------------
  # 7. Panel is balanced: every unit appears in every time period
  # ------------------------------------------------------------------
  all_units_v  <- unique(data[[idname]])
  all_times_v  <- unique(data[[tname]])
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
  # 8. Treatment is absorbing: gname is constant within unit
  # ------------------------------------------------------------------
  ft_by_unit <- tapply(data[[gname]], data[[idname]], function(x) length(unique(x)))
  if (any(ft_by_unit > 1L)) {
    bad <- names(ft_by_unit)[ft_by_unit > 1L]
    stop(sprintf(
      "`%s` (gname) is not constant within unit for %d unit(s) (e.g., %s). ",
      gname, length(bad), bad[1]),
      "Treatment must be absorbing.")
  }

  # ------------------------------------------------------------------
  # 9-10. Control group availability
  # ------------------------------------------------------------------
  # Get time-invariant gname per unit (one row per unit)
  unit_ft <- tapply(data[[gname]], data[[idname]], `[`, 1L)

  if (control_group == "nevertreated") {
    n_never <- sum(is.infinite(unit_ft))
    if (n_never == 0L) {
      stop("No never-treated units found (`gname == Inf`). ",
           "edid() requires never-treated units when `control_group = \"nevertreated\"`.")
    }
  } else {
    # notyettreated
    finite_ft <- unit_ft[is.finite(unit_ft)]
    if (length(finite_ft) == 0L) {
      stop("No finite first-treatment values found; cannot determine last cohort.")
    }
    last_g <- max(finite_ft)
    n_last <- sum(finite_ft == last_g)
    if (n_last == 0L) {
      stop("No units in the last treated cohort found. ",
           "Cannot use `control_group = \"notyettreated\"`.")
    }
  }

  # ------------------------------------------------------------------
  # 11. covariates is deprecated: error with redirect message
  # ------------------------------------------------------------------
  if (!is.null(covariates)) {
    stop(
      "The 'covariates' argument has been replaced by 'xformla'. ",
      "Pass a formula like xformla = ~ X1 + X2."
    )
  }

  # ------------------------------------------------------------------
  # 11b. xformla validation
  # ------------------------------------------------------------------
  if (!is.null(xformla)) {
    if (!inherits(xformla, "formula")) {
      stop("`xformla` must be a one-sided formula (e.g., ~ X1 + X2) or NULL.")
    }
    # Extract RHS variable names (skip ~1)
    rhs_vars <- all.vars(xformla)
    if (length(rhs_vars) > 0L) {
      missing_vars <- setdiff(rhs_vars, names(data))
      if (length(missing_vars) > 0L) {
        stop(sprintf(
          "Variable(s) in `xformla` not found in `data`: %s",
          paste(missing_vars, collapse = ", ")
        ))
      }
      # Check all covariate columns are numeric or factor
      for (v in rhs_vars) {
        if (!is.numeric(data[[v]]) && !is.factor(data[[v]])) {
          stop(sprintf("Covariate column `%s` must be numeric or factor.", v))
        }
      }
    }
  }

  # ------------------------------------------------------------------
  # 12. survey_design != NULL -> stop (stub)
  # ------------------------------------------------------------------
  if (!is.null(survey_design)) {
    stop("survey_design not yet implemented in edid(). ",
         "Pass survey_design = NULL or omit the argument.")
  }

  # ------------------------------------------------------------------
  # 13. alp in (0, 1)
  # ------------------------------------------------------------------
  if (!is.numeric(alp) || length(alp) != 1L || alp <= 0 || alp >= 1) {
    stop("`alp` must be a numeric scalar strictly between 0 and 1.")
  }

  # ------------------------------------------------------------------
  # 14. biters >= 0 integer
  # ------------------------------------------------------------------
  if (!is.numeric(biters) || length(biters) != 1L ||
      biters < 0 || biters != floor(biters)) {
    stop("`biters` must be a non-negative integer.")
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
  # 16. clustervars column checks
  # ------------------------------------------------------------------
  if (!is.null(clustervars)) {
    if (!is.character(clustervars) || length(clustervars) != 1L) {
      stop("`clustervars` must be a character scalar naming a column in `data`, or NULL.")
    }
    if (!clustervars %in% names(data)) {
      stop(sprintf("`clustervars` = \"%s\" is not a column in `data`.", clustervars))
    }
    cl_col <- data[[clustervars]]
    if (anyNA(cl_col)) {
      stop(sprintf("Cluster column `%s` contains NA values.", clustervars))
    }
    # Time-invariant within unit
    cl_by_unit <- tapply(cl_col, data[[idname]], function(x) length(unique(x)))
    if (any(cl_by_unit > 1L)) {
      bad <- names(cl_by_unit)[cl_by_unit > 1L]
      stop(sprintf(
        "Cluster variable `%s` is not time-invariant for %d unit(s) (e.g., %s). ",
        clustervars, length(bad), bad[1]),
        "Cluster variable must be constant within unit.")
    }
  }

  invisible(TRUE)
}
