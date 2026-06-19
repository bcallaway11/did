#' @title trimmer
#'
#' @description A utility function to find observations that appear
#'  to violate support conditions.  This function is not called
#'  anywhere in the code, but it is just useful for debugging some
#'  common issues that users run into.
#'
#' @param g is a particular group (below I pass in 2009)
#' @inheritParams att_gt
#' @param threshold the cutoff for which observations are flagged as
#'  likely violators of the support condition.
#'
#' @return list of ids of observations that likely violate support conditions
#'
#' @export
trimmer <- function(g, tname, idname, gname, xformla, data, control_group="notyettreated", threshold=.999) {
  if (!all(class(data) == "data.frame")) {
    data <- as.data.frame(data)
  }
  validate_finite_numeric_scalar(g, "g")
  validate_column_name(tname, "tname", names(data))
  validate_column_name(idname, "idname", names(data))
  validate_column_name(gname, "gname", names(data))
  validate_xformla(xformla)
  validate_choice_scalar(
    control_group,
    "control_group",
    c("nevertreated", "notyettreated"),
    "control_group must be either 'nevertreated' or 'notyettreated'."
  )
  validate_probability_scalar(threshold, "threshold")

  time.period <- data[[tname]]
  this.data <- data[time.period == (g-1),]
  if (control_group == "notyettreated") {
    # not yet treated
    this.data <- this.data[(this.data[[gname]] >= g) |
                           (this.data[[gname]] == 0), ]
  } else {
    # never treated
    this.data <- this.data[(this.data[[gname]] == g) |
                           (this.data[[gname]] == 0), ]
  }
  this.data$D <- 1 * this.data[[gname]] == g
  this.pscore_reg <- glm(BMisc::toformula("D", BMisc::rhs_vars(xformla)),
                         data=this.data,
                         family=binomial(link="logit"))
  this.pscore <- predict(this.pscore_reg, type="response")
  dropper <- (this.pscore >  threshold) & (this.data$D==1)
  if (sum(dropper) > 0) {
    print("hard to match treated observations: ")
    print(this.data[dropper, idname, drop = FALSE])
    return(this.data[dropper, idname, drop = FALSE])
  }
}

#' @title overlap_logit_fit
#' @description Fit the per-(g,t) overlap-check propensity logit via fastglm's
#'  low-level entry point (`fastglmPure`), skipping the `fastglm()` wrapper's
#'  per-call input coercion and family/deviance bookkeeping. Passing fastglm's own
#'  defaults (`method = 0`, `tol = 1e-8`, `maxit = 100`) makes the fitted values --
#'  and hence the overlap decision (`max(fitted) >= 0.999`) -- bit-identical to the
#'  previous `fastglm::fastglm(x, y, family = binomial())` call, while ~1.4x faster.
#' @param x covariate matrix
#' @param y treatment indicator
#' @return a fastglmPure fit list (with `$fitted.values`)
#' @noRd
overlap_logit_fit <- function(x, y) {
  fastglm::fastglmPure(x = x, y = y, family = stats::binomial(),
                       method = 0L, tol = 1e-8, maxit = 100L)
}

#' @title overlap_check_fail
#' @description Evaluate the per-(g,t) overlap guard, i.e. whether the
#'  preliminary propensity logit gives `max(fitted) >= 0.999`. For an
#'  intercept-only design (`xformla = ~1`, the default) the unweighted logit MLE
#'  fits every unit at `mean(y)`, so the guard reduces to `mean(y) >= 0.999`
#'  with no fit at all -- except within `1e-6` of the 0.999 knife edge, where
#'  the IRLS iterate (tol = 1e-8; observed error up to ~1.5e-9) can land on
#'  either side of the cutoff. There we fall back to the real fit, keeping the
#'  decision bit-identical to always fitting.
#' @param x covariate matrix
#' @param y treatment indicator
#' @param intercept_only whether `x` is a single column of ones
#' @return TRUE if the overlap condition is violated
#' @noRd
overlap_check_fail <- function(x, y, intercept_only = FALSE) {
  if (intercept_only) {
    pbar <- mean(y)
    if (abs(pbar - 0.999) > 1e-6) return(pbar >= 0.999)
    # knife edge: defer to the actual fit below
  }
  max(overlap_logit_fit(x, y)$fitted.values) >= 0.999
}

#' @title rcond_check_fail
#' @description Evaluate the per-(g,t) regression-feasibility guard, i.e.
#'  whether `rcond(crossprod(x[controls, ]))` is below machine epsilon. For an
#'  intercept-only design the control-unit Gram matrix is the scalar count of
#'  control units (n0): `rcond` is exactly 1 when n0 > 0 and exactly 0 when
#'  n0 == 0, so the test reduces to "no control units" with no Gram matrix.
#' @param x covariate matrix
#' @param y treatment indicator (controls are `y == 0`)
#' @param intercept_only whether `x` is a single column of ones
#' @return TRUE if the control-unit design is singular or ill-conditioned
#' @noRd
rcond_check_fail <- function(x, y, intercept_only = FALSE) {
  if (intercept_only) {
    return(!any(y == 0))
  }
  control_covs <- x[y == 0, , drop = FALSE]
  rcond(crossprod(control_covs)) < .Machine$double.eps
}

#' @title Check User-Supplied Names Against Internal Names
#' @description Stop if user arguments reference names that `did` creates and
#' mutates internally.
#' @noRd
check_reserved_did_names <- function(yname, tname, idname, gname, xformla,
                                     weightsname = NULL, clustervars = NULL) {
  reserved_names <- c(".w", ".rowid", ".G", ".C", "post",
                      "asif_never_treated", "treated_first_period")
  xvars <- if (is.null(xformla)) character(0) else all.vars(xformla)
  used_names <- unique(c(yname, tname, idname, gname, weightsname,
                         clustervars, xvars))
  used_names <- used_names[!is.na(used_names) & nzchar(used_names)]
  bad_names <- intersect(used_names, reserved_names)
  if (length(bad_names) > 0L) {
    stop("The following variable name(s) are reserved for internal use by `did`: ",
         paste(bad_names, collapse = ", "), ". Please rename these column(s) ",
         "before calling att_gt().")
  }
}

validate_xformla <- function(xformla) {
  if (!is.null(xformla) && !inherits(xformla, "formula")) {
    stop("xformla must be NULL or a formula.")
  }
  invisible(xformla)
}

validate_character_scalar <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(invisible(x))
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(name, " must be a single non-missing character string.")
  }
  invisible(x)
}

validate_column_name <- function(x, name, data_names, allow_null = FALSE) {
  validate_character_scalar(x, name, allow_null = allow_null)
  if (is.null(x)) return(invisible(x))
  if (!(x %in% data_names)) {
    stop(name, " must be a character scalar and a name of a column from the dataset.")
  }
  invisible(x)
}

validate_column_names <- function(x, name, data_names, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(invisible(x))
  if (!is.character(x) || length(x) == 0L || anyNA(x) || any(!nzchar(x))) {
    stop(name, " must be NULL or a character vector naming column(s) from the dataset.")
  }
  missing <- setdiff(x, data_names)
  if (length(missing) > 0L) {
    stop(name, " contains column name(s) not found in the dataset: ",
         paste(missing, collapse = ", "), ".")
  }
  invisible(x)
}

complete_finite_cases <- function(data, finite_exclude = character(0)) {
  keep <- stats::complete.cases(data)
  if (!length(keep)) return(keep)
  for (nm in names(data)) {
    # Columns in finite_exclude still get the NA/NaN check via complete.cases()
    # above, but skip the is.finite() check below. This preserves Inf in such
    # columns -- notably gname, where Inf is a documented never-treated code
    # ("group status 0 or Inf") that must NOT be dropped as if it were bad data.
    if (nm %in% finite_exclude) next
    x <- data[[nm]]
    if (is.numeric(x)) {
      finite_x <- is.finite(x)
      if (!is.null(dim(finite_x)) && NROW(finite_x) == length(keep)) {
        finite_x <- rowSums(!finite_x) == 0L
      }
      keep <- keep & as.vector(finite_x)
    }
  }
  keep
}

validate_anticipation <- function(anticipation) {
  if (!is.numeric(anticipation)) {
    stop("anticipation must be numeric. Please convert it.")
  }
  if (length(anticipation) != 1L || is.na(anticipation) ||
      !is.finite(anticipation)) {
    stop("anticipation must be a single finite non-missing number. Please check your arguments.")
  }
  if (anticipation < 0 || anticipation != round(anticipation)) {
    stop("anticipation must be a non-negative whole number. Please check your arguments.")
  }
  invisible(anticipation)
}

validate_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be a single logical (TRUE or FALSE).")
  }
  invisible(x)
}

validate_numeric_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be a single non-missing number.")
  }
  invisible(x)
}

validate_finite_numeric_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(name, " must be a single finite non-missing number.")
  }
  invisible(x)
}

validate_finite_numeric_vector <- function(x, name, len) {
  if (!is.numeric(x) || length(x) != len || anyNA(x) || any(!is.finite(x))) {
    stop(name, " must be a numeric vector of length ", len,
         " with only finite non-missing values.")
  }
  invisible(x)
}

validate_positive_numeric_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop(name, " must be a single positive finite number.")
  }
  invisible(x)
}

validate_alp <- function(alp, name = "alp") {
  if (!is.numeric(alp) || length(alp) != 1 || is.na(alp) || alp <= 0 || alp >= 1) {
    stop(name, " must be a single number strictly between 0 and 1.")
  }
  invisible(alp)
}

validate_probability_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0 || x >= 1) {
    stop(name, " must be a single finite number strictly between 0 and 1.")
  }
  invisible(x)
}

validate_positive_whole_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) ||
      !is.finite(x) || x < 1 || x != round(x)) {
    stop(name, " must be a single positive whole number.")
  }
  invisible(x)
}

validate_nonnegative_whole_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) ||
      !is.finite(x) || x < 0 || x != round(x)) {
    stop(name, " must be a single non-negative whole number.")
  }
  invisible(x)
}

validate_choice_scalar <- function(x, name, choices, message = NULL) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !(x %in% choices)) {
    if (is.null(message)) {
      message <- paste0(name, " must be one of: ", paste(choices, collapse = ", "), ".")
    }
    stop(message)
  }
  invisible(x)
}

validate_optional_choice_scalar <- function(x, name, choices, message = NULL) {
  if (is.null(x)) return(invisible(x))
  validate_choice_scalar(x, name, choices, message)
}

validate_optional_numeric_scalar <- function(x, name) {
  if (is.null(x)) return(invisible(x))
  validate_numeric_scalar(x, name)
}

#' @title get_wide_data
#' @description A utility function to convert long data to wide data, i.e., takes a 2 period dataset and turns it into a cross sectional dataset.
#'
#' @param data data.table used in function
#' @param yname name of outcome variable that can change over time
#' @param idname unique id
#' @param tname time period name
#'
#' @return data from first period with .y0 (outcome in first period), .y1 (outcome in second period),
#' and .dy (change in outcomes over time) appended to it
#' @noRd
get_wide_data <- function(data, yname, idname, tname) {
  # check if data is data.table
  if (!is.data.table(data)) {
    stop("data must be a data.table")
  }

  # check that only 2 periods of data are available
  if (length(unique(data[[tname]])) != 2) {
    stop("data must have exactly 2 periods")
  }

  # making computations — use set() to avoid copy triggered by $<-
  y_vals <- data[[yname]]
  set(data, j = ".y1", value = data.table::shift(y_vals, -1))
  set(data, j = ".y0", value = y_vals)
  set(data, j = ".dy", value = data[[".y1"]] - y_vals)

  # Subset to first period's rows
  # Pre-extract to avoid data.table .checkTypos when column name matches variable name
  .time_vals <- data[[tname]]
  first.period <- min(.time_vals)
  data <- data[.time_vals == first.period]

  return(data)
}

#' @title Check balanced panel data
#' @description A utility function to check if your dataset is a balanced panel dataset.
#'
#' Precondition: at most one row per (id_col, time_col) combination -- enforced
#' upstream by validate_args(). Under that uniqueness, every unit has at most one
#' row per period, so the panel is balanced iff nrow == n_ids * n_periods.
#'
#' @param data data.table used in function
#' @param id_col name of id column in the dataset
#' @param time_col name of time column in the dataset
#'
#' @return Boolean indicating if the dataset is balanced or not.
#' @noRd
check_balance <- function(data, id_col, time_col) {

  # Check if every unit is observed in every time period (no grouping needed
  # given the (id, time) uniqueness precondition above)
  is_balanced <- nrow(data) == data.table::uniqueN(data[[id_col]]) * data.table::uniqueN(data[[time_col]])

  return(is_balanced)
}
