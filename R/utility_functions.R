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

  time.period <- data[,tname]
  this.data <- data[time.period == (g-1),]
  if (control_group == "notyettreated") {
    # not yet treated
    this.data <- this.data[(this.data[,gname] >= g) |
                           (this.data[,gname] == 0), ]
  } else {
    # never treated
    this.data <- this.data[(this.data[,gname] == g) |
                           (this.data[,gname] == 0), ]
  }
  this.data$D <- 1*this.data[,gname]==g
  this.pscore_reg <- glm(BMisc::toformula("D", BMisc::rhs_vars(xformla)),
                         data=this.data,
                         family=binomial(link="logit"))
  this.pscore <- predict(this.pscore_reg, type="response")
  dropper <- (this.pscore >  threshold) & (this.data$D==1)
  if (sum(dropper) > 0) {
    print("hard to match treated observations: ")
    print(this.data[dropper,idname])
    return(this.data[dropper,idname])
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
