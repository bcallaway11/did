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
  this.pscore_reg <- glm(BMisc::toformula("D", BMisc::rhs.vars(xformla)),
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
#' @param data data.table used in function
#' @param id_col name of id column in the dataset
#' @param time_col name of time column in the dataset
#'
#' @return Boolean indicating if the dataset is balanced or not.
#' @noRd
check_balance <- function(data, id_col, time_col) {

  # Count the number of observations per unit (idname)
  panel_counts <- data[, .N, by = c(id_col)]

  # Determine the maximum number of time periods for any unit
  max_time_periods <- data.table::uniqueN(data[[time_col]])

  # Check if every unit has the same number of time periods as max_time_periods
  is_balanced <- all(panel_counts$N == max_time_periods)

  return(is_balanced)
}
