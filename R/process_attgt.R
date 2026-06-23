#' @title Process Results from [compute.att_gt()]
#'
#' @param attgt.list list of results from [compute.att_gt()]
#'
#' @return list with elements:
#' \item{group}{which group a set of results belongs to}
#' \item{tt}{which time period a set of results belongs to}
#' \item{att}{the group time average treatment effect}
#'
#' @export
process_attgt <- function(attgt.list) {
  if (!is.list(attgt.list) || length(attgt.list) == 0L) {
    stop("attgt.list must be a non-empty list of group-time result objects.")
  }
  get_cell_value <- function(x, field, allow_na = FALSE) {
    if (!is.list(x)) {
      stop("Each attgt.list element must be a list-like group-time result object.")
    }
    value <- x[[field]]
    if (allow_na && length(value) == 1L && is.na(value)) {
      return(as.numeric(value))
    }
    if (!is.numeric(value) || length(value) != 1L ||
        (!allow_na && is.na(value))) {
      stop("Each attgt.list element must contain a numeric scalar '", field, "'.")
    }
    value
  }
  group <- vapply(attgt.list, get_cell_value, numeric(1), field = "group")
  att   <- vapply(attgt.list, get_cell_value, numeric(1), field = "att", allow_na = TRUE)
  tt    <- vapply(attgt.list, get_cell_value, numeric(1), field = "year")

  list(group=group, att=att, tt=tt)
}
