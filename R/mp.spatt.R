#' @title Group-Time Average Treatment Effects with Multiple Periods
#'
#' @description Deprecated function for computing group-time average treatment effects
#'
#' @param ... extra arguments
#' 
#' @export
mp.spatt <- function(...) {
  .Deprecated(new="att_gt",
              msg="mp.spatt is now deprecated.  Please use att_gt method instead.  It provides the same functionality along with some new improvements."
              )
}
