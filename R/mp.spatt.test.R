#' @title mp.spatt.test
#'
#' @description Deprecated function for integrated moments test for conditional
#'  parallel trends holding in all pre-treatment time periods across all groups
#'
#' @inheritParams mp.spatt
#'
#' @export
mp.spatt.test <- function(...) {
  .Deprecated(new="conditional_did_pretest",
              msg="mp.spatt.test is now deprecated.  Please use conditional_did_pretest method instead.  It provides the same functionality along with some new improvements."
              )
}
