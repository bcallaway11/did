## The idea here is to combine the weighting function with Y and run the previous
## code for computing group-time average treatment effects
#' @title mp.spatt.test
#'
#' @description integrated moments test for conditional common trends holding in all pre-treatment time
#'  periods across all groups
#'
#' @inheritParams mp.spatt
#' @param weightfun A function that takes in two arguments, X and u, to compute
#'  the weighting function for the test.  The default is \code{1*(X <= u)}
#' @param xformlalist A list of formulas for the X variables.  This allows to
#'  test using different specifications for X, if desired
#' @param clustervarlist A list of cluster variables.  This allows to conduct
#'  the test using different levels of clustering, if desired.
#'
#' @examples
#' \dontrun{
#' data(mpdta)
#' mptest <- mp.spatt.test(lemp ~ treat, xformlalist=list(~lpop), data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year", clustervarlist=list(NULL))
#' summary(mptest[[1]])
#' }
#'
#' data(mpdta)
#' mptest <- mp.spatt.test(lemp ~ treat, xformlalist=list(NULL), data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year", clustervarlist=list(NULL))
#' summary(mptest[[1]])
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#'
#' @return list containing test results
#' @export
mp.spatt.test <- function(...) {
  .Deprecated(new="conditional_did_pretest",
              msg="mp.spatt.test is now deprecated.  Please use conditional_did_pretest method instead.  It provides the same functionality along with some new improvements."
              )
}
