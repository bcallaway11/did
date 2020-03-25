#' @title mp.spatt
#'
#' @description \code{mp.spatt} computes the ATT in the case where there are more
#'  than two periods of data and allowing for treatment to occur at different points in time
#'  extending the method of Abadie (2005).  This method relies on once individuals are treated
#'  they remain in the treated state for the duration.
#'
#' @param formla The formula y ~ d where y is the outcome and d is the
#'  treatment indicator (d should be binary)
#' @param xformla A optional one sided formula for additional covariates that
#'  will be adjusted for.  E.g ~ age + education.  Additional covariates can
#'  also be passed by name using the x paramater.
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param aggte boolean for whether or not to compute aggregate treatment effect parameters, default TRUE
#' @param w A vector of weights for each observation (not implemented)
#' @param panel Boolean indicating whether the data is panel or repeated cross
#'  sections
#' @param idname The individual (cross-sectional unit) id name
#' @param first.treat.name The name of the variable in \code{data} that contains the first
#'  period when a particular observation is treated.  This should be a positive
#'  number for all observations in treated groups.  It should be 0 for observations
#'  in the untreated group.
#' @param alp the significance level, default is 0.05
#' @param method The method for estimating the propensity score when covariates
#'  are included
#' @param se Boolean whether or not to compute standard errors
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier boostrap.  If standard errors are clustered, then one
#'  must set \code{bstrap=TRUE}.
#' @param biters The number of boostrap iterations to use.  The default is 100,
#'  and this is only applicable if \code{bstrap=TRUE}.
#' @param clustervars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level.
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability \code{1-alp}.  The default is \code{FALSE}
#'  and the resulting standard errors will be pointwise.
#' @param citers Computing uniform confidence bands requires the bootstrap,
#'  if \code{cband = TRUE}, then this is the number of boostrap iterations
#'  to compute the conidence band.  The default is 100.
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors#' (not implemented)
#' @param pl Boolean for whether or not to use parallel processing
#' @param cores The number of cores to use for parallel processing
#' @param printdetails Boolean for showing detailed results or not
#'
#' @examples
#' data(mpdta)
#'
#' ## with covariates
#' out1 <- mp.spatt(lemp ~ treat, xformla=~lpop, data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year",
#'                 bstrap=FALSE, se=TRUE, cband=FALSE)
#' ## summarize the group-time average treatment effects
#' summary(out1)
#' ## summarize the aggregated treatment effect parameters
#' summary(out1$aggte)
#'
#' ## without any covariates
#' out2 <- mp.spatt(lemp ~ treat, xformla=NULL, data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year",
#'                 bstrap=FALSE, se=TRUE, cband=FALSE)
#' summary(out2)
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#'
#' @return \code{MP} object
#'
#' @export
mp.spatt <- function(...) {
## mp.spatt <- function(formla, xformla=NULL, data, tname,
##                      aggte=TRUE, w=NULL, panel=FALSE,
##                      idname=NULL, first.treat.name, alp=0.05,
##                      method="logit", se=TRUE,
##                      bstrap=FALSE, biters=100, clustervars=NULL,
##                      cband=FALSE, citers=100,
##                      seedvec=NULL, pl=FALSE, cores=2,
##                      printdetails=TRUE) {

  .Deprecated(new="att_gt",
              msg="mp.spatt is now deprecated.  Please use att_gt method instead.  It provides the same functionality along with some new improvements."
              )
}
