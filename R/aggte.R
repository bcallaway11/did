#' @title Aggregate Group-Time Average Treatment Effects
#'
#' @description A function to take group-time average treatment effects
#'  and aggregate them into a smaller number of parameters.  There are
#'  several possible aggregations including "simple", "dynamic", "selective",
#'  and "calendar."
#'
#' @param MP an MP object (i.e., the results of the \code{att_gt} method)
#' @param type Which type of aggregated treatment effect parameter to compute.
#'   The default is "simple" (this just computes a weighted average of all
#'   group-time average treatment effects with weights proportional to group
#'   size).  Other options are "dynamic" (this computes average effects across
#'   different lengths of exposure to the treatment and is similar to an
#'   "event study"; here the overall effect averages the effect of the
#'   treatment across all positive lengths of exposure); "selective" (this
#'   computes average treatment effects across different groups; here
#'   the overall effect averages the effect across different groups); and
#'   "calendar" (this computes average treatment effects across different
#'   time periods; here the overall effect averages the effect across each
#'   time period).
#' @param balance.e If set (and if one computes dynamic effects), it balances
#'  the sample with respect to event time.  For example, if \code{balance.e=2},
#'  \code{aggte} will drop groups that are not exposed to treatment for
#'  at least three periods. (the initial period when \code{e=0} as well as the
#'  next two periods when \code{e=1} and the \code{e=2}).  This ensures that
#'  the composition of groups does not change when event time changes.
#' @param na.rm Logical value if we are to remove missing Values from analyses. Defaults is FALSE.
#'
#' @return AGGTEobj
#' @export
aggte <- function(MP, type="simple", balance.e=NULL, na.rm = FALSE) {
  compute.aggte(MP, type, balance.e, na.rm)
}
