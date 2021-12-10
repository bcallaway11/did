#' @title Aggregate Group-Time Average Treatment Effects
#'
#' @description A function to take group-time average treatment effects
#'  and aggregate them into a smaller number of parameters.  There are
#'  several possible aggregations including "simple", "dynamic", "group",
#'  and "calendar."
#'
#' @param MP an MP object (i.e., the results of the [att_gt()] method)
#' @param type Which type of aggregated treatment effect parameter to compute.
#'   One option is "simple" (this just computes a weighted average of all
#'   group-time average treatment effects with weights proportional to group
#'   size).  Other options are "dynamic" (this computes average effects across
#'   different lengths of exposure to the treatment and is similar to an
#'   "event study"; here the overall effect averages the effect of the
#'   treatment across all positive lengths of exposure); "group" (this
#'   is the default option and
#'   computes average treatment effects across different groups; here
#'   the overall effect averages the effect across different groups); and
#'   "calendar" (this computes average treatment effects across different
#'   time periods; here the overall effect averages the effect across each
#'   time period).
#' @param balance_e If set (and if one computes dynamic effects), it balances
#'  the sample with respect to event time.  For example, if `balance.e=2`,
#'  `aggte` will drop groups that are not exposed to treatment for
#'  at least three periods. (the initial period when `e=0` as well as the
#'  next two periods when `e=1` and the `e=2`).  This ensures that
#'  the composition of groups does not change when event time changes.
#' @param min_e For event studies, this is the smallest event time to compute
#'  dynamic effects for.  By default, `min_e = -Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param max_e For event studies, this is the largest event time to compute
#'  dynamic effects for.  By default, `max_e = Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param na.rm Logical value if we are to remove missing Values from analyses. Defaults is FALSE.
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier bootstrap.  If standard errors are clustered, then one
#'  must set `bstrap=TRUE`. Default is value set in the MP object.  If bstrap is `FALSE`, then analytical
#'  standard errors are reported.
#' @param biters The number of bootstrap iterations to use.  The default is the value set in the MP object,
#'  and this is only applicable if `bstrap=TRUE`.
#'
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability `1-alp`.  In order to compute uniform confidence
#'  bands, `bstrap` must also be set to `TRUE`.  The default is
#'  the value set in the MP object
#' @param alp the significance level, default is value set in the MP object.
#' @param clustervars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level. Default is the variables set in the MP object

#'
#' @return An [`AGGTEobj`] object that holds the results from the
#'  aggregation
#'
#' @section Examples:
#'
#
#' Initial ATT(g,t) estimates from [att_gt()]
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' data(mpdta)
#' out <- att_gt(yname="lemp",
#'                tname="year",
#'                idname="countyreal",
#'                gname="first.treat",
#'                xformla=NULL,
#'                data=mpdta)
#' ```
#'
#' You can aggregate the ATT(g,t) in many ways.
#'
#' **Overall ATT:**
#' ```{r, comment = "#>", collapse = TRUE}
#' aggte(out, type = "simple")
#' ```
#'
#' **Dynamic ATT (Event-Study):**
#' ```{r, comment = "#>", collapse = TRUE}
#' aggte(out, type = "dynamic")
#' ```
#'
#' **ATT for each group:**
#' ```{r, comment = "#>", collapse = TRUE}
#' aggte(out, type = "group")
#' ```
#'
#' **ATT for each calendar year:**
#' ```{r, comment = "#>", collapse = TRUE}
#' aggte(out, type = "calendar")
#' ```
#'
#'
#' @export
aggte <- function(MP,
                  type = "group",
                  balance_e = NULL,
                  min_e = -Inf,
                  max_e = Inf,
                  na.rm = FALSE,
                  bstrap = NULL,
                  biters = NULL,
                  cband = NULL,
                  alp = NULL,
                  clustervars = NULL
                  ) {

  call <- match.call()

  compute.aggte(MP = MP,
                type = type,
                balance_e = balance_e,
                min_e = min_e,
                max_e = max_e,
                na.rm = na.rm,
                bstrap = bstrap,
                biters = biters,
                cband = cband,
                alp = alp,
                clustervars = clustervars,
                call = call)
}
