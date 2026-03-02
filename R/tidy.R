#' tidy results
#' @export
generics::tidy

#' tidy results
#' @export
generics::glance

#' Number of observations used to fit an MP object
#'
#' @importFrom stats nobs
#' @param object a model of class MP produced by the [att_gt()] function
#' @param ... additional arguments (ignored)
#' @return Integer. The number of unique units in the data.
#' @export
nobs.MP <- function(object, ...) {
  object$n
}

#' Number of observations used to fit an AGGTEobj object
#'
#' @importFrom stats nobs
#' @param object a model of class AGGTEobj produced by the [aggte()] function
#' @param ... additional arguments (ignored)
#' @return Integer. The number of unique units in the data.
#' @export
nobs.AGGTEobj <- function(object, ...) {
  if (object$DIDparams$faster_mode) {
    object$DIDparams$id_count
  } else {
    object$DIDparams$n
  }
}

#' Tidy an MP object into a data frame
#'
#' Returns a tidy data frame of group-time average treatment effect estimates
#' from an [att_gt()] result.
#'
#' @importFrom generics tidy
#' @param x a model of class MP produced by the [att_gt()] function
#' @param ... Additional arguments to tidying method.
#'
#' @return A data frame with one row per ATT(g,t) estimate and columns:
#'   \item{term}{ATT(g,t) label}
#'   \item{group}{the treatment cohort g}
#'   \item{time}{the time period t}
#'   \item{estimate}{the ATT(g,t) point estimate}
#'   \item{std.error}{standard error}
#'   \item{statistic}{t-statistic (\code{estimate / std.error})}
#'   \item{p.value}{two-sided pointwise p-value (\code{2*(1-pnorm(|t|))}).
#'     Marginal per-estimate; does \strong{not} account for multiple testing
#'     across ATT(g,t) cells.}
#'   \item{conf.low, conf.high}{simultaneous confidence band limits, using the
#'     bootstrap uniform critical value when \code{bstrap=TRUE} and
#'     \code{cband=TRUE}, otherwise pointwise}
#'   \item{point.conf.low, point.conf.high}{pointwise confidence interval limits
#'     using \code{qnorm(1 - alp/2)}, always}
#'
#' @export
tidy.MP <- function(x, ...) {
  out <- data.frame(
    term      = paste0('ATT(',x$group,",",x$t, ")"),
    group     = x$group,
    time      = x$t,
    estimate  = x$att,
    std.error = x$se,
    statistic = x$att / x$se,
    p.value   = 2 * (1 - stats::pnorm(abs(x$att / x$se))),
    conf.low  = x$att - x$c * x$se,
    conf.high = x$att + x$c * x$se,
    point.conf.low  = x$att - stats::qnorm(1 - x$alp/2) * x$se,
    point.conf.high = x$att + stats::qnorm(1 - x$alp/2) * x$se)
  out
}

#' glance model characteristics from MP objects
#'
#' @importFrom generics glance
#' @param x a model of class MP produced by the [att_gt()] function
#' @param ... other arguments passed to methods
#' @export
glance.MP <- function(x, ...) {
  out <- data.frame(
    nobs          = x$n,
    ngroup        = x$DIDparams$nG,
    ntime         = x$DIDparams$nT,
    control.group = x$DIDparams$control_group,
    est.method    = x$DIDparams$est_method)
  out
}

#' Tidy an AGGTEobj into a data frame
#'
#' Returns a tidy data frame of aggregated treatment effect estimates from an
#' [aggte()] result.
#'
#' @importFrom generics tidy
#' @param x a model of class AGGTEobj produced by the [aggte()] function
#' @param ... Additional arguments to tidying method.
#'
#' @return A data frame whose columns depend on \code{type}:
#'   \item{type}{the aggregation type: \code{"simple"}, \code{"dynamic"},
#'     \code{"group"}, or \code{"calendar"}}
#'   \item{term}{label for each estimate}
#'   \item{estimate}{point estimate}
#'   \item{std.error}{standard error}
#'   \item{statistic}{t-statistic (\code{estimate / std.error})}
#'   \item{p.value}{two-sided pointwise p-value (\code{2*(1-pnorm(|t|))}).
#'     Marginal per-estimate; does \strong{not} account for multiple testing
#'     across event times or groups.}
#'   \item{conf.low, conf.high}{simultaneous confidence band limits.  When
#'     \code{bstrap=TRUE} and \code{cband=TRUE} these use the bootstrap uniform
#'     critical value (\code{crit.val.egt}); otherwise they equal the pointwise
#'     intervals.  For \code{type="simple"} and the overall average row of
#'     \code{type="group"}, a single scalar is returned so simultaneous and
#'     pointwise coincide.}
#'   \item{point.conf.low, point.conf.high}{pointwise confidence interval limits
#'     always using \code{qnorm(1 - alp/2)}.}
#'
#' @details
#' The key distinction between \code{conf.low}/\code{conf.high} and
#' \code{point.conf.low}/\code{point.conf.high} is that the former accounts for
#' multiple testing across all estimates (simultaneous coverage), while the
#' latter provides marginal (per-estimate) coverage only.  Use the simultaneous
#' bands when you want to make joint inferences across all event times or groups.
#'
#' @export
tidy.AGGTEobj<- function(x, ...) {
  if(x$type == "dynamic"){
    out <- data.frame(
      type          = x$type,
      term = paste0('ATT(', x$egt, ")"),
      event.time= x$egt,
      estimate  = x$att.egt,
      std.error = x$se.egt,
      statistic = x$att.egt / x$se.egt,
      p.value   = 2 * (1 - stats::pnorm(abs(x$att.egt / x$se.egt))),
      conf.low  = x$att.egt - x$crit.val.egt * x$se.egt,
      conf.high = x$att.egt + x$crit.val.egt * x$se.egt,
      point.conf.low  = x$att.egt - stats::qnorm(1 - x$DIDparams$alp/2) * x$se.egt,
      point.conf.high = x$att.egt + stats::qnorm(1 - x$DIDparams$alp/2) * x$se.egt)
  }
  if(x$type == "group"){
    out <- data.frame(
      type     = x$type,
      term = c(paste0('ATT(Average)'), paste0('ATT(', x$egt, ")")),
      group    = c('Average', x$egt),
      estimate  = c(x$overall.att, x$att.egt),
      std.error = c(x$overall.se, x$se.egt),
      statistic = c(x$overall.att, x$att.egt) / c(x$overall.se, x$se.egt),
      p.value   = 2 * (1 - stats::pnorm(abs(c(x$overall.att, x$att.egt) / c(x$overall.se, x$se.egt)))),
      conf.low  = c(x$overall.att - stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se, x$att.egt - x$crit.val.egt * x$se.egt),
      conf.high = c(x$overall.att + stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se, x$att.egt + x$crit.val.egt * x$se.egt),
      point.conf.low  = c(x$overall.att - stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se, x$att.egt - stats::qnorm(1 - x$DIDparams$alp/2) * x$se.egt),
      point.conf.high = c(x$overall.att + stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se,x$att.egt + stats::qnorm(1 - x$DIDparams$alp/2) * x$se.egt))
     }
  
  if(x$type == "calendar"){
    out <- data.frame(
      type      = x$type,
      time      = x$egt,
      term = paste0('ATT(', x$egt, ")"),
      estimate  = x$att.egt,
      std.error = x$se.egt,
      statistic = x$att.egt / x$se.egt,
      p.value   = 2 * (1 - stats::pnorm(abs(x$att.egt / x$se.egt))),
      conf.low  = x$att.egt - x$crit.val.egt * x$se.egt,
      conf.high = x$att.egt + x$crit.val.egt * x$se.egt,
      point.conf.low  = x$att.egt - stats::qnorm(1 - x$DIDparams$alp/2) * x$se.egt,
      point.conf.high = x$att.egt + stats::qnorm(1 - x$DIDparams$alp/2) * x$se.egt)
  }

  if(x$type == "simple"){
    out <- data.frame(
      type      = x$type,
      term      = 'ATT(simple average)',
      estimate  = x$overall.att,
      std.error = x$overall.se,
      statistic = x$overall.att / x$overall.se,
      p.value   = 2 * (1 - stats::pnorm(abs(x$overall.att / x$overall.se))),
      conf.low  = x$overall.att - stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se,
      conf.high = x$overall.att + stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se,
      point.conf.low  = x$overall.att - stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se,
      point.conf.high = x$overall.att + stats::qnorm(1 - x$DIDparams$alp/2) * x$overall.se)
  }
  
  out
}

#' glance model characteristics from AGGTEobj objects
#'
#' @importFrom generics glance
#' @param x a model of class AGGTEobj produced by the [aggte()] function
#' @param ... other arguments passed to methods
#' @export
glance.AGGTEobj<- function(x, ...) {
  if(x$DIDparams$faster_mode) {
    out <- data.frame(
      type          = x$type,
      nobs          = x$DIDparams$id_count,
      ngroup        = nrow(x$DIDparams$cohort_counts),
      ntime         = x$DIDparams$time_periods_count,
      control.group = x$DIDparams$control_group,
      est.method    = x$DIDparams$est_method)  
  } else {
    out <- data.frame(
      type          = x$type,
      nobs          = x$DIDparams$n,
      ngroup        = x$DIDparams$nG,
      ntime         = x$DIDparams$nT,
      control.group = x$DIDparams$control_group,
      est.method    = x$DIDparams$est_method)
  }

  out
}
