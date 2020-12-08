#' tidy results
#' @export
generics::tidy

#' tidy results
#' @export
generics::glance

#' tidy results from MP objects
#' 
#' @importFrom generics tidy
#' @param x a model of class MP produced by the `att_gt` function
#' @export
tidy.MP <- function(x, ...) {
  out <- data.frame(
    group     = x$group,
    term      = x$t,
    estimate  = x$att,
    std.error = x$se,
    conf.low  = x$att - x$c * x$se,
    conf.high = x$att + x$c * x$se)
  out
}

#' glance model characteristics from MP objects
#' 
#' @importFrom generics glance
#' @param x a model of class MP produced by the `att_gt` function
#' @export
glance.MP <- function(x, ...) {
  out <- data.frame(
    nobs          = x$n,
    ngroup        = x$DIDparams$nG,
    ntime         = x$DIDparams$nT,
    control.group = x$DIDparams$control_group)
  out
}
