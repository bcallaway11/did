#' @title MP
#'
#' @description Multi-period objects that hold results for group-time average treatment effects
#'
#' @param group which group (defined by period first treated) an group-time average treatment effect is for
#' @param t which time period a group-time average treatment effect is for
#' @param att the group-average treatment effect for group \code{group} and time
#'  period \code{t}
#' @param c simulataneous critical value if one is obtaining simultaneous confidence
#'  bands. Otherwise it reports the critival value based on pointwise normal
#'  approximation.
#' @param V the variance matrix for group-time average treatment effects
#' @param se standard errors for group-time average treatment effects
#' @param inffunc the influence function for estimating group-time average treatment effects
#' @param n the number of observations
#' @param W the Wald statistic for pre-testing the common trends assumption
#' @param Wpval the p-value of the Wald statistic for pre-testing the
#'  common trends assumption
#' @param aggte an aggregate treatment effects object
#' @param alp the significance level, default is 0.05
#' @param DIDparams a \code{\link{DIDparams}} object.  A way to optionally return the parameters
#'  of the call to \code{\link{att_gt}} or \code{\link{conditional_did_pretest}}.
#'
#' @return MP object
#' @export
MP <- function(group, t, att, V, se, c, inffunc, n=NULL, W=NULL, Wpval=NULL, aggte=NULL, alp = 0.05, DIDparams=NULL) {
  out <- list(group=group, t=t, att=att, V=V, se=se, c=c, inffunc=inffunc, n=n, W=W, Wpval=Wpval, aggte=aggte, alp = alp, DIDparams=DIDparams)
  class(out) <- "MP"
  out
}

#' @title summary.MP
#'
#' @description prints a summary of a \code{MP} object
#'
#' @param object an \code{MP} object
#' @param ... extra arguments
#'
#' @export
summary.MP <- function(object, ...) {
  mpobj <- object
  out <- cbind(mpobj$group, mpobj$t, mpobj$att, mpobj$se)#sqrt(diag(mpobj$V)/mpobj$n))
  citation()
  colnames(out) <- c("group", "time", "att","se")
  cat("\n")
  print(knitr::kable(out))
  cat("\n\n")
  if (!is.null(mpobj$Wpval)) {
    cat("P-value for pre-test of parallel trends assumption:  ")
    cat(as.character(mpobj$Wpval))
    cat("\n\n")
  }
}

#' @title print.MP
#'
#' @description prints value of a \code{MP} object
#'
#' @param x a \code{MP} object
#' @param ... extra arguments
#'
#' @export
print.MP <- function(x,...) {
  summary.MP(x,...)
}
