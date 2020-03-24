#' @title AGGTE
#'
#' @description \code{AGGTE} class for aggregate treatment effects
#'
#' @param simple.att simple weighted average of group-time average treatment
#'  effects
#' @param simple.se the standard error for \code{simple.att}
#' @param dynamic.att aggregated group-time average treatment effects when
#'  there are dynamic treatment effects
#' @param dynamic.se the standard error for \code{dynamic.att}
#' @param dynamic.att.e aggregated group-time average treatment effects
#'  when there are dynamic treatment effects for each length of exposure
#'  to treatment
#' @param dynamic.se.e the standard error for \code{dynamic.att.e}
#' @param c.dynamic the (simultaneous) critical value \code{dynamic.att.e}
#' @param group vector of all groups
#' @param times vector of all times
#' @param dyn.inf.func.e influence function for event studies
#' @param simple.att.inf.func influence function for simple average of ATT(g,t)
#' @param dynamic.att.inf.func influence function for time-average of event-study
#' @param e the relative times used in the event-study

## AGGTE <- function(simple.att=NULL, simple.se=NULL,
##                   dynamic.att=NULL, dynamic.se=NULL,
##                   dynamic.att.e=NULL, dynamic.se.e=NULL, c.dynamic=NULL,
##                   dyn.inf.func.e = NULL,
##                   simple.att.inf.func = NULL,
##                   dynamic.att.inf.func = NULL,
##                   group=NULL, times=NULL,
##                   e = NULL) {

##   out <- list(simple.att=simple.att, simple.se=simple.se,
##               dynamic.att=dynamic.att, dynamic.se=dynamic.se,
##               dynamic.att.e=dynamic.att.e, dynamic.se.e=dynamic.se.e, c.dynamic=c.dynamic,
##               # Influence functions
##               dyn.inf.func.e = dyn.inf.func.e,
##               simple.att.inf.func = simple.att.inf.func,
##               dynamic.att.inf.func = dynamic.att.inf.func,

##               group=group,times=times, e = e)
##   class(out) <- "AGGTE"
##   out
## }

AGGTEobj <- function(overall.att=NULL,
                     overall.se=NULL,
                     type="simple",
                     egt=NULL,
                     att.egt=NULL,
                     se.egt=NULL,
                     crit.val.egt=NULL) {
  out <- list(overall.att=overall.att,
              overall.se=overall.se,
              type=type,
              egt=egt,
              att.egt=att.egt,
              se.egt=se.egt,
              crit.val.egt=crit.val.egt)
  class(out)  <- "AGGTEobj"
  out
}


summary.AGGTEobj <- function(object, ...) {
  out1 <- cbind(object$overall.att, object$overall.se)
  citation()
  cat("\n")
  cat("Overall ATT:  ")
  colnames(out1) <- c("att","se")
  cat("\n")
  print(knitr::kable(out1))
  cat("\n\n")

  # handle cases depending on type
  if (object$type %in% c("selective","dynamic","calendar")) {
  
    if (object$type=="dynamic") { c1name <- "event time"; cat("Dynamic Effects:") }
    if (object$type=="selective") { c1name <- "group"; cat("Group Effects:") }
    if (object$type=="calendar") { c1name <- "time"; cat("Time Effects:") } 
    out2 <- cbind(object$egt, object$att.egt, object$se.egt)
    colnames(out2) <- c(c1name, "att", "se")
    cat("\n")
    print(knitr::kable(out2))
  } 
}
