#' @title AGGTEobj
#'
#' @description Objects of this class hold results on aggregated
#'  group-time average treatment effects
#'
#' @title Aggregate Treatment Effect Parameters Object
#'
#' @description An object for holding aggregated treatment effect parameters.
#'
#' @param overall.att The estimated overall ATT
#' @param overall.se Standard error for overall ATT
#' @param type Which type of aggregated treatment effect parameter.
#'  Possibilities are "simple" (the default), "dynamic" (for dynamic effects /
#'  event studies), "selective" (for selective treatment timing / group
#'  specific treatment effects), and "calendar" (for time effects)
#' @param egt Holds the length of exposure (for dynamic effects), the
#'  group (for selective treatment timing), or the time period (for calendar
#'  time effects)
#' @param att.egt The ATT specific to egt
#' @param se.egt The standard error specific to egt
#' @param crit.val.egt A critical value for computing uniform confidence
#'  bands for dynamic effects, selective treatment timing, or time period
#'  effects.
#' @param inf.function The influence function of the chosen aggregated parameters
#'
#' @return an AGGTEobj
#' @export
AGGTEobj <- function(overall.att = NULL,
                     overall.se = NULL,
                     type = "simple",
                     egt = NULL,
                     att.egt = NULL,
                     se.egt = NULL,
                     crit.val.egt = NULL,
                     inf.function = NULL) {

  out <- list(overall.att = overall.att,
              overall.se = overall.se,
              type = type,
              egt = egt,
              att.egt = att.egt,
              se.egt = se.egt,
              crit.val.egt = crit.val.egt,
              inf.function = inf.function)

  class(out)  <- "AGGTEobj"
  out
}

#' @title Summary Aggregate Treatment Effect Parameter Objects
#'
#' @description A function to summarize aggregated treatment effect parameters.
#'
#' @param object an AGGTEobj object
#' @param ... other arguments
#'
#' @export
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
