#' @title AGGTEobj
#'
#' @description Objects of this class hold results on aggregated
#'  group-time average treatment effects
#'
#' @title Aggregate Treatment Effect Parameters Object
#'
#' @description An object for holding aggregated treatment effect parameters.
#'
#' @inheritParams aggte
#' @inheritParams compute.aggte
#' @param overall.att The estimated overall ATT
#' @param overall.se Standard error for overall ATT
#' @param egt Holds the length of exposure (for dynamic effects), the
#'  group (for selective treatment timing), or the time period (for calendar
#'  time effects)
#' @param att.egt The ATT specific to egt
#' @param se.egt The standard error specific to egt
#' @param crit.val.egt A critical value for computing uniform confidence
#'  bands for dynamic effects, selective treatment timing, or time period
#'  effects.
#' @param inf.function The influence function of the chosen aggregated parameters
#' @param DIDparams A DIDparams object
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
                     inf.function = NULL,
                     min_e = NULL,
                     max_e = NULL,
                     balance_e = NULL,
                     call=NULL,
                     DIDparams=NULL) {

  out <- list(overall.att = overall.att,
              overall.se = overall.se,
              type = type,
              egt = egt,
              att.egt = att.egt,
              se.egt = se.egt,
              crit.val.egt = crit.val.egt,
              inf.function = inf.function,
              min_e = min_e,
              max_e = max_e,
              balance_e = balance_e,
              call = call,
              DIDparams = DIDparams)

  class(out)  <- "AGGTEobj"
  out
}

#' @title Summary Aggregate Treatment Effect Parameter Objects
#'
#' @description A function to summarize aggregated treatment effect parameters.
#'
#' @param object an \code{AGGTEobj} object
#' @param ... other arguments
#'
#' @export
summary.AGGTEobj <- function(object, ...) {

  # call
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")

  #citation
  citation()
  cat("\n")

  # overall estimates
  alp <- object$DIDparams$alp
  pointwise_cval <- qnorm(1-alp/2)
  overall_cband_upper <- object$overall.att + pointwise_cval*object$overall.se
  overall_cband_lower <- object$overall.att - pointwise_cval*object$overall.se
  out1 <- cbind.data.frame(object$overall.att, object$overall.se, overall_cband_lower, overall_cband_upper)
  out1 <- round(out1, 4)
  overall_sig <- (overall_cband_upper < 0) | (overall_cband_lower > 0)
  overall_sig[is.na(overall_sig)] <- FALSE
  overall_sig_text <- ifelse(overall_sig, "*", "")
  out1 <- cbind.data.frame(out1, overall_sig_text)
  cat("\n")
  #cat("Overall ATT:  \n")
  if (object$type=="dynamic") cat("Overall summary of ATT\'s based on event-study/dynamic aggregation:  \n")
  if (object$type=="group") cat("Overall summary of ATT\'s based on group/cohort aggregation:  \n")
  if (object$type=="calendar") cat("Overall summary of ATT\'s based on calendar time aggregation:  \n")
  colnames(out1) <- c("ATT","   Std. Error", paste0("    [ ",100*(1-object$DIDparams$alp),"% "), "Conf. Int.]","")
  print(out1, row.names=FALSE)
  cat("\n\n")

  # handle cases depending on type
  if (object$type %in% c("group","dynamic","calendar")) {

    # header
    if (object$type=="dynamic") { c1name <- "Event time"; cat("Dynamic Effects:") }
    if (object$type=="group") { c1name <- "Group"; cat("Group Effects:") }
    if (object$type=="calendar") { c1name <- "Time"; cat("Time Effects:") }

    cat("\n")
    cband_text1a <- paste0(100*(1-object$DIDparams$alp),"% ")
    cband_text1b <- ifelse(object$DIDparams$bstrap,
                           ifelse(object$DIDparams$cband, "Simult. ", "Pointwise "),
                           "Pointwise ")
    cband_text1 <- paste0("[", cband_text1a, cband_text1b)

    cband_lower <- object$att.egt - object$crit.val.egt*object$se.egt
    cband_upper <- object$att.egt + object$crit.val.egt*object$se.egt

    sig <- (cband_upper < 0) | (cband_lower > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")

    out2 <- cbind.data.frame(object$egt, object$att.egt, object$se.egt, cband_lower, cband_upper)
    out2 <- round(out2, 4)
    out2 <- cbind.data.frame(out2, sig_text)

    colnames(out2) <- c(c1name, "Estimate","Std. Error", cband_text1, "Conf. Band]", "")
    print(out2, row.names=FALSE, justify = "centre")
  }
  cat("---\n")
  cat("Signif. codes: `*' confidence band does not cover 0")
  cat("\n\n")

  # set control group text
  control_group <- object$DIDparams$control_group
  control_group_text <- NULL
  if (control_group == "nevertreated") {
    control_group_text <- "Never Treated"
  } else if (control_group == "notyettreated") {
    control_group_text <- "Not Yet Treated"
  }

  if (!is.null(control_group)) {
    cat("Control Group:  ")
    cat(control_group_text)
    cat(",  ")
  }

  # anticipation periods
  cat("Anticipation Periods:  ")
  cat(object$DIDparams$anticipation)
  cat("\n")

  # estimation method text
  est_method <- object$DIDparams$est_method
  if ( is(est_method,"character") ) {
    est_method_text <- est_method
    if (est_method == "dr") {
      est_method_text <- "Doubly Robust"
    } else if (est_method == "ipw") {
      est_method_text <- "Inverse Probability Weighting"
    } else if (est_method == "reg") {
      est_method_text <- "Outcome Regression"
    }

    cat("Estimation Method:  ")
    cat(est_method_text)
    cat("\n")
  }
}

#' @title print.AGGTEobj
#'
#' @description prints value of a \code{AGGTEobj} object
#'
#' @param x a \code{AGGTEobj} object
#' @param ... extra arguments
#'
#' @export
print.AGGTEobj <- function(x,...) {
  summary.AGGTEobj(x,...)
}
