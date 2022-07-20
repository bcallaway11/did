#' @title MP
#'
#' @description Multi-period objects that hold results for group-time average treatment effects
#'
#' @param group which group (defined by period first treated) an group-time average treatment effect is for
#' @param t which time period a group-time average treatment effect is for
#' @param att the group-average treatment effect for group \code{group} and time
#'  period \code{t}
#' @param c simultaneous critical value if one is obtaining simultaneous confidence
#'  bands. Otherwise it reports the critical value based on pointwise normal
#'  approximation.
#' @param V_analytical Analytical estimator for the asymptotic variance-covariance matrix for group-time average treatment effects
#' @param se standard errors for group-time average treatment effects. If bootstrap is set to TRUE, this provides bootstrap-based se.
#' @param inffunc the influence function for estimating group-time average treatment effects
#' @param n the number of unique cross-sectional units (unique values of idname)
#' @param W the Wald statistic for pre-testing the common trends assumption
#' @param Wpval the p-value of the Wald statistic for pre-testing the
#'  common trends assumption
#' @param aggte an aggregate treatment effects object
#' @param alp the significance level, default is 0.05
#' @param DIDparams a [`DIDparams`] object.  A way to optionally return the parameters
#'  of the call to [att_gt()] or [conditional_did_pretest()].
#'
#' @return MP object
#' @export
MP <- function(group, t, att, V_analytical, se, c, inffunc, n=NULL, W=NULL, Wpval=NULL, aggte=NULL, alp = 0.05, DIDparams=NULL) {
  out <- list(group=group, t=t, att=att, V_analytical=V_analytical, se=se, c=c,
  inffunc=inffunc, n=n, W=W, Wpval=Wpval, aggte=aggte, alp = alp,
  DIDparams=DIDparams, call=DIDparams$call)
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

  # call
  cat("\n")
  cat("Call:\n")
  print(mpobj$DIDparams$call)
  cat("\n")

  # citation
  citation()
  cat("\n")

  # group time average treatment effects
  cat("Group-Time Average Treatment Effects:\n")

  cband_text1a <- paste0(100*(1-mpobj$alp),"% ")
  cband_text1b <- ifelse(mpobj$DIDparams$bstrap,
                         ifelse(mpobj$DIDparams$cband, "Simult. ", "Pointwise "),
                         "Pointwise ")
  cband_text1 <- paste0("[", cband_text1a, cband_text1b)

  cband_lower <- mpobj$att - mpobj$c*mpobj$se
  cband_upper <- mpobj$att + mpobj$c*mpobj$se

  sig <- (cband_upper < 0) | (cband_lower > 0)
  sig[is.na(sig)] <- FALSE
  sig_text <- ifelse(sig, "*", "")

  out <- cbind.data.frame(mpobj$group, mpobj$t, mpobj$att, mpobj$se, cband_lower, cband_upper)
  out <- round(out,4)
  out <- cbind.data.frame(out, sig_text)


  colnames(out) <- c("Group", "Time", "ATT(g,t)","Std. Error", cband_text1, "Conf. Band]", "")
  print(out, row.names=FALSE)
  cat("---\n")
  cat("Signif. codes: `*' confidence band does not cover 0")
  cat("\n\n")

  # report pre-test
  if (!is.null(mpobj$Wpval)) {
    cat("P-value for pre-test of parallel trends assumption:  ")
    cat(as.character(mpobj$Wpval))
    cat("\n")
  }



  # set control group text
  control_group <- mpobj$DIDparams$control_group
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
  cat(mpobj$DIDparams$anticipation)
  cat("\n")

  # estimation method text
  est_method <- mpobj$DIDparams$est_method
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
