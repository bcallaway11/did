# edid-aggte.R
# Aggregation function for edid_fit objects, mirroring the aggte() interface.

#' Aggregate edid_fit estimates
#'
#' Provides the same user-facing interface as \code{\link[did]{aggte}} but
#' accepts an \code{edid_fit} object produced by \code{\link{edid}}.
#'
#' @param edid_fit_obj An \code{edid_fit} object returned by \code{edid()}.
#' @param type Character scalar: aggregation type. One of
#'   \code{"simple"} (overall ATT), \code{"dynamic"} (event-study),
#'   \code{"group"} (cohort-level ATT), or \code{"calendar"} (not implemented).
#' @param balance_e Integer or \code{NULL}: if not \code{NULL}, restricts the
#'   dynamic aggregation to relative times in
#'   \eqn{[-\text{balance\_e}, \text{balance\_e}]}.
#' @param min_e Numeric: minimum relative time to include in dynamic output.
#'   Default \code{-Inf}.
#' @param max_e Numeric: maximum relative time to include in dynamic output.
#'   Default \code{Inf}.
#' @param na.rm Logical: whether to remove NA ATT entries before aggregating.
#'   Default \code{FALSE}.
#'
#' @return An S3 object of class \code{c("AGGTEobj_edid", "list")} with fields
#'   matching \code{AGGTEobj} where possible:
#'   \describe{
#'     \item{\code{att.egt}}{Vector of ATT estimates for each index.}
#'     \item{\code{se.egt}}{Vector of standard errors.}
#'     \item{\code{egt}}{Vector of indices (relative time, group, etc.).}
#'     \item{\code{type}}{The aggregation type string.}
#'     \item{\code{overall.att}}{Scalar overall ATT.}
#'     \item{\code{overall.se}}{Scalar overall SE.}
#'     \item{\code{alp}}{Significance level used.}
#'     \item{\code{call}}{The matched call.}
#'   }
#'
#' @seealso \code{\link{edid}}, \code{\link[did]{aggte}}
#' @export
aggte_edid <- function(
  edid_fit_obj,
  type      = c("simple", "dynamic", "group", "calendar"),
  balance_e = NULL,
  min_e     = -Inf,
  max_e     = Inf,
  na.rm     = FALSE
) {
  mc   <- match.call()
  type <- match.arg(type)

  if (!inherits(edid_fit_obj, "edid_fit")) {
    stop("`edid_fit_obj` must be an object of class `edid_fit` returned by edid().")
  }

  alp <- edid_fit_obj$alpha

  # Helper: extract overall ATT + SE from the edid_fit object
  .get_overall <- function(obj) {
    ov  <- obj$overall
    att <- if (!is.null(ov)) ov$att else NA_real_
    se  <- if (!is.null(ov)) ov$se  else NA_real_
    list(att = att, se = se)
  }

  if (type == "calendar") {
    stop("aggte_edid() does not support type = \"calendar\". ",
         "edid() does not compute calendar-time treatment effects.")
  }

  # ------------------------------------------------------------------
  # simple: return overall ATT already computed
  # ------------------------------------------------------------------
  if (type == "simple") {
    ov_info <- .get_overall(edid_fit_obj)
    out <- list(
      att.egt    = ov_info$att,
      se.egt     = ov_info$se,
      egt        = NA_real_,
      type       = type,
      overall.att = ov_info$att,
      overall.se  = ov_info$se,
      alp        = alp,
      call       = mc
    )
    class(out) <- c("AGGTEobj_edid", "list")
    return(out)
  }

  # ------------------------------------------------------------------
  # dynamic: event-study, filter by min_e / max_e / balance_e
  # ------------------------------------------------------------------
  if (type == "dynamic") {
    es_list <- edid_fit_obj$event_study
    if (is.null(es_list) || length(es_list) == 0L) {
      stop("No event-study results in edid_fit_obj. ",
           "Re-run edid() with aggregate = \"event_study\" or \"all\".")
    }

    # Convert to flat vectors
    e_vals  <- vapply(es_list, function(x) x$e,   numeric(1L))
    att_vec <- vapply(es_list, function(x) x$att, numeric(1L))
    se_vec  <- vapply(es_list, function(x) x$se,  numeric(1L))

    # Apply balance_e filter first (symmetric window)
    if (!is.null(balance_e)) {
      keep_bal <- (e_vals >= -balance_e) & (e_vals <= balance_e)
      e_vals   <- e_vals[keep_bal]
      att_vec  <- att_vec[keep_bal]
      se_vec   <- se_vec[keep_bal]
    }

    # Apply min_e / max_e filter
    keep_range <- (e_vals >= min_e) & (e_vals <= max_e)
    e_vals  <- e_vals[keep_range]
    att_vec <- att_vec[keep_range]
    se_vec  <- se_vec[keep_range]

    # Optionally drop NAs
    if (na.rm) {
      keep_na <- !is.na(att_vec)
      e_vals  <- e_vals[keep_na]
      att_vec <- att_vec[keep_na]
      se_vec  <- se_vec[keep_na]
    }

    ov_info <- .get_overall(edid_fit_obj)
    out <- list(
      att.egt     = att_vec,
      se.egt      = se_vec,
      egt         = e_vals,
      type        = type,
      overall.att = ov_info$att,
      overall.se  = ov_info$se,
      alp         = alp,
      call        = mc
    )
    class(out) <- c("AGGTEobj_edid", "list")
    return(out)
  }

  # ------------------------------------------------------------------
  # group: cohort-level ATTs
  # ------------------------------------------------------------------
  if (type == "group") {
    gr_list <- edid_fit_obj$group
    if (is.null(gr_list) || length(gr_list) == 0L) {
      stop("No group results in edid_fit_obj. ",
           "Re-run edid() with aggregate = \"group\" or \"all\".")
    }

    g_vals  <- vapply(gr_list, function(x) x$group, numeric(1L))
    att_vec <- vapply(gr_list, function(x) x$att,   numeric(1L))
    se_vec  <- vapply(gr_list, function(x) x$se,    numeric(1L))

    if (na.rm) {
      keep_na <- !is.na(att_vec)
      g_vals  <- g_vals[keep_na]
      att_vec <- att_vec[keep_na]
      se_vec  <- se_vec[keep_na]
    }

    ov_info <- .get_overall(edid_fit_obj)
    out <- list(
      att.egt     = att_vec,
      se.egt      = se_vec,
      egt         = g_vals,
      type        = type,
      overall.att = ov_info$att,
      overall.se  = ov_info$se,
      alp         = alp,
      call        = mc
    )
    class(out) <- c("AGGTEobj_edid", "list")
    return(out)
  }
}

#' Print method for AGGTEobj_edid objects
#'
#' Prints aggregated treatment effects in a format similar to
#' \code{print.AGGTEobj}.
#'
#' @param x an \code{AGGTEobj_edid} object
#' @param ... additional arguments (currently ignored)
#'
#' @return \code{x} invisibly
#' @export
print.AGGTEobj_edid <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  alp <- x$alp
  pointwise_cval  <- stats::qnorm(1 - alp / 2)

  # Overall ATT summary
  ov_att <- x$overall.att
  ov_se  <- x$overall.se
  if (!is.null(ov_att) && !is.na(ov_att)) {
    ov_lo  <- ov_att - pointwise_cval * ov_se
    ov_hi  <- ov_att + pointwise_cval * ov_se
    ov_sig <- (ov_hi < 0) | (ov_lo > 0)
    if (is.na(ov_sig)) ov_sig <- FALSE
    ov_sig_text <- if (ov_sig) "*" else ""

    if (x$type == "dynamic") {
      cat("Overall summary of ATT's based on event-study/dynamic aggregation:  \n")
    } else if (x$type == "group") {
      cat("Overall summary of ATT's based on group/cohort aggregation:  \n")
    } else {
      cat("Overall ATT:  \n")
    }

    out1 <- cbind.data.frame(round(ov_att, 4), round(ov_se, 4),
                              round(ov_lo, 4), round(ov_hi, 4),
                              ov_sig_text)
    colnames(out1) <- c("ATT", "   Std. Error",
                        paste0("    [ ", 100 * (1 - alp), "% "),
                        "Conf. Int.]", "")
    print(out1, row.names = FALSE)
    cat("\n\n")
  }

  # Per-index table for dynamic / group
  if (x$type %in% c("dynamic", "group")) {
    if (x$type == "dynamic") {
      c1name <- "Event time"
      cat("Dynamic Effects:\n")
    } else {
      c1name <- "Group"
      cat("Group Effects:\n")
    }

    cband_text1 <- paste0("[", 100 * (1 - alp), "% Pointwise ")

    cband_lower <- x$att.egt - pointwise_cval * x$se.egt
    cband_upper <- x$att.egt + pointwise_cval * x$se.egt

    sig <- (cband_upper < 0) | (cband_lower > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")

    out2 <- cbind.data.frame(x$egt, x$att.egt, x$se.egt, cband_lower, cband_upper)
    out2 <- round(out2, 4)
    out2 <- cbind.data.frame(out2, sig_text)
    colnames(out2) <- c(c1name, "Estimate", "Std. Error",
                        cband_text1, "Conf. Band]", "")
    print(out2, row.names = FALSE, justify = "centre")
  }

  cat("---\n")
  cat("Signif. codes: `*' confidence band does not cover 0")
  cat("\n\n")
  cat("Estimation Method:  Efficient DiD (Chen, Sant'Anna & Xie 2025)\n")

  invisible(x)
}

#' Summary method for AGGTEobj_edid objects
#'
#' Delegates to \code{print.AGGTEobj_edid}.
#'
#' @param object an \code{AGGTEobj_edid} object
#' @param ... additional arguments (currently ignored)
#'
#' @return \code{object} invisibly
#' @export
summary.AGGTEobj_edid <- function(object, ...) {
  print.AGGTEobj_edid(object, ...)
  invisible(object)
}
