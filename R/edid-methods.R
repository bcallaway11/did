# edid-methods.R
# S3 methods for the edid_fit class. The aggregations ($overall/$event_study/$group/$calendar/$simple)
# are did::AGGTEobj objects (built via aggte_edid -> did::aggte), so these methods read AGGTEobj fields
# (overall.att/overall.se, att.egt/se.egt/egt) and the per-element influence functions in $inf.function.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Pull the per-element (n x K) and overall (n) influence functions out of a did::AGGTEobj, regardless of
# aggregation type (compute.aggte names them dynamic./selective./calendar.inf.func.* and simple.att).
.edid_agg_if <- function(a) {
  if (is.null(a)) return(list(egt = NULL, overall = NULL))
  inf <- a$inf.function
  list(
    egt     = inf$dynamic.inf.func.e %||% inf$selective.inf.func.g %||% inf$calendar.inf.func.t,
    overall = inf$dynamic.inf.func %||% inf$selective.inf.func %||% inf$calendar.inf.func %||% inf$simple.att
  )
}

#' Print method for edid_fit objects
#'
#' Displays the ATT(g,t) table in the same style as \code{print.MP} / \code{summary.MP}, followed by
#' footer metadata.
#'
#' @param x an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#' @return \code{x} invisibly
#' @export
print.edid_fit <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Group-Time Average Treatment Effects:\n")

  alp <- x$alpha
  cband_text1a <- paste0(100 * (1 - alp), "% ")
  # Simultaneous iff a > qnorm crit was actually applied: the multiplier bootstrap ran (bstrap with the
  # multiplier cband) OR the analytic sup-t band was built (analytic cband with cband = TRUE).
  simult <- (isTRUE(x$bstrap) && identical(x$cband_method, "multiplier")) ||
            (identical(x$cband_method, "analytic") && isTRUE(x$cband))
  cband_text1b <- ifelse(simult, "Simult. ", "Pointwise ")
  cband_text1  <- paste0("[", cband_text1a, cband_text1b)

  att_df <- x$att_gt
  if (!is.null(att_df) && nrow(att_df) > 0L) {
    ci_lower <- att_df$ci_lower
    ci_upper <- att_df$ci_upper
    sig <- (ci_upper < 0) | (ci_lower > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")
    out <- cbind.data.frame(att_df$group, att_df$time, att_df$att, att_df$se, ci_lower, ci_upper)
    out <- round(out, 4)
    out <- cbind.data.frame(out, sig_text)
    colnames(out) <- c("Group", "Time", "ATT(g,t)", "Std. Error", cband_text1, "Conf. Band]", "")
    print(out, row.names = FALSE)
  } else {
    cat("  (no cells)\n")
  }

  cat("---\n")
  cat("Signif. codes: `*' confidence band does not cover 0")
  cat("\n\n")

  cat("Control Group:  "); cat("Never Treated"); cat(",  ")
  cat("Anticipation Periods:  "); cat(x$anticipation); cat("\n")
  cat("Estimation Method:  Efficient DiD (Chen, Sant'Anna & Xie 2025)\n")
  pt_text <- if (x$pt_assumption == "all") "PT-All" else "PT-Post"
  cat("PT Assumption:  "); cat(pt_text); cat("\n")
  invisible(x)
}

#' Summary method for edid_fit objects
#'
#' Prints the ATT(g,t) table (MP style) followed by the requested aggregations, each a
#' \code{did::AGGTEobj} printed with did's own \code{print.AGGTEobj}.
#'
#' @param object an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#' @return \code{object} invisibly
#' @export
summary.edid_fit <- function(object, ...) {
  print.edid_fit(object, ...)
  ov <- object$overall
  if (!is.null(ov)) {
    cat(sprintf("\nOverall ATT (%s):  %s  (SE %s)\n", ov$type %||% "overall",
                .fmt_or_na(ov$overall.att), .fmt_or_na(ov$overall.se)))
  }
  for (nm in c("event_study", "group", "calendar")) {
    a <- object[[nm]]
    if (!is.null(a) && inherits(a, "AGGTEobj")) {
      cat(sprintf("\n--- %s ---\n", nm)); print(a)
    }
  }
  cat("\n")
  invisible(object)
}

# Internal formatting helper
.fmt_or_na <- function(x) {
  if (is.null(x) || !is.finite(x)) "NA" else sprintf("%.4f", x)
}

#' Extract ATT coefficients from an edid_fit object
#'
#' @param object an \code{edid_fit} object
#' @param which character: one of \code{"att_gt"}, \code{"overall"}, \code{"event_study"}, \code{"group"}
#' @param ... additional arguments (ignored)
#' @return named numeric vector of ATT estimates
#' @export
coef.edid_fit <- function(object, which = c("att_gt", "overall", "event_study", "group"), ...) {
  which <- match.arg(which)
  switch(which,
    att_gt = {
      df <- object$att_gt
      stats::setNames(df$att, paste0("ATT(", df$group, ",", df$time, ")"))
    },
    overall = c(overall = object$overall$overall.att),
    event_study = {
      a <- object$event_study
      if (is.null(a)) return(numeric(0L))
      stats::setNames(a$att.egt, paste0("e=", a$egt))
    },
    group = {
      a <- object$group
      if (is.null(a)) return(numeric(0L))
      stats::setNames(a$att.egt, paste0("g=", a$egt))
    }
  )
}

#' Extract variance-covariance matrix from an edid_fit object
#'
#' For \code{which = "att_gt"} returns the cluster-robust (or i.i.d.) covariance of the cell-level
#' ATT(g,t)'s from the stored influence functions. For the aggregations it returns the covariance implied
#' by the corresponding \code{did::AGGTEobj}'s aggregate influence functions (cluster-robust when
#' \code{clustervars} was set).
#'
#' @param object an \code{edid_fit} object
#' @param which character: one of \code{"att_gt"}, \code{"overall"}, \code{"event_study"}, \code{"group"}
#' @param ... additional arguments (ignored)
#' @return square numeric matrix
#' @export
vcov.edid_fit <- function(object, which = c("att_gt", "overall", "event_study", "group"), ...) {
  which <- match.arg(which)
  n  <- object$n
  ci <- object$cluster_indices

  # cluster-robust (or i.i.d.) covariance of the columns of an n x K influence-function matrix
  .cross_mat <- function(M) {
    M <- as.matrix(M)
    if (is.null(ci)) return(crossprod(M) / n^2)
    G <- length(unique(ci))
    if (G <= 1L) return(matrix(NA_real_, ncol(M), ncol(M)))
    CS <- rowsum(M, ci)
    (G / (G - 1)) * crossprod(CS) / n^2
  }

  if (which == "att_gt") {
    df  <- object$att_gt
    nms <- paste0("ATT(", df$group, ",", df$time, ")")
    if (!is.null(object$eif)) {
      v <- .cross_mat(object$eif); dimnames(v) <- list(nms, nms); return(v)
    }
    v <- diag(df$se^2, nrow = nrow(df)); dimnames(v) <- list(nms, nms); return(v)
  }

  if (which == "overall") {
    g <- .edid_agg_if(object$overall)
    if (is.null(g$overall)) return(matrix(NA_real_, 1L, 1L))
    return(matrix(.cross_mat(matrix(g$overall, ncol = 1L)), 1L, 1L,
                  dimnames = list("overall", "overall")))
  }

  a   <- if (which == "event_study") object$event_study else object$group
  g   <- .edid_agg_if(a)
  if (is.null(a) || is.null(g$egt)) return(matrix(NA_real_, 0L, 0L))
  pre <- if (which == "event_study") "e=" else "g="
  nms <- paste0(pre, a$egt)
  v   <- .cross_mat(g$egt); dimnames(v) <- list(nms, nms); v
}

#' Coerce edid_fit to a data.frame
#'
#' @param x an \code{edid_fit} object
#' @param row.names ignored; included for S3 generic consistency
#' @param optional ignored; included for S3 generic consistency
#' @param ... not used
#' @param which character: one of \code{"att_gt"}, \code{"overall"}, \code{"event_study"}, \code{"group"}
#' @return data.frame
#' @export
as.data.frame.edid_fit <- function(x, row.names = NULL, optional = FALSE, ...,
                                    which = c("att_gt", "overall", "event_study", "group")) {
  which <- match.arg(which)
  z <- stats::qnorm(1 - x$alpha / 2)
  switch(which,
    att_gt = x$att_gt,
    overall = {
      a <- x$overall
      data.frame(att = a$overall.att, se = a$overall.se,
                 ci_lower = a$overall.att - z * a$overall.se,
                 ci_upper = a$overall.att + z * a$overall.se, stringsAsFactors = FALSE)
    },
    event_study = {
      a <- x$event_study
      if (is.null(a)) return(data.frame(e = numeric(0L), att = numeric(0L), se = numeric(0L),
                                        ci_lower = numeric(0L), ci_upper = numeric(0L)))
      data.frame(e = a$egt, att = a$att.egt, se = a$se.egt,
                 ci_lower = a$att.egt - z * a$se.egt, ci_upper = a$att.egt + z * a$se.egt,
                 stringsAsFactors = FALSE)
    },
    group = {
      a <- x$group
      if (is.null(a)) return(data.frame(group = numeric(0L), att = numeric(0L), se = numeric(0L),
                                        ci_lower = numeric(0L), ci_upper = numeric(0L)))
      data.frame(group = a$egt, att = a$att.egt, se = a$se.egt,
                 ci_lower = a$att.egt - z * a$se.egt, ci_upper = a$att.egt + z * a$se.egt,
                 stringsAsFactors = FALSE)
    }
  )
}
