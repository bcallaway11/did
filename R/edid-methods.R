# edid-methods.R
# S3 methods for the edid_fit class.

#' Print method for edid_fit objects
#'
#' Calls \code{summary.edid_fit()} and returns the object invisibly.
#'
#' @param x an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#'
#' @return \code{x} invisibly
#' @export
print.edid_fit <- function(x, ...) {
  summary(x, ...)
  invisible(x)
}

#' Summary method for edid_fit objects
#'
#' Prints a structured summary of the EDiD estimation results including
#' metadata, the cell-level ATT table, and aggregated estimates.
#'
#' @param object an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#'
#' @return \code{object} invisibly
#' @export
summary.edid_fit <- function(object, ...) {
  cat("\n=== Efficient Difference-in-Differences (EDiD) ===\n\n")

  # Metadata
  cat(sprintf("PT assumption     : %s\n", object$pt_assumption))
  cat(sprintf("Control group     : %s\n", object$control_group))
  cat(sprintf("Anticipation      : %d\n", object$anticipation))
  cat(sprintf("Units (n)         : %d\n", object$n))
  cat(sprintf("Time periods (T)  : %d\n", object$T_periods))
  cat(sprintf("Treatment cohorts : %s\n",
              paste(object$treatment_groups, collapse = ", ")))
  cat(sprintf("Inference type    : %s\n", object$inference_type))
  if (!is.null(object$cluster)) {
    cat(sprintf("Cluster variable  : %s\n", object$cluster))
  }
  cat(sprintf("Significance (alpha): %.3f\n\n", object$alpha))

  # Cell-level ATT table (compact)
  cat("--- Cell-level ATT(g, t) ---\n")
  if (!is.null(object$att_gt) && nrow(object$att_gt) > 0L) {
    att_print <- object$att_gt
    att_print$pre <- ifelse(att_print$is_pre, "pre", "post")
    # Round numeric columns for display
    fmt_num <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))
    out <- data.frame(
      group    = att_print$group,
      time     = att_print$time,
      type     = att_print$pre,
      att      = fmt_num(att_print$att),
      se       = fmt_num(att_print$se),
      ci_lower = fmt_num(att_print$ci_lower),
      ci_upper = fmt_num(att_print$ci_upper),
      p_value  = fmt_num(att_print$p_value),
      stringsAsFactors = FALSE
    )
    print(out, row.names = FALSE, quote = FALSE)
  } else {
    cat("  (no cells)\n")
  }

  # Overall ATT
  if (!is.null(object$overall)) {
    cat("\n--- Overall ATT ---\n")
    ov <- object$overall
    cat(sprintf("  ATT = %.4f  SE = %s  CI = [%s, %s]  p = %s\n",
                ov$att,
                .fmt_or_na(ov$se),
                .fmt_or_na(ov$ci_lower),
                .fmt_or_na(ov$ci_upper),
                .fmt_or_na(ov$p_value)))
  }

  # Event-study
  if (!is.null(object$event_study) && length(object$event_study) > 0L) {
    cat("\n--- Event-Study ATT(e) ---\n")
    for (es in object$event_study) {
      cat(sprintf("  e = %3g: ATT = %.4f  SE = %s  p = %s\n",
                  es$e,
                  es$att,
                  .fmt_or_na(es$se),
                  .fmt_or_na(es$p_value)))
    }
  }

  # Group
  if (!is.null(object$group) && length(object$group) > 0L) {
    cat("\n--- Group ATT(g) ---\n")
    for (gr in object$group) {
      cat(sprintf("  g = %g: ATT = %.4f  SE = %s  p = %s\n",
                  gr$group,
                  gr$att,
                  .fmt_or_na(gr$se),
                  .fmt_or_na(gr$p_value)))
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
#' @param which character: one of \code{"att_gt"}, \code{"overall"},
#'   \code{"event_study"}, \code{"group"}
#' @param ... additional arguments (ignored)
#'
#' @return named numeric vector of ATT estimates
#' @export
coef.edid_fit <- function(
  object,
  which = c("att_gt", "overall", "event_study", "group"),
  ...
) {
  which <- match.arg(which)
  switch(which,
    att_gt = {
      df  <- object$att_gt
      nms <- paste0("ATT(", df$group, ",", df$time, ")")
      stats::setNames(df$att, nms)
    },
    overall = {
      c(overall = object$overall$att)
    },
    event_study = {
      vals <- vapply(object$event_study, function(x) x$att, numeric(1L))
      es   <- object$event_study
      nms  <- vapply(es, function(x) paste0("e=", x$e), character(1L))
      stats::setNames(vals, nms)
    },
    group = {
      vals <- vapply(object$group, function(x) x$att, numeric(1L))
      nms  <- vapply(object$group, function(x) paste0("g=", x$group), character(1L))
      stats::setNames(vals, nms)
    }
  )
}

#' Extract variance-covariance matrix from an edid_fit object
#'
#' Returns the outer product of aggregated EIF vectors, scaled by \eqn{1/n^2}.
#' When bootstrap inference is used, returns a diagonal matrix of bootstrap
#' variance estimates.
#'
#' @param object an \code{edid_fit} object
#' @param which character: one of \code{"att_gt"}, \code{"overall"},
#'   \code{"event_study"}, \code{"group"}
#' @param ... additional arguments (ignored)
#'
#' @return square numeric matrix
#' @export
vcov.edid_fit <- function(
  object,
  which = c("att_gt", "overall", "event_study", "group"),
  ...
) {
  which <- match.arg(which)
  n     <- object$n

  if (which == "overall") {
    eif_v <- object$overall$eif_agg
    if (is.null(eif_v)) return(matrix(NA_real_, 1L, 1L))
    return(matrix(sum(eif_v^2) / n^2, nrow = 1L, ncol = 1L,
                  dimnames = list("overall", "overall")))
  }

  if (which == "event_study") {
    es_list <- object$event_study
    if (is.null(es_list) || length(es_list) == 0L) return(matrix(NA_real_, 0L, 0L))
    nms  <- vapply(es_list, function(x) paste0("e=", x$e), character(1L))
    K    <- length(es_list)
    vcv  <- matrix(NA_real_, K, K, dimnames = list(nms, nms))
    for (j in seq_len(K)) {
      eif_j <- es_list[[j]]$eif_agg
      if (is.null(eif_j)) next
      for (k in seq_len(K)) {
        eif_k <- es_list[[k]]$eif_agg
        if (is.null(eif_k)) next
        vcv[j, k] <- sum(eif_j * eif_k) / n^2
      }
    }
    return(vcv)
  }

  if (which == "group") {
    gr_list <- object$group
    if (is.null(gr_list) || length(gr_list) == 0L) return(matrix(NA_real_, 0L, 0L))
    nms  <- vapply(gr_list, function(x) paste0("g=", x$group), character(1L))
    K    <- length(gr_list)
    vcv  <- matrix(NA_real_, K, K, dimnames = list(nms, nms))
    for (j in seq_len(K)) {
      eif_j <- gr_list[[j]]$eif_agg
      if (is.null(eif_j)) next
      for (k in seq_len(K)) {
        eif_k <- gr_list[[k]]$eif_agg
        if (is.null(eif_k)) next
        vcv[j, k] <- sum(eif_j * eif_k) / n^2
      }
    }
    return(vcv)
  }

  # att_gt: use cell-level EIFs from eif_matrix if available
  if (!is.null(object$eif)) {
    eif_mat <- object$eif
    K   <- ncol(eif_mat)
    vcv <- (t(eif_mat) %*% eif_mat) / n^2
    df  <- object$att_gt
    nms <- paste0("ATT(", df$group, ",", df$time, ")")
    dimnames(vcv) <- list(nms, nms)
    return(vcv)
  }
  # Fallback: diagonal from stored SEs
  df  <- object$att_gt
  nms <- paste0("ATT(", df$group, ",", df$time, ")")
  K   <- nrow(df)
  vcv <- diag(df$se^2, nrow = K)
  dimnames(vcv) <- list(nms, nms)
  vcv
}

#' Coerce edid_fit to a data.frame
#'
#' @param x an \code{edid_fit} object
#' @param which character: one of \code{"att_gt"}, \code{"overall"},
#'   \code{"event_study"}, \code{"group"}
#' @param ... additional arguments (ignored)
#'
#' @return data.frame
#' @export
as.data.frame.edid_fit <- function(
  x,
  which = c("att_gt", "overall", "event_study", "group"),
  ...
) {
  which <- match.arg(which)
  switch(which,
    att_gt = {
      x$att_gt
    },
    overall = {
      ov <- x$overall
      data.frame(
        att      = ov$att,
        se       = ov$se,
        ci_lower = ov$ci_lower,
        ci_upper = ov$ci_upper,
        t_stat   = ov$t_stat,
        p_value  = ov$p_value,
        stringsAsFactors = FALSE
      )
    },
    event_study = {
      es_list <- x$event_study
      if (is.null(es_list) || length(es_list) == 0L) {
        return(data.frame(e = numeric(0L), att = numeric(0L), se = numeric(0L),
                          ci_lower = numeric(0L), ci_upper = numeric(0L),
                          p_value = numeric(0L)))
      }
      do.call(rbind, lapply(es_list, function(es) {
        data.frame(e = es$e, att = es$att, se = es$se,
                   ci_lower = es$ci_lower, ci_upper = es$ci_upper,
                   t_stat = es$t_stat, p_value = es$p_value,
                   stringsAsFactors = FALSE)
      }))
    },
    group = {
      gr_list <- x$group
      if (is.null(gr_list) || length(gr_list) == 0L) {
        return(data.frame(group = numeric(0L), att = numeric(0L), se = numeric(0L),
                          ci_lower = numeric(0L), ci_upper = numeric(0L),
                          p_value = numeric(0L)))
      }
      do.call(rbind, lapply(gr_list, function(gr) {
        data.frame(group = gr$group, att = gr$att, se = gr$se,
                   ci_lower = gr$ci_lower, ci_upper = gr$ci_upper,
                   t_stat = gr$t_stat, p_value = gr$p_value,
                   stringsAsFactors = FALSE)
      }))
    }
  )
}
