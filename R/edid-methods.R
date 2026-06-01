# edid-methods.R
# S3 methods for the edid_fit class.

#' Print method for edid_fit objects
#'
#' Displays the ATT(g,t) table in the same style as \code{print.MP} /
#' \code{summary.MP}, followed by footer metadata.
#'
#' @param x an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#'
#' @return \code{x} invisibly
#' @export
print.edid_fit <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Group-Time Average Treatment Effects:\n")

  alp <- x$alpha
  # CI label: pointwise for analytical; simult if bootstrap
  cband_text1a <- paste0(100 * (1 - alp), "% ")
  cband_text1b <- ifelse(isTRUE(x$bstrap), "Simult. ", "Pointwise ")
  cband_text1  <- paste0("[", cband_text1a, cband_text1b)

  att_df <- x$att_gt
  if (!is.null(att_df) && nrow(att_df) > 0L) {
    ci_lower <- att_df$ci_lower
    ci_upper <- att_df$ci_upper

    sig <- (ci_upper < 0) | (ci_lower > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")

    out <- cbind.data.frame(
      att_df$group,
      att_df$time,
      att_df$att,
      att_df$se,
      ci_lower,
      ci_upper
    )
    out <- round(out, 4)
    out <- cbind.data.frame(out, sig_text)
    colnames(out) <- c("Group", "Time", "ATT(g,t)", "Std. Error",
                       cband_text1, "Conf. Band]", "")
    print(out, row.names = FALSE)
  } else {
    cat("  (no cells)\n")
  }

  cat("---\n")
  cat("Signif. codes: `*' confidence band does not cover 0")
  cat("\n\n")

  # Control group footer
  cg <- x$control_group
  cg_text <- if (cg == "nevertreated")   "Never Treated"    else
             if (cg == "notyettreated") "Not Yet Treated"  else cg
  cat("Control Group:  ")
  cat(cg_text)
  cat(",  ")
  cat("Anticipation Periods:  ")
  cat(x$anticipation)
  cat("\n")

  cat("Estimation Method:  Efficient DiD (Chen, Sant'Anna & Xie 2025)\n")

  pt_text <- if (x$pt_assumption == "all") "PT-All" else "PT-Post"
  cat("PT Assumption:  ")
  cat(pt_text)
  cat("\n")

  invisible(x)
}

#' Summary method for edid_fit objects
#'
#' Prints a structured summary of the EDiD estimation results including
#' the ATT(g,t) table (in MP style) plus aggregated overall, event-study,
#' and group estimates.
#'
#' @param object an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#'
#' @return \code{object} invisibly
#' @export
summary.edid_fit <- function(object, ...) {
  # First, print the ATT(g,t) table in MP format
  print.edid_fit(object, ...)

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
#' Returns the IID asymptotic covariance: the outer product of the (per-unit)
#' aggregated influence-function vectors, scaled by \eqn{1/n^2}.
#'
#' @note This currently returns the IID covariance regardless of how \code{object}
#'   was fit. When \code{clustervars} or \code{bstrap = TRUE} was used, the SEs
#'   reported by \code{summary()} apply a cluster-robust / bootstrap adjustment that
#'   \code{vcov()} does not yet mirror, so the two can differ in those cases.
#'   Cluster- and bootstrap-aware \code{vcov} is a planned addition.
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
  ci    <- object$cluster_indices   # NULL unless clustervars was set; aligned with the per-unit EIFs

  # Cross-product of two per-unit EIF vectors: cluster-robust (sum within clusters first, with the
  # (G/(G-1)) finite-cluster factor) when ci is non-NULL, IID otherwise. This is the SAME estimator as
  # safe_inference_edid(), so diag(vcov(.)) reproduces the SE^2 that summary() reports -- including when
  # the clustering level is coarser than the unit id and bstrap = FALSE.
  .cross <- function(ej, ek) {
    if (is.null(ci)) return(sum(ej * ek) / n^2)
    G <- length(unique(ci))
    if (G <= 1L) return(NA_real_)
    (G / (G - 1)) * sum(drop(rowsum(ej, ci)) * drop(rowsum(ek, ci))) / n^2
  }

  if (which == "overall") {
    eif_v <- object$overall$eif_agg
    if (is.null(eif_v)) return(matrix(NA_real_, 1L, 1L))
    return(matrix(.cross(eif_v, eif_v), nrow = 1L, ncol = 1L,
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
        vcv[j, k] <- .cross(eif_j, eif_k)
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
        vcv[j, k] <- .cross(eif_j, eif_k)
      }
    }
    return(vcv)
  }

  # att_gt: cluster-robust crossproduct of the cell-level EIF matrix (n x K) when available.
  if (!is.null(object$eif)) {
    eif_mat <- object$eif
    df  <- object$att_gt
    nms <- paste0("ATT(", df$group, ",", df$time, ")")
    if (is.null(ci)) {
      vcv <- crossprod(eif_mat) / n^2
    } else {
      G  <- length(unique(ci))
      vcv <- if (G <= 1L) {
        matrix(NA_real_, ncol(eif_mat), ncol(eif_mat))
      } else {
        CS <- rowsum(eif_mat, ci)                 # G x K cluster sums
        (G / (G - 1)) * crossprod(CS) / n^2
      }
    }
    dimnames(vcv) <- list(nms, nms)
    return(vcv)
  }
  # Fallback: diagonal from stored SEs (which already carry the cluster-robust adjustment).
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
#' @param row.names ignored; included for S3 generic consistency
#' @param optional ignored; included for S3 generic consistency
#' @param ... not used; reserved for \code{which} (see Details)
#' @param which character: one of \code{"att_gt"}, \code{"overall"},
#'   \code{"event_study"}, \code{"group"}
#'
#' @return data.frame
#' @export
as.data.frame.edid_fit <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  which = c("att_gt", "overall", "event_study", "group")
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
