# edid-aggte.R
# Aggregation function for edid_fit objects, mirroring the aggte() interface.

#' Aggregate edid_fit estimates
#'
#' Provides the same user-facing interface as \code{\link[did]{aggte}} but
#' accepts an \code{edid_fit} object produced by \code{\link{edid}}.
#'
#' @param edid_fit_obj An \code{edid_fit} object returned by \code{edid()}.
#' @param type Character scalar: aggregation type, mirroring \code{did::aggte()}. One of
#'   \code{"simple"} (cohort-share-weighted average over all post-treatment cells),
#'   \code{"dynamic"} (event-study: average of \eqn{ES(e)} over \eqn{e \ge 0}),
#'   \code{"group"} (cohort-level overall ATTs), or \code{"calendar"} (equal-weight average over
#'   post-treatment calendar periods).
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

  # Helper: extract the "simple" overall ATT + SE from the edid_fit object. The headline
  # edid()$overall is the dynamic (event-study) average, so the simple cohort-share aggregate is
  # stored in $overall_simple; fall back to $overall when only the simple aggregate was computed.
  .get_overall <- function(obj) {
    ov  <- if (!is.null(obj$overall_simple)) obj$overall_simple else obj$overall
    att <- if (!is.null(ov)) ov$att else NA_real_
    se  <- if (!is.null(ov)) ov$se  else NA_real_
    list(att = att, se = se)
  }

  # ------------------------------------------------------------------
  # calendar: per-calendar-time effects + equal-weight overall (e>=0 post periods)
  # ------------------------------------------------------------------
  if (type == "calendar") {
    cal_list <- edid_fit_obj$calendar
    if (is.null(cal_list) || length(cal_list) == 0L) {
      stop("No calendar results in edid_fit_obj. ",
           "Re-run edid() with aggregate = \"calendar\" or \"all\".")
    }
    orig_idx <- seq_along(cal_list)
    t_vals   <- vapply(cal_list, function(x) x$time, numeric(1L))
    att_vec  <- vapply(cal_list, function(x) x$att,  numeric(1L))
    se_vec   <- vapply(cal_list, function(x) x$se,   numeric(1L))
    if (na.rm) {
      keep <- !is.na(att_vec)
      orig_idx <- orig_idx[keep]; t_vals <- t_vals[keep]
      att_vec  <- att_vec[keep];  se_vec <- se_vec[keep]
    }
    # Calendar overall: equal-weight average over post calendar periods (matching
    # did::aggte(type = "calendar")), with cluster-robust SE from the averaged per-period EIFs.
    fin <- is.finite(att_vec)
    cal_overall_att <- NA_real_; cal_overall_se <- NA_real_
    if (any(fin)) {
      cal_overall_att <- mean(att_vec[fin])
      eifs <- lapply(orig_idx[fin], function(ii) cal_list[[ii]]$eif_agg)
      if (!any(vapply(eifs, is.null, logical(1L)))) {
        agg_eif <- Reduce(`+`, eifs) / length(eifs)
        cal_overall_se <- safe_inference_edid(agg_eif, edid_fit_obj$cluster_indices,
                                              alp, cal_overall_att)$se
      }
    }
    out <- list(
      att.egt     = att_vec,
      se.egt      = se_vec,
      egt         = t_vals,
      type        = type,
      overall.att = cal_overall_att,
      overall.se  = cal_overall_se,
      alp         = alp,
      call        = mc
    )
    class(out) <- c("AGGTEobj_edid", "list")
    return(out)
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

    # Convert to flat vectors, tracking original es_list indices
    orig_idx <- seq_along(es_list)
    e_vals   <- vapply(es_list, function(x) x$e,   numeric(1L))
    att_vec  <- vapply(es_list, function(x) x$att, numeric(1L))
    se_vec   <- vapply(es_list, function(x) x$se,  numeric(1L))

    # Apply balance_e filter first (symmetric window)
    if (!is.null(balance_e)) {
      keep_bal <- (e_vals >= -balance_e) & (e_vals <= balance_e)
      orig_idx <- orig_idx[keep_bal]
      e_vals   <- e_vals[keep_bal]
      att_vec  <- att_vec[keep_bal]
      se_vec   <- se_vec[keep_bal]
    }

    # Apply min_e / max_e filter
    keep_range <- (e_vals >= min_e) & (e_vals <= max_e)
    orig_idx <- orig_idx[keep_range]
    e_vals   <- e_vals[keep_range]
    att_vec  <- att_vec[keep_range]
    se_vec   <- se_vec[keep_range]

    # Optionally drop NAs
    if (na.rm) {
      keep_na <- !is.na(att_vec)
      orig_idx <- orig_idx[keep_na]
      e_vals   <- e_vals[keep_na]
      att_vec  <- att_vec[keep_na]
      se_vec   <- se_vec[keep_na]
    }

    # Dynamic overall: equal-weight average of post-treatment ES(e),
    # matching did::aggte(type = "dynamic").
    post_mask <- e_vals >= 0 & is.finite(att_vec)
    dyn_overall_att <- NA_real_
    dyn_overall_se  <- NA_real_

    if (any(post_mask)) {
      dyn_overall_att <- mean(att_vec[post_mask])

      # Use orig_idx to access the correct es_list elements
      post_orig_idx <- orig_idx[post_mask]
      post_eifs <- lapply(post_orig_idx, function(ii) es_list[[ii]]$eif_agg)
      if (!any(vapply(post_eifs, is.null, logical(1L)))) {
        n_units <- edid_fit_obj$n
        agg_eif <- Reduce(`+`, post_eifs) / length(post_eifs)
        # cluster-robust (matches the internal aggregate_*_edid path); reduces to the
        # i.i.d. sqrt(sum(eif^2)/n^2) when cluster_indices is NULL.
        dyn_overall_se <- safe_inference_edid(agg_eif, edid_fit_obj$cluster_indices,
                                              edid_fit_obj$alpha, dyn_overall_att)$se
      }
    }

    out <- list(
      att.egt     = att_vec,
      se.egt      = se_vec,
      egt         = e_vals,
      type        = type,
      overall.att = dyn_overall_att,
      overall.se  = dyn_overall_se,
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

    orig_idx <- seq_along(gr_list)
    g_vals   <- vapply(gr_list, function(x) x$group, numeric(1L))
    att_vec  <- vapply(gr_list, function(x) x$att,   numeric(1L))
    se_vec   <- vapply(gr_list, function(x) x$se,    numeric(1L))

    if (na.rm) {
      keep_na  <- !is.na(att_vec)
      orig_idx <- orig_idx[keep_na]
      g_vals   <- g_vals[keep_na]
      att_vec  <- att_vec[keep_na]
      se_vec   <- se_vec[keep_na]
    }

    # Group overall: cohort-share-weighted average of per-group ATTs
    # with WIF correction for estimated weights, matching
    # did::aggte(type = "group") in compute.aggte.R lines 300-325.
    valid_mask <- is.finite(att_vec)
    grp_overall_att <- NA_real_
    grp_overall_se  <- NA_real_

    if (any(valid_mask)) {
      g_valid      <- g_vals[valid_mask]
      att_valid    <- att_vec[valid_mask]
      valid_orig   <- orig_idx[valid_mask]
      n_valid      <- length(g_valid)

      # Cohort-share weights
      cf <- edid_fit_obj$cohort_fractions
      if (!is.null(cf)) {
        pgg <- vapply(g_valid, function(gv) {
          val <- cf[[as.character(gv)]]
          if (is.null(val)) 1.0 else val
        }, numeric(1L))
      } else {
        pgg <- rep(1, n_valid)
      }
      pg_norm <- pgg / sum(pgg)

      grp_overall_att <- sum(pg_norm * att_valid)

      # SE with WIF correction (matching compute.aggte.R)
      grp_eifs <- lapply(valid_orig, function(ii) gr_list[[ii]]$eif_agg)
      has_eifs <- !any(vapply(grp_eifs, is.null, logical(1L)))
      G_vec <- edid_fit_obj$unit_cohorts

      if (has_eifs && !is.null(G_vec)) {
        n_units <- edid_fit_obj$n
        eif_mat <- do.call(cbind, grp_eifs)

        # WIF: influence function of the estimated pi_g weights
        # Following compute.aggte.R wif() at lines 604-621
        S <- sum(pgg)
        wif_mat <- sapply(seq_len(n_valid), function(k) {
          (1 * (G_vec == g_valid[k]) - pgg[k]) / S
        })
        wif_denom <- rowSums(sapply(seq_len(n_valid), function(k) {
          1 * (G_vec == g_valid[k]) - pgg[k]
        }))
        wif_mat <- wif_mat - wif_denom %*% t(pgg / S^2)

        # Overall IF = weighted sum of per-group IFs + WIF * att
        agg_eif <- drop(eif_mat %*% pg_norm) +
                   drop(wif_mat %*% att_valid)
        grp_overall_se <- safe_inference_edid(agg_eif, edid_fit_obj$cluster_indices,
                                              edid_fit_obj$alpha, grp_overall_att)$se
      } else if (has_eifs) {
        # Fallback: no WIF (if unit_cohorts not available)
        n_units <- edid_fit_obj$n
        eif_mat <- do.call(cbind, grp_eifs)
        agg_eif <- drop(eif_mat %*% pg_norm)
        grp_overall_se <- safe_inference_edid(agg_eif, edid_fit_obj$cluster_indices,
                                              edid_fit_obj$alpha, grp_overall_att)$se
      }
    }

    out <- list(
      att.egt     = att_vec,
      se.egt      = se_vec,
      egt         = g_vals,
      type        = type,
      overall.att = grp_overall_att,
      overall.se  = grp_overall_se,
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
