# edid-frontier.R
# Reported-parameter robustness frontier for edid fits
# (Section 5.3, Theorem 5.2 of Chen, Sant'Anna & Xie 2025).

#' Robustness frontier for reported event-study contrasts
#'
#' Implements the reported-parameter robustness frontier of Theorem 5.2 in
#' Chen, Sant'Anna & Xie (2025). For a scalar event-study summary
#' \eqn{\theta} (a single \eqn{ES(e)} or the average \eqn{ES_{avg}}), let
#' \eqn{\widehat\theta_R} be the efficient (PT-All) estimate,
#' \eqn{\widehat\theta_U} the conservative (PT-Post) estimate, and
#' \eqn{\xi = \psi_U - \psi_R} the estimator-difference influence function.
#' The scalar Hausman diagnostic of eqn (5.5) is
#' \deqn{H_{\theta,n} = n(\widehat\theta_U - \widehat\theta_R)^2 / \widehat{D},
#'   \qquad \widehat{D} = \widehat{E}[\xi^2],}
#' which is asymptotically \eqn{\chi^2_1} under PT-All. For a tolerance
#' \eqn{\tau > 0} on the acceptable variance inflation, the set of estimates in
#' the affine class \eqn{\widehat\theta_R + \lambda(\widehat\theta_U -
#' \widehat\theta_R)} whose first-order variance does not exceed
#' \eqn{(1 + \tau^2) V_R / n} is the frontier interval of eqn (5.6):
#' \deqn{\widehat\theta_R \pm \tau \sqrt{H_{\theta,n}}\;
#'   \widehat{se}(\widehat\theta_R).}
#' The frontier quantifies how far the reported estimate can move if the
#' researcher relaxes the stronger PT-All restrictions while paying a
#' transparent precision cost; it does not bound movement under arbitrary
#' violations of parallel trends (for that, see honest-confidence-interval
#' approaches, which are complementary).
#'
#' @inheritParams edid_hausman
#' @param parameter Which scalar summaries to report: any of
#'   \code{"event_study"} (each post-treatment \eqn{ES(e)}) and
#'   \code{"overall"} (\eqn{ES_{avg}}). Default: both.
#' @param tau Numeric vector of positive tolerance parameters. Default
#'   \code{c(0.25, 0.5, 1)}, the grid recommended in Section 5.3 of the paper.
#'
#' @details
#' \eqn{\widehat{D}} is estimated from the per-unit influence-function
#' difference (cluster-robust when the fits carry cluster assignments), which
#' makes it nonnegative in finite samples. Theorem 5.2 requires \eqn{D > 0}:
#' when the two estimators coincide for a coordinate (a just-identified
#' contrast, \eqn{\xi \approx 0}), the statistic is 0/0 and the frontier
#' degenerates to the point estimate; the implementation guards this case and
#' reports \eqn{H = 0} (so the frontier radius is exactly zero) instead of
#' \code{NaN}. P-values use \code{pchisq(lower.tail = FALSE)} to avoid
#' underflow for large statistics.
#'
#' @return An object of class \code{edid_frontier} whose \code{$table} is a
#'   data.frame with one row per (parameter, tau): \code{parameter}, \code{e},
#'   \code{theta_R}, \code{theta_U}, \code{se_R}, \code{H} (the eqn (5.5)
#'   statistic), \code{sqrt_H}, \code{p_value}, \code{tau}, \code{radius}
#'   (\eqn{= \tau \sqrt{H}\, \widehat{se}_R}), \code{frontier_low},
#'   \code{frontier_high}, and the efficient confidence limits \code{ci_low},
#'   \code{ci_high} (pointwise, at the restricted fit's \code{alp}).
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025). Efficient
#'   Difference-in-Differences and Event Study Estimators. Section 5.3,
#'   Theorem 5.2 and Remark 5.3. \cr
#'   Andrews, I., Chen, J., & Tecchio, J. (2025). Overidentification and
#'   Misspecification-Robust Inference. arXiv:2508.13076. \cr
#'   Hausman, J. A. (1978). Specification Tests in Econometrics.
#'   \emph{Econometrica}, 46(6), 1251-1271.
#'
#' @seealso \code{\link{edid}}, \code{\link{edid_hausman}},
#'   \code{\link{edid_adaptive}}
#'
#' @examples
#' \donttest{
#' df <- data.frame(
#'   id   = rep(1:120, each = 6),
#'   time = rep(1:6, 120),
#'   g    = rep(sample(c(3, 5, Inf), 120, replace = TRUE), each = 6)
#' )
#' df$y <- rnorm(120)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
#'   rnorm(nrow(df), 0, 0.5)
#' fit_R <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
#'               aggregate = "event_study", cband = FALSE)
#' fit_U <- edid(df, "y", "id", "time", "g", pt_assumption = "post",
#'               aggregate = "event_study", cband = FALSE)
#' edid_frontier(fit_U, fit_R)
#' }
#'
#' @export
edid_frontier <- function(fit_unrestricted, fit_restricted,
                          parameter = c("event_study", "overall"),
                          tau = c(0.25, 0.5, 1),
                          e_set = NULL, data = NULL) {
  parameter <- match.arg(parameter, several.ok = TRUE)
  if (!is.numeric(tau) || length(tau) == 0L || any(!is.finite(tau)) || any(tau <= 0)) {
    stop("`tau` must be a vector of positive finite tolerances.", call. = FALSE)
  }
  tau <- sort(unique(as.numeric(tau)))
  .edid_toolkit_check_fits(fit_unrestricted, fit_restricted)

  # Refit both legs in the EFFICIENT plug-in configuration: the robustness frontier
  # theta_R +/- tau * sqrt(H) * se(theta_R) is the Andrews-Chen-Tecchio (2025) range result, stated for the
  # efficient inverse-variance variance (Prop 5.2: se(theta_R) is the EFFICIENT estimator's SE), not a
  # misspecification-robust SE. Point estimates are unchanged; `data` is recovered from the call when NULL.
  data <- .edid_recover_data(fit_restricted, data, parent.frame())
  fit_unrestricted <- .edid_plugin_refit(fit_unrestricted, data, parent.frame())
  fit_restricted   <- .edid_plugin_refit(fit_restricted,   data, parent.frame())

  n     <- fit_restricted$n
  ci    <- fit_restricted$cluster_indices
  n_eff <- .edid_overid_n_eff(fit_restricted)  # passed to the scalar Hausman for signature parity (no-op in 1-D)
  alpha <- fit_restricted$alpha %||% 0.05
  z     <- stats::qnorm(1 - alpha / 2)

  # Full-J (worst-case scope) radius, reported ALONGSIDE the directed Hausman radius (F1). The directed
  # radius theta_R +/- tau*sqrt(H_theta)*se is the actual movement from relaxing PT-All to PT-Post; the
  # full-J radius theta_R +/- tau*sqrt(J)*se is the Andrews-Chen-Tecchio worst-case scope over ALL
  # admissible reweightings (J = the omnibus over-identification statistic for the cells feeding the
  # reported scalar). A `fragile` flag marks a width driven by over-id rank-deficiency (Q -> n_eff)
  # rather than genuine over-identifying signal -- the do-what-is-right caveat: do not read a wide
  # full-J frontier as real sensitivity when the joint contrast is poorly identified.
  ov_par <- c(if ("overall" %in% parameter) "overall", if ("event_study" %in% parameter) "event_study")
  ovJ <- tryCatch(suppressWarnings(suppressMessages(
           edid_overid(fit_restricted, data = data, parameter = ov_par, e_set = e_set))),
           error = function(e) NULL)
  .lookupJ <- function(label, e) {
    if (is.null(ovJ) || is.null(ovJ$table)) return(list(J = NA_real_, df = NA_integer_, Qn = NA_real_))
    tb  <- ovJ$table
    row <- if (identical(label, "ES_avg")) tb[tb$parameter == "overall", , drop = FALSE]
           else tb[tb$parameter == "event_study" & tb$e == e, , drop = FALSE]
    if (nrow(row) != 1L) return(list(J = NA_real_, df = NA_integer_, Qn = NA_real_))
    list(J = row$J_statistic, df = row$df, Qn = row$nominal_df)
  }

  # Assemble the scalar coordinates: each ES(e) over e_set, then ES_avg.
  coords <- list()
  if ("event_study" %in% parameter) {
    e_set <- .edid_shared_e_set(fit_unrestricted, fit_restricted, e_set)
    pU <- .edid_param_ifs(fit_unrestricted, "event_study", e_set)
    pR <- .edid_param_ifs(fit_restricted,  "event_study", e_set)
    for (j in seq_along(e_set)) {
      coords[[length(coords) + 1L]] <- list(
        label = sprintf("ES(%g)", e_set[j]), e = e_set[j],
        theta_U = pU$est[j], theta_R = pR$est[j],
        psi_U = pU$IF[, j], psi_R = pR$IF[, j])
    }
  }
  if ("overall" %in% parameter) {
    oU <- .edid_param_ifs(fit_unrestricted, "overall")
    oR <- .edid_param_ifs(fit_restricted,  "overall")
    coords[[length(coords) + 1L]] <- list(
      label = "ES_avg", e = NA_real_,
      theta_U = oU$est, theta_R = oR$est,
      psi_U = oU$IF[, 1L], psi_R = oR$IF[, 1L])
  }

  rows <- vector("list", length(coords) * length(tau))
  r <- 0L
  for (co in coords) {
    V_R  <- as.numeric(n * cluster_cov_edid(matrix(co$psi_R, ncol = 1L), ci, n))
    se_R <- sqrt(V_R / n)
    d    <- co$theta_U - co$theta_R
    sc   <- .edid_scalar_hausman(d, co$psi_U - co$psi_R, n, ci, v_scale = V_R, n_eff = n_eff)
    fj   <- .lookupJ(co$label, co$e)
    # Fragility keyed to the EFFECTIVE over-id rank (df_full), not the nominal Q-p: a design with small
    # genuine rank but large nominal redundancy is benign (the low-rank thesis), so nominal over-fires.
    # A rank-deficient scope (uncomputable joint J: NA statistic with a finite, here zero, df) is fragile.
    fragile <- (is.na(fj$J) && is.finite(fj$df)) ||
      (is.finite(fj$df) && is.finite(n_eff) && n_eff > 0 && fj$df > 0.25 * n_eff)
    for (tt in tau) {
      radius   <- tt * sqrt(sc$H) * se_R                                  # directed Hausman radius
      radiusJ  <- if (is.finite(fj$J)) tt * sqrt(fj$J) * se_R else NA_real_  # full-J worst-case radius
      r <- r + 1L
      rows[[r]] <- data.frame(
        parameter = co$label, e = co$e,
        theta_R = co$theta_R, theta_U = co$theta_U,
        se_R = se_R, H = sc$H, sqrt_H = sqrt(sc$H), p_value = sc$p_value,
        tau = tt, radius = radius,
        frontier_low  = co$theta_R - radius,
        frontier_high = co$theta_R + radius,
        # Full-J (worst-case scope) frontier reported alongside (F1):
        J_full = fj$J, df_full = fj$df, radius_fullJ = radiusJ,
        frontier_fullJ_low  = if (is.finite(radiusJ)) co$theta_R - radiusJ else NA_real_,
        frontier_fullJ_high = if (is.finite(radiusJ)) co$theta_R + radiusJ else NA_real_,
        fragile = fragile,
        ci_low  = co$theta_R - z * se_R,
        ci_high = co$theta_R + z * se_R,
        stringsAsFactors = FALSE)
    }
  }
  table <- do.call(rbind, rows)
  rownames(table) <- NULL

  out <- list(table = table, tau = tau, e_set = e_set, n = n,
              clustered = !is.null(ci), alpha = alpha)
  class(out) <- c("edid_frontier", "list")
  out
}

#' @describeIn edid_frontier Print method.
#' @param x an \code{edid_frontier} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_frontier <- function(x, digits = 4, ...) {
  cat("\nRobustness frontier for reported event-study contrasts\n")
  cat("(Chen, Sant'Anna & Xie 2025, Theorem 5.2, eqn 5.6)\n")
  cat(sprintf("  tau grid: {%s}%s\n", paste(x$tau, collapse = ", "),
              if (isTRUE(x$clustered)) "; cluster-robust" else ""))
  cat("  directed  : radius      = tau * sqrt(H)      * se(theta_R)  [sensitivity to relaxing PT-All]\n")
  cat("  worst-case: radius_fullJ = tau * sqrt(J_full) * se(theta_R)  [scope over all admissible reweightings]\n")
  if (!is.null(x$table$fragile) && any(isTRUE(x$table$fragile)))
    cat("  NOTE: `fragile = TRUE` rows -- the full-J width is driven by over-id rank-deficiency (Q ~ n_eff),\n        not genuine signal; read the per-cell edid_overid()$cells / edid_sargan instead.\n")
  cat("\n")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1L))
  tab[num] <- lapply(tab[num], function(z) signif(z, digits))
  print(tab, row.names = FALSE)
  invisible(x)
}
