# edid-sargan.R
# Incremental Sargan moment-selection procedure for edid fits
# (Section 5.1 of Chen, Sant'Anna & Xie 2025, building on Chen & Santos 2018).

# Holm-Bonferroni step-down: reject p_(l) if p_(l) < alpha / (L + 1 - l),
# stopping at the first non-rejection (Holm 1979). Returns the per-input
# threshold (aligned to the original order) and the rejection indicator.
.edid_holm <- function(p, alpha) {
  L   <- length(p)
  ord <- order(p)
  thr_sorted <- alpha / (L + 1 - seq_len(L))
  rejected   <- logical(L)
  for (i in seq_len(L)) {
    if (p[ord[i]] < thr_sorted[i]) rejected[ord[i]] <- TRUE else break
  }
  threshold <- numeric(L)
  threshold[ord] <- thr_sorted
  list(threshold = threshold, rejected = rejected)
}

# Refit edid() with a restricted moment set: same data and options as the
# original call, but pt_assumption = "all" with the supplied `moment_set`, no
# aggregation at fit time (aggte_edid is called afterwards), pointwise/no
# bands, no bootstrap. The influence-function convention of the refits is set
# by `inference`:
#   "match_fit"   (default): copy the fitted object's EFFECTIVE misspec_robust /
#                 estimation_effect / higher_order flags and its bs_df into every
#                 refit, so xi = psi_base - psi_aug is a difference of influence
#                 functions in the SAME convention as the fit. (The edid_* session
#                 options -- e.g. edid_omega_method -- apply to the refits as they
#                 stand at call time, so keep them as at fit time.)
#   "plugin_fast": the cheap plug-in configuration (all three channels off). The
#                 incremental statistic is a quadratic form in Var-hat(xi) built
#                 from EIF differences and does not require efficiency, so under
#                 correct specification the chi-square null distribution is the
#                 same asymptotically -- but the omitted weight-estimation and
#                 first-step channels DO move Var-hat(xi) in finite samples, so
#                 p-values can differ from the fit's own inference convention.
.edid_refit_moment_set <- function(fit, data, moment_set, envir = parent.frame(),
                                   inference = "match_fit") {
  # Estimation arguments come from the fit's stored snapshot ($args), never from
  # re-evaluating the call in the caller's environment (see .edid_refit_args):
  # a caller variable mutated after fitting (e.g. a reassigned xformla) must not
  # silently change the refit, and wrapper-built calls carry `..N` promises that
  # cannot be re-evaluated at all. `envir` is used only by the legacy fallback
  # for fits that predate the snapshot.
  args <- .edid_refit_args(fit, envir)
  args[["cband_method"]] <- NULL          # let the analytic default apply (bstrap is off below)
  args$data              <- data
  args$pt_assumption     <- "all"
  args$moment_set        <- moment_set
  args$aggregate         <- "none"
  args$cband             <- FALSE
  args$bstrap            <- FALSE
  if (identical(inference, "match_fit")) {
    # The fit stores the EFFECTIVE flags (after edid()'s master-switch / applicability downgrades), so the
    # refits reproduce the convention that actually produced the fit's influence functions.
    args$misspec_robust    <- isTRUE(fit$misspec_robust)
    args$estimation_effect <- isTRUE(fit$estimation_effect)
    args$higher_order      <- isTRUE(fit$higher_order)
    if (!is.null(fit$bs_df)) args$bs_df <- fit$bs_df
  } else {
    args$misspec_robust    <- FALSE
    args$estimation_effect <- FALSE
    args$higher_order      <- FALSE
  }
  do.call(edid, args)
}

#' Incremental Sargan moment-selection procedure for edid
#'
#' Implements the incremental Sargan procedure of Section 5.1 in Chen,
#' Sant'Anna & Xie (2025), building on Chen & Santos (2018). Let
#' \eqn{\mathcal{M}} be the just-identified PT-Post base moment set: for each
#' target cohort \eqn{g}, the single restriction \eqn{(g' = g,
#' t_{pre} = g - 1)} (the most recent pre-treatment period under
#' \code{anticipation}). For each candidate pair \eqn{(g', t_{pre})} that
#' supplies an additional PT-All restriction beyond the base, the augmented set
#' \eqn{\mathcal{M}_{g', t_{pre}}} extends \eqn{\mathcal{M}} by that single
#' restriction wherever it is a valid pair for a target cohort. For each
#' candidate, a Hausman-type statistic (the eqn (5.3) form, with the positive
#' semi-definite influence-function-difference covariance) compares the
#' post-treatment event-study vector under \eqn{\mathcal{M}_{g', t_{pre}}}
#' against the base, with degrees of freedom equal to the rank of the
#' IF-difference covariance (generically, the number of event-study
#' coefficients the added restriction moves). The resulting p-values are then
#' screened by the Holm-Bonferroni step-down procedure at familywise level
#' \code{alpha}: ordering \eqn{p_{(1)} \le \cdots \le p_{(L)}}, reject
#' \eqn{p_{(\ell)}} if \eqn{p_{(\ell)} < \alpha / (L + 1 - \ell)}, stopping at
#' the first non-rejection. Rejected candidates are moment restrictions the
#' data reject; the non-rejected candidates form the admissible extension of
#' the base set.
#'
#' @param fit_restricted An \code{edid_fit}, normally from
#'   \code{edid(..., pt_assumption = "all")}. The fit supplies the design
#'   (cohorts, periods, anticipation) and the estimation options for the
#'   internal refits; the candidate pairs are always enumerated under PT-All.
#' @param data The panel data used to estimate \code{fit_restricted}, or
#'   \code{NULL} (default), in which case the data expression stored in the
#'   fit's call is re-evaluated in the caller's environment (the
#'   \code{update()} idiom). Supply \code{data} explicitly when the original
#'   object is no longer reachable by that name. All other estimation
#'   arguments are taken from the fit's stored argument snapshot
#'   (\code{fit_restricted$args}), never re-evaluated from the call, and the
#'   refits verify that \code{data} reproduces the fitted sample (same
#'   \code{n} and unit ids).
#' @param alpha Familywise error rate for the Holm-Bonferroni step-down.
#'   Default \code{0.05}.
#' @param e_set Numeric vector of post-treatment event times over which the
#'   event-study comparison is computed, or \code{NULL} (default: all finite
#'   post-treatment event times of the base fit).
#' @param inference Influence-function convention for the internal refits.
#'   \code{"match_fit"} (default) copies the fitted object's \emph{effective}
#'   \code{misspec_robust}, \code{estimation_effect}, and \code{higher_order}
#'   flags and its \code{bs_df} into every refit (base and augmented), so the
#'   test statistic is built from influence functions in the same convention
#'   as \code{fit_restricted}; the refits also run under the session's
#'   \code{edid_*} options (\code{edid_omega_method}, ...), so keep those as at
#'   fit time. \code{"plugin_fast"} uses the cheap plug-in configuration
#'   (all three channels off): substantially faster when the fit carries the
#'   \code{misspec_robust} channels (each refit then skips the weight-estimation
#'   and first-step corrections), and asymptotically valid under correct
#'   specification -- the statistic's chi-square null distribution is the same
#'   in the limit -- but its finite-sample \eqn{\widehat{Var}(\xi)} omits the
#'   weight-estimation and first-step contributions, so p-values can differ
#'   from \code{"match_fit"}.
#'
#' @details
#' Each candidate requires one refit of \code{edid()} with the internal
#' \code{moment_set} restriction (so \eqn{L + 1} fits in total), with no bands
#' and no bootstrap; the influence-function convention of the refits follows
#' \code{inference} (see above). The
#' statistic is a quadratic form in the variance of the influence-function
#' difference, which is valid without efficiency of either estimator. Both the
#' base and augmented estimators are rebuilt from the same per-pair objects,
#' so \eqn{\xi = \psi_{\mathcal{M}} - \psi_{\mathcal{M}_{g',t_{pre}}}} is a
#' clean influence-function difference. The IF-difference covariance is
#' typically rank-deficient here (adding one restriction moves only the event
#' times the affected cohorts feed), so the statistic is the pseudoinverse
#' quadratic form with \eqn{df = \mathrm{rank}}; see the Andrews (1987) caveat
#' in \code{\link{edid_hausman}}.
#'
#' @return An object of class \code{edid_sargan}: a list with elements
#'   \code{table} (one row per candidate: \code{gp}, \code{tpre},
#'   \code{H_statistic}, \code{df}, \code{p_value}, \code{holm_threshold},
#'   \code{rejected}), \code{base} (the PT-Post base moment set as a
#'   \code{(g, gp, tpre)} data.frame), \code{admissible} (candidates not
#'   rejected), \code{alpha}, \code{L}, \code{e_set}, \code{n},
#'   \code{inference} (the refit convention used). Returns
#'   \code{NULL} (with a message) when the model is just-identified (no
#'   candidate restrictions).
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025). Efficient
#'   Difference-in-Differences and Event Study Estimators. Section 5.1. \cr
#'   Chen, X., & Santos, A. (2018). Overidentification in Regular Models.
#'   \emph{Econometrica}, 86(5), 1771-1817. \cr
#'   Holm, S. (1979). A Simple Sequentially Rejective Multiple Test Procedure.
#'   \emph{Scandinavian Journal of Statistics}, 6(2), 65-70.
#'
#' @seealso \code{\link{edid}} (the \code{moment_set} argument),
#'   \code{\link{edid_hausman}}
#'
#' @examples
#' \donttest{
#' df <- data.frame(
#'   id   = rep(1:150, each = 6),
#'   time = rep(1:6, 150),
#'   g    = rep(sample(c(3, 5, Inf), 150, replace = TRUE), each = 6)
#' )
#' df$y <- rnorm(150)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
#'   rnorm(nrow(df), 0, 0.5)
#' fit <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
#'             aggregate = "event_study", cband = FALSE)
#' edid_sargan(fit, data = df)
#' }
#'
#' @export
edid_sargan <- function(fit_restricted, data = NULL, alpha = 0.05, e_set = NULL,
                        inference = c("match_fit", "plugin_fast")) {
  if (!inherits(fit_restricted, "edid_fit")) {
    stop("`fit_restricted` must be an `edid_fit` object returned by edid().", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a numeric scalar in (0, 1).", call. = FALSE)
  }
  inference <- match.arg(inference)
  if (is.null(data)) {
    data <- tryCatch(as.data.frame(eval(fit_restricted$call$data, envir = parent.frame())),
                     error = function(e) NULL)
    if (is.null(data) || !is.data.frame(data) || nrow(data) == 0L) {
      stop("Could not recover the estimation data from the fit's call; pass `data` explicitly.",
           call. = FALSE)
    }
  }

  fit <- fit_restricted
  tg  <- sort(fit$treatment_groups[is.finite(fit$treatment_groups) & fit$treatment_groups != 0])
  tp  <- sort(fit$time_periods)
  p1  <- min(tp)
  ant <- fit$anticipation %||% 0L

  # Full PT-All pair enumeration per target cohort (the same enumeration edid()
  # itself uses), the PT-Post base pair per cohort, and the candidate
  # additional restrictions in the paper's order (by g', then t_pre).
  full_pairs <- stats::setNames(lapply(tg, function(g) {
    enumerate_valid_pairs_edid(target_g = g, treatment_groups = tg, time_periods = tp,
                               period_1 = p1, pt_assumption = "all", anticipation = ant)
  }), as.character(tg))

  base_rows <- vector("list", length(tg))
  for (k in seq_along(tg)) {
    g  <- tg[k]
    pg <- full_pairs[[as.character(g)]]
    self <- pg[is.finite(pg$gp) & pg$gp == g, , drop = FALSE]
    if (nrow(self) == 0L) next                      # cohort with no pre period: no base moment
    base_rows[[k]] <- data.frame(g = g, gp = g, tpre = max(self$tpre))
  }
  base_ms <- do.call(rbind, base_rows)
  if (is.null(base_ms) || nrow(base_ms) == 0L) {
    stop("No PT-Post base moments are available (no cohort has a usable pre-treatment period).",
         call. = FALSE)
  }

  extra_rows <- vector("list", length(tg))
  for (k in seq_along(tg)) {
    g   <- tg[k]
    pg  <- full_pairs[[as.character(g)]]
    if (nrow(pg) == 0L) next
    bg  <- base_ms[base_ms$g == g, , drop = FALSE]
    if (nrow(bg) == 1L) {
      is_base <- pg$gp == bg$gp & pg$tpre == bg$tpre
      pg <- pg[!is_base, , drop = FALSE]
    }
    if (nrow(pg) > 0L) extra_rows[[k]] <- pg
  }
  extra <- unique(do.call(rbind, extra_rows))
  if (is.null(extra) || nrow(extra) == 0L) {
    message("No additional moment restrictions to test (the model is just-identified).")
    return(invisible(NULL))
  }
  extra <- extra[order(extra$gp, extra$tpre), , drop = FALSE]
  rownames(extra) <- NULL
  L <- nrow(extra)

  # Base estimator M: every cell uses only its PT-Post base pair. All test fits
  # (base and augmented) are rebuilt in the same cheap configuration so the
  # influence-function differences are clean.
  caller_env <- parent.frame()
  fit_base <- .edid_refit_moment_set(fit, data, base_ms, envir = caller_env, inference = inference)
  if (!identical(fit_base$n, fit$n) || !identical(fit_base$all_units, fit$all_units)) {
    stop("The data used for the refits does not match the fitted sample (n or unit ids ",
         "differ from `fit_restricted`); pass the original estimation data via `data`.",
         call. = FALSE)
  }
  pB <- .edid_param_ifs(fit_base, "event_study", e_set)
  e_set <- pB$e
  n  <- fit_base$n
  ci <- fit_base$cluster_indices
  # Absolute variance scale of the base estimator's coordinates, for the
  # degenerate-contrast guard in .edid_if_diff_quadform (see edid-hausman.R).
  vB <- diag(as.matrix(n * cluster_cov_edid(pB$IF, ci, n)))

  results <- data.frame(
    gp = extra$gp, tpre = extra$tpre,
    H_statistic = numeric(L), df = integer(L), p_value = numeric(L),
    holm_threshold = numeric(L), rejected = logical(L)
  )

  for (l in seq_len(L)) {
    gp_l <- extra$gp[l]; tp_l <- extra$tpre[l]
    # M_{g',tpre}: base + the candidate restriction wherever it is a valid
    # (non-base) pair for the target cohort.
    add_rows <- vector("list", length(tg))
    for (k in seq_along(tg)) {
      g  <- tg[k]
      pg <- full_pairs[[as.character(g)]]
      if (any(pg$gp == gp_l & pg$tpre == tp_l)) {
        add_rows[[k]] <- data.frame(g = g, gp = gp_l, tpre = tp_l)
      }
    }
    ms_l <- unique(rbind(base_ms, do.call(rbind, add_rows)))
    fit_aug <- .edid_refit_moment_set(fit, data, ms_l, envir = caller_env, inference = inference)
    pA <- .edid_param_ifs(fit_aug, "event_study", e_set)

    d  <- pA$est - pB$est
    xi <- pB$IF - pA$IF
    vA <- diag(as.matrix(n * cluster_cov_edid(pA$IF, ci, n)))
    qf <- .edid_if_diff_quadform(d, xi, n, ci, v_scale = max(vB, vA))
    results$H_statistic[l] <- qf$statistic
    results$df[l]          <- qf$df
    results$p_value[l]     <- qf$p_value
  }

  holm <- .edid_holm(results$p_value, alpha)
  results$holm_threshold <- holm$threshold
  results$rejected       <- holm$rejected

  out <- list(
    table      = results,
    base       = base_ms,
    admissible = results[!results$rejected, c("gp", "tpre"), drop = FALSE],
    alpha      = alpha,
    L          = L,
    e_set      = e_set,
    n          = n,
    clustered  = !is.null(ci),
    inference  = inference
  )
  class(out) <- c("edid_sargan", "list")
  out
}

#' @describeIn edid_sargan Print method.
#' @param x an \code{edid_sargan} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_sargan <- function(x, digits = 4, ...) {
  cat("\nIncremental Sargan procedure (Chen, Sant'Anna & Xie 2025, Section 5.1)\n")
  cat(sprintf("  Base set: PT-Post pairs (g' = g, t_pre = g - 1); %d candidate restriction(s)%s\n",
              x$L, if (isTRUE(x$clustered)) "; cluster-robust" else ""))
  cat(sprintf("  Event-study comparison over E = {%s}; Holm-Bonferroni FWER alpha = %s\n",
              paste(x$e_set, collapse = ", "), format(x$alpha)))
  if (identical(x$inference, "plugin_fast")) {
    cat("  Refit convention: plug-in influence functions only; excludes weight-estimation\n")
    cat("  and first-step corrections (inference = 'plugin_fast').\n")
  }
  cat("\n")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1L))
  tab[num] <- lapply(tab[num], function(z) signif(z, digits))
  print(tab, row.names = FALSE)
  n_rej <- sum(x$table$rejected)
  cat(sprintf("\n%d of %d candidate restriction(s) rejected; %d admissible beyond the base set.\n",
              n_rej, x$L, x$L - n_rej))
  invisible(x)
}
