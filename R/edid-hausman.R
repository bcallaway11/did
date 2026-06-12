# edid-hausman.R
# Hausman-type specification test of PT-All vs PT-Post for edid fits
# (Section 5.1, Theorem 5.1 of Chen, Sant'Anna & Xie 2025), plus the shared
# internals used by the Section-5 toolkit (edid_sargan, edid_frontier,
# edid_adaptive).

# ---------------------------------------------------------------------------
# Shared toolkit internals
# ---------------------------------------------------------------------------

# Validate that two edid fits are comparable: same sample (n, unit order, cohort
# assignment = the data fingerprint available on the fit), same design (periods,
# cohorts, anticipation), and same clustering. `require_pt = TRUE` additionally
# warns unless (unrestricted, restricted) = (PT-Post, PT-All), the pairing the
# Section-5 results are stated for.
.edid_toolkit_check_fits <- function(fit_unrestricted, fit_restricted, require_pt = TRUE) {
  if (!inherits(fit_unrestricted, "edid_fit") || !inherits(fit_restricted, "edid_fit")) {
    stop("`fit_unrestricted` and `fit_restricted` must be `edid_fit` objects returned by edid().",
         call. = FALSE)
  }
  fu <- fit_unrestricted; fr <- fit_restricted
  if (!identical(fu$n, fr$n)) {
    stop("The two fits have different sample sizes (n = ", fu$n, " vs ", fr$n,
         "); they must be estimated on the same data.", call. = FALSE)
  }
  if (!is.null(fu$all_units) && !is.null(fr$all_units) && !identical(fu$all_units, fr$all_units)) {
    stop("The two fits carry different unit identifiers / unit order; the per-unit influence ",
         "functions cannot be differenced. Fit both estimators on the same data.", call. = FALSE)
  }
  if (!is.null(fu$unit_cohorts) && !is.null(fr$unit_cohorts) &&
      !isTRUE(all.equal(fu$unit_cohorts, fr$unit_cohorts))) {
    stop("The two fits assign units to different cohorts (data fingerprint mismatch); ",
         "they must be estimated on the same data.", call. = FALSE)
  }
  if (!isTRUE(all.equal(fu$time_periods, fr$time_periods)) ||
      !isTRUE(all.equal(fu$treatment_groups, fr$treatment_groups))) {
    stop("The two fits have different time periods or treatment cohorts.", call. = FALSE)
  }
  if (!identical(fu$anticipation, fr$anticipation)) {
    stop("The two fits use different `anticipation` values.", call. = FALSE)
  }
  cu <- fu$cluster_indices; cr <- fr$cluster_indices
  if (is.null(cu) != is.null(cr) || (!is.null(cu) && !identical(cu, cr))) {
    stop("The two fits use different cluster assignments (`clustervars`); the estimator-",
         "difference covariance must be computed under a common clustering.", call. = FALSE)
  }
  if (require_pt &&
      !(identical(fu$pt_assumption, "post") && identical(fr$pt_assumption, "all"))) {
    warning("Expected `fit_unrestricted` from edid(pt_assumption = \"post\") and ",
            "`fit_restricted` from edid(pt_assumption = \"all\") (the conservative vs efficient ",
            "pairing of Section 5 of Chen, Sant'Anna & Xie 2025). Got pt_assumption = \"",
            fu$pt_assumption, "\" vs \"", fr$pt_assumption, "\"; results are only meaningful if ",
            "the restricted fit imposes strictly more moment restrictions.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Leg-health guards for the Section-5 toolkit (broken-leg + few-cluster)
# ---------------------------------------------------------------------------

# Number of clusters feeding a fit's cluster-robust covariance (G; n when unclustered,
# i.e. unit-level / iid). The few-cluster guard compares this to EDID_FEWCLUSTER_MIN.
.edid_n_clusters <- function(fit) {
  ci <- fit$cluster_indices
  if (is.null(ci)) (fit$n %||% length(unique(fit$all_units))) else length(unique(ci))
}

# Is a fit a numerically degenerate leg? A non-rejection (or a point estimate) built on
# such a leg is HOLLOW -- the restricted/efficient fit carried extreme propensity ratios,
# a non-credible weight channel, or cross-cohort hedges that carry the estimand, or its
# reported SEs are non-finite. Reads the fit's own $diagnostics read-out (set by edid())
# plus the aggregation SEs the test will difference. Returns a list with a logical
# `broken` and a character vector of `reasons` (empty when healthy). Older fits without
# $diagnostics degrade gracefully to the SE-finiteness check only.
.edid_leg_health <- function(fit, label) {
  reasons <- character(0L)
  d <- fit$diagnostics
  if (!is.null(d)) {
    if (isTRUE(d$unstable)) {
      if ((d$n_extreme_ratio %||% 0L) > 0L)
        reasons <- c(reasons, sprintf("extreme propensity ratios in %d estimation step(s)", d$n_extreme_ratio))
      if ((d$n_psi_unstable %||% 0L) > 0L)
        reasons <- c(reasons, sprintf("weight-estimation channel not a credible IF in %d cell(s)", d$n_psi_unstable))
      if (isTRUE(d$net_hedge_flag))
        reasons <- c(reasons, sprintf("cross-cohort hedges carry the estimand (net hedge mass %.2f ~ gross %.2f)",
                                      d$net_hedge_mass %||% NA_real_, d$gross_hedge_mass %||% NA_real_))
    }
  }
  # Non-finite reported SEs on the headline (overall) aggregation: a hard breakage signal
  # available even without $diagnostics.
  ov <- fit$overall
  se_ov <- if (!is.null(ov)) suppressWarnings(ov$overall.se) else NULL
  if (!is.null(se_ov) && length(se_ov) && any(!is.finite(se_ov))) {
    reasons <- c(reasons, "non-finite reported standard error on the overall aggregation")
  }
  list(broken = length(reasons) > 0L, reasons = unique(reasons), label = label)
}

# Dynamic (event-study) AGGTEobj of a fit: reuse the stored one, else aggregate
# on the fly with the same machinery edid() uses (aggte_edid -> did::aggte).
.edid_dynamic_aggte <- function(fit) {
  a <- fit$event_study
  if (is.null(a)) a <- aggte_edid(fit, type = "dynamic", na.rm = TRUE)
  if (is.null(a) || !inherits(a, "AGGTEobj")) {
    stop("Could not construct the event-study aggregation for an edid fit.", call. = FALSE)
  }
  a
}

# Per-unit influence functions + estimates of the requested parameter vector.
# parameter = "event_study": the ES(e) vector over e_set (default: all finite
# post-treatment e's); parameter = "overall": the scalar ES_avg (the average of
# ES(e) over e >= 0, i.e. the dynamic AGGTEobj's overall). The IFs are the same
# aggregation IFs that vcov.edid_fit() reads (did::aggte's per-element
# influence functions, including the cohort-share weight-estimation correction).
.edid_param_ifs <- function(fit, parameter, e_set = NULL) {
  a <- .edid_dynamic_aggte(fit)
  g <- .edid_agg_if(a)
  if (identical(parameter, "overall")) {
    if (is.null(a$overall.att) || is.null(g$overall)) {
      stop("The dynamic aggregation does not expose the overall ES_avg influence function.",
           call. = FALSE)
    }
    return(list(e = NA_real_, est = a$overall.att,
                IF = matrix(as.numeric(g$overall), ncol = 1L)))
  }
  if (is.null(g$egt) || is.null(a$egt)) {
    stop("The dynamic aggregation does not expose per-e influence functions.", call. = FALSE)
  }
  e_all <- a$egt
  post  <- e_all[e_all >= 0 & is.finite(a$att.egt)]
  if (is.null(e_set)) {
    e_set <- post
  } else {
    bad <- setdiff(e_set, post)
    if (length(bad) > 0L) {
      stop("`e_set` contains event times not available as finite post-treatment ES(e) ",
           "coordinates in this fit: ", paste(bad, collapse = ", "), call. = FALSE)
    }
  }
  e_set <- sort(unique(e_set))
  if (length(e_set) == 0L) stop("No post-treatment event times available.", call. = FALSE)
  idx <- match(e_set, e_all)
  list(e = e_set, est = as.numeric(a$att.egt[idx]),
       IF = as.matrix(g$egt)[, idx, drop = FALSE])
}

# Default e_set for a two-fit comparison: the intersection of both fits'
# finite post-treatment event times.
.edid_shared_e_set <- function(fit_unrestricted, fit_restricted, e_set = NULL) {
  if (!is.null(e_set)) return(sort(unique(e_set)))
  eu <- .edid_param_ifs(fit_unrestricted, "event_study")$e
  er <- .edid_param_ifs(fit_restricted,  "event_study")$e
  shared <- intersect(eu, er)
  if (length(shared) == 0L) {
    stop("The two fits share no finite post-treatment event times.", call. = FALSE)
  }
  sort(shared)
}

# Rank-aware quadratic form for the IF-difference Hausman statistic:
#   H = n * d' D^+ d,   D = n * Var-hat(xi)  (cluster-robust when clustered),
# with df = rank(D) by eigenvalue thresholding and the Moore-Penrose
# pseudoinverse on the rank-deficient branch. The full-rank branch uses the
# exact inverse (numerically identical to solve()). Estimating D from the
# per-unit IF differences xi_i = psi_U,i - psi_R,i makes it positive
# semi-definite in finite samples (the footnote to eqn (5.3) of Chen,
# Sant'Anna & Xie 2025). Caveat (Andrews 1987): the chi^2(rank) limit on the
# rank-deficient branch additionally requires the estimated rank to be
# consistent for the true rank; the eigenvalue threshold is the standard
# practical device but is not a formal guarantee.
#
# `v_scale` is the absolute variance scale of the constituent estimators (the
# largest per-coordinate asymptotic variance of either fit; 1 if unknown). It
# feeds the degenerate-contrast guard, the joint-path mirror of the scalar
# .edid_scalar_hausman() guard: when the two estimators coincide (e.g. both
# fits pinned to the same just-identified moments by the thin-cohort guard),
# xi is pure float dust (differences ~1e-17, D entries ~1e-33) and the
# RELATIVE eigenvalue threshold below would still "find" rank in that noise,
# returning an arbitrary large H with a tiny p-value. A D that is negligible
# on the ABSOLUTE scale of the estimators' own variances carries no testable
# contrast: report H = 0, df = 0, p = 1 with degenerate = TRUE instead.
.edid_if_diff_quadform <- function(d, xi, n, cluster_indices, v_scale = 1) {
  xi <- as.matrix(xi)
  D  <- n * cluster_cov_edid(xi, cluster_indices, n)     # = E_n[xi xi'] when iid
  if (any(!is.finite(D))) {
    return(list(statistic = NA_real_, df = NA_integer_, p_value = NA_real_, D = D,
                degenerate = NA))
  }
  eps_D <- .Machine$double.eps^0.5
  if (max(abs(D)) <= eps_D * max(v_scale, 1)) {          # degenerate contrast: estimators coincide
    return(list(statistic = 0, df = 0L, p_value = 1, D = D, degenerate = TRUE))
  }
  ev  <- eigen(D, symmetric = TRUE)
  mx  <- max(ev$values, 0)
  tol <- mx * sqrt(.Machine$double.eps)
  pos <- ev$values > tol
  rk  <- sum(pos)
  if (rk == 0L) {                                        # degenerate D: no power, report p = 1
    return(list(statistic = 0, df = 0L, p_value = 1, D = D, degenerate = TRUE))
  }
  if (rk == length(d)) {
    H <- as.numeric(n * crossprod(d, solve(D, d)))       # full rank: exact inverse
  } else {
    V <- ev$vectors[, pos, drop = FALSE]
    H <- as.numeric(n * crossprod(d, V %*% (crossprod(V, d) / ev$values[pos])))
  }
  list(statistic = H, df = rk,
       p_value = stats::pchisq(H, df = rk, lower.tail = FALSE), D = D,
       degenerate = FALSE)
}

# Scalar Hausman component H = n d^2 / D with the degenerate-D guard of
# Theorem 5.2 (D > 0 is required; xi ~ 0 makes the statistic 0/0, so report
# H = 0 / p = 1 instead of NaN). The guard is relative to the parameter's
# asymptotic variance scale.
.edid_scalar_hausman <- function(d, xi_vec, n, cluster_indices, v_scale = 1) {
  D <- as.numeric(n * cluster_cov_edid(matrix(xi_vec, ncol = 1L), cluster_indices, n))
  eps_D <- .Machine$double.eps^0.5
  if (!is.finite(D) || D <= eps_D * max(v_scale, 1)) {
    return(list(D = D, H = 0, p_value = 1, degenerate = TRUE))
  }
  H <- n * d^2 / D
  list(D = D, H = H, p_value = stats::pchisq(H, df = 1, lower.tail = FALSE),
       degenerate = FALSE)
}

# ---------------------------------------------------------------------------
# edid_hausman
# ---------------------------------------------------------------------------

#' Hausman test of PT-All against PT-Post for edid fits
#'
#' Implements the Hausman-type specification test of Theorem 5.1 in Chen,
#' Sant'Anna & Xie (2025): it compares the efficient event-study estimator
#' \eqn{\widehat{ES}} (consistent and semiparametrically efficient under
#' PT-All) with the conservative just-identified estimator
#' \eqn{\widecheck{ES}} of eqns (5.1)-(5.2) (consistent under PT-Post alone),
#' via the statistic of eqn (5.3),
#' \deqn{\widehat{H} = n\,(\widehat{ES} - \widecheck{ES})'\,\widehat{D}^{-1}\,
#'   (\widehat{ES} - \widecheck{ES}),}
#' where \eqn{\widehat{D}} is estimated from the per-unit difference of the two
#' estimators' influence functions, \eqn{\xi_i = \psi_{U,i} - \psi_{R,i}} --- the
#' positive semi-definite rendering noted in the footnote to eqn (5.3). Under
#' PT-All, \eqn{\widehat{H} \overset{d}{\to} \chi^2(|\mathcal{E}|)}; rejection
#' is evidence against the additional moment restrictions that PT-All imposes
#' beyond PT-Post.
#'
#' @param fit_unrestricted An \code{edid_fit} from
#'   \code{edid(..., pt_assumption = "post")}: the conservative just-identified
#'   estimator (the paper's staggered just-identification corollary),
#'   consistent under PT-Post alone. Its event-study aggregation (via
#'   \code{did::aggte}) includes the cohort-share weight-estimation
#'   influence-function correction, matching the conservative estimator used in
#'   the paper's empirical application.
#' @param fit_restricted An \code{edid_fit} from
#'   \code{edid(..., pt_assumption = "all")}: the efficient estimator under
#'   PT-All. Both fits must be estimated on the same data with the same
#'   clustering.
#' @param parameter \code{"event_study"} (default) for the joint test over the
#'   post-treatment event-study coefficients \eqn{ES(e), e \in \mathcal{E}}, or
#'   \code{"overall"} for the scalar test on \eqn{ES_{\mathrm{avg}}} (the
#'   average of \eqn{ES(e)} over \eqn{e \ge 0}).
#' @param e_set Numeric vector of post-treatment event times defining
#'   \eqn{\mathcal{E}}, or \code{NULL} (default: the intersection of the two
#'   fits' finite post-treatment event times). Ignored for
#'   \code{parameter = "overall"}.
#'
#' @details
#' The joint statistic uses \eqn{df = \mathrm{rank}(\widehat{D})} by eigenvalue
#' thresholding with a Moore-Penrose pseudoinverse on the rank-deficient branch
#' (the generically full-rank case reproduces the exact-inverse statistic with
#' \eqn{df = |\mathcal{E}|}). Following Andrews (1987), the
#' \eqn{\chi^2(\mathrm{rank})} limit under rank deficiency additionally
#' requires the estimated rank to be consistent; the threshold is the standard
#' practical device, not a formal guarantee. The covariance \eqn{\widehat{D}}
#' is cluster-robust when the fits carry cluster assignments.
#'
#' The returned object also reports the scalar per-coordinate statistics
#' \eqn{H_{\theta,n} = n(\widehat\theta_U - \widehat\theta_R)^2/\widehat{D}}
#' of eqn (5.5) for each \eqn{ES(e)} and for \eqn{ES_{\mathrm{avg}}}, with a
#' degenerate-\eqn{\widehat{D}} guard (coordinates where the two estimators
#' coincide report \eqn{H = 0}, \eqn{p = 1}). The joint statistic carries the
#' same guard: when the two estimators coincide on every coordinate --- e.g.
#' both fits pinned to the same just-identified moments by the thin-cohort
#' guard, so \eqn{\widehat{D}} is numerical noise relative to the estimators'
#' own variances --- the joint contrast is degenerate and is reported as
#' \eqn{H = 0}, \eqn{df = 0}, \eqn{p = 1} (with a message and
#' \code{degenerate = TRUE}) rather than ranking the noise.
#'
#' @return An object of class \code{edid_hausman}: a list with elements
#'   \code{statistic}, \code{df}, \code{p_value} (the joint test),
#'   \code{degenerate} (\code{TRUE} when the joint contrast was degenerate and
#'   the \eqn{H = 0}, \eqn{df = 0}, \eqn{p = 1} guard applied), \code{d}
#'   (the estimate difference vector, unrestricted minus restricted), \code{D}
#'   (the estimated asymptotic covariance of \eqn{\sqrt{n}\,d}), \code{scalar}
#'   (data.frame of per-coordinate eqn (5.5) statistics, including an
#'   \code{ES_avg} row), \code{parameter}, \code{e_set}, \code{n},
#'   \code{clustered}, plus the round-3 sanity guards:
#'   \code{leg_unstable} (\code{TRUE} when a constituent fit is numerically
#'   degenerate -- extreme propensity ratios, a non-credible weight channel,
#'   cross-cohort hedges carrying the estimand, or non-finite reported SEs --
#'   so a non-rejection is \emph{hollow}; a loud warning is also emitted),
#'   \code{leg_reasons} (the per-leg breakage descriptions; empty when healthy),
#'   \code{few_clusters} (\code{TRUE} when the fits carry fewer than 5 clusters,
#'   so the cluster-robust statistic is unreliable -- a few-cluster artifact,
#'   not PT evidence), and \code{n_clusters}.
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025). Efficient
#'   Difference-in-Differences and Event Study Estimators. Section 5.1,
#'   Theorem 5.1. \cr
#'   Hausman, J. A. (1978). Specification Tests in Econometrics.
#'   \emph{Econometrica}, 46(6), 1251-1271. \cr
#'   Andrews, D. W. K. (1987). Asymptotic Results for Generalized Wald Tests.
#'   \emph{Econometric Theory}, 3(3), 348-358.
#'
#' @seealso \code{\link{edid}}, \code{\link{edid_sargan}},
#'   \code{\link{edid_frontier}}, \code{\link{edid_adaptive}}
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
#' edid_hausman(fit_U, fit_R)
#' }
#'
#' @export
edid_hausman <- function(fit_unrestricted, fit_restricted,
                         parameter = c("event_study", "overall"),
                         e_set = NULL) {
  parameter <- match.arg(parameter)
  .edid_toolkit_check_fits(fit_unrestricted, fit_restricted)

  n  <- fit_restricted$n
  ci <- fit_restricted$cluster_indices

  # Broken-leg sanity guard (a hollow non-rejection footgun, Nguyen/Bailey-GB/ACA gate
  # evidence): if EITHER leg is a numerically degenerate fit -- extreme propensity ratios,
  # a non-credible weight channel, cross-cohort hedges carrying the estimand, or non-finite
  # SEs -- the Hausman contrast inherits the breakage and a non-rejection means nothing.
  # Warn loudly, naming the unhealthy leg(s) and reason(s), and flag the result object
  # ($leg_unstable / $leg_reasons) so a hollow p ~ 0.3-0.95 is not mistaken for a pass.
  health_U <- .edid_leg_health(fit_unrestricted, "fit_unrestricted (PT-Post)")
  health_R <- .edid_leg_health(fit_restricted,  "fit_restricted (PT-All)")
  leg_unstable <- health_U$broken || health_R$broken
  leg_reasons  <- character(0L)
  if (leg_unstable) {
    msgs <- c(if (health_R$broken) sprintf("%s: %s", health_R$label, paste(health_R$reasons, collapse = "; ")),
              if (health_U$broken) sprintf("%s: %s", health_U$label, paste(health_U$reasons, collapse = "; ")))
    leg_reasons <- msgs
    warning(paste0(
      "edid_hausman: a constituent fit is numerically degenerate, so this test is HOLLOW -- ",
      "a non-rejection here carries no evidence (the contrast inherits the broken leg's ",
      "instability). ", paste(msgs, collapse = " | "),
      ". Repair the fit (e.g. moment_set = \"own\" to drop cross-cohort pairs, weight_scheme = ",
      "\"averaged\", or a low-dimensional covariate index) before reading the p-value; the result ",
      "is returned with $leg_unstable = TRUE."), call. = FALSE)
  }

  if (parameter == "event_study") {
    e_set <- .edid_shared_e_set(fit_unrestricted, fit_restricted, e_set)
    pU <- .edid_param_ifs(fit_unrestricted, "event_study", e_set)
    pR <- .edid_param_ifs(fit_restricted,  "event_study", e_set)
  } else {
    pU <- .edid_param_ifs(fit_unrestricted, "overall")
    pR <- .edid_param_ifs(fit_restricted,  "overall")
    e_set <- NULL
  }

  d  <- pU$est - pR$est                  # unrestricted minus restricted
  xi <- pU$IF - pR$IF                    # per-unit IF difference (n x |E|)

  # Absolute variance scale of the two estimators (largest per-coordinate
  # asymptotic variance of either fit), for the degenerate-contrast guard.
  v_scale <- max(diag(as.matrix(n * cluster_cov_edid(pU$IF, ci, n))),
                 diag(as.matrix(n * cluster_cov_edid(pR$IF, ci, n))))
  joint <- .edid_if_diff_quadform(d, xi, n, ci, v_scale = v_scale)
  if (isTRUE(joint$degenerate)) {
    message("edid_hausman: the two estimators coincide (the IF-difference covariance is at ",
            "numerical-noise scale relative to the estimators' own variances), so the joint ",
            "contrast is degenerate and is reported as H = 0, df = 0, p = 1. This is expected ",
            "when the fits carry no over-identifying content to test -- e.g. when the ",
            "thin-cohort guard has pinned every cell to its just-identified moment.")
  }

  # Few-cluster guard: with G < EDID_FEWCLUSTER_MIN clusters the cluster-robust D is too
  # noisy / degenerate for the chi-square reference (the ACA gate's 2-3-state cohorts give
  # H ~ 350, p ~ 1e-73 -- a few-cluster artifact, not PT evidence). Flag the statistic as
  # unreliable (message + field); it is still returned.
  n_clusters <- .edid_n_clusters(fit_restricted)
  few_clusters <- !is.null(ci) && n_clusters < EDID_FEWCLUSTER_MIN
  if (few_clusters) {
    warning(sprintf(paste0(
      "edid_hausman: only %d cluster(s) at the clustering level (clustervars). The cluster-robust ",
      "statistic is UNRELIABLE below %d clusters -- the cluster-level covariance is too noisy / ",
      "degenerate for the chi-square reference, so the p-value is not trustworthy (flagged ",
      "$few_clusters = TRUE). Report unit-level (unclustered) toolkit statistics, or aggregate to ",
      "a coarser level with enough clusters."), n_clusters, EDID_FEWCLUSTER_MIN), call. = FALSE)
  }

  # Scalar eqn (5.5) statistics: each ES(e) plus ES_avg (always included).
  oU <- .edid_param_ifs(fit_unrestricted, "overall")
  oR <- .edid_param_ifs(fit_restricted,  "overall")
  lab_e <- if (parameter == "event_study") pU$e else numeric(0L)
  rows  <- vector("list", length(lab_e) + 1L)
  for (j in seq_along(lab_e)) {
    vR <- as.numeric(n * cluster_cov_edid(pR$IF[, j, drop = FALSE], ci, n))
    sc <- .edid_scalar_hausman(d[j], xi[, j], n, ci, v_scale = vR)
    rows[[j]] <- data.frame(
      parameter = sprintf("ES(%g)", lab_e[j]), e = lab_e[j],
      theta_U = pU$est[j], theta_R = pR$est[j], difference = d[j],
      D = sc$D, H = sc$H, p_value = sc$p_value, stringsAsFactors = FALSE)
  }
  d_ov  <- oU$est - oR$est
  xi_ov <- oU$IF[, 1L] - oR$IF[, 1L]
  vR_ov <- as.numeric(n * cluster_cov_edid(oR$IF, ci, n))
  sc_ov <- .edid_scalar_hausman(d_ov, xi_ov, n, ci, v_scale = vR_ov)
  rows[[length(rows)]] <- data.frame(
    parameter = "ES_avg", e = NA_real_,
    theta_U = oU$est, theta_R = oR$est, difference = d_ov,
    D = sc_ov$D, H = sc_ov$H, p_value = sc_ov$p_value, stringsAsFactors = FALSE)
  scalar <- do.call(rbind, rows)
  rownames(scalar) <- NULL

  out <- list(
    statistic  = joint$statistic,
    df         = joint$df,
    p_value    = joint$p_value,
    degenerate = isTRUE(joint$degenerate),
    d          = stats::setNames(d, if (parameter == "event_study") sprintf("e=%g", pU$e) else "overall"),
    D          = joint$D,
    scalar     = scalar,
    parameter  = parameter,
    e_set      = e_set,
    n          = n,
    clustered  = !is.null(ci),
    leg_unstable = isTRUE(leg_unstable),   # a constituent fit is numerically degenerate -> hollow test
    leg_reasons  = leg_reasons,            # per-leg breakage descriptions (empty when healthy)
    few_clusters = isTRUE(few_clusters),   # G < EDID_FEWCLUSTER_MIN -> cluster-robust stat unreliable
    n_clusters   = n_clusters,
    alpha      = fit_restricted$alpha %||% 0.05
  )
  class(out) <- c("edid_hausman", "list")
  out
}

#' @describeIn edid_hausman Print method.
#' @param x an \code{edid_hausman} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_hausman <- function(x, digits = 4, ...) {
  cat("\nHausman test of PT-All vs PT-Post (Chen, Sant'Anna & Xie 2025, Theorem 5.1)\n")
  param_lab <- if (identical(x$parameter, "event_study")) {
    sprintf("event study, E = {%s}", paste(x$e_set, collapse = ", "))
  } else "overall ES_avg"
  cat(sprintf("  Parameter: %s%s\n", param_lab,
              if (isTRUE(x$clustered)) " (cluster-robust)" else ""))
  cat(sprintf("  H = %s on %s df (rank of D-hat), p-value = %s\n",
              format(x$statistic, digits = digits), format(x$df),
              format.pval(x$p_value, digits = digits)))
  if (isTRUE(x$degenerate)) {
    cat("  NOTE: degenerate contrast -- the two estimators coincide, so there is no\n")
    cat("  over-identifying content to test (H = 0, df = 0, p = 1 by construction).\n")
  }
  if (isTRUE(x$leg_unstable)) {
    cat("  WARNING: HOLLOW TEST -- a constituent fit is numerically degenerate, so a\n")
    cat("  non-rejection carries no evidence. Reason(s):\n")
    for (r in x$leg_reasons) cat("    - ", r, "\n", sep = "")
  }
  if (isTRUE(x$few_clusters)) {
    cat(sprintf("  WARNING: only %d cluster(s) -- the cluster-robust statistic is unreliable\n", x$n_clusters))
    cat("  below ", EDID_FEWCLUSTER_MIN, " clusters (few-cluster artifact, not PT evidence).\n", sep = "")
  }
  cat("\nPer-parameter scalar statistics (eqn 5.5):\n")
  tab <- x$scalar
  num <- vapply(tab, is.numeric, logical(1L))
  tab[num] <- lapply(tab[num], function(z) signif(z, digits))
  print(tab, row.names = FALSE)
  cat("\nH0: parallel trends holds across all groups and pre-treatment periods (PT-All).\n")
  invisible(x)
}
