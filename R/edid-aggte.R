# edid-aggte.R
# Aggregation for edid_fit objects: a thin wrapper that builds a did MP (as_MP_edid) and delegates to
# did::aggte(). edid aggregation therefore follows the published CS2021 aggregation definitions exactly
# and inherits did's AGGTEobj methods (print / summary / tidy / ggdid). This was verified numerically
# equivalent to the previous edid-native aggregation to machine precision (att ~4e-16, se ~6e-17).

#' Aggregate edid_fit estimates
#'
#' Aggregates the group-time ATT(g,t) estimates from an \code{edid_fit} object using the same interface
#' and definitions as \code{\link[did]{aggte}}. Internally it builds a \code{did::MP} object from the
#' edid estimates and their influence functions (\code{\link{as_MP_edid}}) and calls
#' \code{did::aggte()}, so the result is a standard \code{did::AGGTEobj}.
#'
#' @param edid_fit_obj An \code{edid_fit} object returned by \code{\link{edid}}.
#' @param type Character scalar, mirroring \code{did::aggte()}: \code{"simple"} (cohort-share-weighted
#'   average over post-treatment cells), \code{"dynamic"} (event-study: average of \eqn{ES(e)} over
#'   \eqn{e \ge 0}), \code{"group"} (cohort-level overalls), or \code{"calendar"} (calendar-period overalls).
#' @param balance_e Integer or \code{NULL}: if not \code{NULL}, balances the cohort composition of
#'   the dynamic aggregation (as in \code{did::aggte}): cohorts observed for fewer than
#'   \code{balance_e} post-treatment periods are dropped, and event times
#'   \eqn{e \in [\text{balance\_e} - (T_{\max} - T_{\min}),\ \text{balance\_e}]} are reported, so
#'   every reported \eqn{e} averages over the same set of cohorts.
#' @param min_e,max_e Numeric: minimum/maximum relative time to include in dynamic output.
#' @param na.rm Logical: drop \code{NA} ATT entries before aggregating. Default \code{FALSE}.
#' @param seed Integer or \code{NULL}: RNG seed for the multiplier-bootstrap path (\code{cband_method =
#'   "multiplier"}). Defaults to the seed stored on the fit, so standalone \code{aggte_edid()} bootstrap
#'   results are reproducible; the caller's RNG stream is restored on exit. Ignored on the analytic path.
#'
#' @return A \code{did::AGGTEobj} (as returned by \code{\link[did]{aggte}}), so the did \code{print},
#'   \code{summary}, and \code{tidy} methods apply.
#'
#' @seealso \code{\link{edid}}, \code{\link{as_MP_edid}}, \code{\link[did]{aggte}}
#' @export
aggte_edid <- function(
  edid_fit_obj,
  type      = c("simple", "dynamic", "group", "calendar"),
  balance_e = NULL,
  min_e     = -Inf,
  max_e     = Inf,
  na.rm     = FALSE,
  seed      = NULL
) {
  mc   <- match.call()
  type <- match.arg(type)
  if (!inherits(edid_fit_obj, "edid_fit")) {
    stop("`edid_fit_obj` must be an object of class `edid_fit` returned by edid().")
  }
  # Feasibility of balance_e: did::aggte keeps only cohorts observed for >= balance_e
  # post-treatment periods; past the longest available window no cohort qualifies and
  # compute.aggte fails with an opaque subscript error.
  if (!is.null(balance_e) && identical(type, "dynamic")) {
    g_fin <- edid_fit_obj$treatment_groups[is.finite(edid_fit_obj$treatment_groups) &
                                           edid_fit_obj$treatment_groups != 0]
    if (length(g_fin) > 0L) {
      max_e <- max(edid_fit_obj$time_periods) - min(g_fin)
      if (balance_e > max_e) {
        stop(sprintf(paste0(
          "`balance_e` = %s exceeds the longest available post-treatment window: the earliest ",
          "cohort (g = %s) is observed for at most e = %s post-treatment periods. ",
          "Choose balance_e <= %s."),
          format(balance_e), format(min(g_fin)), format(max_e), format(max_e)), call. = FALSE)
      }
    }
  }
  # Inference path follows how the fit was produced:
  #  - cband_method = "analytic" (default): aggregate analytically (bstrap = FALSE), then replace the
  #    simultaneous critical value with the MOPM sup-t crit from the cluster-robust covariance of the
  #    aggregate influence functions (no bootstrap; the only path that can carry the higher-order term).
  #  - cband_method = "multiplier": the did multiplier bootstrap (legacy) when edid(bstrap = TRUE).
  use_analytic <- identical(edid_fit_obj$cband_method, "analytic")
  do_boot      <- (!use_analytic && isTRUE(edid_fit_obj$bstrap))
  # Reproducible multiplier bootstrap: seed the RNG (default to the fit's seed) before aggte() -> mboot(),
  # which draws from the global stream. Save/restore .Random.seed so the caller's RNG stream is undisturbed.
  if (do_boot) {
    if (is.null(seed)) seed <- edid_fit_obj$seed
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
      }
      set.seed(seed)
    }
  }
  boot_cband <- do_boot && isTRUE(edid_fit_obj$cband)
  a <- aggte(as_MP_edid(edid_fit_obj, bstrap = do_boot, cband = boot_cband),
             type = type, balance_e = balance_e,
             min_e = min_e, max_e = max_e, na.rm = na.rm,
             bstrap = do_boot, cband = boot_cband)
  if (use_analytic) {
    # Closure that replays THIS aggregation on an att(g,t) vector perturbed at one cell, returning the
    # aggregate per-element att.egt. The aggregation weights do not depend on the att VALUES, so finite-
    # differencing this map recovers the constant cell -> aggregate linear map A (used by the higher-order
    # Wick refinement to map Sigma_quad to aggregate scale). Same aggregation arguments as the call above.
    reaggregate <- function(att_vec) {
      f2 <- edid_fit_obj
      f2$att_gt$att <- att_vec
      aa <- aggte(as_MP_edid(f2, bstrap = FALSE, cband = FALSE), type = type, balance_e = balance_e,
                  min_e = min_e, max_e = max_e, na.rm = na.rm, bstrap = FALSE)
      aa$att.egt %||% aa$overall.att
    }
    a <- .edid_analytic_cband_agg(a, edid_fit_obj, reaggregate)
  }
  a$call <- mc
  a
}

# Replace a did::AGGTEobj's simultaneous critical value (crit.val.egt) with the analytic MOPM sup-t crit
# computed from the cluster-robust covariance of the aggregate per-element influence functions. SEs are
# left untouched (already analytic from aggte(bstrap = FALSE)); only the uniform-band crit changes. Simple
# (single overall) aggregations have no per-element vector, so nothing changes there. When the fit carries
# the higher-order ("Wick") refinement, the aggregate covariance becomes Sigma1_agg + A Sigma_quad A',
# with A the (finite-difference-recovered) cell -> aggregate linear map; the crit then matches the SEs that
# the cell-level path inflated, keeping the aggregate band higher-order-aware too.
.edid_analytic_cband_agg <- function(a, fit, reaggregate = NULL) {
  # Align the clustered finite-sample convention with the cell-level SEs: edid's cell SEs and
  # vcov() apply the G/(G-1) finite-cluster factor (cluster_cov_edid), while did's getSE() does
  # not, so without this the aggregate SEs (and the uniform bands built on them) contradict the
  # cell-level convention of the same fit by sqrt(G/(G-1)). Applied before the higher-order
  # increment below, which is already in the corrected convention (it comes from sigma_quad's
  # cluster-robust V).
  if (!is.null(fit$cluster_indices)) {
    G_cl <- length(unique(fit$cluster_indices))
    if (G_cl > 1L) {
      cl_fac <- sqrt(G_cl / (G_cl - 1))
      if (!is.null(a$se.egt))     a$se.egt     <- a$se.egt * cl_fac
      if (!is.null(a$overall.se)) a$overall.se <- a$overall.se * cl_fac
    }
  }
  g <- .edid_agg_if(a)
  # Total second-order covariance increment: the higher-order ("Wick") Sigma_quad (covariate path), plus the
  # diagonal no-covariate weight-estimation increment (estimation_effect on a no-covariate fit) -- the cell
  # SEs already carry the latter, so the aggregate SEs/bands must add A Sigma A' for cell/aggregate
  # consistency. NULL on fits with neither (the entire classic path), keeping those byte-identical.
  Sigma_so <- .edid_secondorder_sigma(fit)
  if (is.null(g$egt) && !is.null(g$overall) && !is.null(Sigma_so) &&
      !is.null(reaggregate) && !is.null(fit$cells) && !is.null(a$overall.se)) {
    A <- .edid_recover_agg_map(a, fit, reaggregate)
    if (!is.null(A) && nrow(A) == 1L && ncol(A) == nrow(Sigma_so)) {
      inc <- drop(A %*% Sigma_so %*% t(A))
      if (is.finite(inc)) a$overall.se <- sqrt(max(a$overall.se^2 + inc, 0))
    }
  }
  if (!is.null(g$egt) && is.matrix(g$egt) && ncol(g$egt) >= 1L) {
    Sig <- cluster_cov_edid(g$egt, fit$cluster_indices, fit$n)
    if (!is.null(Sigma_so) && !is.null(reaggregate) && !is.null(fit$cells)) {
      A <- .edid_recover_agg_map(a, fit, reaggregate)      # n_agg x K, constant cell -> aggregate weights
      if (!is.null(A)) {
        HO  <- A %*% Sigma_so %*% t(A)        # aggregate-scale second-order covariance increment
        Sig <- Sig + HO
        # Make the aggregate SEs higher-order-aware too (not just the sup-t crit): the cell SEs already include
        # diag(Sigma_quad), so the event-study SEs must include diag(A Sigma_quad A') for cell/aggregate
        # consistency. ADD only this increment to did::aggte's analytic se.egt -- rather than replacing it with
        # sqrt(diag(Sig)), whose cluster_cov_edid finite-cluster convention can differ from did's getSE -- so the
        # non-higher_order path is byte-identical and higher_order only inflates.
        he_inc <- diag(HO)
        if (length(he_inc) == length(a$se.egt))
          a$se.egt <- sqrt(pmax(a$se.egt^2 + he_inc, 0))
        # The overall summary SE too: the overall IF is an EXACT linear combo of the per-element IFs (overall =
        # weighted mean of the dynamic effects). Recover the weights w (overall = egt %*% w) by a RANK-SAFE
        # least-squares solve and add only the higher-order term w'(A Sigma_quad A')w to the existing analytic
        # overall.se (same byte-identity property). A rank-deficient or inconsistent system (collinear
        # event-study influence columns) warns and skips the increment -- never silently (the previous plain
        # solve() swallowed the error and dropped the increment without a trace).
        if (!is.null(g$overall) && !is.null(a$overall.se) && is.matrix(g$egt) && ncol(g$egt) == nrow(HO)) {
          w <- .edid_recover_overall_weights(g$egt, g$overall)
          if (!is.null(w) && length(w) == nrow(HO))
            a$overall.se <- sqrt(max(a$overall.se^2 + drop(crossprod(w, HO %*% w)), 0))
        }
      }
    }
    if (isTRUE(fit$cband)) {
      a$crit.val.egt <- supt_crit_edid(Sig, alp = fit$alpha %||% 0.05, seed = fit$seed)
      # Record that a simultaneous crit is in force so summary.AGGTEobj labels the band
      # correctly (it would otherwise print "Pointwise" because bstrap/cband are FALSE on
      # the analytic path).
      a$DIDparams$cband <- TRUE
    }
  }
  a
}

# Total K x K second-order covariance increment carried by a fit: the higher-order ("Wick") Sigma_quad
# (covariate path, gated on fit$higher_order) plus the no-covariate weight-estimation diagonal
# (estimation_effect on a no-covariate fit; nocov_ee_sigma_edid). Each piece reuses the matrix cached on
# the fit when present and recomputes from fit$cells otherwise (hand-built fits / standalone aggte_edid
# calls; same inputs => bit-identical). Returns NULL when the fit carries neither, so classic fits take
# the unchanged byte-identical path.
.edid_secondorder_sigma <- function(fit) {
  out <- NULL
  if (isTRUE(fit$higher_order) && !is.null(fit$cells)) {
    out <- if (!is.null(fit$sigma_quad)) fit$sigma_quad
           else sigma_quad_edid(fit$cells, fit$cluster_indices, fit$n)
  }
  ee <- if (!is.null(fit$sigma_nocov_ee)) fit$sigma_nocov_ee
        else if (!is.null(fit$cells)) nocov_ee_sigma_edid(fit$cells)
        else NULL
  if (!is.null(ee)) out <- if (is.null(out)) ee else out + ee
  out
}

# Rank-safe recovery of the overall-aggregate weights w solving overall = egt %*% w, by a DIRECT SVD
# least squares on egt (not the normal equations, whose condition number is cond(egt)^2 and whose solve
# can leave a spuriously large residual on strongly correlated influence columns; no new dependencies).
# The overall influence function is an exact linear combination of the per-element event-study influence
# functions, so the system is consistent in exact arithmetic; but the egt columns can be COLLINEAR
# (e.g. duplicated cells feeding one event time, or a degenerate 2-period design), where the previous
# plain solve(crossprod(egt), .) threw and the increment was SILENTLY skipped. Returns w when the
# (full-rank) solution reproduces `overall`; on a rank-deficient system or a non-negligible recovery
# residual ||egt %*% w - overall|| it warns and returns NULL (the caller then SKIPS the higher-order
# overall.se increment -- a deliberate, audible skip: with collinear columns w is not unique and the
# quadratic increment w' HO w is not identified from the least-squares fit alone).
.edid_recover_overall_weights <- function(egt, overall) {
  .skip <- function(why) {
    warning(paste0(
      "higher_order: could not recover the overall-aggregate weights from the per-element influence ",
      "columns (", why, "); the higher-order increment to the OVERALL SE of this aggregation is skipped ",
      "(the per-element SEs and the uniform band keep their increment). This is expected for the 'group' ",
      "aggregation, whose overall influence function carries the estimated cohort-share weights and is ",
      "genuinely outside the column span."), call. = FALSE)
    NULL
  }
  sv <- tryCatch(svd(egt), error = function(e) NULL)
  if (is.null(sv) || !all(is.finite(sv$d))) return(.skip("SVD of the influence columns failed"))
  d_max <- max(sv$d)
  if (!is.finite(d_max) || d_max <= 0) return(.skip("the influence columns are zero"))
  tol      <- d_max * max(dim(egt)) * .Machine$double.eps
  pos      <- sv$d > tol
  rank_def <- any(!pos)
  d_inv    <- ifelse(pos, 1 / sv$d, 0)
  w        <- drop(sv$v %*% (d_inv * drop(crossprod(sv$u, overall))))
  resid    <- sqrt(sum((egt %*% w - overall)^2))
  scale_o  <- sqrt(sum(overall^2))
  if (rank_def) return(.skip("the columns are collinear: the weight vector is not unique"))
  if (!all(is.finite(w)) || resid > 1e-6 * max(scale_o, .Machine$double.eps))
    return(.skip("the recovery residual is non-negligible: the overall IF is not in the column span"))
  w
}

# Recover the constant cell -> aggregate linear map A (n_agg x K) by finite-differencing the aggregation:
# perturb each cell's att(g,t) by eps, re-aggregate, and read (att.egt(perturbed) - att.egt) / eps as that
# cell's column. The aggregation weights do not depend on the att values, so A is exact (the map is linear);
# NA / non-finite columns (NA cells, or cells the aggregation drops) are set to 0. Returns NULL if the
# baseline aggregate egt is unavailable.
.edid_recover_agg_map <- function(a, fit, reaggregate, eps = 1e-4) {
  base <- a$att.egt %||% a$overall.att
  if (is.null(base) || !length(base)) return(NULL)
  att0 <- fit$att_gt$att
  K    <- length(att0)
  A    <- matrix(0, length(base), K)
  for (k in seq_len(K)) {
    if (!is.finite(att0[k])) next                          # NA cell: contributes nothing to any aggregate
    att_p <- att0; att_p[k] <- att0[k] + eps
    col <- tryCatch((reaggregate(att_p) - base) / eps, error = function(e) NULL)
    if (!is.null(col) && length(col) == length(base) && all(is.finite(col))) A[, k] <- col
  }
  A
}
