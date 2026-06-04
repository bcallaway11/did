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
#' @param balance_e Integer or \code{NULL}: if not \code{NULL}, restricts the dynamic aggregation to
#'   relative times in \eqn{[-\text{balance\_e}, \text{balance\_e}]}.
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
  a <- aggte(as_MP_edid(edid_fit_obj), type = type, balance_e = balance_e,
             min_e = min_e, max_e = max_e, na.rm = na.rm,
             bstrap = do_boot)
  if (use_analytic) {
    # Closure that replays THIS aggregation on an att(g,t) vector perturbed at one cell, returning the
    # aggregate per-element att.egt. The aggregation weights do not depend on the att VALUES, so finite-
    # differencing this map recovers the constant cell -> aggregate linear map A (used by the higher-order
    # Wick refinement to map Sigma_quad to aggregate scale). Same aggregation arguments as the call above.
    reaggregate <- function(att_vec) {
      f2 <- edid_fit_obj
      f2$att_gt$att <- att_vec
      aa <- aggte(as_MP_edid(f2), type = type, balance_e = balance_e,
                  min_e = min_e, max_e = max_e, na.rm = na.rm, bstrap = FALSE)
      aa$att.egt
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
  g <- .edid_agg_if(a)
  if (!is.null(g$egt) && is.matrix(g$egt) && ncol(g$egt) >= 1L) {
    Sig <- cluster_cov_edid(g$egt, fit$cluster_indices, fit$n)
    if (isTRUE(fit$higher_order) && !is.null(reaggregate) && !is.null(fit$cells)) {
      Sigma_quad <- sigma_quad_edid(fit$cells, fit$cluster_indices, fit$n)
      A <- .edid_recover_agg_map(a, fit, reaggregate)      # n_agg x K, constant cell -> aggregate weights
      if (!is.null(A)) {
        HO  <- A %*% Sigma_quad %*% t(A)      # aggregate-scale higher-order ("Wick") covariance increment
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
        # weighted mean of the dynamic effects). Recover the weights w (overall = egt %*% w) and add only the
        # higher-order term w'(A Sigma_quad A')w to the existing analytic overall.se (same byte-identity property).
        if (!is.null(g$overall) && !is.null(a$overall.se) && is.matrix(g$egt) && ncol(g$egt) == nrow(HO)) {
          w <- tryCatch(drop(solve(crossprod(g$egt), crossprod(g$egt, g$overall))), error = function(e) NULL)
          if (!is.null(w) && length(w) == nrow(HO))
            a$overall.se <- sqrt(max(a$overall.se^2 + drop(crossprod(w, HO %*% w)), 0))
        }
      }
    }
    if (isTRUE(fit$cband)) {
      a$crit.val.egt <- supt_crit_edid(Sig, alp = fit$alpha %||% 0.05, seed = fit$seed)
    }
  }
  a
}

# Recover the constant cell -> aggregate linear map A (n_agg x K) by finite-differencing the aggregation:
# perturb each cell's att(g,t) by eps, re-aggregate, and read (att.egt(perturbed) - att.egt) / eps as that
# cell's column. The aggregation weights do not depend on the att values, so A is exact (the map is linear);
# NA / non-finite columns (NA cells, or cells the aggregation drops) are set to 0. Returns NULL if the
# baseline aggregate egt is unavailable.
.edid_recover_agg_map <- function(a, fit, reaggregate, eps = 1e-4) {
  base <- a$att.egt
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
