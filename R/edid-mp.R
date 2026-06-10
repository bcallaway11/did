# did-compatible MP construction for edid.
# Lets the did aggregation/inference ecosystem (aggte, tidy, ggdid, summary) operate on edid output,
# so edid does not need its own parallel aggregation. Builds a did::MP object from an edid fit.

#' Build a \code{did::MP} object from an \code{edid} fit
#'
#' Constructs the same \code{MP} object that \code{att_gt()} returns, populated with edid's
#' group-time estimates and their influence functions, so that \code{did::aggte()} (and the rest of the
#' did ecosystem) can aggregate edid output unchanged. edid() always stores the influence functions.
#'
#' @param fit an \code{edid_fit} object returned by \code{\link{edid}}.
#' @param bstrap,biters,clustervars,cband optional overrides; default to the fit's effective
#'   settings. Clustered or bootstrap inference in \code{aggte()} then follows the did conventions.
#' @return a \code{did::MP} object (\code{group}, \code{t}, \code{att}, \code{inffunc}, \code{DIDparams}, ...).
#' @export
as_MP_edid <- function(fit, bstrap = NULL, biters = NULL, clustervars = NULL, cband = NULL) {
  if (!inherits(fit, "edid_fit")) stop("as_MP_edid() expects an 'edid_fit' object.")
  if (is.null(biters)) biters <- if (!is.null(fit$biters)) as.integer(fit$biters) else 1000L
  if (is.null(fit$eif)) {
    stop("as_MP_edid(): the fit does not contain influence functions ($eif); edid() stores them by default.")
  }
  agt <- fit$att_gt
  n   <- fit$n
  inffunc <- as.matrix(fit$eif)
  if (nrow(inffunc) != n || ncol(inffunc) != nrow(agt)) {
    stop(sprintf("as_MP_edid(): eif is %dx%d but expected %dx%d (n x n_cells).",
                 nrow(inffunc), ncol(inffunc), n, nrow(agt)))
  }

  # Time-invariant per-unit data sufficient for compute.aggte()'s group-probability weights:
  # one row per unit at the first period, with never-treated coded 0 (att_gt convention) and unit
  # sampling weight .w = 1 (edid does not use sampling weights).
  g_unit <- fit$unit_cohorts
  g_unit[!is.finite(g_unit)] <- 0
  period_1 <- min(fit$time_periods)
  tinv <- data.frame(fit$all_units, period_1, g_unit, 1)
  names(tinv) <- c(fit$idname, fit$tname, fit$gname, ".w")

  glist <- sort(fit$treatment_groups[is.finite(fit$treatment_groups) & fit$treatment_groups != 0])
  if (is.null(bstrap))      bstrap      <- isTRUE(fit$bstrap) && identical(fit$cband_method, "multiplier")
  if (is.null(cband))       cband       <- isTRUE(fit$cband) && isTRUE(bstrap)
  if (is.null(clustervars)) clustervars <- fit$clustervars
  if (!is.null(clustervars) && is.null(fit$cluster_indices)) {
    stop("as_MP_edid(): `clustervars` requested but the fit carries no cluster assignments; ",
         "refit edid() with `clustervars` to enable clustered aggregation.")
  }
  # cluster column (EIF-aligned), stored under a reserved name: writing it under the caller's
  # column name overwrites the cohort/id column of tinv whenever clustervars coincides with
  # gname/idname, silently corrupting compute.aggte()'s group shares. mboot() resolves the
  # column through DIDparams$clustervars, so the reserved name is self-consistent.
  if (!is.null(clustervars) && !is.null(fit$cluster_indices)) {
    tinv[[".edid_cluster"]] <- fit$cluster_indices
    clustervars <- ".edid_cluster"
  }

  dp <- list(
    yname = NULL, tname = fit$tname, idname = fit$idname, gname = fit$gname,
    data = tinv, panel = TRUE, faster_mode = FALSE,
    tlist = sort(fit$time_periods), glist = glist,
    nG = length(glist), nT = length(sort(fit$time_periods)), est_method = "edid",
    control_group = "nevertreated", anticipation = fit$anticipation,
    bstrap = bstrap, biters = biters, alp = fit$alpha, cband = cband,
    clustervars = clustervars, cluster_vector = fit$cluster_indices, n = n
  )

  mp <- list(
    group = agt$group, t = agt$time, att = agt$att,
    V_analytical = NULL, se = agt$se, c = stats::qnorm(1 - fit$alpha / 2),
    inffunc = inffunc, n = n, W = NULL, Wpval = NULL,
    aggte = NULL, alp = fit$alpha, DIDparams = dp
  )
  class(mp) <- "MP"
  mp
}
