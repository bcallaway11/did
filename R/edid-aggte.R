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
  na.rm     = FALSE
) {
  mc   <- match.call()
  type <- match.arg(type)
  if (!inherits(edid_fit_obj, "edid_fit")) {
    stop("`edid_fit_obj` must be an object of class `edid_fit` returned by edid().")
  }
  # bstrap follows how the fit was produced: edid(bstrap = TRUE) -> the aggregations use the did
  # multiplier bootstrap (simultaneous bands); otherwise analytical.
  a <- aggte(as_MP_edid(edid_fit_obj), type = type, balance_e = balance_e,
             min_e = min_e, max_e = max_e, na.rm = na.rm, bstrap = isTRUE(edid_fit_obj$bstrap))
  a$call <- mc
  a
}
