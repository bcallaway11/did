# edid-pairs.R
# Enumerate the set of valid comparison pairs H_gt for a target cell (g, t).

#' Enumerate valid comparison pairs for a target (g, t) cell
#'
#' Constructs the set \eqn{H_{gt}} of valid \code{(gp, tpre)} pairs used to
#' form identifying DiD moments for cohort \code{target_g} at time \code{target_t}.
#'
#' Under \strong{PT-Post}: returns exactly one pair \code{(Inf, tpre)} with \code{tpre}
#' the most recent observed period strictly before \code{target_g - anticipation}
#' (\code{= target_g - 1 - anticipation} on a unit-spaced grid), or a 0-row data.frame
#' if no observed period precedes the (anticipation-adjusted) treatment onset. When
#' \code{tpre == period_1} the pair is the standard 2x2 DiD moment.
#'
#' Under \strong{PT-All}: iterates over treated cohorts \code{gp} only (the
#' never-treated group is the time control inside every moment, not a comparison
#' cohort). For \code{gp == target_g}: valid \code{tpre} are all periods strictly
#' less than \code{gp - anticipation}, including \code{period_1} (this is the
#' degenerate CS DiD moment whose comparison-cohort EIF term is identically zero).
#' For \code{gp != target_g}: valid \code{tpre} are periods strictly between
#' \code{period_1} and \code{gp - anticipation} (exclusive on both ends).
#' Returns a 0-row data.frame if no valid pairs exist (e.g., single cohort with
#' only one pre-period equal to \code{period_1}).
#'
#' When \code{moment_set} is supplied (advanced; see \code{\link{edid}}), the
#' enumerated pairs for \code{target_g} are intersected with the user-supplied
#' \code{(gp, tpre)} rows for that cohort: rows of \code{moment_set} that are
#' not part of the enumeration are silently ignored (the mechanism can only
#' restrict, never extend, the set of valid identifying moments). With
#' \code{moment_set = NULL} (default) the enumeration is unchanged.
#'
#' @param target_g scalar: treatment cohort being estimated
#' @param treatment_groups sorted numeric vector of all finite cohort values
#' @param time_periods sorted numeric vector of all time periods in the panel
#' @param period_1 scalar: universal first period
#' @param pt_assumption character: \code{"all"} or \code{"post"}
#' @param anticipation integer >= 0
#' @param never_treated_val value used to represent the never-treated cohort
#'   (default \code{Inf})
#' @param moment_set \code{NULL} (default: no restriction) or a data.frame with
#'   columns \code{g}, \code{gp}, \code{tpre} restricting the enumerated pairs
#'   per target cohort (intersection semantics).
#'
#' @return data.frame with columns \code{gp} (comparison cohort) and
#'   \code{tpre} (pre-period). May have 0 rows.
#' @keywords internal
enumerate_valid_pairs_edid <- function(
  target_g,
  treatment_groups,
  time_periods,
  period_1,
  pt_assumption,
  anticipation     = 0L,
  never_treated_val = Inf,
  moment_set       = NULL
) {
  empty <- data.frame(gp = numeric(0L), tpre = numeric(0L))

  # Restrict the enumerated pairs to the user-supplied moment set for this target
  # cohort (intersection; rows not in the enumeration are ignored). NULL = no-op,
  # keeping the default path byte-identical.
  .restrict <- function(pairs) {
    if (is.null(moment_set) || nrow(pairs) == 0L) return(pairs)
    ms_g <- moment_set[moment_set$g == target_g, , drop = FALSE]
    keep <- paste(pairs$gp, pairs$tpre) %in% paste(ms_g$gp, ms_g$tpre)
    out  <- pairs[keep, , drop = FALSE]
    rownames(out) <- NULL
    out
  }

  if (pt_assumption == "post") {
    # -----------------------------------------------------------------------
    # PT-Post: exactly one pair (Inf, baseline), baseline = the last observed period strictly before the
    # effective onset g - anticipation (the most recent clean pre-treatment period). On an integer-spaced
    # grid this equals the period at or before g-1-anticipation; the strict "< g - anticipation" form is
    # also exact on non-integer grids (e.g. periods {1, 1.5, 2, 3}, g = 2 -> baseline 1.5, where the
    # "<= g-1-anticipation" arithmetic would skip back to 1). Under irregular spacing the previous
    # observed period is used rather than dropping the cohort. When baseline == period_1 this is the
    # standard 2x2 DiD.
    # -----------------------------------------------------------------------
    pre_periods <- time_periods[time_periods < target_g - anticipation]
    if (!length(pre_periods)) return(empty)
    return(.restrict(data.frame(gp = never_treated_val, tpre = max(pre_periods))))
  }

  # -------------------------------------------------------------------------
  # PT-All: loop over treated cohorts only (never-treated is NOT a comparison
  # cohort; it appears only as the time control E[Y_inf(t)-Y_inf(tpre)] inside
  # each moment).
  #
  # For g' == target_g: valid tpre = {s : s < eff_start(g')}
  #   -- INCLUDES period_1 (degenerate CS DiD moment; comparison EIF = 0)
  # For g' != target_g: valid tpre = {s : period_1 < s < eff_start(g')}
  #   -- EXCLUDES period_1 (non-degenerate moments only)
  # -------------------------------------------------------------------------
  out_gp   <- numeric(0L)
  out_tpre <- numeric(0L)

  for (gp in treatment_groups) {
    eff_start <- gp - anticipation
    if (gp == target_g) {
      # Self-pair: include period_1
      valid_tpre <- time_periods[time_periods < eff_start]
    } else {
      # Cross-pair: exclude period_1
      valid_tpre <- time_periods[
        time_periods > period_1 & time_periods < eff_start
      ]
    }
    if (length(valid_tpre) > 0L) {
      out_gp   <- c(out_gp,   rep(gp, length(valid_tpre)))
      out_tpre <- c(out_tpre, valid_tpre)
    }
  }

  if (length(out_gp) == 0L) return(empty)
  .restrict(data.frame(gp = out_gp, tpre = out_tpre, stringsAsFactors = FALSE))
}
