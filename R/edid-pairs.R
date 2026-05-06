# edid-pairs.R
# Enumerate the set of valid comparison pairs H_gt for a target cell (g, t).

#' Enumerate valid comparison pairs for a target (g, t) cell
#'
#' Constructs the set \eqn{H_{gt}} of valid \code{(gp, tpre)} pairs used to
#' form identifying DiD moments for cohort \code{target_g} at time \code{target_t}.
#'
#' Under \strong{PT-Post}: returns exactly one pair \code{(Inf, target_g - 1 - anticipation)},
#' or a 0-row data.frame if that pre-period does not exist in \code{time_periods} or
#' equals \code{period_1}.
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
#' @param target_g scalar: treatment cohort being estimated
#' @param treatment_groups sorted numeric vector of all finite cohort values
#' @param time_periods sorted numeric vector of all time periods in the panel
#' @param period_1 scalar: universal first period
#' @param pt_assumption character: \code{"all"} or \code{"post"}
#' @param anticipation integer >= 0
#' @param never_treated_val value used to represent the never-treated cohort
#'   (default \code{Inf})
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
  never_treated_val = Inf
) {
  empty <- data.frame(gp = numeric(0L), tpre = numeric(0L))

  if (pt_assumption == "post") {
    # -----------------------------------------------------------------------
    # PT-Post: exactly one pair (Inf, g - 1 - anticipation)
    # -----------------------------------------------------------------------
    tpre_val <- target_g - 1L - anticipation
    if (!tpre_val %in% time_periods) return(empty)
    if (tpre_val == period_1)         return(empty)
    return(data.frame(gp = never_treated_val, tpre = tpre_val))
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
  data.frame(gp = out_gp, tpre = out_tpre, stringsAsFactors = FALSE)
}
