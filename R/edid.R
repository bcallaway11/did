#' Efficient Difference-in-Differences Estimator
#'
#' Estimates group-time average treatment effects \eqn{ATT(g, t)} for staggered
#' adoption designs using the Efficient DiD (EDiD) estimator of Chen, Sant'Anna
#' & Xie (2025). The estimator combines all valid DiD identifying moments for
#' each \eqn{(g, t)} cell with optimal inverse-covariance weights to achieve
#' minimum asymptotic variance.
#'
#' @param data A \code{data.frame}, \code{data.table}, or tibble in long format
#'   (one row per unit-time observation).
#' @param outcome Character scalar: name of the outcome column (must be numeric
#'   with no missing or non-finite values).
#' @param unit Character scalar: name of the unit identifier column.
#' @param time Character scalar: name of the time period column (numeric).
#' @param first_treat Character scalar: name of the column recording each unit's
#'   first treatment period. Never-treated units should have \code{Inf} (or
#'   \code{NA} is not accepted --- use \code{Inf} explicitly).
#' @param covariates Character vector of covariate column names, or \code{NULL}
#'   (default). \strong{Currently not implemented}: passing non-NULL triggers an
#'   error.
#' @param pt_assumption Parallel-trends assumption regime. One of:
#'   \describe{
#'     \item{\code{"all"}}{PT-All: parallel trends holds for all pre-treatment
#'       periods (default). Uses all valid \eqn{(g', t_{pre})} pairs.}
#'     \item{\code{"post"}}{PT-Post: parallel trends holds only for the period
#'       immediately before treatment. Each cell uses a single DiD moment.}
#'   }
#' @param alpha Significance level for confidence intervals. Default \code{0.05}.
#' @param cluster Character scalar naming a time-invariant cluster variable in
#'   \code{data}, or \code{NULL} for no clustering (default). When supplied,
#'   cluster-robust standard errors are computed via the sandwich EIF formula.
#' @param control_group Control group definition. One of:
#'   \describe{
#'     \item{\code{"never_treated"}}{Use never-treated units (default).}
#'     \item{\code{"last_cohort"}}{Use the last-treated cohort as
#'       pseudo-controls (relabeled as never-treated internally).}
#'   }
#' @param n_bootstrap Non-negative integer: number of multiplier bootstrap draws.
#'   \code{0} (default) returns analytical standard errors only.
#' @param bootstrap_weights Distribution for multiplier weights. One of
#'   \code{"rademacher"} (default), \code{"mammen"}, or \code{"webb"}.
#' @param seed Integer seed for reproducibility of bootstrap draws, or
#'   \code{NULL} (default, no seed set).
#' @param anticipation Non-negative integer: number of anticipation periods.
#'   Default \code{0L}. The effective treatment start for cohort \eqn{g} is
#'   \eqn{g - \text{anticipation}}.
#' @param aggregate Which aggregations to compute. One or more of
#'   \code{"all"} (default), \code{"overall"}, \code{"event_study"},
#'   \code{"group"}, or \code{"none"}.
#' @param balance_e Integer or \code{NULL}: if not \code{NULL}, restricts the
#'   event-study aggregation to relative times in
#'   \eqn{[-\text{balance\_e}, \text{balance\_e}]}.
#' @param survey_design Always \code{NULL}. Survey designs are not yet
#'   implemented; passing a non-NULL value triggers an error.
#' @param store_eif Logical: if \code{TRUE}, store the full \eqn{n \times K}
#'   EIF matrix in \code{edid_fit$eif}. Default \code{FALSE}. The EIF is
#'   always computed internally when \code{n_bootstrap > 0}.
#'
#' @return An object of class \code{edid_fit} (a list) with elements:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{att_gt}}{data.frame of cell-level estimates (group, time,
#'       att, se, ci_lower, ci_upper, t_stat, p_value, is_pre).}
#'     \item{\code{overall}}{List: overall ATT with SE and CI.}
#'     \item{\code{event_study}}{List of per-relative-time ATTs.}
#'     \item{\code{group}}{List of per-cohort ATTs.}
#'     \item{\code{eif}}{EIF matrix or \code{NULL}.}
#'     \item{\code{bootstrap}}{Bootstrap results or \code{NULL}.}
#'   }
#'
#' @references Chen, L., Sant'Anna, P. H. C., & Xie, Y. (2025).
#'   \emph{Efficient Difference-in-Differences}. Working paper.
#'
#' @examples
#' # Simulate a simple balanced panel with staggered adoption
#' set.seed(42)
#' n_units <- 100
#' n_periods <- 6
#' unit_ids  <- rep(1:n_units, each = n_periods)
#' time_ids  <- rep(1:n_periods, times = n_units)
#' # Assign cohorts: 1/3 treated in period 3, 1/3 in period 5, 1/3 never
#' cohort_assign <- rep(
#'   c(3, 5, Inf),
#'   times = c(ceiling(n_units / 3),
#'             ceiling(n_units / 3),
#'             n_units - 2 * ceiling(n_units / 3))
#' )[1:n_units]
#' first_treat_vec <- cohort_assign[unit_ids]
#' # Generate outcomes: ATT = 1 for treated post-treatment
#' treat_effect <- as.numeric(time_ids >= first_treat_vec)
#' y_vals <- 0.5 * time_ids + treat_effect + rnorm(n_units * n_periods, sd = 0.5)
#' panel_df <- data.frame(
#'   id          = unit_ids,
#'   period      = time_ids,
#'   y           = y_vals,
#'   first_treat = first_treat_vec
#' )
#' # Fit EDiD (no-covariate, PT-All, analytical SE)
#' fit <- edid(
#'   data        = panel_df,
#'   outcome     = "y",
#'   unit        = "id",
#'   time        = "period",
#'   first_treat = "first_treat",
#'   pt_assumption = "all"
#' )
#' # View overall ATT
#' fit$overall$att
#' # Extract cell-level estimates
#' head(fit$att_gt)
#'
#' @export
edid <- function(
  data,
  outcome,
  unit,
  time,
  first_treat,
  covariates        = NULL,
  pt_assumption     = c("all", "post"),
  alpha             = 0.05,
  cluster           = NULL,
  control_group     = c("never_treated", "last_cohort"),
  n_bootstrap       = 0L,
  bootstrap_weights = c("rademacher", "mammen", "webb"),
  seed              = NULL,
  anticipation      = 0L,
  aggregate         = c("all", "overall", "event_study", "group", "none"),
  balance_e         = NULL,
  survey_design     = NULL,
  store_eif         = FALSE
) {
  mc <- match.call()

  # ------------------------------------------------------------------
  # Argument matching
  # ------------------------------------------------------------------
  pt_assumption     <- match.arg(pt_assumption)
  control_group     <- match.arg(control_group)
  bootstrap_weights <- match.arg(bootstrap_weights)
  aggregate         <- match.arg(aggregate, several.ok = TRUE)
  # When "all" is present it subsumes the others
  if ("all" %in% aggregate) aggregate <- "all"

  n_bootstrap  <- as.integer(n_bootstrap)
  anticipation <- as.integer(anticipation)

  # ------------------------------------------------------------------
  # Validation
  # ------------------------------------------------------------------
  validate_edid_inputs(
    data          = data,
    outcome       = outcome,
    unit          = unit,
    time          = time,
    first_treat   = first_treat,
    covariates    = covariates,
    pt_assumption = pt_assumption,
    alpha         = alpha,
    cluster       = cluster,
    control_group = control_group,
    n_bootstrap   = n_bootstrap,
    anticipation  = anticipation,
    survey_design = survey_design
  )

  # ------------------------------------------------------------------
  # Panel preparation
  # ------------------------------------------------------------------
  panel_obj <- prepare_edid_panel(
    data          = data,
    outcome       = outcome,
    unit          = unit,
    time          = time,
    first_treat   = first_treat,
    covariates    = covariates,
    cluster       = cluster,
    control_group = control_group,
    anticipation  = anticipation
  )

  # ------------------------------------------------------------------
  # Cell estimation
  # EIF is always needed for aggregated SE computation, not just bootstrap.
  # ------------------------------------------------------------------
  do_any_agg   <- !("none" %in% aggregate)
  need_eif_for_boot <- (n_bootstrap > 0L)
  # need_eif: TRUE whenever we need aggregated inference OR bootstrap
  need_eif_internal <- do_any_agg || need_eif_for_boot

  fit_result <- fit_edid_cells(
    panel_obj     = panel_obj,
    pt_assumption = pt_assumption,
    alpha         = alpha,
    store_eif     = store_eif,
    covariates    = covariates,
    need_eif      = need_eif_internal
  )

  cells      <- fit_result$cells
  eif_matrix <- fit_result$eif_matrix
  cell_index <- fit_result$cell_index

  # ------------------------------------------------------------------
  # Convenience att_gt table
  # ------------------------------------------------------------------
  att_gt_df <- data.frame(
    group    = vapply(cells, function(x) x$group,  numeric(1L)),
    time     = vapply(cells, function(x) x$time,   numeric(1L)),
    att      = vapply(cells, function(x) if (is.null(x$att))  NA_real_ else x$att,  numeric(1L)),
    se       = vapply(cells, function(x) if (is.null(x$se))   NA_real_ else x$se,   numeric(1L)),
    ci_lower = vapply(cells, function(x) if (is.null(x$ci_lower)) NA_real_ else x$ci_lower, numeric(1L)),
    ci_upper = vapply(cells, function(x) if (is.null(x$ci_upper)) NA_real_ else x$ci_upper, numeric(1L)),
    t_stat   = vapply(cells, function(x) if (is.null(x$t_stat))  NA_real_ else x$t_stat,  numeric(1L)),
    p_value  = vapply(cells, function(x) if (is.null(x$p_value)) NA_real_ else x$p_value, numeric(1L)),
    n_pairs  = vapply(cells, function(x) if (is.null(x$n_pairs)) 0L else x$n_pairs, integer(1L)),
    is_pre   = vapply(cells, function(x) x$is_pre, logical(1L)),
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------------
  # Aggregation
  # ------------------------------------------------------------------
  do_overall    <- aggregate %in% c("all", "overall")
  do_event_study <- aggregate %in% c("all", "event_study")
  do_group      <- aggregate %in% c("all", "group")

  overall_res    <- NULL
  event_study_res <- NULL
  group_res      <- NULL

  if (do_overall) {
    overall_res <- aggregate_overall_edid(cells, eif_matrix, cell_index, panel_obj, alpha)
  }
  if (do_event_study) {
    event_study_res <- aggregate_event_study_edid(
      cells, eif_matrix, cell_index, panel_obj, alpha, balance_e
    )
  }
  if (do_group) {
    group_res <- aggregate_group_edid(cells, eif_matrix, cell_index, panel_obj, alpha)
  }

  # ------------------------------------------------------------------
  # Bootstrap
  # ------------------------------------------------------------------
  bootstrap_res <- NULL
  if (n_bootstrap > 0L) {
    if (is.null(eif_matrix)) {
      warning("EIF matrix is NULL; bootstrap cannot be run. ",
              "This should not happen --- please report this issue.")
    } else {
      boot_agg <- if ("all" %in% aggregate || identical(aggregate, "all")) "all" else
                  paste(intersect(aggregate, c("overall", "event_study", "group")), collapse = ",")
      bootstrap_res <- run_multiplier_bootstrap_edid(
        cells             = cells,
        eif_matrix        = eif_matrix,
        cell_index        = cell_index,
        panel_obj         = panel_obj,
        n_bootstrap       = n_bootstrap,
        bootstrap_weights = bootstrap_weights,
        seed              = seed,
        aggregate         = "all",
        balance_e         = balance_e,
        alpha             = alpha
      )
      class(bootstrap_res) <- c("edid_bootstrap", "list")

      # Overwrite SEs/CIs with bootstrap versions
      if (do_overall && !is.null(overall_res) && !is.null(bootstrap_res$overall_draws)) {
        bs_ov <- compute_bootstrap_stats_edid(bootstrap_res$overall_draws, overall_res$att, alpha)
        overall_res$se       <- bs_ov$se_boot
        overall_res$ci_lower <- bs_ov$ci_lower
        overall_res$ci_upper <- bs_ov$ci_upper
        overall_res$p_value  <- bs_ov$p_value_boot
        overall_res$t_stat   <- if (!is.na(bs_ov$se_boot) && bs_ov$se_boot > 0) {
          overall_res$att / bs_ov$se_boot
        } else NA_real_
      }

      if (do_event_study && !is.null(event_study_res) &&
          !is.null(bootstrap_res$event_study_draws)) {
        for (e_nm in names(event_study_res)) {
          draws <- bootstrap_res$event_study_draws[[e_nm]]
          if (is.null(draws)) next
          bs_es <- compute_bootstrap_stats_edid(draws, event_study_res[[e_nm]]$att, alpha)
          event_study_res[[e_nm]]$se       <- bs_es$se_boot
          event_study_res[[e_nm]]$ci_lower <- bs_es$ci_lower
          event_study_res[[e_nm]]$ci_upper <- bs_es$ci_upper
          event_study_res[[e_nm]]$p_value  <- bs_es$p_value_boot
        }
      }

      if (do_group && !is.null(group_res) && !is.null(bootstrap_res$group_draws)) {
        for (g_nm in names(group_res)) {
          draws <- bootstrap_res$group_draws[[g_nm]]
          if (is.null(draws)) next
          bs_gr <- compute_bootstrap_stats_edid(draws, group_res[[g_nm]]$att, alpha)
          group_res[[g_nm]]$se       <- bs_gr$se_boot
          group_res[[g_nm]]$ci_lower <- bs_gr$ci_lower
          group_res[[g_nm]]$ci_upper <- bs_gr$ci_upper
          group_res[[g_nm]]$p_value  <- bs_gr$p_value_boot
        }
        # Also update cell-level SEs from bootstrap if we have per-cell draws
        # (not stored at cell level -- only aggregate-level bootstrap is implemented)
      }
    }
  }

  # ------------------------------------------------------------------
  # EIF matrix storage (only if user requested it)
  # ------------------------------------------------------------------
  eif_export <- if (store_eif) eif_matrix else NULL

  # ------------------------------------------------------------------
  # Construct edid_fit S3 object
  # ------------------------------------------------------------------
  edid_fit <- list(
    call             = mc,
    pt_assumption    = pt_assumption,
    control_group    = control_group,
    alpha            = alpha,
    n                = panel_obj$n,
    T_periods        = panel_obj$T_periods,
    treatment_groups = panel_obj$treatment_groups,
    anticipation     = panel_obj$anticipation,
    inference_type   = if (n_bootstrap > 0L) "bootstrap" else "analytical",
    cluster          = cluster,
    cells            = cells,
    att_gt           = att_gt_df,
    overall          = overall_res,
    event_study      = event_study_res,
    group            = group_res,
    eif              = eif_export,
    bootstrap        = bootstrap_res,
    n_bootstrap      = n_bootstrap,
    bootstrap_weights = bootstrap_weights
  )

  class(edid_fit) <- c("edid_fit", "list")
  edid_fit
}
