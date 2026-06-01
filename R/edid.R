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
#' @param yname Character scalar: name of the outcome column (must be numeric
#'   with no missing or non-finite values).
#' @param idname Character scalar: name of the unit identifier column.
#' @param tname Character scalar: name of the time period column (numeric).
#' @param gname Character scalar: name of the column recording each unit's
#'   first treatment period. Never-treated units should have \code{Inf} or
#'   \code{0} (the \code{att_gt()} convention). \code{0} is automatically
#'   converted to \code{Inf} internally.
#' @param xformla A one-sided formula specifying covariates to condition on,
#'   e.g., \code{~ X1 + X2}. Default \code{NULL} (equivalent to \code{~1},
#'   no covariates). When \code{NULL} or \code{~1}, the efficient no-covariate
#'   path is used. \strong{Note}: The \code{covariates} argument is deprecated
#'   and will error if non-NULL; use \code{xformla} instead.
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
#' @param alp Significance level for confidence intervals. Default \code{0.05}.
#' @param clustervars Character scalar naming a time-invariant cluster variable
#'   in \code{data}, or \code{NULL} for no clustering (default). When supplied,
#'   cluster-robust standard errors are computed via the sandwich EIF formula.
#'   Note: edid() currently supports only a single cluster variable internally.
#' @param control_group Control group definition. One of:
#'   \describe{
#'     \item{\code{"nevertreated"}}{Use never-treated units (default).}
#'     \item{\code{"notyettreated"}}{Use the last-treated cohort as
#'       pseudo-controls (relabeled as never-treated internally).}
#'   }
#' @param bstrap Logical: whether to use multiplier bootstrap inference.
#'   Default \code{FALSE} (analytical standard errors). When \code{TRUE},
#'   \code{biters} bootstrap draws are used.
#' @param biters Positive integer: number of multiplier bootstrap iterations.
#'   Default \code{1000L}. Only used when \code{bstrap = TRUE}.
#' @param bootstrap_weights Distribution for multiplier weights. One of
#'   \code{"rademacher"} (default), \code{"mammen"}, or \code{"webb"}.
#' @param seed Integer seed for reproducibility of bootstrap draws, or
#'   \code{NULL} (default, no seed set).
#' @param anticipation Non-negative integer: number of anticipation periods.
#'   Default \code{0L}. The effective treatment start for cohort \eqn{g} is
#'   \eqn{g - \text{anticipation}}.
#' @param aggregate Which aggregations to compute. One or more of
#'   \code{"all"} (default), \code{"overall"}, \code{"event_study"},
#'   \code{"group"}, \code{"calendar"}, or \code{"none"}. \code{"event_study"}
#'   reports the cohort-share-weighted event-study parameters \eqn{ES(e)};
#'   \code{"group"} averages \eqn{ATT(g,t)} within each cohort; \code{"calendar"}
#'   averages \eqn{ATT(g,t)} across the cohorts treated by each calendar period;
#'   \code{"overall"} returns the simple cohort-share aggregate over all
#'   post-treatment cells. \code{"all"} computes every aggregation.
#' @param balance_e Integer or \code{NULL}: if not \code{NULL}, restricts the
#'   event-study aggregation to relative times in
#'   \eqn{[-\text{balance\_e}, \text{balance\_e}]}.
#' @param survey_design Always \code{NULL}. Survey designs are not yet
#'   implemented; passing a non-NULL value triggers an error.
#' @param store_eif Logical: if \code{TRUE}, store the full \eqn{n \times K}
#'   EIF matrix in \code{edid_fit$eif}. Default \code{FALSE}. The EIF is
#'   always computed internally when \code{bstrap = TRUE}.
#' @param weights How the per-pair generated-outcome moments are combined in the
#'   covariate path. \code{"efficient"} (default) uses the semiparametric-efficient
#'   pointwise weights \eqn{w(X_i)=\Omega^*(X_i)^{-1}\mathbf 1/(\mathbf 1'\Omega^*(X_i)^{-1}\mathbf 1)},
#'   estimated by kernel and stabilized by two finite-sample regularizations: data-driven shrinkage of
#'   \eqn{\hat\Omega^*(X_i)} toward the pooled \eqn{\bar\Omega^*} (intensity \eqn{\hat\lambda\to0}) and a
#'   relative eigenvalue floor that vanishes with the sample size. Both are asymptotically inactive, so this
#'   feasible estimator is asymptotically equivalent to the efficient estimator and attains the efficiency
#'   bound in the limit (it is not exactly bound-attaining in finite samples).
#'   The constant-weight alternatives remain consistent for \eqn{ATT(g,t)} (any weights summing to one
#'   identify the estimand, with no rate condition on the weights) but do not attain the bound:
#'   \code{"averaged"} inverts the covariate-averaged conditional covariance \eqn{\bar\Omega^*};
#'   \code{"gmm"} inverts the unconditional moment covariance \eqn{\hat S}; \code{"uniform"} assigns
#'   equal weight \eqn{1/H} to the \eqn{H} non-collinear moments.
#'
#' @return An object of class \code{edid_fit} (a list) with elements:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{att_gt}}{data.frame of cell-level estimates (group, time,
#'       att, se, ci_lower, ci_upper, t_stat, p_value, is_pre).}
#'     \item{\code{overall}}{List: the HEADLINE overall ATT, with SE and CI. By default this is the
#'       DYNAMIC event-study average over relative times \eqn{e \ge 0} (the paper's main object,
#'       \code{= aggte_edid(type = "dynamic")$overall.att}); it falls back to the simple aggregate
#'       when no event study is computed.}
#'     \item{\code{overall_simple}}{List: the did-"simple" overall ATT (cohort-share-weighted average
#'       over all post cells). Equals \code{aggte_edid(type = "simple")$overall.att}.}
#'     \item{\code{overall_group}}{List: the cohort-share-weighted group overall ATT. Present when
#'       \code{group} is computed. Equals \code{aggte_edid(type = "group")$overall.att}.}
#'     \item{\code{overall_calendar}}{List: the calendar-time overall ATT (equal-weight average over
#'       post-treatment calendar periods), with SE and CI. Present when \code{calendar} is computed.
#'       Equals \code{aggte_edid(type = "calendar")$overall.att}.}
#'     \item{\code{event_study}}{List of per-relative-time event-study estimates \eqn{ES(e)}, one per
#'       event time \eqn{e}.}
#'     \item{\code{group}}{List of per-cohort overall ATTs (within-cohort averages of \eqn{ATT(g,t)}).}
#'     \item{\code{calendar}}{List of per-calendar-period averages of \eqn{ATT(g,t)} across the cohorts
#'       treated by that period, or \code{NULL} when not requested.}
#'     \item{\code{eif}}{The \eqn{n \times K} efficient-influence-function matrix, or \code{NULL} when
#'       \code{store_eif = FALSE}.}
#'     \item{\code{bootstrap}}{Bootstrap results or \code{NULL}.}
#'     \item{\code{bstrap}}{Logical: whether bootstrap inference was used.}
#'   }
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025).
#'   \emph{Efficient Difference-in-Differences and Event Study Estimators}.
#'   Working paper.
#'
#' @keywords models
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
#'   data          = panel_df,
#'   yname         = "y",
#'   idname        = "id",
#'   tname         = "period",
#'   gname         = "first_treat",
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
  yname,
  idname,
  tname,
  gname,
  xformla           = NULL,
  covariates        = NULL,
  pt_assumption     = c("all", "post"),
  alp               = 0.05,
  clustervars       = NULL,
  control_group     = c("nevertreated", "notyettreated"),
  bstrap            = FALSE,
  biters            = 1000L,
  bootstrap_weights = c("rademacher", "mammen", "webb"),
  seed              = NULL,
  anticipation      = 0L,
  aggregate         = c("all", "overall", "event_study", "group", "calendar", "none"),
  balance_e         = NULL,
  survey_design     = NULL,
  store_eif         = FALSE,
  weights           = c("efficient", "averaged", "gmm", "uniform")
) {
  weight_method <- match.arg(weights)
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

  anticipation <- as.integer(anticipation)

  # ------------------------------------------------------------------
  # Bootstrap: derive internal n_bootstrap from bstrap + biters
  # ------------------------------------------------------------------
  n_bootstrap_internal <- if (bstrap) as.integer(biters) else 0L

  # ------------------------------------------------------------------
  # Accept G=0 (att_gt convention) or G=Inf (edid native) for never-treated
  # Convert 0 -> Inf internally, matching att_gt's internal transformation
  data <- as.data.frame(data)
  # Only the numeric att_gt convention uses G=0 for never-treated. Guard with
  # is.numeric so a factor/character gname is NOT silently coerced to integer codes
  # here (which would relabel cohorts); it reaches validate_edid_inputs and errors.
  if (is.numeric(data[[gname]])) {
    zero_nt <- is.finite(data[[gname]]) & data[[gname]] == 0
    if (any(zero_nt)) {
      data[[gname]] <- ifelse(zero_nt, Inf, data[[gname]])
    }
  }

  # ------------------------------------------------------------------
  # Validation
  # ------------------------------------------------------------------
  validate_edid_inputs(
    data          = data,
    yname         = yname,
    idname        = idname,
    tname         = tname,
    gname         = gname,
    xformla       = xformla,
    covariates    = covariates,
    pt_assumption = pt_assumption,
    alp           = alp,
    clustervars   = clustervars,
    control_group = control_group,
    biters        = n_bootstrap_internal,
    anticipation  = anticipation,
    survey_design = survey_design
  )

  # ------------------------------------------------------------------
  # Panel preparation
  # ------------------------------------------------------------------
  panel_obj <- prepare_edid_panel(
    data          = data,
    yname         = yname,
    idname        = idname,
    tname         = tname,
    gname         = gname,
    xformla       = xformla,
    covariates    = covariates,
    clustervars   = clustervars,
    control_group = control_group,
    anticipation  = anticipation
  )

  # ------------------------------------------------------------------
  # Cell estimation
  # EIF is always needed for aggregated SE computation, not just bootstrap.
  # ------------------------------------------------------------------
  do_any_agg   <- !("none" %in% aggregate)
  need_eif_for_boot <- (n_bootstrap_internal > 0L)
  # need_eif: TRUE whenever we need aggregated inference OR bootstrap
  need_eif_internal <- do_any_agg || need_eif_for_boot

  fit_result <- fit_edid_cells(
    panel_obj     = panel_obj,
    pt_assumption = pt_assumption,
    alpha         = alp,
    store_eif     = store_eif,
    xformla       = xformla,
    need_eif      = need_eif_internal,
    seed          = seed,
    weight_method = weight_method
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
  do_calendar   <- aggregate %in% c("all", "calendar")

  overall_res    <- NULL
  event_study_res <- NULL
  group_res      <- NULL
  calendar_res   <- NULL

  if (do_overall) {
    overall_res <- aggregate_overall_edid(cells, eif_matrix, cell_index, panel_obj, alp)
  }
  if (do_event_study) {
    event_study_res <- aggregate_event_study_edid(
      cells, eif_matrix, cell_index, panel_obj, alp, balance_e
    )
  }
  if (do_group) {
    group_res <- aggregate_group_edid(cells, eif_matrix, cell_index, panel_obj, alp)
  }
  if (do_calendar) {
    calendar_res <- aggregate_calendar_edid(cells, eif_matrix, cell_index, panel_obj, alp)
  }

  # ------------------------------------------------------------------
  # Bootstrap
  # ------------------------------------------------------------------
  bootstrap_res <- NULL
  if (n_bootstrap_internal > 0L) {
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
        n_bootstrap       = n_bootstrap_internal,
        bootstrap_weights = bootstrap_weights,
        seed              = seed,
        aggregate         = "all",
        balance_e         = balance_e,
        alpha             = alp
      )
      class(bootstrap_res) <- c("edid_bootstrap", "list")

      # Overwrite SEs/CIs with bootstrap versions
      if (do_overall && !is.null(overall_res) && !is.null(bootstrap_res$overall_b)) {
        bs_ov <- compute_bootstrap_stats_edid(bootstrap_res$overall_b, overall_res$att, alp)
        overall_res$se       <- bs_ov$se_boot
        overall_res$ci_lower <- bs_ov$ci_lower
        overall_res$ci_upper <- bs_ov$ci_upper
        overall_res$p_value  <- bs_ov$p_value_boot
        overall_res$t_stat   <- if (!is.na(bs_ov$se_boot) && bs_ov$se_boot > 0) {
          overall_res$att / bs_ov$se_boot
        } else NA_real_
      }

      if (do_event_study && !is.null(event_study_res) &&
          !is.null(bootstrap_res$event_study_b)) {
        for (e_nm in names(event_study_res)) {
          draws <- bootstrap_res$event_study_b[[e_nm]]
          if (is.null(draws)) next
          bs_es <- compute_bootstrap_stats_edid(draws, event_study_res[[e_nm]]$att, alp)
          event_study_res[[e_nm]]$se       <- bs_es$se_boot
          event_study_res[[e_nm]]$ci_lower <- bs_es$ci_lower
          event_study_res[[e_nm]]$ci_upper <- bs_es$ci_upper
          event_study_res[[e_nm]]$p_value  <- bs_es$p_value_boot
        }
      }

      if (do_group && !is.null(group_res) && !is.null(bootstrap_res$group_b)) {
        for (g_nm in names(group_res)) {
          draws <- bootstrap_res$group_b[[g_nm]]
          if (is.null(draws)) next
          bs_gr <- compute_bootstrap_stats_edid(draws, group_res[[g_nm]]$att, alp)
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
    alpha            = alp,
    n                = panel_obj$n,
    T_periods        = panel_obj$T_periods,
    treatment_groups = panel_obj$treatment_groups,
    cohort_fractions = panel_obj$cohort_fractions,
    unit_cohorts     = panel_obj$unit_cohorts,
    anticipation     = panel_obj$anticipation,
    inference_type   = if (n_bootstrap_internal > 0L) "bootstrap" else "analytical",
    clustervars      = clustervars,
    cluster_indices  = panel_obj$cluster_indices,  # for cluster-robust re-aggregation in aggte_edid
    xformla          = xformla,
    bstrap           = bstrap,
    cells            = cells,
    att_gt           = att_gt_df,
    overall          = overall_res,
    event_study      = event_study_res,
    group            = group_res,
    calendar         = calendar_res,
    eif              = eif_export,
    bootstrap        = bootstrap_res,
    n_bootstrap      = n_bootstrap_internal,
    bootstrap_weights = bootstrap_weights
  )

  class(edid_fit) <- c("edid_fit", "list")

  # Headline `$overall` = the DYNAMIC event-study average over relative times e >= 0 -- the paper's
  # main object (ES_avg). The did-"simple" cohort-share aggregate over all post cells is kept as
  # `$overall_simple`, and `$overall_group` is the cohort-share-weighted group average. Each matches
  # the corresponding aggte_edid(type=)$overall.att (computed via aggte_edid -- single source of
  # truth, so the WIF + cluster-robust SE flow through). The per-relative-time `$event_study` and
  # per-cohort `$group` lists are UNCHANGED. When no event study is computed, `$overall` falls back
  # to the simple aggregate.
  edid_fit$overall_simple <- overall_res
  zq <- stats::qnorm(1 - alp / 2)
  to_ov <- function(a) {
    att <- a$overall.att; se <- a$overall.se
    tstat <- if (is.finite(se) && se > 0) att / se else NA_real_
    list(att = att, se = se, ci_lower = att - zq * se, ci_upper = att + zq * se,
         t_stat = tstat,
         p_value = if (is.finite(tstat)) 2 * stats::pnorm(-abs(tstat)) else NA_real_,
         type = a$type)
  }
  if (do_event_study && !is.null(event_study_res) && length(event_study_res) > 0L) {
    dyn <- tryCatch(aggte_edid(edid_fit, type = "dynamic", na.rm = TRUE), error = function(e) NULL)
    if (!is.null(dyn)) edid_fit$overall <- to_ov(dyn)   # headline = dynamic ES_avg (e >= 0)
  }
  if (do_group && !is.null(group_res) && length(group_res) > 0L) {
    grp <- tryCatch(aggte_edid(edid_fit, type = "group"), error = function(e) NULL)
    if (!is.null(grp)) edid_fit$overall_group <- to_ov(grp)
  }
  if (do_calendar && !is.null(calendar_res) && length(calendar_res) > 0L) {
    cal <- tryCatch(aggte_edid(edid_fit, type = "calendar", na.rm = TRUE), error = function(e) NULL)
    if (!is.null(cal)) edid_fit$overall_calendar <- to_ov(cal)
  }

  edid_fit
}
