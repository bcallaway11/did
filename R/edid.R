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
#' @param correct_first_step Logical (default \code{FALSE}). If \code{TRUE}, the influence function
#'   is augmented with the first-step nuisance-estimation correction of Ackerberg, Chen and Hahn (2012)
#'   for the sieve nuisances (conditional means and propensity ratios) entering the doubly-robust moment.
#'   The influence-function moments are Neyman orthogonal, so this correction is asymptotically negligible
#'   under correct specification (it leaves the variance bound unchanged in the limit); it provides
#'   finite-sample robustness when a first-step nuisance is misspecified, where the doubly-robust point
#'   estimate remains consistent. It is a practical (numerical-derivative) form of the two-step variance
#'   estimator and is supported only for the default plug-in nuisances (covariate path).
#'   \strong{Scope:} the correction is for the sieve nuisances (m, r) that enter the generated outcome,
#'   computed with the estimated efficient weights held FIXED. It does \emph{not} correct the
#'   weight-estimation channel (the kernel \eqn{\Omega^*}, its Ledoit-Wolf shrinkage, and the eigenvalue
#'   floor that map to \eqn{w(X)}); that channel is asymptotically negligible separately but is not part of
#'   this correction. With the rich default sieve the correction is empirically small; its value is
#'   robustness when a nuisance is genuinely misspecified.
#'
#' @references Ackerberg, D., Chen, X., and Hahn, J. (2012). A Practical Asymptotic Variance Estimator
#'   for Two-Step Semiparametric Estimators. \emph{Review of Economics and Statistics}, 94(2), 481-498.
#'
#' @return An object of class \code{edid_fit} (a list) with elements:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{att_gt}}{data.frame of cell-level estimates (group, time,
#'       att, se, ci_lower, ci_upper, t_stat, p_value, is_pre).}
#'     \item{\code{overall}}{A \code{did::AGGTEobj}: the HEADLINE aggregation -- the dynamic event-study
#'       average over relative times \eqn{e \ge 0} (the paper's main object) when an event study is
#'       requested, otherwise the cohort-share "simple" aggregate.}
#'     \item{\code{simple}}{A \code{did::AGGTEobj} for the cohort-share-weighted average over all
#'       post-treatment cells (\code{= aggte_edid(type = "simple")}); present when \code{overall}/\code{all}
#'       is requested.}
#'     \item{\code{event_study}}{A \code{did::AGGTEobj} for the event study \eqn{ES(e)}: per relative time
#'       (\code{att.egt}/\code{egt}) plus the dynamic overall.}
#'     \item{\code{group}}{A \code{did::AGGTEobj} for the per-cohort overall ATTs.}
#'     \item{\code{calendar}}{A \code{did::AGGTEobj} for the per-calendar-period averages of
#'       \eqn{ATT(g,t)}, or \code{NULL} when not requested.}
#'     \item{\code{eif}}{The \eqn{n \times K} efficient-influence-function matrix (always stored).}
#'     \item{\code{bstrap}}{Logical: whether bootstrap inference was used. Under \code{bstrap = TRUE} the
#'       cell SEs and the aggregations use the did multiplier bootstrap (\code{\link[did]{mboot}} /
#'       \code{\link[did]{aggte}}).}
#'   }
#'   The aggregation slots are standard \code{did::AGGTEobj} objects, so \code{summary}, \code{tidy}, and
#'   \code{ggdid} work on them directly.
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
  seed              = NULL,
  anticipation      = 0L,
  aggregate         = c("all", "overall", "event_study", "group", "calendar", "none"),
  balance_e         = NULL,
  survey_design     = NULL,
  weights           = c("efficient", "averaged", "gmm", "uniform"),
  correct_first_step = FALSE
) {
  weight_method <- match.arg(weights)
  mc <- match.call()

  # ------------------------------------------------------------------
  # Argument matching
  # ------------------------------------------------------------------
  pt_assumption     <- match.arg(pt_assumption)
  control_group     <- match.arg(control_group)
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
    store_eif     = TRUE,                # edid always retains the EIF (used by the aggregations)
    xformla       = xformla,
    need_eif      = need_eif_internal,
    seed          = seed,
    weight_method = weight_method,
    correct_first_step = isTRUE(correct_first_step)
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

  # Aggregations are computed below, AFTER the edid_fit object exists, via aggte_edid() -- i.e. through
  # did::aggte() on the edid MP. So $overall/$event_study/$group/$calendar/$simple are standard
  # did::AGGTEobj objects and inherit did's print/summary/tidy methods (full att_gt-compatibility).
  overall_res <- event_study_res <- group_res <- calendar_res <- NULL

  # ------------------------------------------------------------------
  # Bootstrap. The multiplier bootstrap runs through did::mboot on the cell influence functions for the
  # cell-level SEs + simultaneous critical value; the aggregations bootstrap through aggte_edid() ->
  # did::aggte(bstrap = TRUE) below (so $overall/$event_study/... carry bootstrap SEs and uniform bands).
  # Reproducible via `seed`.
  # ------------------------------------------------------------------
  if (bstrap) {
    if (!is.null(seed)) set.seed(seed)
    bdp <- list(idname = idname, tname = tname, clustervars = clustervars,
                biters = as.integer(biters), alp = alp, panel = TRUE, faster_mode = FALSE,
                true_repeated_cross_sections = FALSE, allow_unbalanced_panel = FALSE)
    if (!is.null(clustervars)) {
      # time-invariant cluster data, one row per unit in the influence-function (all_units) order
      bdp$data <- stats::setNames(
        data.frame(panel_obj$all_units, min(panel_obj$time_periods), panel_obj$cluster_indices),
        c(idname, tname, clustervars))
    }
    bb <- mboot(eif_matrix, bdp, pl = FALSE)
    ok <- is.finite(bb$se)
    att_gt_df$se[ok]       <- bb$se[ok]
    att_gt_df$ci_lower[ok] <- att_gt_df$att[ok] - bb$crit.val * bb$se[ok]
    att_gt_df$ci_upper[ok] <- att_gt_df$att[ok] + bb$crit.val * bb$se[ok]
    att_gt_df$t_stat[ok]   <- att_gt_df$att[ok] / bb$se[ok]
    att_gt_df$p_value[ok]  <- 2 * stats::pnorm(-abs(att_gt_df$t_stat[ok]))
  }

  # ------------------------------------------------------------------
  # EIF matrix storage. edid always retains the influence functions (as att_gt() always returns
  # $inffunc): aggte_edid()/as_MP_edid() build the did MP from them.
  # ------------------------------------------------------------------
  eif_export <- eif_matrix

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
    all_units        = panel_obj$all_units,        # metadata for did-compatible MP construction (as_MP_edid)
    idname           = idname,
    tname            = tname,
    gname            = gname,
    time_periods     = panel_obj$time_periods,
    panel            = TRUE,
    anticipation     = panel_obj$anticipation,
    inference_type   = if (n_bootstrap_internal > 0L) "bootstrap" else "analytical",
    correct_first_step = isTRUE(correct_first_step),
    clustervars      = clustervars,
    cluster_indices  = panel_obj$cluster_indices,  # for cluster-robust re-aggregation in aggte_edid
    xformla          = xformla,
    bstrap           = bstrap,
    biters           = as.integer(biters),         # used by aggte_edid()/as_MP_edid() for bstrap = TRUE
    cells            = cells,
    att_gt           = att_gt_df,
    overall          = overall_res,
    event_study      = event_study_res,
    group            = group_res,
    calendar         = calendar_res,
    eif              = eif_export
  )

  class(edid_fit) <- c("edid_fit", "list")

  # Compute the requested aggregations as did::AGGTEobj objects via aggte_edid() (= did::aggte on the
  # edid MP), so they inherit did's print/summary/tidy methods (full att_gt-compatibility). The dynamic
  # AGGTEobj carries BOTH the headline event-study average (overall.att) and the per-relative-time ES(e)
  # (att.egt / egt); pre-treatment leads are dropped via na.rm. `$overall` is the headline AGGTEobj:
  # the dynamic event-study average when an event study is requested, else the cohort-share "simple"
  # aggregate. `$event_study`, `$group`, `$calendar`, `$simple` are the corresponding AGGTEobj objects.
  .agg <- function(ty) tryCatch(aggte_edid(edid_fit, type = ty, na.rm = TRUE), error = function(e) NULL)
  if (do_event_study) edid_fit$event_study <- .agg("dynamic")
  if (do_group)       edid_fit$group       <- .agg("group")
  if (do_calendar)    edid_fit$calendar    <- .agg("calendar")
  if (do_overall)     edid_fit$simple      <- .agg("simple")
  edid_fit$overall <- if (!is.null(edid_fit$event_study)) edid_fit$event_study else edid_fit$simple

  edid_fit
}
