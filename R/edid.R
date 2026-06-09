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
#' @param bstrap Logical: whether to use multiplier bootstrap inference.
#'   Default \code{FALSE} (analytical standard errors). When \code{TRUE},
#'   \code{biters} bootstrap draws are used.
#' @param biters Positive integer: number of multiplier bootstrap iterations.
#'   Default \code{1000L}. Only used when \code{bstrap = TRUE}.
#' @param cband Logical: whether to report simultaneous (uniform) confidence bands across the cells and the
#'   event-study / group coefficients. Default \code{TRUE}; \code{FALSE} gives pointwise bands.
#' @param cband_method Character: how the simultaneous critical value is computed. \code{"analytic"}
#'   (default) is the Montiel Olea & Plagborg-Moller sup-t critical value from the analytic, cluster-robust
#'   coefficient covariance -- no bootstrap needed, and the only method compatible with \code{higher_order}.
#'   \code{"multiplier"} uses the did multiplier bootstrap (\code{\link[did]{mboot}}) when
#'   \code{bstrap = TRUE} and reproduces the prior behavior exactly. With very few clusters or very small
#'   samples the multiplier bootstrap can be the safer choice.
#' @param higher_order Logical (default \code{FALSE}). If \code{TRUE}, adds the higher-order ("Wick")
#'   nuisance-estimation variance refinement: the degenerate second-order U-statistic contribution from
#'   estimating the first-step sieve nuisances (\eqn{m}, \eqn{r}) is added to the analytic coefficient
#'   covariance, so BOTH the reported cell standard errors (\eqn{\sqrt{\mathrm{diag}(\Sigma_1 +
#'   \Sigma_{quad})}}) and the sup-t critical value come from the same higher-order-aware covariance.
#'   Because \eqn{\Sigma_{quad}} is positive semi-definite, the cell SEs are never below the plug-in SEs.
#'   The refinement requires the analytic sup-t path (a degenerate-U term cannot be carried by the
#'   multiplier bootstrap, so \code{cband_method = "multiplier"} is coerced to \code{"analytic"} with a
#'   warning) and a covariate formula (with no covariates the nuisances have no first-step coefficients and
#'   the term is exactly zero, so \code{xformla = NULL} errors). It is asymptotically negligible under
#'   correct specification; its value is finite-sample honesty in covariate-rich designs.
#' @param misspec_robust Logical (default \code{TRUE}). Master switch for misspecification-robust standard
#'   errors: when \code{TRUE}, the reported SE accounts for \emph{every} applicable estimation effect --- the
#'   weight-estimation channel (described below), the first-step nuisance ACH correction
#'   (\code{estimation_effect}), and the higher-order ("Wick") nuisance term (\code{higher_order}) --- each
#'   enabled only where it applies and silently skipped where it does not (no covariates; multiplier path;
#'   \code{weight_scheme = "uniform"}), so default calls do not warn. An explicitly-set \code{estimation_effect}
#'   or \code{higher_order} overrides that piece, and \code{misspec_robust = FALSE} reverts to the plug-in
#'   efficient-IF SE. The weight-estimation channel \eqn{\psi_\Omega} is the first-step estimation effect of the
#'   efficient weights \eqn{w(X) = \Omega^{-1}\mathbf{1}/(\mathbf{1}'\Omega^{-1}\mathbf{1})} (the sibling of
#'   \code{estimation_effect}'s nuisance correction that it explicitly leaves out). It yields standard errors
#'   that are robust to misspecification of the weighting model: under correct specification \eqn{\psi_\Omega}
#'   is first-order zero (it vanishes at the \eqn{\sqrt{n}} rate, so the SE converges to the efficient SE),
#'   while under misspecification it accounts for the resulting estimand drift to the weighted pseudo-true
#'   \eqn{\theta_w}. Because \eqn{\psi_\Omega} is a genuine per-unit influence function it is folded into the
#'   EIF, so the cell SEs, every aggregation, the clustered covariance, and the sup-t bands all inherit it. The
#'   reported variance is that of the augmented influence function \eqn{\mathrm{Var}(\mathrm{eif} + \psi_\Omega)}:
#'   it equals the plug-in variance under correct specification (\eqn{\psi_\Omega \to 0}) and corrects it under
#'   misspecification -- which may move a standard error up \emph{or} down, since the plug-in SE is then
#'   inconsistent (unlike \code{higher_order}, whose positive semi-definite \eqn{\Sigma_{quad}} only inflates).
#'   It composes additively with \code{estimation_effect} (which corrects the
#'   nuisance channel) and \code{higher_order}, and -- unlike \code{higher_order} -- it does \emph{not} coerce
#'   \code{cband_method} (a real influence function is carried by the multiplier bootstrap). Supported for the
#'   covariate path with \code{weight_scheme} in \code{c("efficient", "averaged", "gmm")} and plug-in nuisances (cross-fitted
#'   nuisances, \code{K > 1}, are not supported and error). For \code{"gmm"} the weight inverts the unconditional
#'   sample covariance \eqn{C}, a second moment that (unlike the linear ATT moment) is not protected by Neyman
#'   orthogonality, so the channel includes an Ackerberg-Chen-Hahn correction for the first-step (\eqn{r}, \eqn{m})
#'   nuisance estimation entering \eqn{C}. It is \emph{not} available for \code{"uniform"} (fixed weights have no
#'   estimation channel; warns and falls back to the plug-in SE). In weak-overlap / small-\eqn{n} cells where the
#'   efficient pointwise shrinkage \eqn{\lambda} is non-negligible the efficient channel omits an \eqn{O(\lambda)}
#'   term and the SE there is a (warned) approximation to the weight-channel variance.
#' @param trim_level Numeric (default \code{200}). Overlap-trimming threshold (covariate path only). A
#'   comparison observation is dropped from an \eqn{(g,t)} cell's moment and efficient weights when its
#'   estimated propensity ratio \eqn{r(X)} OR inverse propensity \eqn{1/p_{g'}(X)} has absolute value
#'   \eqn{\ge} \code{trim_level} -- mirroring DRDID's \code{trim.level = 0.995} (a control IPW-weight cap of
#'   \eqn{\approx 200}). The observation still contributes to nuisance estimation; only its outcome-side
#'   weight is zeroed. This guards against severe lack of overlap and redefines the target to the overlap
#'   sub-population (as in DRDID). \code{trim_level = Inf} disables trimming. No effect on the no-covariate path.
#' @param seed Integer seed for reproducibility of the bootstrap draws / the analytic sup-t simulation, or
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
#' @param weight_scheme How the per-pair generated-outcome moments are combined in the
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
#' @param estimation_effect Logical (default \code{FALSE}). If \code{TRUE}, the influence function
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
#' @section Advanced options (set via \code{options()}):
#' These global options expose escape hatches and tuning knobs for the covariate path. All have safe
#' defaults; they are intended for diagnostics, reproducibility studies, and large-\eqn{n} scaling. Except
#' where noted, they change only the reported standard errors / bands, not the point estimate \eqn{ATT(g,t)}.
#' \describe{
#'   \item{\code{edid_mc_cores}}{Integer (default \code{1L}). Number of forked workers for the
#'     embarrassingly-parallel \eqn{(g,t)} cell loop and the per-cohort nuisance prebuild, via
#'     \code{\link[parallel]{mclapply}}. A value \code{> 1} gives a wall-clock speed-up on multi-core
#'     machines and is numerically \emph{identical} to the serial path (the cells are independent). It is
#'     fork-based, so it has no effect on Windows (leave at \code{1L}); peak memory grows roughly linearly
#'     in the number of workers.}
#'   \item{\code{edid_omega_method}}{How the conditional covariance \eqn{\Omega^*(X)} is built.
#'     \code{"kernel"} (default) is the fast BLAS Nadaraya-Watson build; \code{"kernel_orig"} is the exact
#'     original per-pair build (a reference that agrees with \code{"kernel"} to roughly \code{1e-13});
#'     \code{"sieve"} is an \eqn{O(np)} series build that avoids the \eqn{n \times n} kernel matrix and so
#'     scales past its memory wall at large \eqn{n}. \strong{Note:} the sieve uses a different smoother (it
#'     changes the point estimate slightly) and is not wired for the \code{misspec_robust} weight-estimation
#'     channel, so pair it with \code{misspec_robust = FALSE}.}
#'   \item{\code{edid_pd_blend}}{Logical (default \code{FALSE}). When \code{TRUE}, a per-unit
#'     \eqn{\Omega^*(X_i)} that is genuinely indefinite is blended toward the pooled \eqn{\bar\Omega^*} by the
#'     minimum amount that restores positive-definiteness (closed form via Weyl's inequality), instead of
#'     relying on the eigenvalue floor alone. Useful for the sieve at small \eqn{n} / large \eqn{H}; it never
#'     fires on the well-conditioned default kernel. Changes the variance where it activates.}
#'   \item{\code{edid_hessian}}{\code{"analytic"} (default) uses the exact closed-form per-cell Hessian for
#'     the \code{higher_order} term; \code{"fd"} forces the finite-difference fallback (slower and less
#'     accurate, kept as an oracle).}
#'   \item{\code{edid_ach}}{\code{"analytic"} (default) uses the exact closed-form Ackerberg-Chen-Hahn
#'     first-step correction; \code{"fd"} forces the finite-difference oracle (validation only).}
#'   \item{\code{edid_shrink_lambda}}{Numeric, or \code{NA} (default) for the data-driven Ledoit-Wolf
#'     shrinkage intensity of the pointwise \eqn{\hat\Omega^*(X_i)} toward \eqn{\bar\Omega^*}. \code{0}
#'     disables shrinkage; a value in \eqn{[0,1]} fixes the intensity.}
#'   \item{\code{edid_eig_tol}}{Numeric, or \code{NA} (default) for the rate-based relative eigenvalue floor
#'     \eqn{n^{-a}}. A positive value sets the floor directly (condition-number cap \eqn{= 1/}\code{tol}).}
#' }
#' Other \code{edid_*} options are internal development / diagnostic hooks (e.g. \code{edid_fixed_weights},
#' \code{edid_fixed_wpw}, \code{edid_store_psiomega}, \code{edid_psiomega_fd}) and are unsupported.
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
#'     \item{\code{bstrap}}{Logical: whether the multiplier bootstrap was requested. \code{bstrap = TRUE}
#'       with \code{cband_method} left at its default selects the multiplier bootstrap, so the cell SEs and
#'       the aggregations use the did multiplier bootstrap (\code{\link[did]{mboot}} / \code{\link[did]{aggte}});
#'       under an explicit \code{cband_method = "analytic"} (or \code{higher_order = TRUE}) inference is
#'       analytic regardless of \code{bstrap}.}
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
  bstrap            = FALSE,
  biters            = 1000L,
  seed              = NULL,
  anticipation      = 0L,
  aggregate         = c("all", "overall", "event_study", "group", "calendar", "none"),
  balance_e         = NULL,
  survey_design     = NULL,
  weight_scheme     = c("efficient", "averaged", "gmm", "uniform"),
  estimation_effect = FALSE,
  cband             = TRUE,
  cband_method      = c("analytic", "multiplier"),
  higher_order      = FALSE,
  misspec_robust    = TRUE,
  trim_level        = 200
) {
  weight_method <- match.arg(weight_scheme)
  cband_method_explicit <- !missing(cband_method)   # was cband_method passed, or left at its default?
  ee_explicit   <- !missing(estimation_effect)      # did the user set the fine-grained flags explicitly?
  ho_explicit   <- !missing(higher_order)
  mr_explicit   <- !missing(misspec_robust)
  cband_method  <- match.arg(cband_method)
  estimation_effect <- isTRUE(estimation_effect)
  higher_order  <- isTRUE(higher_order)
  misspec_robust <- isTRUE(misspec_robust)
  has_cov       <- !(is.null(xformla) || identical(deparse(xformla, width.cutoff = 500L), "~1"))
  mc <- match.call()

  # ------------------------------------------------------------------
  # Higher-order ("Wick") variance refinement: validation / coercion
  # ------------------------------------------------------------------
  # The refinement puts the degenerate second-order U-statistic nuisance-estimation variance INTO the
  # coefficient covariance that the analytic sup-t crit and SEs read from. It is gated to the analytic
  # cband for two structural reasons:
  #   (1) the multiplier bootstrap resamples the (first-order) influence functions and cannot carry a
  #       degenerate-U higher-order term, so it is coerced to the analytic path (with a warning);
  #   (2) it needs the first-step sieve coefficients -- with no covariates the nuisances are unconditional
  #       means with no coefficients and Sigma_quad is identically zero, so xformla is required (error).
  # Only an EXPLICIT higher_order = TRUE validates/coerces here; the misspec_robust master switch enables the
  # higher-order term only where it already applies (covariates + analytic), so it never errors/coerces.
  if (higher_order && ho_explicit) {
    if (cband_method != "analytic") {
      warning("higher_order = TRUE requires cband_method = 'analytic' (the multiplier bootstrap cannot ",
              "carry the degenerate-U higher-order term); coercing cband_method to 'analytic'.",
              call. = FALSE)
      cband_method <- "analytic"
    }
    if (!has_cov) {
      stop("higher_order = TRUE requires a covariate formula (xformla): with no covariates the sieve ",
           "nuisances are unconditional means with no first-step coefficients, so the higher-order ",
           "variance is exactly zero. Supply xformla, or use higher_order = FALSE.", call. = FALSE)
    }
  }

  # ------------------------------------------------------------------
  # Backward-compatible bootstrap entry point
  # ------------------------------------------------------------------
  # bstrap = TRUE is the legacy switch for the multiplier bootstrap. Under the new default
  # cband_method = "analytic" a bare bstrap = TRUE would otherwise be silently ignored (no bootstrap runs
  # at the cell OR aggregation level). So when the user requests bstrap = TRUE WITHOUT explicitly choosing a
  # cband_method -- and is not in the higher_order path, which requires the analytic covariance -- select the
  # multiplier bootstrap, as the bstrap documentation promises. An explicit cband_method always wins.
  if (isTRUE(bstrap) && !cband_method_explicit && !(higher_order && ho_explicit)) {
    cband_method <- "multiplier"
  }

  # ------------------------------------------------------------------
  # misspec_robust master switch (default TRUE)
  # ------------------------------------------------------------------
  # misspec_robust = TRUE makes the reported SE account for EVERY applicable estimation effect: the
  # weight-estimation channel, the first-step nuisance ACH correction (estimation_effect), and the
  # higher-order ("Wick") nuisance term. Each is enabled only where it applies and silently skipped where it
  # does not (no covariates; multiplier path; weight_scheme = "uniform"), so default calls never warn. An
  # explicitly-set fine-grained flag overrides the bundle; misspec_robust = FALSE reverts to the plug-in
  # efficient-IF SE (honoring any individually-set estimation_effect / higher_order).
  if (misspec_robust) {
    if (!ee_explicit) estimation_effect <- has_cov
    if (!ho_explicit) higher_order      <- has_cov && cband_method == "analytic"
    if (!mr_explicit) misspec_robust    <- has_cov && weight_method != "uniform"
  }

  # ------------------------------------------------------------------
  # Argument matching
  # ------------------------------------------------------------------
  pt_assumption     <- match.arg(pt_assumption)
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
    estimation_effect = isTRUE(estimation_effect),
    higher_order  = higher_order,
    misspec_robust = misspec_robust,
    trim_level    = trim_level
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
  # Higher-order ("Wick") covariance Sigma_quad -- computed ONCE over all cells
  # ------------------------------------------------------------------
  # When higher_order is on, the SAME Sigma_quad is consumed by the cell-level band below (subset [ok, ok])
  # AND by every aggregation (aggte_edid() -> .edid_analytic_cband_agg, up to 4 times for aggregate = "all").
  # Sigma_quad[k, j] depends only on cells k and j, so sigma_quad_edid(cells)[ok, ok] is bit-identical to
  # sigma_quad_edid(cells[ok]); computing the full matrix once removes the redundant aggregate recomputes.
  sigma_quad_full <- if (isTRUE(higher_order))
    sigma_quad_edid(cells, panel_obj$cluster_indices, panel_obj$n) else NULL

  # ------------------------------------------------------------------
  # Aggregation
  # ------------------------------------------------------------------
  do_overall    <- any(aggregate %in% c("all", "overall"))
  do_event_study <- any(aggregate %in% c("all", "event_study"))
  do_group      <- any(aggregate %in% c("all", "group"))
  do_calendar   <- any(aggregate %in% c("all", "calendar"))

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
  # Multiplier-bootstrap cell SEs + simultaneous critical value (cband_method = "multiplier", the legacy
  # path). Untouched so cband_method = "multiplier" reproduces the previous behavior exactly.
  if (bstrap && cband_method == "multiplier") {
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

  # Analytic simultaneous (sup-t) uniform bands for the cell ATT(g,t) vector (cband_method = "analytic",
  # the default). Montiel Olea-Plagborg-Moller critical value from the cluster-robust analytic covariance
  # of the EIFs -- no bootstrap needed. sqrt(diag(Sigma)) equals safe_inference_edid()'s SE, so the
  # reported SEs are unchanged; only the band crit changes (pointwise z when cband = FALSE).
  if (cband_method == "analytic" && !is.null(eif_matrix)) {
    ok <- is.finite(att_gt_df$se) & att_gt_df$se > 0
    if (any(ok)) {
      Sig <- cluster_cov_edid(eif_matrix[, ok, drop = FALSE], panel_obj$cluster_indices, panel_obj$n)
      # Higher-order refinement: add the degenerate-U "Wick" covariance Sigma_quad to the first-order
      # Sigma1 so BOTH the reported SE (sqrt(diag(Sigma_HO))) and the sup-t crit come from the SAME
      # higher-order-aware covariance (no first-order-vs-higher-order splice). Sigma_quad >= 0 on the
      # diagonal, so the inflated SE is never below the plug-in SE.
      if (isTRUE(higher_order)) {
        # Reuse the once-computed full Sigma_quad; [ok, ok] == sigma_quad_edid(cells[ok]) exactly (per-(k,j)
        # entries depend only on cells k and j, so dropping the non-`ok` cells leaves the survivors unchanged).
        Sig <- Sig + sigma_quad_full[ok, ok, drop = FALSE]
      }
      bnd <- analytic_bands_edid(att_gt_df$att[ok], Sig, alp = alp, cband = isTRUE(cband), seed = seed)
      att_gt_df$se[ok]       <- bnd$se
      att_gt_df$ci_lower[ok] <- bnd$ci_lower
      att_gt_df$ci_upper[ok] <- bnd$ci_upper
      # Keep t_stat / p_value consistent with the (possibly higher-order-inflated) analytic SE. Under the
      # non-higher_order path bnd$se equals the plug-in SE, so these recompute to byte-identical values.
      att_gt_df$t_stat[ok]   <- att_gt_df$att[ok] / bnd$se
      att_gt_df$p_value[ok]  <- 2 * stats::pnorm(-abs(att_gt_df$t_stat[ok]))
    }
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
    inference_type   = if (n_bootstrap_internal > 0L && cband_method == "multiplier") "bootstrap" else "analytical",
    estimation_effect = isTRUE(estimation_effect),
    clustervars      = clustervars,
    cluster_indices  = panel_obj$cluster_indices,  # for cluster-robust re-aggregation in aggte_edid
    xformla          = xformla,
    bstrap           = bstrap,
    cband            = isTRUE(cband),
    cband_method     = cband_method,               # "analytic" (default) or "multiplier"; used by aggte_edid()
    higher_order     = higher_order,               # opt-in higher-order ("Wick") variance refinement
    misspec_robust   = fit_result$misspec_robust,  # EFFECTIVE flag (FALSE if the guards downgraded gmm/uniform/no-cov)
    seed             = seed,                        # for reproducible analytic sup-t crit in the aggregations
    biters           = as.integer(biters),         # used by aggte_edid()/as_MP_edid() for bstrap = TRUE
    cells            = cells,
    sigma_quad       = sigma_quad_full,            # higher-order ("Wick") K x K covariance (NULL unless higher_order); reused by aggte_edid()
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
  .agg <- function(ty) tryCatch(aggte_edid(edid_fit, type = ty, balance_e = balance_e, na.rm = TRUE), error = function(e) NULL)
  if (do_event_study) edid_fit$event_study <- .agg("dynamic")
  if (do_group)       edid_fit$group       <- .agg("group")
  if (do_calendar)    edid_fit$calendar    <- .agg("calendar")
  if (do_overall)     edid_fit$simple      <- .agg("simple")
  edid_fit$overall <- if (!is.null(edid_fit$event_study)) edid_fit$event_study else edid_fit$simple

  edid_fit
}
