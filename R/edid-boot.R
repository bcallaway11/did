# edid-boot.R
# Finite-sample bootstrap inference tools for edid fits:
#
#   edid_refit_bootstrap()        -- nuisance-REFITTING nonparametric cluster bootstrap: resample units
#                                    (or whole clusters), re-run the FULL edid() pipeline per draw.
#   edid_perturbation_bootstrap() -- sieve-coefficient PERTURBATION bootstrap: re-inject the estimated
#                                    nuisance-coefficient variability WITHOUT refitting anything.
#
# Both are standalone post-fit tools; nothing here changes edid() defaults (the analytic SE remains the
# default). Calibration provenance (the Chen-Sant'Anna-Xie efficient-DiD inference study, paper repo
# `Efficient_DiD_Claude/edid_inference_tests/`):
#   docs/coverage_all_estimands_findings.md  -- the nuisance-refitting bootstrap is the small-n remedy:
#     ~0.92-0.95 coverage on the weak-overlap long-horizon cells where the analytic SE covers ~0.85-0.92
#     at n = 500 (Tables B); reference harness scripts/edid_ach_aggregation_sim.R (boot_estimands).
#   docs/perturbation_bootstrap_findings.md  -- the no-refit sieve-coefficient perturbation recovers ~88%
#     of the refit bootstrap's small-n coverage gain at matmul cost (within ~1pp coverage at n = 500,
#     essentially exact by n >= 1000) for BOTH the uniform and the production efficient weight schemes;
#     reference harnesses scripts/perturbation_bootstrap.R (uniform, INDEP variant) and
#     scripts/perturbation_efficient.R (efficient, W held fixed).

# ---------------------------------------------------------------------------
# Shared internals
# ---------------------------------------------------------------------------

# Recover the estimation data from the fit's stored call when `data = NULL`
# (the update() idiom used by edid_sargan()).
.edid_boot_recover_data <- function(fit, data, envir) {
  if (is.null(data)) {
    data <- tryCatch(as.data.frame(eval(fit$call$data, envir = envir)),
                     error = function(e) NULL)
    if (is.null(data) || !is.data.frame(data) || nrow(data) == 0L) {
      stop("Could not recover the estimation data from the fit's call; pass `data` explicitly.",
           call. = FALSE)
    }
  }
  as.data.frame(data)
}

# Cheap per-draw refit configuration: exact point estimates + the requested
# aggregations, nothing else. The fit's weight_scheme / xformla / pt_assumption /
# anticipation / trim_level / moment_set / clustervars are preserved through the
# fit's stored argument snapshot (.edid_refit_args, edid-utils.R; legacy fits
# fall back to call re-evaluation); bands, the multiplier bootstrap, and ALL analytic
# estimation-effect corrections are switched off: the bootstrap only consumes
# per-draw POINT estimates, and the nonparametric resample already carries the
# weight-estimation and first-step nuisance-estimation channels (the very
# channels misspec_robust / estimation_effect / higher_order approximate
# analytically), so computing those SE corrections per draw would only slow each
# refit without changing the draw distribution.
.edid_boot_cheap_args <- function(fit, envir, aggregate) {
  args <- .edid_refit_args(fit, envir)
  args[["cband_method"]] <- NULL          # analytic default applies (bstrap is off below)
  args$aggregate         <- aggregate
  args$cband             <- FALSE
  args$bstrap            <- FALSE
  args$misspec_robust    <- FALSE
  args$estimation_effect <- FALSE
  args$higher_order      <- FALSE
  args$seed              <- NULL          # the per-draw refit is RNG-free in this configuration
  args$cores             <- 1L            # parallelism lives at the draw level
  args
}

.edid_boot_check_count <- function(value, name, min = 1L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < min || value != floor(value)) {
    stop(sprintf("`%s` must be an integer scalar >= %d.", name, min), call. = FALSE)
  }
  as.integer(value)
}

# Per-call base seed. With `seed` supplied, draws use the deterministic per-draw
# seeds seed + b (b = 1..B), so results are identical for any `cores` and any
# draw scheduling. With seed = NULL a base seed is taken from the current RNG
# stream (results differ across calls, but cores-invariance still holds within
# a call). seed + B must stay a valid integer for set.seed(), hence the
# B-aware range check / sampling headroom.
.edid_boot_seed_base <- function(seed, B) {
  if (is.null(seed)) {
    return(sample.int(.Machine$integer.max - B - 1L, 1L))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed != floor(seed) || abs(seed) > .Machine$integer.max - B - 1L) {
    stop(sprintf("`seed` must be NULL or an integer scalar with |seed| <= %d (so that seed + B is a valid seed).",
                 .Machine$integer.max - B - 1L), call. = FALSE)
  }
  as.integer(seed)
}

# Snapshot the caller's RNG state and return a restorer. Call AFTER the base
# seed is drawn: with seed = NULL the base-seed draw advances the caller's
# stream (so back-to-back unseeded calls differ), while the per-draw
# set.seed() stomping inside the draw loop is undone on exit -- the same
# save/restore discipline as aggte_edid()'s seeded multiplier bootstrap.
.edid_boot_rng_guard <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    return(function() assign(".Random.seed", old_seed, envir = .GlobalEnv))
  }
  function() NULL
}

# Dispatch draws serially or via fork-based parallel::mclapply (not on Windows,
# matching edid()'s own `cores` semantics). Each draw re-seeds itself
# (set.seed(seed_base + b)), so the two paths are numerically identical -- and
# any draw whose forked worker dies without delivering a result (mclapply
# returns NULL; the classic cause is a fork-unsafe multithreaded BLAS such as
# macOS Accelerate) is simply recomputed serially, byte-identical to what the
# fork would have produced. Forked-worker death is therefore a lost speed-up,
# never a lost / silently-NA draw.
.edid_boot_lapply <- function(X, FUN, cores, label) {
  if (cores > 1L && .Platform$OS.type != "windows") {
    res <- suppressWarnings(parallel::mclapply(X, FUN, mc.cores = cores))
    bad <- vapply(res, is.null, logical(1L))
    if (any(bad)) {
      warning(sprintf(paste0(
        "%s: %d of %d forked draws did not deliver a result (typically a fork-unsafe ",
        "multithreaded BLAS, e.g. macOS Accelerate); they were recomputed serially, so the ",
        "results are exact and identical to cores = 1. For an actual parallel speed-up use a ",
        "fork-safe (e.g. single-threaded) BLAS."), label, sum(bad), length(X)), call. = FALSE)
      res[bad] <- lapply(X[bad], FUN)
    }
    return(res)
  }
  lapply(X, FUN)
}

# Nonparametric panel resample at the independence level: units when the fit is
# unclustered, whole clusters when clustered. Duplicate draws are re-indexed to
# DISTINCT unit (and cluster) ids so the rebuilt panel is a valid panel of the
# original size -- the harness construction (paper repo,
# edid_inference_tests/scripts/edid_ach_aggregation_sim.R, boot_estimands()).
.edid_boot_resample <- function(data, idname, clustervar = NULL) {
  if (is.null(clustervar)) {
    byid <- split(seq_len(nrow(data)), data[[idname]])
    n_id <- length(byid)
    samp <- sample.int(n_id, n_id, replace = TRUE)
    rows <- unlist(byid[samp], use.names = FALSE)
    out  <- data[rows, , drop = FALSE]
    # re-index: the b-th drawn unit becomes unit b, so duplicate draws are distinct units
    out[[idname]] <- rep.int(seq_along(samp), lengths(byid)[samp])
    rownames(out) <- NULL
    return(out)
  }
  # Clustered: a cluster's units move TOGETHER (the resampling unit is the
  # cluster). Each sampled occurrence becomes a new distinct cluster id, and its
  # units get globally distinct new unit ids, so a cluster drawn twice enters as
  # two separate clusters.
  bycl <- split(seq_len(nrow(data)), data[[clustervar]])
  n_cl <- length(bycl)
  samp <- sample.int(n_cl, n_cl, replace = TRUE)
  rows <- unlist(bycl[samp], use.names = FALSE)
  out  <- data[rows, , drop = FALSE]
  out[[clustervar]] <- rep.int(seq_along(samp), lengths(bycl)[samp])
  key  <- paste(out[[clustervar]], out[[idname]], sep = "\r")
  out[[idname]] <- match(key, unique(key))
  rownames(out) <- NULL
  out
}

# Column-wise bootstrap summary: SE = sd over draws (the harness convention),
# the symmetric normal-quantile CI est +/- z * se_boot (the convention under
# which the refit bootstrap's coverage was validated -- the harness computes
# coverage as |att - true| <= qnorm(1 - alpha/2) * sd(draws)), and the
# equal-tailed percentile CI of the draws as a secondary report.
.edid_boot_stats <- function(est, draws, alpha) {
  draws <- as.matrix(draws)
  z     <- stats::qnorm(1 - alpha / 2)
  n_ok  <- colSums(is.finite(draws))
  se    <- rep(NA_real_, ncol(draws))
  pct_l <- rep(NA_real_, ncol(draws))
  pct_u <- rep(NA_real_, ncol(draws))
  for (j in seq_len(ncol(draws))) {
    v <- draws[, j]
    v <- v[is.finite(v)]
    if (length(v) >= 2L) {
      se[j]  <- stats::sd(v)
      qq     <- stats::quantile(v, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
      pct_l[j] <- qq[1L]
      pct_u[j] <- qq[2L]
    }
  }
  data.frame(
    se_boot   = se,
    n_boot    = n_ok,
    ci_lower  = est - z * se,
    ci_upper  = est + z * se,
    pct_lower = pct_l,
    pct_upper = pct_u
  )
}

# Aggregation slot / did::aggte type lookup shared by both tools.
.edid_boot_agg_slot <- c(event_study = "event_study", overall = "simple",
                         group = "group", calendar = "calendar")
.edid_boot_agg_type <- c(event_study = "dynamic", overall = "simple",
                         group = "group", calendar = "calendar")

# Named per-coefficient vector of an AGGTEobj: c(e<egt> coefficients, overall).
.edid_boot_aggte_vec <- function(a) {
  v <- numeric(0L)
  if (!is.null(a$egt) && length(a$egt)) {
    v <- stats::setNames(as.numeric(a$att.egt), paste0("e", a$egt))
  }
  c(v, overall = as.numeric(a$overall.att %||% NA_real_))
}

# Baseline AGGTEobj for an aggregation type: the fit's stored slot when present,
# else aggregate on the fly with the same machinery edid() uses.
.edid_boot_baseline_aggte <- function(fit, agg_name, balance_e) {
  a <- fit[[.edid_boot_agg_slot[[agg_name]]]]
  if (is.null(a)) {
    ty <- .edid_boot_agg_type[[agg_name]]
    a  <- aggte_edid(fit, type = ty,
                     balance_e = if (identical(ty, "dynamic")) balance_e else NULL,
                     na.rm = TRUE)
  }
  if (!inherits(a, "AGGTEobj")) {
    stop(sprintf("Could not construct the '%s' aggregation for the fit.", agg_name), call. = FALSE)
  }
  a
}

# ---------------------------------------------------------------------------
# Tool 1: nuisance-refitting cluster bootstrap
# ---------------------------------------------------------------------------

#' Nuisance-refitting cluster bootstrap for edid
#'
#' Finite-sample bootstrap inference for a fitted \code{\link{edid}} model by
#' the nonparametric cluster bootstrap that \emph{re-estimates everything} per
#' draw: units (or, for clustered fits, whole clusters) are resampled with
#' replacement, the panel is rebuilt (duplicate draws are re-indexed to
#' distinct units/clusters), and the full \code{edid()} pipeline -- first-step
#' sieve nuisances, the conditional-covariance weights, overlap trimming, every
#' \eqn{ATT(g,t)} cell, and the requested aggregations -- is re-run on each
#' resampled panel with the fit's own configuration (\code{weight_scheme},
#' \code{xformla}, \code{pt_assumption}, \code{anticipation},
#' \code{trim_level}, \code{moment_set}, clustering).
#'
#' @param fit An \code{edid_fit} returned by \code{\link{edid}}.
#' @param data The panel data used to estimate \code{fit}, or \code{NULL}
#'   (default), in which case the data expression stored in the fit's call is
#'   re-evaluated in the caller's environment (the \code{update()} idiom).
#'   Supply \code{data} explicitly when the original object is no longer
#'   reachable by that name.
#' @param B Number of bootstrap draws. Default \code{199L}. Each draw is a FULL
#'   \code{edid()} refit, so the cost is roughly \code{B} times the original
#'   fit (in its cheap configuration; see Details) -- budget accordingly,
#'   especially with \code{weight_scheme = "efficient"} at large \eqn{n}.
#' @param seed Integer base seed, or \code{NULL}. Draw \code{b} uses the
#'   deterministic per-draw seed \code{seed + b}, so results are reproducible
#'   and identical for any \code{cores} value and any draw scheduling. With
#'   \code{NULL}, a base seed is taken from the current RNG stream.
#' @param cores Number of forked workers for the draw loop
#'   (\code{\link[parallel]{mclapply}}); fork-based, so no effect on Windows.
#'   Numerically identical to \code{cores = 1L}: every draw re-seeds itself, and
#'   a draw whose forked worker dies without delivering (e.g. under a
#'   fork-unsafe multithreaded BLAS such as macOS Accelerate) is recomputed
#'   serially with the same per-draw seed, with a warning.
#' @param agg Which aggregations to bootstrap alongside the \eqn{ATT(g,t)}
#'   cells: any subset of \code{c("event_study", "overall", "group",
#'   "calendar")} (default: all four). \code{"overall"} is the cohort-share
#'   "simple" aggregate; \code{"event_study"} includes the dynamic overall
#'   (the average of \eqn{ES(e)} over \eqn{e \ge 0}).
#'
#' @details
#' \strong{When to use it.} The analytic (efficient-influence-function)
#' standard error that \code{edid()} reports is asymptotically valid, but in
#' small samples it can under-cover the weak-overlap \emph{long-horizon} cells
#' (a cohort evaluated several periods after treatment) and the event-study /
#' overall aggregates that load on them: in the calibration study the analytic
#' intervals covered ~0.85-0.92 on those cells at \eqn{n = 500} while this
#' bootstrap covered ~0.92-0.95. The resample re-estimates the nuisances and
#' the weights on every draw, so it captures the finite-sample
#' nuisance-estimation and weight-estimation variability nonparametrically --
#' the channels a plug-in (or multiplier-bootstrap) SE misses, because those
#' resample plug-in influence functions with the first step held fixed. It is
#' the recommended small-\eqn{n} inference for long-horizon estimands; for
#' large \eqn{n} the analytic SE is calibrated and \eqn{B} full refits buy
#' little.
#'
#' \strong{Per-draw configuration.} Each refit runs \code{edid()} with
#' \code{cband = FALSE}, \code{bstrap = FALSE}, and
#' \code{misspec_robust = estimation_effect = higher_order = FALSE}: the
#' bootstrap consumes only per-draw \emph{point estimates}, and the resample
#' itself already carries the estimation-effect channels those analytic SE
#' corrections approximate, so computing them per draw would only slow each
#' refit without changing the draw distribution. Point estimates are identical
#' to the fit's configuration. Per-draw warnings are suppressed (a resample of
#' a small cohort routinely triggers the small-cohort fallbacks); a draw that
#' errors entirely is recorded in \code{n_failed} and skipped, with a warning
#' when more than 5\% of draws fail. Draws in which a particular cell or
#' aggregate is unavailable (e.g. a small cohort absent from the resample)
#' enter that coordinate as \code{NA} and are dropped coordinate-wise; the
#' per-coordinate draw count is reported as \code{n_boot}.
#'
#' \strong{Confidence intervals.} \code{ci_lower} / \code{ci_upper} are the
#' symmetric normal-quantile intervals \eqn{\widehat{att} \pm z_{1-\alpha/2}\,
#' \widehat{se}_{boot}} with \eqn{\widehat{se}_{boot}} the standard deviation
#' over draws -- the convention under which the procedure's coverage was
#' validated. Equal-tailed percentile intervals of the draws are also reported
#' (\code{pct_lower} / \code{pct_upper}). The level is the fit's \code{alp}.
#'
#' @return An object of class \code{edid_refit_bootstrap}: a list with
#'   \describe{
#'     \item{\code{att_gt}}{data.frame, one row per \eqn{ATT(g,t)} cell:
#'       \code{group}, \code{time}, \code{att} (original estimate),
#'       \code{se_analytic} (the fit's reported SE), \code{se_boot},
#'       \code{n_boot}, \code{ci_lower}, \code{ci_upper} (symmetric
#'       normal-quantile, bootstrap SE), \code{pct_lower}, \code{pct_upper}
#'       (percentile), \code{is_pre}.}
#'     \item{\code{aggregates}}{named list (one element per requested
#'       \code{agg}) of data.frames with the same bootstrap columns over the
#'       aggregation coefficients (rows \code{e<egt>} plus \code{overall}),
#'       with a \code{parameter} label column in place of
#'       (\code{group}, \code{time}) and no \code{is_pre}.}
#'     \item{\code{B}, \code{n_failed}, \code{failed_messages}}{draw counts and
#'       up to 5 distinct error messages from failed draws.}
#'     \item{\code{resample}, \code{n_resample_units}}{\code{"cluster"} or
#'       \code{"unit"}, and the number of resampled blocks.}
#'     \item{\code{alpha}, \code{seed}, \code{call}}{inference level, the base
#'       seed used, and the matched call.}
#'   }
#'
#' @seealso \code{\link{edid}}, \code{\link{edid_perturbation_bootstrap}} (a
#'   no-refit approximation at a fraction of the cost).
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025).
#'   \emph{Efficient Difference-in-Differences and Event Study Estimators}.
#'   Working paper.
#'
#' @examples
#' \donttest{
#' set.seed(20260610)
#' df <- data.frame(
#'   id   = rep(1:150, each = 4),
#'   time = rep(1:4, 150),
#'   g    = rep(sample(c(2, 3, Inf), 150, replace = TRUE), each = 4),
#'   x1   = rep(rnorm(150), each = 4)
#' )
#' df$y <- rnorm(150)[df$id] + 0.2 * df$time + 0.3 * df$x1 * df$time +
#'   1 * (df$time >= df$g) + rnorm(nrow(df), 0, 0.5)
#' fit <- edid(df, "y", "id", "time", "g", xformla = ~ x1,
#'             weight_scheme = "uniform", aggregate = "event_study",
#'             cband = FALSE)
#' edid_refit_bootstrap(fit, data = df, B = 49L, seed = 1L)
#' }
#'
#' @export
edid_refit_bootstrap <- function(fit, data = NULL, B = 199L, seed = NULL, cores = 1L,
                                 agg = c("event_study", "overall", "group", "calendar")) {
  if (!inherits(fit, "edid_fit")) {
    stop("`fit` must be an `edid_fit` object returned by edid().", call. = FALSE)
  }
  B     <- .edid_boot_check_count(B, "B", min = 2L)
  cores <- .edid_boot_check_count(cores, "cores", min = 1L)
  agg   <- match.arg(agg, several.ok = TRUE)
  mc    <- match.call()

  caller <- parent.frame()
  data   <- .edid_boot_recover_data(fit, data, caller)
  args   <- .edid_boot_cheap_args(fit, caller, aggregate = agg)

  idname     <- fit$idname
  clustervar <- fit$clustervars
  alpha      <- fit$alpha %||% 0.05
  balance_e  <- args$balance_e

  # ---- baseline (original-fit) estimates the draws are centered on ----------
  agt       <- fit$att_gt
  cell_keys <- paste(agt$group, agt$time, sep = "_")
  base_objs <- stats::setNames(
    lapply(agg, function(a) .edid_boot_baseline_aggte(fit, a, balance_e)), agg)
  base_aggs <- lapply(base_objs, .edid_boot_aggte_vec)

  seed_base <- .edid_boot_seed_base(seed, B)
  restore_rng <- .edid_boot_rng_guard()
  on.exit(restore_rng(), add = TRUE)

  # ---- the draws -------------------------------------------------------------
  one_draw <- function(b) {
    set.seed(seed_base + b)
    dfb <- .edid_boot_resample(data, idname, clustervar)
    fb  <- tryCatch(
      suppressMessages(suppressWarnings(do.call(edid, c(list(data = dfb), args)))),
      error = function(e) e)
    if (inherits(fb, "error")) return(list(error = conditionMessage(fb)))
    out_cells <- stats::setNames(fb$att_gt$att, paste(fb$att_gt$group, fb$att_gt$time, sep = "_"))
    out_aggs  <- lapply(agg, function(a) {
      obj <- fb[[.edid_boot_agg_slot[[a]]]]
      if (is.null(obj)) return(NULL)
      .edid_boot_aggte_vec(obj)
    })
    names(out_aggs) <- agg
    list(cells = out_cells, aggs = out_aggs)
  }
  res <- .edid_boot_lapply(seq_len(B), one_draw, cores, label = "edid_refit_bootstrap")

  # ---- align draws on the original fit's coordinates -------------------------
  # A draw fails either inside one_draw (tryCatch -> list(error = msg)) or, under
  # mclapply, in the fork itself (NULL / an atomic "try-error" -- never index those with `$`).
  failed   <- vapply(res, function(r) {
    is.null(r) || inherits(r, "try-error") || (is.list(r) && !is.null(r$error))
  }, logical(1L))
  n_failed <- sum(failed)
  fail_msg <- unique(vapply(res[failed], function(r) {
    if (is.list(r) && !is.null(r$error)) r$error
    else if (is.null(r)) "forked worker delivered no result"
    else trimws(as.character(r)[1L])
  }, character(1L)))
  if (n_failed > 0.05 * B) {
    warning(sprintf(paste0(
      "edid_refit_bootstrap: %d of %d bootstrap draws (%.1f%%) failed entirely and were skipped ",
      "(first error: %s). Bootstrap SEs are computed from the remaining draws; with this many ",
      "failures (typically a cohort or comparison group too small to survive resampling) treat ",
      "them with caution."), n_failed, B, 100 * n_failed / B, fail_msg[1L]), call. = FALSE)
  }

  cell_draws <- matrix(NA_real_, nrow = B, ncol = length(cell_keys),
                       dimnames = list(NULL, cell_keys))
  agg_draws  <- stats::setNames(lapply(agg, function(a) {
    matrix(NA_real_, nrow = B, ncol = length(base_aggs[[a]]),
           dimnames = list(NULL, names(base_aggs[[a]])))
  }), agg)
  for (b in seq_len(B)) {
    if (failed[b]) next
    rb <- res[[b]]
    mi <- match(cell_keys, names(rb$cells))
    cell_draws[b, !is.na(mi)] <- rb$cells[mi[!is.na(mi)]]
    for (a in agg) {
      va <- rb$aggs[[a]]
      if (is.null(va)) next
      ma <- match(colnames(agg_draws[[a]]), names(va))
      agg_draws[[a]][b, !is.na(ma)] <- va[ma[!is.na(ma)]]
    }
  }

  # ---- assemble ---------------------------------------------------------------
  att_tab <- cbind(
    data.frame(group = agt$group, time = agt$time, att = agt$att, se_analytic = agt$se),
    .edid_boot_stats(agt$att, cell_draws, alpha),
    data.frame(is_pre = agt$is_pre)
  )
  rownames(att_tab) <- NULL

  agg_tabs <- stats::setNames(lapply(agg, function(a) {
    est <- base_aggs[[a]]
    ao  <- base_objs[[a]]
    se_an <- c(if (!is.null(ao$egt) && length(ao$egt)) as.numeric(ao$se.egt) else numeric(0L),
               as.numeric(ao$overall.se %||% NA_real_))
    tab <- cbind(
      data.frame(parameter = names(est), att = as.numeric(est), se_analytic = se_an),
      .edid_boot_stats(as.numeric(est), agg_draws[[a]], alpha)
    )
    rownames(tab) <- NULL
    tab
  }), agg)

  out <- list(
    att_gt            = att_tab,
    aggregates        = agg_tabs,
    B                 = B,
    n_failed          = n_failed,
    failed_messages   = utils::head(fail_msg, 5L),
    resample          = if (is.null(clustervar)) "unit" else "cluster",
    n_resample_units  = if (is.null(clustervar)) length(unique(data[[idname]]))
                        else length(unique(data[[clustervar]])),
    alpha             = alpha,
    seed              = seed_base,
    call              = mc
  )
  class(out) <- c("edid_refit_bootstrap", "list")
  out
}

#' @describeIn edid_refit_bootstrap Print method.
#' @param x an \code{edid_refit_bootstrap} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_refit_bootstrap <- function(x, digits = 4, ...) {
  cat("\nNuisance-refitting cluster bootstrap for edid\n")
  cat(sprintf("  B = %d full edid() refits (%d failed); resampling: %s blocks (%d)\n",
              x$B, x$n_failed, x$resample, x$n_resample_units))
  cat(sprintf("  CI: att +/- z * se_boot at level %.3g (percentile CIs also reported)\n\n",
              1 - x$alpha))
  fmt <- function(tab) {
    num <- vapply(tab, is.numeric, logical(1L))
    tab[num] <- lapply(tab[num], function(z) signif(z, digits))
    tab
  }
  cat("Group-time ATT(g,t):\n")
  print(fmt(x$att_gt), row.names = FALSE)
  for (a in names(x$aggregates)) {
    cat(sprintf("\nAggregation '%s':\n", a))
    print(fmt(x$aggregates[[a]]), row.names = FALSE)
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Tool 2: sieve-coefficient perturbation bootstrap (no refit)
# ---------------------------------------------------------------------------

# Matrix square root of a coefficient covariance: Cholesky with a tiny jitter,
# eigen square root (negative eigenvalues clipped at 0) as the fallback --
# ported from the validation harnesses (perturbation_bootstrap.R::sqrt_cov).
.edid_boot_sqrt_cov <- function(V) {
  V <- (V + t(V)) / 2
  L <- tryCatch(t(chol(V + diag(1e-10, ncol(V)))), error = function(e) NULL)
  if (!is.null(L)) return(L)
  ed <- eigen(V, symmetric = TRUE)
  ed$vectors %*% diag(sqrt(pmax(ed$values, 0)), length(ed$values))
}

#' Sieve-coefficient perturbation bootstrap for edid (no refit)
#'
#' A cheap, no-refit finite-sample variance correction for a fitted
#' \code{\link{edid}} model. The first-step sieve nuisances (the conditional
#' means \eqn{m} and propensity ratios \eqn{r}) are estimated, and in small
#' samples that estimation injects higher-order variability that the plug-in
#' efficient-influence-function SE misses. This tool re-creates that
#' variability \emph{without re-solving anything}: for each first-step nuisance
#' \eqn{k} it forms the sieve-coefficient sandwich covariance
#' \eqn{\widehat{V}_{\theta,k} = n^{-2} H_k^{-1}\,
#' (\mathrm{score}_k'\mathrm{score}_k)\, H_k^{-1}} from the stored M-estimator
#' pieces, draws \eqn{\theta_k^* \sim N(\hat\theta_k, \widehat{V}_{\theta,k})}
#' \emph{independently across distinct coefficient blocks} (the validated
#' "INDEP" variant; the joint draw is dominated by it) -- nuisance entries that
#' share ONE underlying fitted coefficient vector, like the
#' \code{ratio_method = "coherent"} ratios and inverse propensities that all
#' read the same joint multinomial-logit system, are dedup'd on the aux's
#' \code{coef_id} and share a single draw per replication, mapped into each
#' entry through its own chain-rule Jacobian (independent draws there would
#' break the exact cross-entry coupling of the shared system) -- recomputes the
#' doubly-robust generated outcomes nonlinearly at the perturbed predictions
#' \eqn{\hat\nu + B_k(\theta_k^* - \hat\theta_k)} with the weights \eqn{W} and
#' the overlap-trim masks held FIXED at the original fit, and reads off the
#' perturbed estimate per draw. By Neyman orthogonality the first-order term is
#' \eqn{\approx 0}, so the draw variance \eqn{Var_b(att^*)} estimates the
#' higher-order nuisance-estimation variance, and the reported combined SE is
#' \deqn{\widehat{se}_{comb} = \sqrt{\widehat{se}_{plug}^2 + Var_b(att^*)}.}
#'
#' @param fit An \code{edid_fit} from \code{edid()} with a covariate formula
#'   and \code{weight_scheme} \code{"efficient"} or \code{"uniform"} (the two
#'   schemes the construction was validated for; other schemes error -- use
#'   \code{\link{edid_refit_bootstrap}} there). Without covariates there are no
#'   first-step sieve coefficients to perturb and the function errors.
#' @param data The panel data used to estimate \code{fit}, or \code{NULL}
#'   (default: re-evaluate the data expression stored in the fit's call in the
#'   caller's environment).
#' @param B Number of perturbation draws. Default \code{499L}. Each draw costs
#'   a few matrix products per cell (no nuisance/weight re-solve), so large
#'   \code{B} is cheap.
#' @param seed Integer base seed, or \code{NULL}. Draw \code{b} re-seeds with
#'   \code{seed + b} and consumes the nuisance draws in a fixed order, so
#'   results are reproducible and identical for any \code{cores} value.
#' @param cores Number of forked workers for the draw loop (fork-based; no
#'   effect on Windows). Numerically identical to \code{cores = 1L}; draws a
#'   forked worker fails to deliver are recomputed serially with the same
#'   per-draw seed (see \code{\link{edid_refit_bootstrap}}).
#' @param agg Which aggregations to report alongside the cells: any subset of
#'   \code{c("event_study", "overall", "group", "calendar")} (default all).
#'   Aggregates are linear in the cells with weights held fixed at the original
#'   fit (the cohort-share weight-estimation effect is not perturbed),
#'   matching the validation harness.
#'
#' @details
#' \strong{Calibration provenance.} In the Chen-Sant'Anna-Xie inference study
#' this construction recovers ~88\% of the nuisance-refitting bootstrap's
#' small-sample coverage improvement on the hardest weak-overlap long-horizon
#' cell at \eqn{n = 500} (coverage ~0.85 plug-in, ~0.91 perturbation, ~0.92
#' refit bootstrap), and is essentially equivalent to the refit bootstrap from
#' \eqn{n \gtrsim 1000}, at a tiny fraction of its cost -- for both supported
#' weight schemes. Use it as the cheap default small-sample check;
#' \code{\link{edid_refit_bootstrap}} remains the most reliable choice at the
#' smallest sample sizes (it also captures the weight-estimation channel and
#' the heavy resampling tail).
#'
#' \strong{What is (and is not) recomputed.} The tool re-derives the first-step
#' nuisances, weights, and plug-in influence functions from \code{data} with
#' the package's own internal estimators (plug-in regime, \code{K = 1}) under
#' the fit's configuration, and verifies that the re-derived cell estimates
#' reproduce \code{fit$att_gt$att} exactly; a mismatch (changed data, or
#' \code{edid_omega_method} / shrinkage / eigen-floor options differing from
#' fit time) is an error, not a silent miscalibration. Nuisances whose sieve
#' fit fell back to a constant (tiny cohorts) carry no estimated coefficients
#' and are left unperturbed, mirroring the higher-order machinery.
#'
#' \strong{Confidence intervals.} Only the Wald interval
#' \eqn{\widehat{att} \pm z_{1-\alpha/2}\,\widehat{se}_{comb}} is reported (at
#' the fit's \code{alp}): the perturbation draws simulate the
#' nuisance-estimation \emph{channel}, not the full sampling distribution of
#' the estimator, so percentile intervals of the draws would be meaningless --
#' and the study found the small-sample coverage gap to be one of scale, with
#' percentile refinements unnecessary. The combined SE pairs the perturbation
#' variance with the \emph{plug-in} EIF SE (reported as \code{se_plug};
#' cluster-robust when the fit is clustered). Do not stack it on top of the
#' \code{higher_order} ("Wick") analytic refinement -- that term is the
#' quadratic approximation of the same channel this tool simulates.
#'
#' @return An object of class \code{edid_perturbation_bootstrap}: a list with
#'   \describe{
#'     \item{\code{att_gt}}{data.frame, one row per cell: \code{group},
#'       \code{time}, \code{att}, \code{se_analytic} (the fit's reported SE),
#'       \code{se_plug} (plug-in EIF SE), \code{se_pert} (sd of the
#'       perturbation draws), \code{se_combined}, \code{n_pert},
#'       \code{ci_lower}, \code{ci_upper} (Wald, combined SE), \code{is_pre}.}
#'     \item{\code{aggregates}}{named list of data.frames for the requested
#'       aggregations (rows \code{e<egt>} plus \code{overall}; columns
#'       \code{parameter}, \code{att}, \code{se_plug}, \code{se_pert},
#'       \code{se_combined}, \code{n_pert}, \code{ci_lower}, \code{ci_upper}),
#'       built from the fixed linear cell-to-aggregate map of the original
#'       fit.}
#'     \item{\code{B}, \code{n_failed}, \code{n_nuisances},
#'       \code{weight_scheme}, \code{alpha}, \code{seed}, \code{call}}{draw
#'       counts, the number of perturbed nuisance functions, and metadata.}
#'   }
#'
#' @seealso \code{\link{edid}}, \code{\link{edid_refit_bootstrap}} (the
#'   gold-standard refitting bootstrap this tool approximates).
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025).
#'   \emph{Efficient Difference-in-Differences and Event Study Estimators}.
#'   Working paper. \cr
#'   Ackerberg, D., Chen, X., and Hahn, J. (2012). A Practical Asymptotic
#'   Variance Estimator for Two-Step Semiparametric Estimators. \emph{Review of
#'   Economics and Statistics}, 94(2), 481-498.
#'
#' @examples
#' \donttest{
#' set.seed(20260610)
#' df <- data.frame(
#'   id   = rep(1:150, each = 4),
#'   time = rep(1:4, 150),
#'   g    = rep(sample(c(2, 3, Inf), 150, replace = TRUE), each = 4),
#'   x1   = rep(rnorm(150), each = 4)
#' )
#' df$y <- rnorm(150)[df$id] + 0.2 * df$time + 0.3 * df$x1 * df$time +
#'   1 * (df$time >= df$g) + rnorm(nrow(df), 0, 0.5)
#' fit <- edid(df, "y", "id", "time", "g", xformla = ~ x1,
#'             weight_scheme = "uniform", aggregate = "event_study",
#'             cband = FALSE)
#' edid_perturbation_bootstrap(fit, data = df, B = 199L, seed = 1L)
#' }
#'
#' @export
edid_perturbation_bootstrap <- function(fit, data = NULL, B = 499L, seed = NULL, cores = 1L,
                                        agg = c("event_study", "overall", "group", "calendar")) {
  if (!inherits(fit, "edid_fit")) {
    stop("`fit` must be an `edid_fit` object returned by edid().", call. = FALSE)
  }
  B     <- .edid_boot_check_count(B, "B", min = 2L)
  cores <- .edid_boot_check_count(cores, "cores", min = 1L)
  agg   <- match.arg(agg, several.ok = TRUE)
  mc    <- match.call()

  caller <- parent.frame()
  data   <- .edid_boot_recover_data(fit, data, caller)
  args   <- .edid_refit_args(fit, caller)

  # Fit configuration from the stored argument snapshot (the %||% defaults cover
  # legacy fits recovered through the call-re-evaluation fallback).
  ws <- args$weight_scheme %||% "efficient"
  ws <- match.arg(ws, c("efficient", "averaged", "gmm", "uniform"))
  if (!ws %in% c("efficient", "uniform")) {
    stop(sprintf(paste0(
      "edid_perturbation_bootstrap() supports weight_scheme = 'efficient' or 'uniform' (the two ",
      "schemes the perturbation construction was validated for); this fit uses '%s'. Use ",
      "edid_refit_bootstrap(), which re-estimates everything and covers any scheme."), ws),
      call. = FALSE)
  }
  xformla <- fit$xformla
  has_cov <- !is.null(xformla) && inherits(xformla, "formula") && length(all.vars(xformla)) > 0L
  if (!has_cov) {
    stop(paste0(
      "edid_perturbation_bootstrap() requires a covariate fit (xformla): with no covariates the ",
      "first-step nuisances are unconditional means with no sieve coefficients to perturb, and ",
      "there is no nuisance-estimation channel for this tool to correct."), call. = FALSE)
  }
  trim_level <- args$trim_level %||% 200
  # Sieve df of the fit (integer or "ic"). The rebuild below must use the SAME first-step basis
  # as the fit; under "ic" the per-fit selection is deterministic given the (verified-identical)
  # original data, so re-running it reproduces the fit's selected dimensions exactly.
  bs_df_fit  <- args$bs_df %||% 4L
  # Cross-cohort ratio construction of the fit (the rebuild must match it exactly; the
  # exactness guard below would otherwise reject). Snapshot fallback = edid()'s default.
  ratio_method_fit <- fit$ratio_method %||% args$ratio_method %||% "coherent"
  yname  <- args$yname
  if (is.null(yname)) stop("Could not recover `yname` from the fit's call.", call. = FALSE)
  alpha      <- fit$alpha %||% 0.05
  balance_e  <- args$balance_e
  pt         <- fit$pt_assumption

  # ---- rebuild the panel exactly as edid() does ------------------------------
  if (is.numeric(data[[fit$gname]])) {
    zero_nt <- is.finite(data[[fit$gname]]) & data[[fit$gname]] == 0
    if (any(zero_nt)) data[[fit$gname]] <- ifelse(zero_nt, Inf, data[[fit$gname]])
  }
  panel <- prepare_edid_panel(
    data = data, yname = yname, idname = fit$idname, tname = fit$tname, gname = fit$gname,
    xformla = xformla, clustervars = fit$clustervars, anticipation = fit$anticipation)
  if (is.null(panel$covariate_matrix)) {
    stop("The covariate matrix could not be rebuilt from `data` and the fit's xformla.", call. = FALSE)
  }
  if (!identical(panel$n, fit$n) ||
      !isTRUE(all.equal(panel$unit_cohorts, fit$unit_cohorts))) {
    stop("`data` does not reproduce the panel the fit was estimated on (n or cohort assignment ",
         "differs); pass the original estimation data.", call. = FALSE)
  }
  n  <- panel$n
  tg <- panel$treatment_groups
  tp <- panel$time_periods
  p1 <- panel$period_1

  # ---- smoother + kernel hoist (mirrors fit_edid_cells) -----------------------
  omega_method <- getOption("edid_omega_method", "kernel")
  omega_fun <- switch(omega_method,
    sieve       = compute_omega_star_sieve_edid,
    kernel_orig = compute_omega_star_cov_edid,
    compute_omega_star_kernel_fast_edid)
  kern_bw <- NULL; kern_K <- NULL
  if (ws == "efficient" && !identical(omega_method, "sieve")) {
    kk <- build_kernel_weights_edid(panel$covariate_matrix)
    kern_bw <- kk$bw; kern_K <- kk$K
    .ks <- rowSums(kern_K); .ksq <- rowSums(kern_K^2)
    attr(kern_K, "edid_m_eff") <- stats::median(.ks^2 / pmax(.ksq, .Machine$double.eps))
  }

  # ---- first-step nuisances, plug-in regime (mirrors fit_edid_cells' .gbuild) -
  fold_id <- rep(1L, n)
  # Thin-cohort guard threshold of the ORIGINAL fit (legacy fits without the field map to 2,
  # the inert legacy threshold), so the rebuilt pair sets reproduce the fit's exactly --
  # otherwise the exactness guard below would reject a guarded fit.
  mpu_fit <- fit$min_pair_units %||% 2L
  cohort_sizes_boot <- stats::setNames(
    vapply(tg, function(gg) sum(panel$unit_cohorts == gg), numeric(1L)),
    as.character(tg))
  gcache  <- stats::setNames(lapply(tg, function(g) {
    pairs_g <- enumerate_valid_pairs_edid(
      target_g = g, treatment_groups = tg, time_periods = tp, period_1 = p1,
      pt_assumption = pt, anticipation = panel$anticipation, moment_set = fit$moment_set)
    pairs_g <- apply_thin_cohort_guard_edid(g, pairs_g, cohort_sizes_boot, mpu_fit, pt)$pairs
    gb <- list(pairs = pairs_g, pfn = NULL, prop_ratios = NULL, r_aux = NULL,
               inv_p = NULL, trim_keep = NULL)
    if (nrow(pairs_g) == 0L) return(gb)
    pfn <- pairs_g
    self_cmp <- is.finite(pfn$gp) & (pfn$gp == g); if (any(self_cmp)) pfn$gp[self_cmp] <- Inf
    cross_pairs <- pairs_g[is.finite(pairs_g$gp) & pairs_g$gp != g, , drop = FALSE]
    if (nrow(cross_pairs) > 0L)
      pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cross_pairs$tpre))))
    gb$pfn <- pfn
    pr <- suppressWarnings(estimate_all_propensity_ratios(
      panel_obj = panel, g = g, pairs = pfn, bs_df = bs_df_fit, K_folds = 1L,
      fold_id = fold_id, return_aux = TRUE, ratio_method = ratio_method_fit))
    gb$prop_ratios <- pr$predictions
    gb$r_aux       <- pr$aux
    gb$inv_p <- suppressWarnings(estimate_all_inverse_propensities(
      panel_obj = panel, g = g, pairs = pairs_g, bs_df = bs_df_fit, K_folds = 1L, fold_id = fold_id,
      ratio_method = ratio_method_fit))
    gb$trim_keep <- build_trim_keep_edid(gb$prop_ratios, gb$inv_p, trim_level, n)
    gb
  }), as.character(tg))

  # Global conditional-mean cache over the distinct (gp, period) combos any cell
  # uses (mirrors fit_edid_cells' mcache; plug-in fits are deterministic).
  iter_periods <- tp[tp != p1]
  combo_set <- unique(do.call(rbind, lapply(tg, function(g) {
    pfn <- gcache[[as.character(g)]]$pfn; if (is.null(pfn)) return(NULL)
    rbind(expand.grid(gp = unique(pfn$gp), period = iter_periods, KEEP.OUT.ATTRS = FALSE),
          data.frame(gp = pfn$gp, period = pfn$tpre))
  })))
  if (is.null(combo_set) || nrow(combo_set) == 0L) {
    stop("No estimable cells to perturb (no cohort has a valid comparison pair).", call. = FALSE)
  }
  mcache_pred <- list(); mcache_aux <- list()
  for (i in seq_len(nrow(combo_set))) {
    mm <- suppressWarnings(estimate_all_conditional_means(
      panel_obj = panel,
      pairs = data.frame(gp = combo_set$gp[i], tpre = combo_set$period[i]),
      t_val = combo_set$period[i], bs_df = bs_df_fit, K_folds = 1L, fold_id = fold_id,
      return_aux = TRUE))
    mcache_pred <- c(mcache_pred, mm$predictions)
    mcache_aux  <- c(mcache_aux,  mm$aux)
  }

  # ---- per-cell rebuild: generated outcomes, FIXED weights, plug-in EIF -------
  agt      <- fit$att_gt
  K_cells  <- nrow(agt)
  att0     <- rep(NA_real_, K_cells)
  se_plug  <- rep(NA_real_, K_cells)
  eif_plug <- matrix(NA_real_, n, K_cells)
  cell_info <- vector("list", K_cells)
  for (k in seq_len(K_cells)) {
    if (!is.finite(agt$att[k])) next                     # NA cell in the fit: nothing to perturb
    g  <- agt$group[k]; t <- agt$time[k]
    gk <- as.character(g)
    gb <- gcache[[gk]]
    if (is.null(gb) || is.null(gb$pairs) || nrow(gb$pairs) == 0L) next
    go_res <- compute_generated_outcomes_cov_edid(
      panel, g, t, gb$pairs, gb$prop_ratios, mcache_pred, pt,
      trim_keep = gb$trim_keep, return_trim_info = TRUE)
    go <- go_res$gen_out
    # Mirror fit_edid_cells: DEAD pairs (own overlap mask retains no treated mass) are dropped from the
    # pair set BEFORE weights, so the rebuilt cell reproduces the fit's surviving moment stack exactly
    # (the exactness guard below would otherwise fail under a binding trim). The surviving pairs are
    # stored in cell_info so every perturbation draw rebuilds the SAME moment set.
    pairs_k <- gb$pairs
    keep_k  <- go_res$keep
    mkept_k <- go_res$m_kept
    if (!is.null(go_res$dead) && any(go_res$dead)) {
      alive   <- !go_res$dead
      pairs_k <- pairs_k[alive, , drop = FALSE]
      rownames(pairs_k) <- NULL
      go      <- go[, alive, drop = FALSE]
      if (!is.null(keep_k)) { keep_k <- keep_k[, alive, drop = FALSE]; mkept_k <- mkept_k[alive] }
      if (nrow(pairs_k) == 0L) next
    }
    if (anyNA(go)) next
    if (ws == "efficient") {
      kp_cache <- new.env(parent = emptyenv())
      omega_arr <- omega_fun(panel, g, t, pairs_k, gb$prop_ratios, mcache_pred,
                             gb$inv_p, bw = kern_bw, K_mat = kern_K,
                             return_pointwise = TRUE, kp_cache = kp_cache,
                             keep = if (!is.null(keep_k) && ncol(keep_k) > 0L) keep_k[, 1L] else NULL)
      W  <- compute_pointwise_weights_edid(omega_arr, d = ncol(panel$covariate_matrix))
      wY <- rowSums(go * W)
    } else {
      H  <- nrow(pairs_k)
      W  <- rep(1 / H, H)
      wY <- drop(go %*% W)
    }
    att_k <- mean(wY)
    eif_k <- compute_eif_cov_edid(panel, go, W, att_k, g, keep_k, mkept_k)
    att0[k]      <- att_k
    se_plug[k]   <- safe_inference_edid(eif_k, panel$cluster_indices, alpha, att_k)$se
    eif_plug[, k] <- eif_k
    cell_info[[k]] <- list(g = g, t = t, gkey = gk, W = W, pairs = pairs_k)
  }

  # Exactness guard: the re-derived cells must reproduce the fit (same data, same
  # edid_* options as at fit time), otherwise the perturbation would be centered
  # on a different estimator. Catches changed data and changed smoother /
  # shrinkage / eigen-floor options.
  chk <- is.finite(att0) & is.finite(agt$att)
  scale_att <- 1 + max(abs(agt$att[chk]), 0, na.rm = TRUE)
  if (!identical(which(is.finite(att0)), which(is.finite(agt$att))) ||
      (any(chk) && max(abs(att0[chk] - agt$att[chk])) > 1e-6 * scale_att)) {
    stop(paste0(
      "edid_perturbation_bootstrap() could not reproduce the fit's cell estimates from `data`: ",
      "either the data changed, or the edid_* options in force (edid_omega_method, ",
      "edid_shrink_lambda, edid_eig_tol, ...) differ from fit time. Restore them and retry."),
      call. = FALSE)
  }
  used <- which(vapply(cell_info, Negate(is.null), logical(1L)))
  if (length(used) == 0L) {
    stop("No estimable cells to perturb (all cells are NA).", call. = FALSE)
  }

  # ---- the perturbable nuisances ----------------------------------------------
  # One entry per DISTINCT first-step nuisance function: the propensity ratios
  # r_{g,g'} (per target cohort) and the conditional means m_{g',s} (shared
  # across cohorts). V_theta = H^-1 (S'S) H^-1 / n^2 is the sieve-coefficient
  # sandwich from the stored M-estimator pieces (same objects and normalization
  # as the validation harnesses' fit_m/fit_r). Fallback fits (no coefficients)
  # are skipped, mirroring edid_nuisance_blocks().
  #
  # COEFFICIENT-BLOCK identity (cid). Independent per-target fits draw one
  # Gaussian coefficient vector each (the validated "INDEP" variant) -- their cid
  # defaults to the entry's own name, preserving the legacy draw stream bit for
  # bit. Entries that SHARE one underlying fitted coefficient vector -- the
  # ratio_method = "coherent" entries, whose r_{g,g'} all read the SAME joint
  # multinomial-logit system (aux field coef_id) -- must share ONE draw per
  # replication: independent draws per entry would break the perfect cross-entry
  # coupling (e.g. dr_{g,g'} and dr_{g',g} move reciprocally through the shared
  # blocks) and mis-state the perturbation variance of every aggregate. The
  # draws are therefore dedup'd on cid: one N(0, V_theta) coefficient draw per
  # DISTINCT cid, mapped into each entry through ITS OWN chain-rule Jacobian B.
  infos <- list()
  for (g in tg) {
    gb <- gcache[[as.character(g)]]
    for (key in names(gb$r_aux)) {
      a <- gb$r_aux[[key]]
      if (is.null(a) || isTRUE(a$is_fallback) || is.null(a$B_test)) next
      infos[[paste0("r:", g, ":", key)]] <- list(
        type = "r", gkey = as.character(g), key = key, B = a$B_test,
        cid = a$coef_id %||% paste0("r:", g, ":", key),
        V = (a$H_inv %*% crossprod(a$score_mat) %*% a$H_inv) / n^2)
    }
  }
  for (key in names(mcache_aux)) {
    a <- mcache_aux[[key]]
    if (is.null(a) || isTRUE(a$is_fallback) || is.null(a$B_test)) next
    infos[[paste0("m:", key)]] <- list(
      type = "m", gkey = NA_character_, key = key, B = a$B_test,
      cid = a$coef_id %||% paste0("m:", key),
      V = (a$H_inv %*% crossprod(a$score_mat) %*% a$H_inv) / n^2)
  }
  if (length(infos) == 0L) {
    stop(paste0(
      "All first-step nuisance fits fell back to constants (no estimated sieve coefficients to ",
      "perturb); the perturbation bootstrap is not informative here. Use edid_refit_bootstrap()."),
      call. = FALSE)
  }
  info_names <- names(infos)                       # fixed draw order => cores-invariant RNG
  # One Cholesky/eigen square root per DISTINCT coefficient block (entries sharing a cid carry
  # bitwise-identical score_mat/H_inv, hence the same V; computed once at the first encounter).
  Lmap <- list()
  for (nm in info_names) {
    cid <- infos[[nm]]$cid
    if (is.null(Lmap[[cid]])) Lmap[[cid]] <- .edid_boot_sqrt_cov(infos[[nm]]$V)
    infos[[nm]]$V <- NULL                          # V no longer needed; keep the draw objects lean
  }

  # ---- the draws ----------------------------------------------------------------
  seed_base <- .edid_boot_seed_base(seed, B)
  restore_rng <- .edid_boot_rng_guard()
  on.exit(restore_rng(), add = TRUE)
  one_draw <- function(b) {
    set.seed(seed_base + b)
    # one independent Gaussian coefficient draw per DISTINCT coefficient block (cid; "INDEP"
    # across genuinely independent fits): the SAME draw enters every cell -- and every
    # ENTRY -- that consumes that block, so both the cross-cell correlation needed by the
    # aggregations and the cross-entry coupling of shared-system (coherent) nuisances are
    # inherited coherently. Draws happen at the first encounter of each cid in the fixed
    # info_names order, so fits with no shared blocks reproduce the legacy stream exactly.
    pr_pert <- lapply(gcache, function(gb) gb$prop_ratios)
    cm_pert <- mcache_pred
    delta <- list()                                # cid -> the replication's coefficient draw
    for (nm in info_names) {
      ii <- infos[[nm]]
      if (is.null(delta[[ii$cid]])) {
        Lc <- Lmap[[ii$cid]]
        delta[[ii$cid]] <- as.vector(Lc %*% stats::rnorm(ncol(Lc)))
      }
      shift <- as.vector(ii$B %*% delta[[ii$cid]])
      if (ii$type == "r") {
        pr_pert[[ii$gkey]][[ii$key]] <- pr_pert[[ii$gkey]][[ii$key]] + shift
      } else {
        cm_pert[[ii$key]] <- cm_pert[[ii$key]] + shift
      }
    }
    out <- rep(NA_real_, K_cells)
    for (k in used) {
      ci_k <- cell_info[[k]]
      gb   <- gcache[[ci_k$gkey]]
      # W, the overlap-trim masks, AND the surviving (dead-pair-dropped) pair set are FROZEN at the
      # original fit; only the nuisance predictions move (the harness construction:
      # att* = mean(W o go(theta*))). ci_k$pairs is the surviving set, so the rebuilt columns line up
      # with the frozen W under a binding trim.
      gob <- compute_generated_outcomes_cov_edid(
        panel, ci_k$g, ci_k$t, ci_k$pairs, pr_pert[[ci_k$gkey]], cm_pert, pt,
        trim_keep = gb$trim_keep)
      out[k] <- if (is.matrix(ci_k$W)) mean(rowSums(gob * ci_k$W)) else mean(drop(gob %*% ci_k$W))
    }
    out
  }
  res <- .edid_boot_lapply(seq_len(B), function(b) tryCatch(one_draw(b), error = function(e) e),
                           cores, label = "edid_perturbation_bootstrap")

  failed   <- vapply(res, function(r) {
    is.null(r) || inherits(r, "try-error") || inherits(r, "error")
  }, logical(1L))
  n_failed <- sum(failed)
  if (n_failed > 0.05 * B) {
    warning(sprintf(
      "edid_perturbation_bootstrap: %d of %d perturbation draws (%.1f%%) failed and were skipped.",
      n_failed, B, 100 * n_failed / B), call. = FALSE)
  }
  draws <- matrix(NA_real_, nrow = B, ncol = K_cells)
  for (b in seq_len(B)) if (!failed[b]) draws[b, ] <- res[[b]]

  # ---- per-cell combined SE + Wald CI -------------------------------------------
  z <- stats::qnorm(1 - alpha / 2)
  se_pert <- rep(NA_real_, K_cells)
  n_pert  <- colSums(is.finite(draws))
  for (k in used) {
    v <- draws[, k]; v <- v[is.finite(v)]
    if (length(v) >= 2L) se_pert[k] <- stats::sd(v)
  }
  se_comb <- sqrt(se_plug^2 + ifelse(is.na(se_pert), 0, se_pert)^2)
  se_comb[!is.finite(se_plug)] <- NA_real_
  att_tab <- data.frame(
    group = agt$group, time = agt$time, att = agt$att,
    se_analytic = agt$se, se_plug = se_plug, se_pert = se_pert,
    se_combined = se_comb, n_pert = n_pert,
    ci_lower = agt$att - z * se_comb, ci_upper = agt$att + z * se_comb,
    is_pre = agt$is_pre)
  rownames(att_tab) <- NULL

  # ---- aggregations: fixed linear map over the cells ------------------------------
  # The aggregations are exactly linear in the att(g,t) vector with weights that
  # do not depend on the att values; recover the map by finite differences of
  # did::aggte (the same device aggte_edid's higher-order path uses), then push
  # the perturbation draws and the plug-in EIFs through it. Cohort-share weights
  # are held fixed (matching the validation harness); the refit bootstrap is the
  # tool that re-estimates them.
  reagg <- function(ty, att_vec) {
    f2 <- fit
    f2$att_gt$att <- att_vec
    aa <- aggte(as_MP_edid(f2, bstrap = FALSE, cband = FALSE), type = ty,
                balance_e = if (identical(ty, "dynamic")) balance_e else NULL,
                na.rm = TRUE, bstrap = FALSE)
    v <- numeric(0L)
    if (!is.null(aa$egt) && length(aa$egt)) {
      v <- stats::setNames(as.numeric(aa$att.egt), paste0("e", aa$egt))
    }
    c(v, overall = as.numeric(aa$overall.att %||% NA_real_))
  }
  draws_z <- draws; draws_z[, !is.finite(att0)] <- 0     # structurally-NA cells carry weight 0
  eif_z   <- eif_plug; eif_z[, !is.finite(att0)] <- 0
  agg_tabs <- stats::setNames(vector("list", length(agg)), agg)
  for (a in agg) {
    ty     <- .edid_boot_agg_type[[a]]
    base_v <- tryCatch(reagg(ty, att0), error = function(e) NULL)
    if (is.null(base_v) || !length(base_v)) next
    A <- matrix(0, length(base_v), K_cells)
    for (k in which(is.finite(att0))) {
      att_p <- att0; att_p[k] <- att0[k] + 1e-4
      col <- tryCatch((reagg(ty, att_p) - base_v) / 1e-4, error = function(e) NULL)
      if (!is.null(col) && length(col) == length(base_v) && all(is.finite(col))) A[, k] <- col
    }
    # zero-intercept linearity check: the map must reproduce the baseline exactly
    fin <- is.finite(base_v)
    lin <- drop(A %*% ifelse(is.finite(att0), att0, 0))
    if (any(fin) && max(abs(lin[fin] - base_v[fin])) > 1e-6 * (1 + max(abs(base_v[fin])))) {
      warning(sprintf(paste0(
        "edid_perturbation_bootstrap: the '%s' aggregation is not an exact linear map of the ",
        "cells here; skipping it (cells are unaffected)."), a), call. = FALSE)
      next
    }
    agg_draws  <- draws_z %*% t(A)                       # B x n_agg, NA rows propagate
    if_agg     <- eif_z %*% t(A)                         # n x n_agg plug-in IFs (fixed weights)
    plug_a     <- sqrt(pmax(diag(cluster_cov_edid(if_agg, panel$cluster_indices, n)), 0))
    pert_a     <- rep(NA_real_, length(base_v))
    npert_a    <- colSums(is.finite(agg_draws))
    for (j in seq_along(base_v)) {
      v <- agg_draws[, j]; v <- v[is.finite(v)]
      if (length(v) >= 2L) pert_a[j] <- stats::sd(v)
    }
    comb_a <- sqrt(plug_a^2 + ifelse(is.na(pert_a), 0, pert_a)^2)
    tab <- data.frame(
      parameter = names(base_v), att = as.numeric(base_v),
      se_plug = plug_a, se_pert = pert_a, se_combined = comb_a, n_pert = npert_a,
      ci_lower = as.numeric(base_v) - z * comb_a,
      ci_upper = as.numeric(base_v) + z * comb_a)
    rownames(tab) <- NULL
    agg_tabs[[a]] <- tab
  }
  agg_tabs <- Filter(Negate(is.null), agg_tabs)

  out <- list(
    att_gt        = att_tab,
    aggregates    = agg_tabs,
    B             = B,
    n_failed      = n_failed,
    n_nuisances   = length(infos),
    weight_scheme = ws,
    alpha         = alpha,
    seed          = seed_base,
    call          = mc
  )
  class(out) <- c("edid_perturbation_bootstrap", "list")
  out
}

#' @describeIn edid_perturbation_bootstrap Print method.
#' @param x an \code{edid_perturbation_bootstrap} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_perturbation_bootstrap <- function(x, digits = 4, ...) {
  cat("\nSieve-coefficient perturbation bootstrap for edid (no nuisance refit)\n")
  cat(sprintf("  B = %d coefficient draws (%d failed) over %d first-step nuisance(s); weights ('%s') held fixed\n",
              x$B, x$n_failed, x$n_nuisances, x$weight_scheme))
  cat(sprintf("  SE: sqrt(se_plug^2 + Var_b(att*)); CI: att +/- z * se_combined at level %.3g\n\n",
              1 - x$alpha))
  fmt <- function(tab) {
    num <- vapply(tab, is.numeric, logical(1L))
    tab[num] <- lapply(tab[num], function(z) signif(z, digits))
    tab
  }
  cat("Group-time ATT(g,t):\n")
  print(fmt(x$att_gt), row.names = FALSE)
  for (a in names(x$aggregates)) {
    cat(sprintf("\nAggregation '%s':\n", a))
    print(fmt(x$aggregates[[a]]), row.names = FALSE)
  }
  invisible(x)
}
