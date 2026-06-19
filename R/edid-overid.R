# edid-overid.R
# Omnibus over-identification J for edid fits: the joint test that ALL admissible
# elementary identifying moments agree, over the cells feeding a reported object.
# The reference df is the RANK of the contrast covariance (typically far below the
# nominal Q - p: the over-identification is intrinsically low-rank / model-level).
# This is the model-level over-identification statistic of
# Andrews, Chen & Tecchio (2025) (the "report the J" object), distinct from the
# 2-leg PT-All-vs-PT-Post Hausman test (edid_hausman, df <= |E|) and the
# one-at-a-time incremental Sargan procedure (edid_sargan). Built on the same
# refit + IF-difference-quadratic-form machinery as those (see edid-sargan.R,
# edid-hausman.R), so it inherits the AHT effective-df F + eigen-ridge size
# correction of .edid_if_diff_quadform().

# Assemble the joint over-identification statistic for a set of target cells.
# `cell_elems` is the per-cell list of elementary estimators (each list(att,
# ifvec, is_base)). For each cell with >= 2 elementary moments we form q-1
# reference-anchored contrasts (elem_j - elem_ref, j != ref); the reference is
# the PT-Post self base pair when present (is_base), else the first elementary.
# The Wald quadratic form n d' D^+ d is invariant to the reference choice (the
# contrast space is fixed), so the statistic tests mutual agreement of all
# elementary moments. The NOMINAL contrast count is sum_cell (q_cell - 1) = Q - p,
# but the reference df is rank(D) <= Q - p -- typically strictly smaller, because
# the contrasts are linearly dependent (shared trend restrictions): the genuine
# low-rank over-id. nominal_df is reported alongside df for transparency.
.edid_overid_assemble <- function(cell_keys, cell_elems, n, ci, v_scale, n_eff, rel_tol = NULL) {
  d <- numeric(0L)
  xi_cols <- vector("list", 0L)
  n_mom <- 0L; n_par <- 0L
  for (key in cell_keys) {
    els <- cell_elems[[key]]
    if (is.null(els) || length(els) < 2L) next          # just-identified cell: no over-id content
    n_mom <- n_mom + length(els)
    n_par <- n_par + 1L
    ref_idx <- which(vapply(els, function(e) isTRUE(e$is_base), logical(1L)))
    ref_idx <- if (length(ref_idx)) ref_idx[1L] else 1L
    ref <- els[[ref_idx]]
    for (j in seq_along(els)) {
      if (j == ref_idx) next
      d <- c(d, els[[j]]$att - ref$att)
      xi_cols[[length(xi_cols) + 1L]] <- els[[j]]$ifvec - ref$ifvec
    }
  }
  if (!length(d)) {
    return(list(statistic = 0, df = 0L, p_value = 1, n_moments = n_mom,
                n_params = n_par, nominal_df = n_mom - n_par, degenerate = TRUE))
  }
  xi <- do.call(cbind, xi_cols)
  # Relative rank floor. Default (NULL) -> "auto" = the EFFECTIVE-RANK floor r_bare/n_eff computed in
  # the engine (keyed to the genuine bare rank, not the nominal stacked dimension): a true no-op for the
  # clean low-rank no-covariate spectrum and the spectral-gap cut for the covariate decaying spectrum.
  # A numeric rel_tol overrides; 0 reproduces the bare numerical cut.
  rt <- if (is.null(rel_tol)) "auto" else rel_tol
  qf <- .edid_if_diff_quadform(d, xi, n, ci, v_scale = v_scale, n_eff = n_eff, rel_tol = rt)
  list(statistic = qf$statistic, df = qf$df, p_value = qf$p_value,
       n_moments = n_mom, n_params = n_par, nominal_df = n_mom - n_par,
       degenerate = isTRUE(qf$degenerate))
}

#' Omnibus over-identification test for edid fits (the joint J)
#'
#' Computes the omnibus over-identification statistic for an efficient edid fit:
#' the joint test that all admissible elementary identifying moments agree, over
#' the cells feeding a reported object. \eqn{Q} is the total number of elementary
#' identifying moments and \eqn{p} the number of over-identified cells, so
#' \eqn{Q - p = \sum_{(g,t)} (q_{g,t} - 1)} is the \emph{nominal} over-identification
#' count. The reference degrees of freedom, however, is the \emph{rank} of the
#' contrast covariance, which for DiD is typically \emph{far below} \eqn{Q - p}
#' (the elementary moments share the never-treated control and the comparison
#' trend restrictions, so the over-identification is intrinsically low-rank --
#' a single model-level object; see \strong{Details}). Using \eqn{\chi^2(Q - p)}
#' instead would severely under-reject.
#'
#' This is the model-level over-identification statistic emphasised by Andrews,
#' Chen & Tecchio (2025) (the "report the J" object): unlike \code{\link{edid_hausman}}
#' (a 2-leg PT-All-vs-PT-Post contrast, \eqn{df \le |\mathcal{E}|}) and
#' \code{\link{edid_sargan}} (one added restriction at a time), it tests the
#' \emph{full} set of over-identifying restrictions jointly.
#'
#' @details
#' For each target cohort \eqn{g} the admissible pairs \eqn{(g', t_{pre})} are
#' enumerated under PT-All (the same enumeration \code{edid()} uses); each single
#' pair is refit as a just-identified estimator via the internal
#' \code{moment_set} mechanism, in the efficient plug-in configuration (all three
#' estimation-effect channels off: the over-identification contrast lives on the
#' efficient inverse-variance covariance, not a misspecification-robust variance;
#' see \code{\link{edid_sargan}}). For each cell \eqn{(g,t)} the elementary
#' estimators are contrasted against a reference (the self-pair \eqn{g'=g} -- the
#' never-treated-comparison moment at the most recent baseline -- when present,
#' else the first elementary), giving \eqn{q_{g,t} - 1} contrasts; the Wald
#' statistic is invariant to the reference. The contrasts are stacked across
#' the cells in scope and the statistic is the IF-difference quadratic form
#' \eqn{n\, d' \widehat{D}^{+} d} with \eqn{df = \mathrm{rank}(\widehat{D})},
#' carrying the AHT effective-df F reference of \code{\link{edid_hausman}}.
#' The per-cell anchoring is valid for the over-identification \emph{test} (the
#' Wald form is invariant to the reference); the per-cell vs full-system
#' distinction matters only for the robustness-frontier identity, not here.
#'
#' Total refits: \eqn{\sum_g q_g} just-identified fits (no bands, no bootstrap).
#'
#' \strong{Validation scope.} Monte Carlo size/power studies (2026-06-19) validated nominal size with
#' strong power on the no-covariate path (i.i.d., AR(1), clustered), the covariate path (nominal size,
#' \code{mean(stat)/df}\eqn{\approx 1}, high power across \eqn{n} and 1--2 covariates), the
#' observation-weighted path, AND genuinely clustered designs with real within-cluster correlation
#' (clustered + covariate and clustered + no-covariate: the floor recovers the true rank, nominal size,
#' power \eqn{\approx 0.95}). The DiD over-identification is intrinsically \emph{low-rank} (the elementary
#' moments share the never-treated time control and the comparison cohorts' trend restrictions), so the
#' effective degrees of freedom is the \emph{rank} of the contrast covariance, \emph{not} the naive
#' moment-count \eqn{Q-p}; the statistic is therefore a single model-level object (per-cell, per-horizon,
#' and overall coincide up to the cells in scope). The covariate-adjusted contrast covariance has a
#' \emph{decaying} eigenvalue spectrum (vs the exact zeros of the no-covariate case), so the rank is set
#' by the relative floor \code{rel_tol} (default \code{"auto"} = \eqn{r_{\mathrm{bare}}/n_{\mathrm{eff}}},
#' keyed to the genuine bare rank). This is a \strong{division of labor} with the AHT F: the floor sets
#' the RANK (which spectral directions are genuine over-id content vs the covariate noise tail; denominator
#' \eqn{n_{\mathrm{eff}}}, the spectral noise scale), while the AHT F handles the few-cluster sampling
#' RELIABILITY of the survivors (\eqn{m = G_{\mathrm{eff}}-1}). The floor is \emph{motivated by} the
#' random-matrix noise scale but justified empirically (it lands in the spectral gap; verified by the
#' per-direction calibration and the size/power MC, including the genuinely clustered designs); it is a
#' true no-op for the clean no-covariate spectrum (even at small \eqn{n_{\mathrm{eff}}}). The enumerated
#' moment set replicates \code{edid()}'s thin-cohort guard so the over-identification is tested over
#' exactly the FITTED moments.
#'
#' \strong{Rank-deficiency / saturation.} A genuinely just-identified design returns \code{NULL} with a
#' message. The joint \eqn{J} is reported as \code{NA} (\code{rank_deficient = TRUE}, never a misleading
#' \eqn{p}) in two regimes, with a fragile-regime \code{warning} routing to the per-cell breakdown
#' (\code{$cells}) / \code{\link{edid_sargan}}: (i) the relative floor drops every direction
#' (\eqn{r_{\mathrm{bare}}\ge n_{\mathrm{eff}}}); and (ii) \strong{cluster-rank saturation} -- a defensive
#' guard for genuinely \emph{few-cluster} fits, when the bare numerical rank reaches the cluster-robust
#' ceiling (\eqn{r_{\mathrm{bare}} = G_{\mathrm{eff}}-1}, \eqn{G_{\mathrm{eff}}} = number of clusters; a
#' centered cluster sandwich has rank \eqn{\le G_{\mathrm{eff}}-1}), so the over-identifying dimension
#' meets/exceeds the cluster budget and the joint is not reliably estimable. This fires ONLY under coarse
#' clustering with a large over-id; it does \emph{not} fire under the \strong{default unit-level
#' clustering} (\eqn{G_{\mathrm{eff}} = n}, thousands of pieces), where the fixed low-rank over-id is
#' comfortably supported. (E.g. Bailey-Goodman-Bacon, clustered at the county = unit level per the original
#' paper, has \eqn{G_{\mathrm{eff}}\approx 3059 \gg} its structural rank 83 and computes a normal joint J,
#' \eqn{df = 15}, \eqn{p\approx 5.6\times 10^{-5}} -- the over-id \emph{rejects}.) \strong{Scope:} like
#' \code{edid()}, this assumes a balanced-panel / fixed-unit structure; it is not validated for repeated
#' cross-sections.
#'
#' @param fit An \code{edid_fit}, normally from
#'   \code{edid(..., pt_assumption = "all")}. Supplies the design (cohorts,
#'   periods, anticipation) and the estimation options for the internal refits.
#' @param data The panel data used to estimate \code{fit}, or \code{NULL}
#'   (default), in which case the data expression in the fit's call is
#'   re-evaluated in the caller's environment (the \code{update()} idiom). Pass
#'   \code{data} explicitly when the original object is no longer reachable.
#' @param parameter Which scopes to report, any of \code{"overall"} (all treated
#'   cells; the \eqn{ES_{avg}} footprint), \code{"event_study"} (one J per
#'   post-treatment horizon \eqn{e}, over cells \eqn{(g, g+e)}), and
#'   \code{"att_gt"} (one J per over-identified cell; a localizer). Default:
#'   \code{c("overall", "event_study")}.
#' @param e_set Numeric vector of post-treatment event times for
#'   \code{"event_study"}, or \code{NULL} (default: all finite post-treatment
#'   horizons present in the fit).
#' @param rel_tol Relative eigenvalue floor for the contrast-covariance rank
#'   determination. \code{NULL} (default) uses \code{"auto"} = the effective-rank
#'   floor \eqn{r_{\mathrm{bare}}/n_{\mathrm{eff}}} (\eqn{r_{\mathrm{bare}}} = the
#'   genuine numerical rank), which recovers the effective over-identification
#'   rank on the covariate path and is a true no-op for the no-covariate spectrum.
#'   A numeric value (e.g. \code{0.01}) overrides it; \code{0} reproduces the bare
#'   numerical-rank cut used by \code{\link{edid_hausman}}/\code{\link{edid_sargan}}.
#'
#' @return An object of class \code{edid_overid}: a list with \code{table} (one
#'   row per requested scope level: \code{parameter}, \code{e},
#'   \code{J_statistic}, \code{df}, \code{p_value}, \code{n_moments},
#'   \code{n_params}, \code{nominal_df}), \code{cells} (per-cell breakdown when
#'   requested or always computed for transparency), \code{n}, \code{clustered},
#'   and \code{rank_deficient} (\code{TRUE} when at least one scope's joint J is
#'   \code{NA} -- the floor dropped all directions or the cluster-rank saturated;
#'   read \code{$cells} / \code{\link{edid_sargan}} there).
#'
#' @references Andrews, I., Chen, J., & Tecchio, O. (2025). The purpose of an
#'   estimator is what it does: Misspecification, estimands, and
#'   over-identification. arXiv:2508.13076. \cr
#'   Chen, X., & Santos, A. (2018). Overidentification in Regular Models.
#'   \emph{Econometrica}, 86(5), 1771-1817. \cr
#'   Hansen, L. P. (1982). Large Sample Properties of Generalized Method of
#'   Moments Estimators. \emph{Econometrica}, 50(4), 1029-1054.
#'
#' @seealso \code{\link{edid_hausman}}, \code{\link{edid_sargan}},
#'   \code{\link{edid_frontier}}
#'
#' @examples
#' \donttest{
#' df <- data.frame(
#'   id   = rep(1:200, each = 5),
#'   time = rep(1:5, 200),
#'   g    = rep(sample(c(3, 4, 5, Inf), 200, replace = TRUE), each = 5)
#' )
#' df$y <- rnorm(200)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
#'   rnorm(nrow(df), 0, 0.5)
#' fit <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
#'             aggregate = "event_study", cband = FALSE)
#' edid_overid(fit, data = df)
#' }
#'
#' @export
edid_overid <- function(fit, data = NULL,
                        parameter = c("overall", "event_study", "att_gt"),
                        e_set = NULL, rel_tol = NULL) {
  parameter <- match.arg(parameter, several.ok = TRUE)
  if (!inherits(fit, "edid_fit")) {
    stop("`fit` must be an `edid_fit` object returned by edid().", call. = FALSE)
  }
  # Both paths are size-validated by MC (2026-06-19): no-covariate (iid / AR(1) / clustered, via the
  # AHT effective-df F) and covariate (the rel_tol = "auto" = r_bare/n_eff effective-rank floor recovers
  # the effective over-id rank; nominal size and strong power across n in {150,300} x {1,2} covariates,
  # and on genuinely clustered designs with within-cluster correlation -- round-3b).
  # A fragile-regime guard (below, after the contrast dimension is known) warns only when the
  # over-id dimension approaches n_eff, where ANY rank determination is unreliable (the project's
  # over-id rank-deficiency regime), not merely because covariates are present.
  if (is.null(data)) {
    data <- tryCatch(as.data.frame(eval(fit$call$data, envir = parent.frame())),
                     error = function(e) NULL)
    if (is.null(data) || !is.data.frame(data) || nrow(data) == 0L) {
      stop("Could not recover the estimation data from the fit's call; pass `data` explicitly.",
           call. = FALSE)
    }
  }

  # att-reproduction safety: confirm the recovered/supplied data reproduces the FIT's point estimates
  # (the hardened .edid_plugin_refit guard), not merely n + unit ids -- a wrong-but-same-shape data set
  # (e.g. a data-symbol collision, or shuffled outcomes) would otherwise yield a silently wrong J. The
  # plug-in cache is keyed on the fit fingerprint + nrow(data), NOT data content, so a cache hit would
  # BYPASS the guard; force a cache MISS for this one validation call so the reproduction check runs.
  .pc <- options(edid_plugin_cache = FALSE)
  invisible(tryCatch(.edid_plugin_refit(fit, data), finally = options(.pc)))

  tg  <- sort(fit$treatment_groups[is.finite(fit$treatment_groups) & fit$treatment_groups != 0])
  tp  <- sort(fit$time_periods); p1 <- min(tp); ant <- fit$anticipation %||% 0L

  # Cohort sizes (units per finite cohort) to replicate edid()'s THIN-COHORT GUARD on the enumeration,
  # so the over-identification moment set matches the FITTED moment set (apply_thin_cohort_guard_edid in
  # .gbuild). Without this, edid_overid would test cross-pairs the fit excised / pin-overridden moments a
  # thin cohort the fit reduced to just-identified, reporting a spurious over-id "pass".
  .idn <- fit$args$idname; .gn <- fit$args$gname
  cohort_sizes <- if (!is.null(.idn) && !is.null(.gn) && all(c(.idn, .gn) %in% names(data))) {
    .u <- !duplicated(data[[.idn]]); .tab <- table(data[[.gn]][.u])
    .cs <- stats::setNames(as.numeric(.tab), names(.tab))
    .cs[is.finite(suppressWarnings(as.numeric(names(.cs))))]
  } else NULL
  mpu <- fit$min_pair_units %||% 5L

  full_pairs <- stats::setNames(lapply(tg, function(g) {
    pr <- enumerate_valid_pairs_edid(target_g = g, treatment_groups = tg, time_periods = tp,
                                     period_1 = p1, pt_assumption = "all", anticipation = ant)
    if (!is.null(cohort_sizes) && nrow(pr) > 0L)
      pr <- apply_thin_cohort_guard_edid(g, pr, cohort_sizes, mpu, "all")$pairs
    pr
  }), as.character(tg))

  base_pair <- stats::setNames(lapply(tg, function(g) {
    pg <- full_pairs[[as.character(g)]]
    self <- pg[is.finite(pg$gp) & pg$gp == g, , drop = FALSE]
    if (nrow(self) == 0L) return(NULL)
    list(gp = g, tpre = max(self$tpre))
  }), as.character(tg))

  caller_env <- parent.frame()
  n   <- fit$n
  ci  <- fit$cluster_indices
  n_eff <- .edid_overid_n_eff(fit)
  # Absolute variance scale of the cell estimators, for the degenerate-contrast guard.
  v_scale <- {
    se <- fit$att_gt$se
    if (is.null(se) || !any(is.finite(se))) 1 else n * max(se[is.finite(se)]^2)
  }

  ck <- function(g, t) paste0(g, "|", t)
  cell_elems <- list()
  checked_sample <- FALSE

  for (g in tg) {
    pg <- full_pairs[[as.character(g)]]
    if (nrow(pg) == 0L) next
    bp <- base_pair[[as.character(g)]]
    for (r in seq_len(nrow(pg))) {
      gp_r <- pg$gp[r]; tp_r <- pg$tpre[r]
      ms <- data.frame(g = g, gp = gp_r, tpre = tp_r)
      rf <- suppressWarnings(.edid_refit_moment_set(fit, data, ms, envir = caller_env))
      if (!checked_sample) {
        if (!identical(rf$n, n) || !identical(rf$all_units, fit$all_units)) {
          stop("The data used for the refits does not match the fitted sample (n or unit ids ",
               "differ from `fit`); pass the original estimation data via `data`.", call. = FALSE)
        }
        checked_sample <- TRUE
      }
      agt <- rf$att_gt; eif <- rf$eif
      if (is.null(agt) || is.null(eif)) next
      is_base <- !is.null(bp) && gp_r == bp$gp && tp_r == bp$tpre
      rows <- which(agt$group == g & is.finite(agt$att))
      for (rr in rows) {
        ifvec <- eif[, rr]
        if (any(!is.finite(ifvec))) next
        key <- ck(g, agt$time[rr])
        cell_elems[[key]] <- c(cell_elems[[key]],
          list(list(g = g, t = agt$time[rr], att = agt$att[rr],
                    ifvec = ifvec, is_base = is_base)))
      }
    }
  }

  # Catalogue all cell keys actually observed, with their (g, t, e).
  keys <- names(cell_elems)
  # No over-identifying content if there are no cells OR every cell has a single elementary moment
  # (just-identified -- e.g. PT-Post, a single-pre-period design, or after the thin-cohort guard pinned
  # every cohort). Return NULL + message rather than a misleading all-zero / p = 1 table.
  if (!length(keys) || all(vapply(cell_elems, length, integer(1L)) < 2L)) {
    message("No over-identifying content: every contributing cell is just-identified (Q = p).")
    return(invisible(NULL))
  }
  meta <- do.call(rbind, lapply(keys, function(k) {
    e1 <- cell_elems[[k]][[1L]]
    data.frame(key = k, g = e1$g, t = e1$t, e = e1$t - e1$g, q = length(cell_elems[[k]]),
               stringsAsFactors = FALSE)
  }))

  # ----- per-cell breakdown (always computed; cheap, and the localizer) -----
  cell_rows <- vector("list", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    res <- .edid_overid_assemble(meta$key[i], cell_elems, n, ci, v_scale, n_eff, rel_tol = rel_tol)
    cell_rows[[i]] <- data.frame(
      g = meta$g[i], t = meta$t[i], e = meta$e[i], q = meta$q[i],
      J_statistic = res$statistic, df = res$df, p_value = res$p_value,
      stringsAsFactors = FALSE)
  }
  cells <- do.call(rbind, cell_rows)
  cells <- cells[order(cells$g, cells$t), , drop = FALSE]
  rownames(cells) <- NULL

  # ----- scoped statistics -----
  rows <- list()
  if ("overall" %in% parameter) {
    # ES_avg footprint: post-treatment cells (e >= 0); pre-treatment cells are placebo.
    keys_post <- meta$key[meta$e >= 0]
    res <- .edid_overid_assemble(keys_post, cell_elems, n, ci, v_scale, n_eff, rel_tol = rel_tol)
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = "overall", e = NA_real_,
      J_statistic = res$statistic, df = res$df, p_value = res$p_value,
      n_moments = res$n_moments, n_params = res$n_params, nominal_df = res$nominal_df,
      stringsAsFactors = FALSE)
  }
  if ("event_study" %in% parameter) {
    es <- sort(unique(meta$e[meta$e >= 0]))
    if (!is.null(e_set)) es <- intersect(es, e_set)
    for (e in es) {
      keys_e <- meta$key[meta$e == e]
      res <- .edid_overid_assemble(keys_e, cell_elems, n, ci, v_scale, n_eff, rel_tol = rel_tol)
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = "event_study", e = e,
        J_statistic = res$statistic, df = res$df, p_value = res$p_value,
        n_moments = res$n_moments, n_params = res$n_params, nominal_df = res$nominal_df,
        stringsAsFactors = FALSE)
    }
  }
  table <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(table)) rownames(table) <- NULL

  # Fragile-regime warning, keyed to the EFFECTIVE rank (not the nominal Q-p, which over-fires on benign
  # high-redundancy designs -- the low-rank thesis). Two cases: (i) a scope is rank-deficient (NA J --
  # the engine could not compute the joint statistic, the genuinely extreme regime); (ii) the effective
  # over-id rank consumes a large share of n_eff (df/n_eff > 0.25), where the joint reference is fragile.
  # In both, prefer the localized reads ($cells / edid_sargan).
  if (!is.null(table)) {
    rd <- any(is.na(table$J_statistic))
    hi <- is.finite(n_eff) && n_eff > 0 && any(is.finite(table$df) & table$df > 0.25 * n_eff)
    if (rd) {
      warning("edid_overid: the joint statistic is rank-deficient / uncomputable for at least one scope ",
              "(NA reported) -- the over-identifying dimension exceeds what the effective sample resolves. ",
              "Read the per-cell breakdown ($cells) or edid_sargan instead of the joint J.", call. = FALSE)
    } else if (hi) {
      warning("edid_overid: the effective over-identification rank consumes a large share of the effective ",
              "sample size (df/n_eff > 0.25); the joint chi-square reference is fragile here. Prefer the ",
              "per-cell breakdown ($cells) or edid_sargan.", call. = FALSE)
    }
  }

  # rank_deficient: at least one scope's joint J is uncomputable (NA) -- the over-identifying dimension
  # exceeds what the effective sample / cluster budget resolves (the saturation / extreme regime). Keyed
  # to NA J, the same signal as the fragile warning above; lets callers (e.g. edid_frontier) and tests
  # detect the regime programmatically and route to the per-cell breakdown / edid_sargan.
  rank_deficient <- !is.null(table) && any(is.na(table$J_statistic))
  out <- list(table = table, cells = cells, n = n, clustered = !is.null(ci),
              parameter = parameter, e_set = e_set, rank_deficient = rank_deficient)
  class(out) <- c("edid_overid", "list")
  out
}

#' @describeIn edid_overid Print method.
#' @param x an \code{edid_overid} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_overid <- function(x, digits = 4, ...) {
  cat("\nOmnibus over-identification test (joint J; df = rank of contrast covariance)\n")
  cat("(Andrews, Chen & Tecchio 2025; the model-level 'report the J' object)\n")
  cat(sprintf("  All admissible elementary moments tested jointly%s\n",
              if (isTRUE(x$clustered)) "; cluster-robust" else ""))
  cat("  Refit convention: efficient plug-in influence functions (over-id contrast on the\n")
  cat("  efficient inverse-variance covariance); df = rank(D-hat), AHT effective-df F.\n\n")
  if (!is.null(x$table)) {
    tab <- x$table
    num <- vapply(tab, is.numeric, logical(1L))
    tab[num] <- lapply(tab[num], function(z) signif(z, digits))
    print(tab, row.names = FALSE)
  } else {
    cat("  (no scoped statistics requested)\n")
  }
  invisible(x)
}
