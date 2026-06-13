# edid-methods.R
# S3 methods for the edid_fit class. The aggregations ($overall/$event_study/$group/$calendar/$simple)
# are did::AGGTEobj objects (built via aggte_edid -> did::aggte), so these methods read AGGTEobj fields
# (overall.att/overall.se, att.egt/se.egt/egt) and the per-element influence functions in $inf.function.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Pull the per-element (n x K) and overall (n) influence functions out of a did::AGGTEobj, regardless of
# aggregation type (compute.aggte names them dynamic./selective./calendar.inf.func.* and simple.att).
.edid_agg_if <- function(a) {
  if (is.null(a)) return(list(egt = NULL, overall = NULL))
  inf <- a$inf.function
  list(
    egt     = inf$dynamic.inf.func.e %||% inf$selective.inf.func.g %||% inf$calendar.inf.func.t,
    overall = inf$dynamic.inf.func %||% inf$selective.inf.func %||% inf$calendar.inf.func %||% inf$simple.att
  )
}

.edid_sigma_quad <- function(fit) {
  if (!isTRUE(fit$higher_order)) return(NULL)
  K <- nrow(fit$att_gt)
  Sig <- fit$sigma_quad
  if (is.null(Sig) && !is.null(fit$cells)) {
    Sig <- sigma_quad_edid(fit$cells, fit$cluster_indices, fit$n)
  }
  if (is.null(Sig)) return(NULL)
  Sig <- as.matrix(Sig)
  if (length(dim(Sig)) != 2L || any(dim(Sig) != c(K, K))) return(NULL)
  Sig
}

.edid_agg_na_rm <- function(a) {
  if (!is.null(a$call) && !is.null(a$call$na.rm)) {
    val <- tryCatch(eval(a$call$na.rm), error = function(e) NULL)
    if (is.logical(val) && length(val) == 1L && !is.na(val)) return(val)
  }
  TRUE
}

.edid_agg_reaggregate <- function(fit, a) {
  type <- a$type %||% NULL
  if (is.null(type) || !(type %in% c("simple", "dynamic", "group", "calendar"))) return(NULL)
  balance_e <- a$balance_e %||% NULL
  min_e     <- a$min_e %||% -Inf
  max_e     <- a$max_e %||% Inf
  na.rm     <- .edid_agg_na_rm(a)
  function(att_vec) {
    f2 <- fit
    f2$att_gt$att <- att_vec
    aa <- aggte(as_MP_edid(f2, bstrap = FALSE, cband = FALSE), type = type, balance_e = balance_e,
                min_e = min_e, max_e = max_e, na.rm = na.rm, bstrap = FALSE)
    aa$att.egt %||% aa$overall.att
  }
}

.edid_agg_higher_order_cov <- function(a, fit) {
  if (!isTRUE(fit$higher_order) || is.null(a)) return(NULL)
  g <- .edid_agg_if(a)
  if (is.null(g$egt) && is.null(g$overall)) return(NULL)
  Sig_quad <- .edid_sigma_quad(fit)
  reaggregate <- .edid_agg_reaggregate(fit, a)
  if (is.null(Sig_quad) || is.null(reaggregate)) return(NULL)
  A <- .edid_recover_agg_map(a, fit, reaggregate)
  if (is.null(A) || ncol(A) != nrow(Sig_quad)) return(NULL)
  A %*% Sig_quad %*% t(A)
}

# ---------------------------------------------------------------------------
# Stability diagnostics (the $diagnostics field of an edid_fit)
# ---------------------------------------------------------------------------

# Net / gross cross-cohort hedge mass over the POST cells, the cheap red flag of a
# poisoned over-identified fit (Nguyen / Bailey-GB / ACA gate evidence). For each post
# cell the cross-cohort pairs are the finite-gp pairs with gp != group (the never-treated
# anchor gp = Inf and the own-cohort self pairs are excluded); gross = mean over cells of
# sum|w_cross|, net = mean over cells of |sum w_cross|. A healthy efficient fit hedges
# (gross negative mass offsets gross positive, so net << gross); a broken fit has the
# cross-cohort "hedges" carrying the estimand (net ~= gross, gross negative mass ~ 0).
# Returns NA when no post cell carries a cross-cohort pair (no over-identification to
# hedge -- e.g. PT-Post, single-cohort, or moment_set = "own").
.edid_net_hedge_mass <- function(cells) {
  if (is.null(cells) || length(cells) == 0L) return(list(net = NA_real_, gross = NA_real_, n_cells = 0L))
  net <- numeric(0L); gross <- numeric(0L)
  for (cc in cells) {
    if (isTRUE(cc$is_pre)) next                       # post cells only
    w  <- cc$weights; pr <- cc$pairs
    if (is.null(w) || length(w) == 0L || is.null(pr) || nrow(pr) != length(w)) next
    cross <- is.finite(pr$gp) & pr$gp != cc$group     # cross-cohort control-variate pairs
    if (!any(cross)) next
    wc <- as.numeric(w[cross])
    net   <- c(net,   abs(sum(wc)))
    gross <- c(gross, sum(abs(wc)))
  }
  if (length(net) == 0L) return(list(net = NA_real_, gross = NA_real_, n_cells = 0L))
  list(net = mean(net), gross = mean(gross), n_cells = length(net))
}

# Assemble the $diagnostics object from the raw per-fit counts (threaded out of
# fit_edid_cells) plus the completed cells. Pure read-out: no moment, weight, or estimate
# is touched. The booleans are the toolkit's machine-readable red flags; `$unstable` is
# the single summary the broken-leg Hausman guard keys on.
.edid_build_diagnostics <- function(raw, cells, pt_assumption, weight_scheme, min_pair_units) {
  if (is.null(raw)) raw <- list()
  hedge <- .edid_net_hedge_mass(cells)
  n_extreme <- as.integer(raw$n_extreme_ratio %||% 0L)
  n_psi     <- as.integer(raw$n_psi_unstable %||% 0L)
  n_drop    <- as.integer(raw$n_pairs_dropped %||% 0L)
  n_full    <- as.integer(raw$n_fulltrim %||% 0L)
  # Over-identified efficient covariate fit whose cross-cohort hedges carry the estimand:
  # net hedge mass at/above the calibrated flag AND essentially equal to gross (no
  # offsetting negative mass). Only meaningful where hedging is possible (n_cells > 0).
  net_flag <- isTRUE(is.finite(hedge$net) && hedge$n_cells > 0L &&
                     hedge$net >= EDID_NET_HEDGE_FLAG &&
                     hedge$net >= 0.95 * (hedge$gross %||% Inf))
  # "Unstable leg": the conditions under which a non-rejection / point estimate from this
  # fit is not trustworthy -- extreme propensity ratios entered, the weight channel was
  # not a credible IF, or the cross-cohort hedges carry the estimand. (Dead pairs / full
  # trims alone redefine the estimand but are reported separately; they do not by
  # themselves flag the fit as numerically broken.)
  unstable <- (n_extreme > 0L) || (n_psi > 0L) || net_flag
  list(
    n_extreme_ratio      = n_extreme,
    n_psi_unstable       = n_psi,
    n_pairs_dropped      = n_drop,
    n_fulltrim           = n_full,
    net_hedge_mass       = hedge$net,
    gross_hedge_mass     = hedge$gross,
    net_hedge_flag       = net_flag,
    min_finite_cohort    = raw$min_finite_cohort %||% NA_integer_,
    small_cohorts        = raw$small_cohorts,         # finite cohorts in [min_pair_units, comfort); NULL if none
    cohort_sizes         = raw$cohort_sizes,
    use_cov_path         = isTRUE(raw$use_cov_path),
    unstable             = isTRUE(unstable)
  )
}

#' Print method for edid_fit objects
#'
#' Displays the ATT(g,t) table in the same style as \code{print.MP} / \code{summary.MP}, followed by
#' footer metadata.
#'
#' @param x an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#' @return \code{x} invisibly
#' @export
print.edid_fit <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Group-Time Average Treatment Effects:\n")

  alp <- x$alpha
  cband_text1a <- paste0(100 * (1 - alp), "% ")
  # Simultaneous iff a > qnorm crit was actually applied: cband = TRUE with either the multiplier
  # bootstrap or the analytic sup-t construction. bstrap = TRUE with cband = FALSE keeps pointwise CIs.
  simult <- isTRUE(x$cband) &&
            ((isTRUE(x$bstrap) && identical(x$cband_method, "multiplier")) ||
             identical(x$cband_method, "analytic"))
  cband_text1b <- ifelse(simult, "Simult. ", "Pointwise ")
  cband_text1  <- paste0("[", cband_text1a, cband_text1b)

  att_df <- x$att_gt
  if (!is.null(att_df) && nrow(att_df) > 0L) {
    ci_lower <- att_df$ci_lower
    ci_upper <- att_df$ci_upper
    sig <- (ci_upper < 0) | (ci_lower > 0)
    sig[is.na(sig)] <- FALSE
    sig_text <- ifelse(sig, "*", "")
    out <- cbind.data.frame(att_df$group, att_df$time, att_df$att, att_df$se, ci_lower, ci_upper)
    out <- round(out, 4)
    out <- cbind.data.frame(out, sig_text)
    colnames(out) <- c("Group", "Time", "ATT(g,t)", "Std. Error", cband_text1, "Conf. Band]", "")
    print(out, row.names = FALSE)
  } else {
    cat("  (no cells)\n")
  }

  cat("---\n")
  cat("Signif. codes: `*' confidence band does not cover 0")
  cat("\n\n")

  cat("Control Group:  "); cat("Never Treated"); cat(",  ")
  cat("Anticipation Periods:  "); cat(x$anticipation); cat("\n")
  cat("Estimation Method:  Efficient DiD (Chen, Sant'Anna & Xie 2025)\n")
  pt_text <- if (x$pt_assumption == "all") "PT-All" else "PT-Post"
  cat("PT Assumption:  "); cat(pt_text); cat("\n")
  .edid_thin_radar_note(x)
  invisible(x)
}

# Thin-cohort radar note (fix 2): printed -- not warned -- when the over-identified
# covariate fit has finite cohorts in the [min_pair_units, comfort) band, where the
# efficient-weight machinery can be unreliable (the Nguyen 14/33-unit gate failure). The
# data is in $diagnostics$small_cohorts; this only formats it for print/summary. Shown
# only on the covariate PT-All path (moot otherwise); no-op when there are no such cohorts.
.edid_thin_radar_note <- function(x) {
  d <- x$diagnostics
  if (is.null(d) || is.null(d$small_cohorts) || nrow(d$small_cohorts) == 0L) return(invisible())
  has_cov <- !is.null(x$xformla) && inherits(x$xformla, "formula") && length(all.vars(x$xformla)) > 0L
  if (!identical(x$pt_assumption, "all") || !has_cov) return(invisible())
  sc <- d$small_cohorts
  cat(sprintf(paste0("Thin-cohort radar:  cohort(s) %s (%d-%d units) are above the hard guard but ",
                     "below a comfortable size;\n  on this covariate PT-All path the over-identified ",
                     "efficient weights can be unreliable for cohorts this small\n  (see ",
                     "$diagnostics; consider weight_scheme = \"averaged\" or moment_set = \"own\").\n"),
              paste(format(sc$cohort, trim = TRUE, scientific = FALSE), collapse = ", "),
              min(sc$n_units), max(sc$n_units)))
  invisible()
}

#' Summary method for edid_fit objects
#'
#' Prints the ATT(g,t) table (MP style) followed by the requested aggregations, each a
#' \code{did::AGGTEobj} printed with did's own \code{print.AGGTEobj}.
#'
#' @param object an \code{edid_fit} object
#' @param ... additional arguments (currently ignored)
#' @return \code{object} invisibly
#' @export
summary.edid_fit <- function(object, ...) {
  print.edid_fit(object, ...)
  ov <- object$overall
  if (!is.null(ov)) {
    cat(sprintf("\nOverall ATT (%s):  %s  (SE %s)\n", ov$type %||% "overall",
                .fmt_or_na(ov$overall.att), .fmt_or_na(ov$overall.se)))
  }
  for (nm in c("event_study", "group", "calendar")) {
    a <- object[[nm]]
    if (!is.null(a) && inherits(a, "AGGTEobj")) {
      cat(sprintf("\n--- %s ---\n", nm)); print(a)
    }
  }
  cat("\n")
  invisible(object)
}

# Internal formatting helper
.fmt_or_na <- function(x) {
  if (is.null(x) || !is.finite(x)) "NA" else sprintf("%.4f", x)
}

#' Extract ATT coefficients from an edid_fit object
#'
#' @param object an \code{edid_fit} object
#' @param which character: one of \code{"att_gt"}, \code{"overall"}, \code{"event_study"}, \code{"group"}
#' @param ... additional arguments (ignored)
#' @return named numeric vector of ATT estimates
#' @export
coef.edid_fit <- function(object, which = c("att_gt", "overall", "event_study", "group"), ...) {
  which <- match.arg(which)
  switch(which,
    att_gt = {
      df <- object$att_gt
      stats::setNames(df$att, paste0("ATT(", df$group, ",", df$time, ")"))
    },
    overall = {
      if (is.null(object$overall)) return(numeric(0L))
      c(overall = object$overall$overall.att)
    },
    event_study = {
      a <- object$event_study
      if (is.null(a)) return(numeric(0L))
      stats::setNames(a$att.egt, paste0("e=", a$egt))
    },
    group = {
      a <- object$group
      if (is.null(a)) return(numeric(0L))
      stats::setNames(a$att.egt, paste0("g=", a$egt))
    }
  )
}

#' Extract variance-covariance matrix from an edid_fit object
#'
#' For \code{which = "att_gt"} returns the cluster-robust (or i.i.d.) covariance of the cell-level
#' ATT(g,t)'s from the stored influence functions. For the aggregations it returns the covariance implied
#' by the corresponding \code{did::AGGTEobj}'s aggregate influence functions (cluster-robust when
#' \code{clustervars} was set).
#'
#' @param object an \code{edid_fit} object
#' @param which character: one of \code{"att_gt"}, \code{"overall"}, \code{"event_study"}, \code{"group"}
#' @param ... additional arguments (ignored)
#' @return square numeric matrix
#' @export
vcov.edid_fit <- function(object, which = c("att_gt", "overall", "event_study", "group"), ...) {
  which <- match.arg(which)
  n  <- object$n
  ci <- object$cluster_indices

  # cluster-robust (or i.i.d.) covariance of the columns of an n x K influence-function matrix
  .cross_mat <- function(M) {
    M <- as.matrix(M)
    if (is.null(ci)) return(crossprod(M) / n^2)
    G <- length(unique(ci))
    if (G <= 1L) return(matrix(NA_real_, ncol(M), ncol(M)))
    CS <- rowsum(M, ci)
    (G / (G - 1)) * crossprod(CS) / n^2
  }

  if (which == "att_gt") {
    df  <- object$att_gt
    nms <- paste0("ATT(", df$group, ",", df$time, ")")
    if (!is.null(object$eif)) {
      v <- .cross_mat(object$eif)
      Sig_quad <- .edid_sigma_quad(object)
      if (!is.null(Sig_quad) && all(dim(Sig_quad) == dim(v))) v <- v + Sig_quad
      dimnames(v) <- list(nms, nms); return(v)
    }
    v <- diag(df$se^2, nrow = nrow(df)); dimnames(v) <- list(nms, nms); return(v)
  }

  if (which == "overall") {
    g <- .edid_agg_if(object$overall)
    if (is.null(g$overall)) return(matrix(NA_real_, 1L, 1L))
    v <- matrix(.cross_mat(matrix(g$overall, ncol = 1L)), 1L, 1L,
                dimnames = list("overall", "overall"))
    HO <- .edid_agg_higher_order_cov(object$overall, object)
    if (!is.null(HO) && is.null(g$egt) && all(dim(HO) == c(1L, 1L))) {
      v[1L, 1L] <- v[1L, 1L] + HO[1L, 1L]
    } else if (!is.null(HO) && !is.null(g$egt) && is.matrix(g$egt) && ncol(g$egt) == nrow(HO)) {
      w <- tryCatch(drop(solve(crossprod(g$egt), crossprod(g$egt, g$overall))),
                    error = function(e) NULL)
      if (!is.null(w) && length(w) == nrow(HO)) v[1L, 1L] <- v[1L, 1L] + drop(crossprod(w, HO %*% w))
    }
    return(v)
  }

  a   <- if (which == "event_study") object$event_study else object$group
  g   <- .edid_agg_if(a)
  if (is.null(a) || is.null(g$egt)) return(matrix(NA_real_, 0L, 0L))
  pre <- if (which == "event_study") "e=" else "g="
  nms <- paste0(pre, a$egt)
  v   <- .cross_mat(g$egt)
  HO  <- .edid_agg_higher_order_cov(a, object)
  if (!is.null(HO) && all(dim(HO) == dim(v))) v <- v + HO
  dimnames(v) <- list(nms, nms); v
}

#' Coerce edid_fit to a data.frame
#'
#' @param x an \code{edid_fit} object
#' @param row.names ignored; included for S3 generic consistency
#' @param optional ignored; included for S3 generic consistency
#' @param ... not used
#' @param which character: one of \code{"att_gt"}, \code{"overall"}, \code{"event_study"}, \code{"group"}
#' @return data.frame
#' @export
as.data.frame.edid_fit <- function(x, row.names = NULL, optional = FALSE, ...,
                                    which = c("att_gt", "overall", "event_study", "group")) {
  which <- match.arg(which)
  z <- stats::qnorm(1 - x$alpha / 2)
  # Per-element aggregation CIs reuse the stored crit (sup-t when cband = TRUE), so the
  # column semantics match `which = "att_gt"` (which returns the stored bands) instead of
  # silently downgrading the aggregations to pointwise.
  .agg_crit <- function(a) {
    cv <- a$crit.val.egt
    if (length(cv) == 1L && is.finite(cv)) cv else z
  }
  switch(which,
    att_gt = x$att_gt,
    overall = {
      a <- x$overall
      if (is.null(a)) return(data.frame(att = numeric(0L), se = numeric(0L),
                                        ci_lower = numeric(0L), ci_upper = numeric(0L)))
      data.frame(att = a$overall.att, se = a$overall.se,
                 ci_lower = a$overall.att - z * a$overall.se,
                 ci_upper = a$overall.att + z * a$overall.se, stringsAsFactors = FALSE)
    },
    event_study = {
      a <- x$event_study
      if (is.null(a)) return(data.frame(e = numeric(0L), att = numeric(0L), se = numeric(0L),
                                        ci_lower = numeric(0L), ci_upper = numeric(0L)))
      cv <- .agg_crit(a)
      data.frame(e = a$egt, att = a$att.egt, se = a$se.egt,
                 ci_lower = a$att.egt - cv * a$se.egt, ci_upper = a$att.egt + cv * a$se.egt,
                 stringsAsFactors = FALSE)
    },
    group = {
      a <- x$group
      if (is.null(a)) return(data.frame(group = numeric(0L), att = numeric(0L), se = numeric(0L),
                                        ci_lower = numeric(0L), ci_upper = numeric(0L)))
      cv <- .agg_crit(a)
      data.frame(group = a$egt, att = a$att.egt, se = a$se.egt,
                 ci_lower = a$att.egt - cv * a$se.egt, ci_upper = a$att.egt + cv * a$se.egt,
                 stringsAsFactors = FALSE)
    }
  )
}
