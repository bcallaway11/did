# edid-weights.R
# Weight-decomposition accessor + heatmap for edid() fits: the paper's signature
# diagnostic for how each (g', t_pre) identifying moment is leveraged per cell.

utils::globalVariables(c("tpre_f", "gp_f", "weight", "cell_lab"))

#' Extract the per-pair efficiency weights from an \code{edid} fit
#'
#' Returns the weight that each identifying DiD moment --- a comparison-cohort /
#' pre-treatment-period pair \eqn{(g', t_{pre})} --- receives in each
#' \eqn{ATT(g,t)} cell, in tidy (long) form. This is the weight-decomposition
#' diagnostic of Chen, Sant'Anna & Xie (2025): the efficient estimand weights
#' the generated outcomes by
#' \eqn{w(X) = \Omega_{gt}^*(X)^{-1}\mathbf 1 / (\mathbf 1'\Omega_{gt}^*(X)^{-1}\mathbf 1)},
#' and the paper recommends \emph{plotting the expected value of these weights}
#' "so that one can have a better understanding of how each pre-treatment
#' period and comparison group is leveraged for efficiency considerations"
#' (Section 4; the weight heatmaps in the paper's simulations and empirical
#' application are exactly this object). See \code{\link{edid_weight_plot}} for
#' the companion heatmap.
#'
#' @section What the weight column contains:
#' For \code{weight_scheme = "efficient"} (the default) on the covariate path,
#' the reported weight for pair \eqn{j} is the \emph{mean pointwise} weight
#' \eqn{\mathbb{E}_n[\hat w_j(X_i)]} --- the sample analogue of the expected
#' weight the paper recommends plotting (the per-unit weights \eqn{\hat w(X_i)}
#' themselves vary with \eqn{X_i} and are not stored). For the constant-weight
#' schemes (\code{"averaged"}, \code{"gmm"}, \code{"uniform"}) and for the
#' no-covariate path, the reported weight \emph{is} the constant weight vector
#' used by the estimator. In every case the weights of a cell sum to one (up to
#' floating point), since each per-unit weight vector sums to one by
#' construction.
#'
#' @section Negative weights are legitimate:
#' Do not be alarmed by negative entries. As the paper's Remark on negative
#' weights (their \code{rem:negative-weights}) explains, any covariate-specific
#' weights summing to one identify \eqn{ATT(g,t)} here: the DiD model is
#' \emph{overidentified} and the conditional ATT is \emph{homogeneous} across
#' baseline periods and comparison groups, so --- unlike the negative-weight
#' pathologies of two-way fixed-effects estimators --- non-convex weights do not
#' threaten the causal interpretation. Negative weights simply indicate that a
#' moment is being used to difference out noise correlated with other moments.
#'
#' @param fit An \code{edid_fit} object returned by \code{\link{edid}}.
#' @param type Character: the level at which weights are reported. Currently
#'   only \code{"att_gt"} (one row per cell-pair) is available.
#'
#' @return A data.frame with one row per (cell, pair):
#'   \describe{
#'     \item{\code{group}, \code{time}}{the \eqn{ATT(g,t)} cell.}
#'     \item{\code{gp}}{comparison cohort \eqn{g'} (\code{Inf} = never-treated).}
#'     \item{\code{tpre}}{pre-treatment baseline period \eqn{t_{pre}}.}
#'     \item{\code{weight}}{the pair's weight (see the section above).}
#'     \item{\code{n_pairs}}{number of pairs in the cell (repeated per row).}
#'     \item{\code{condition_num}}{condition number of the cell's (averaged)
#'       \eqn{\Omega^*} (repeated per row); may be \code{NA} when its
#'       computation was skipped or failed (e.g. cheap paths for
#'       \code{"uniform"}/\code{"gmm"} weights, or a degenerate covariance).}
#'   }
#'   Cells with no stored weights (no valid pairs, or unidentified after full
#'   overlap trimming) are excluded; they are recorded in the attribute
#'   \code{"na_cells"} (a data.frame with columns \code{group}, \code{time}).
#'
#' @seealso \code{\link{edid_weight_plot}}, \code{\link{edid}}.
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025).
#'   \emph{Efficient Difference-in-Differences and Event Study Estimators}.
#'   Working paper.
#'
#' @examples
#' set.seed(7)
#' n <- 60; Tt <- 5
#' df <- data.frame(id = rep(1:n, each = Tt), time = rep(1:Tt, n))
#' df$g <- rep(sample(c(3, 4, Inf), n, replace = TRUE), each = Tt)
#' df$y <- rnorm(n)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
#'   rnorm(n * Tt, 0, 0.5)
#' fit <- edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE)
#' w <- edid_weights(fit)
#' head(w)
#' # weights sum to one within each cell:
#' tapply(w$weight, paste(w$group, w$time), sum)
#'
#' @export
edid_weights <- function(fit, type = c("att_gt")) {
  if (!inherits(fit, "edid_fit")) {
    stop("`fit` must be an `edid_fit` object returned by edid().", call. = FALSE)
  }
  type  <- match.arg(type)
  cells <- fit$cells
  if (is.null(cells) || length(cells) == 0L) {
    stop("This fit carries no stored cells ($cells); the per-pair weights are not available. ",
         "Refit with edid(), which stores them by default.", call. = FALSE)
  }

  rows     <- vector("list", length(cells))
  na_group <- numeric(0L)
  na_time  <- numeric(0L)

  for (k in seq_along(cells)) {
    cc <- cells[[k]]
    if (is.null(cc$weights) || length(cc$weights) == 0L) {
      # NA cell: no valid pairs, or unidentified at this trim_level (weights never formed).
      na_group <- c(na_group, cc$group)
      na_time  <- c(na_time, cc$time)
      next
    }
    pr <- cc$pairs
    if (is.null(pr) || nrow(pr) != length(cc$weights)) {
      stop(sprintf(paste0(
        "Cell (g=%g, t=%g) does not carry its (gp, tpre) pair keys; the fit predates the ",
        "weight-labeling support. Refit with the current edid() to use edid_weights()."),
        cc$group, cc$time), call. = FALSE)
    }
    rows[[k]] <- data.frame(
      group         = cc$group,
      time          = cc$time,
      gp            = pr$gp,
      tpre          = pr$tpre,
      weight        = as.numeric(cc$weights),
      n_pairs       = cc$n_pairs,
      condition_num = if (is.null(cc$condition_num)) NA_real_ else cc$condition_num,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(group = numeric(0L), time = numeric(0L), gp = numeric(0L),
                      tpre = numeric(0L), weight = numeric(0L), n_pairs = integer(0L),
                      condition_num = numeric(0L), stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  attr(out, "na_cells") <- data.frame(group = na_group, time = na_time,
                                      stringsAsFactors = FALSE)
  out
}

#' Heatmap of the per-pair efficiency weights (paper-style weight decomposition)
#'
#' Plots the weight that each identifying \eqn{(g', t_{pre})} moment receives in
#' each \eqn{ATT(g,t)} cell as a heatmap, in the style of the weight-decomposition
#' figures of Chen, Sant'Anna & Xie (2025): pre-treatment baseline period on the
#' horizontal axis, comparison cohort \eqn{g'} on the vertical axis, one facet per
#' \eqn{(g,t)} cell, and a diverging fill centered at zero (with symmetric limits)
#' so that legitimately negative weights are immediately visible (see
#' \code{\link{edid_weights}} on why negative weights are not a concern here).
#' For \code{weight_scheme = "efficient"} the fill is the mean pointwise weight
#' \eqn{\mathbb{E}_n[\hat w(X_i)]} --- the expected-weight object the paper
#' recommends plotting; for the constant schemes it is the constant weight.
#'
#' @param fit An \code{edid_fit} object returned by \code{\link{edid}}.
#' @param cells \code{NULL} (default: plot every cell with stored weights), or a
#'   data.frame/matrix with columns \code{group} and \code{time} (or two unnamed
#'   columns in that order) selecting the cells to plot.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{edid_weights}} for the underlying tidy data.
#'
#' @examples
#' set.seed(7)
#' n <- 60; Tt <- 5
#' df <- data.frame(id = rep(1:n, each = Tt), time = rep(1:Tt, n))
#' df$g <- rep(sample(c(3, 4, Inf), n, replace = TRUE), each = Tt)
#' df$y <- rnorm(n)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
#'   rnorm(n * Tt, 0, 0.5)
#' fit <- edid(df, "y", "id", "time", "g", aggregate = "none", cband = FALSE)
#' edid_weight_plot(fit)
#' # a single cell:
#' edid_weight_plot(fit, cells = data.frame(group = 3, time = 4))
#'
#' @export
edid_weight_plot <- function(fit, cells = NULL) {
  w <- edid_weights(fit)   # validates `fit` and errors informatively when cells are missing
  if (nrow(w) == 0L) {
    stop("No cell in this fit has stored weights (all cells are NA); nothing to plot.",
         call. = FALSE)
  }

  if (!is.null(cells)) {
    cells <- as.data.frame(cells)
    if (!all(c("group", "time") %in% names(cells))) {
      if (ncol(cells) >= 2L) names(cells)[1:2] <- c("group", "time")
      else stop("`cells` must have columns `group` and `time` (or two columns in that order).",
                call. = FALSE)
    }
    keep <- paste(w$group, w$time) %in% paste(cells$group, cells$time)
    if (!any(keep)) {
      stop("None of the requested `cells` matches a cell with stored weights in this fit.",
           call. = FALSE)
    }
    w <- w[keep, , drop = FALSE]
  }

  # Ordered discrete axes: tpre ascending; comparison cohorts ascending with the
  # never-treated (Inf) row last. Facets ordered by (g, t).
  tp_lev <- sort(unique(w$tpre))
  gp_lev <- sort(unique(w$gp))                       # numeric sort puts Inf last
  w$tpre_f   <- factor(w$tpre, levels = tp_lev)
  w$gp_f     <- factor(paste0("g' = ", w$gp), levels = paste0("g' = ", gp_lev))
  cell_keys  <- unique(w[order(w$group, w$time), c("group", "time"), drop = FALSE])
  lab_levels <- sprintf("ATT(g = %s, t = %s)", cell_keys$group, cell_keys$time)
  w$cell_lab <- factor(sprintf("ATT(g = %s, t = %s)", w$group, w$time), levels = lab_levels)

  # Symmetric fill limits so 0 is exactly the midpoint color and negative weights
  # are as visible as positive ones (diverging palette centered at 0).
  m <- max(abs(w$weight), na.rm = TRUE)
  if (!is.finite(m) || m <= 0) m <- 1

  ggplot2::ggplot(w, ggplot2::aes(x = tpre_f, y = gp_f, fill = weight)) +
    ggplot2::geom_tile(color = "grey85") +
    ggplot2::facet_wrap(~cell_lab) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-m, m), name = "E[w(X)]") +
    ggplot2::labs(
      x = expression("pre-treatment baseline period " * t[pre]),
      y = "comparison cohort g'",
      title = "EDiD weight decomposition",
      subtitle = "Mean weight on each (g', t_pre) identifying moment, per ATT(g,t) cell"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
