# edid-nocov.R
# No-covariate path for the EDiD estimator:
#   compute_omega_star_nocov_edid()
#   compute_efficient_weights_edid()
#   compute_generated_outcomes_nocov_edid()
#   compute_eif_nocov_edid()

# ---------------------------------------------------------------------------
# Helper: get column index from panel_obj
# ---------------------------------------------------------------------------
.col <- function(panel_obj, period_val) {
  panel_obj$period_to_col[[as.character(period_val)]]
}

# ---------------------------------------------------------------------------
# Omega* covariance matrix (H x H)
# ---------------------------------------------------------------------------

#' Compute the Omega* covariance matrix for the no-covariate EDiD path
#'
#' Builds the \eqn{H \times H} sample covariance matrix of the identifying
#' moments for cell \code{(target_g, target_t)}.
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param pt_assumption \code{"all"} or \code{"post"}
#'
#' @return numeric matrix H x H
#' @keywords internal
compute_omega_star_nocov_edid <- function(
  target_g, target_t, pairs, panel_obj, pt_assumption
) {
  H   <- nrow(pairs)
  n   <- panel_obj$n
  ow  <- panel_obj$outcome_wide

  mask_g   <- panel_obj$cohort_masks[[as.character(target_g)]]
  mask_inf <- panel_obj$never_treated_mask
  n_g      <- sum(mask_g)
  n_inf    <- sum(mask_inf)

  col_t  <- .col(panel_obj, target_t)
  col_1  <- .col(panel_obj, panel_obj$period_1)

  if (pt_assumption == "post") {
    # PT-Post: 1x1 matrix = var of standard DiD moment
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)

    delta_g   <- ow[mask_g,   col_t] - ow[mask_g,   col_base]
    delta_inf <- ow[mask_inf, col_t] - ow[mask_inf, col_base]

    omega <- matrix(
      cov_nn_edid(delta_g, delta_g) / n_g +
        cov_nn_edid(delta_inf, delta_inf) / n_inf,
      nrow = 1L, ncol = 1L
    )
    return(omega)
  }

  # ---------------------------------------------------------------------------
  # PT-All: H x H matrix, entry-by-entry
  # ---------------------------------------------------------------------------
  # Pre-compute treated-group change (same for all j, k)
  delta_g_t_1 <- ow[mask_g, col_t] - ow[mask_g, col_1]

  # Pre-compute never-treated changes for each unique tpre
  unique_tpre <- unique(pairs$tpre)
  delta_inf_cache <- vector("list", length(unique_tpre))
  names(delta_inf_cache) <- as.character(unique_tpre)
  for (tp in unique_tpre) {
    col_pre <- .col(panel_obj, tp)
    delta_inf_cache[[as.character(tp)]] <-
      ow[mask_inf, col_t] - ow[mask_inf, col_pre]
  }

  # Pre-compute comparison-cohort changes for each unique (gp, tpre)
  unique_gp_tpre <- unique(pairs[, c("gp", "tpre")])
  delta_gp_cache <- list()
  for (rr in seq_len(nrow(unique_gp_tpre))) {
    gp_val  <- unique_gp_tpre$gp[rr]
    tp_val  <- unique_gp_tpre$tpre[rr]
    key     <- paste0(gp_val, "_", tp_val)
    mask_gp <- panel_obj$cohort_masks[[as.character(gp_val)]]
    col_pre <- .col(panel_obj, tp_val)
    delta_gp_cache[[key]] <- ow[mask_gp, col_pre] - ow[mask_gp, col_1]
  }

  omega <- matrix(0, nrow = H, ncol = H)

  for (j in seq_len(H)) {
    gp_j   <- pairs$gp[j]
    tpre_j <- pairs$tpre[j]
    key_j  <- paste0(gp_j, "_", tpre_j)
    n_gp_j <- sum(panel_obj$cohort_masks[[as.character(gp_j)]])
    delta_inf_j <- delta_inf_cache[[as.character(tpre_j)]]
    delta_gp_j  <- delta_gp_cache[[key_j]]

    for (k in seq_len(H)) {
      if (k < j) {
        omega[j, k] <- omega[k, j]  # symmetric
        next
      }
      gp_k   <- pairs$gp[k]
      tpre_k <- pairs$tpre[k]
      key_k  <- paste0(gp_k, "_", tpre_k)
      n_gp_k <- sum(panel_obj$cohort_masks[[as.character(gp_k)]])
      delta_inf_k <- delta_inf_cache[[as.character(tpre_k)]]
      delta_gp_k  <- delta_gp_cache[[key_k]]

      # Term A: treated group variance (always present; same for all j, k)
      term_a <- cov_nn_edid(delta_g_t_1, delta_g_t_1) / n_g

      # Term B: never-treated cross-covariance
      term_b <- cov_nn_edid(delta_inf_j, delta_inf_k) / n_inf

      # Term C_j: non-zero only if gp_j == target_g
      term_cj <- 0
      if (is.finite(gp_j) && gp_j == target_g) {
        term_cj <- cov_nn_edid(delta_g_t_1, delta_gp_j) / n_g
      }

      # Term C_k: non-zero only if gp_k == target_g
      term_ck <- 0
      if (is.finite(gp_k) && gp_k == target_g) {
        term_ck <- cov_nn_edid(delta_g_t_1, delta_gp_k) / n_g
      }

      # Term D: non-zero only if gp_j == gp_k
      term_d <- 0
      if (gp_j == gp_k) {  # works for both finite and Inf
        term_d <- cov_nn_edid(delta_gp_j, delta_gp_k) / n_gp_j
      }

      omega[j, k] <- term_a + term_b - term_cj - term_ck + term_d
    }
  }

  omega
}

# ---------------------------------------------------------------------------
# Pole-target Ledoit-Wolf shrinkage of Omega* (no-covariate path; nocov_shrink)
# ---------------------------------------------------------------------------

# Validate the experimental no-covariate shrinkage target configuration.
# The default is the existing i.i.d. pole (rho = 0).  AR(1) is deliberately
# opt-in until simulation/application gates justify any automatic choice.
.validate_nocov_ar1_rho_edid <- function(rho, name = "rho") {
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop(sprintf("`%s` must be a finite numeric scalar in (-1, 1).", name), call. = FALSE)
  }
  as.numeric(rho)
}

#' @keywords internal
nocov_shrink_target_edid <- function() {
  target <- getOption("edid_nocov_shrink_target", "iid")
  if (!is.character(target) || length(target) != 1L || is.na(target)) {
    stop("`options(edid_nocov_shrink_target=...)` must be either 'iid' or 'ar1'.",
         call. = FALSE)
  }
  target <- tolower(target)
  if (!target %in% c("iid", "ar1")) {
    stop("`options(edid_nocov_shrink_target=...)` must be either 'iid' or 'ar1'.",
         call. = FALSE)
  }
  rho <- 0
  if (identical(target, "ar1")) {
    rho <- .validate_nocov_ar1_rho_edid(
      getOption("edid_nocov_ar1_rho", NA_real_),
      "options(edid_nocov_ar1_rho)"
    )
  }
  list(target = target, rho = rho)
}

#' Serial-correlation target structure matrix for a no-covariate cell's moment covariance
#'
#' Builds the \eqn{H \times H} matrix \eqn{S_\rho} such that under shocks
#' \eqn{\varepsilon_{i,t}} with variance \eqn{\sigma^2} and serial correlation
#' \eqn{Corr(\varepsilon_{i,a}, \varepsilon_{i,b}) = \rho^{|a-b|}} (plus
#' arbitrary unit effects and deterministic period effects, which difference
#' out), the population covariance of the cell's identifying moments is exactly
#' \eqn{\sigma^2 S_\rho} at the sample cohort sizes. The default \code{rho = 0}
#' is the original i.i.d.-shock pole, so
#' \deqn{Cov(\varepsilon_a-\varepsilon_b, \varepsilon_c-\varepsilon_d)/\sigma^2
#'   = 1\{a=c\} - 1\{a=d\} - 1\{b=c\} + 1\{b=d\}.}
#' The builder is the term-by-term mirror of
#' \code{compute_omega_star_nocov_edid()} (PT-All branch), replacing every
#' empirical covariance \code{cov_nn_edid(delta_a_b, delta_c_d)} by the
#' corresponding AR(1) kernel. Entries depend only on the pair set, group sizes
#' (shares), and \code{rho}. For \code{rho = 0} they reduce exactly to the
#' paper's closed-form i.i.d. pole covariance (the imputation/network algebra).
#' All edge cases (\code{tpre == period_1} degenerate self pairs, shared base
#' periods) are handled by the kernel mechanically, exactly as the empirical
#' builder handles them through zero/overlapping difference vectors.
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#'   (PT-All enumeration: \code{gp} finite)
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param rho AR(1) serial-correlation target in \code{(-1, 1)}. The default
#'   \code{0} is the original i.i.d. pole.
#'
#' @return numeric matrix H x H (unit-\eqn{\sigma^2} target covariance)
#' @keywords internal
compute_pole_structure_nocov_edid <- function(target_g, target_t, pairs, panel_obj, rho = 0) {
  rho <- .validate_nocov_ar1_rho_edid(rho)
  H  <- nrow(pairs)
  t1 <- panel_obj$period_1
  n_g   <- sum(panel_obj$cohort_masks[[as.character(target_g)]])
  n_inf <- sum(panel_obj$never_treated_mask)
  n_gp  <- vapply(pairs$gp, function(gp)
    sum(panel_obj$cohort_masks[[as.character(gp)]]), numeric(1L))

  times_all <- unique(as.numeric(c(target_t, t1, pairs$tpre)))
  if (rho < 0) {
    gaps <- abs(outer(times_all, times_all, "-"))
    if (any(abs(gaps - round(gaps)) > 1e-8)) {
      stop("negative `rho` requires integer-spaced time periods.", call. = FALSE)
    }
  }

  corr <- function(x, y) {
    lag <- abs(as.numeric(x) - as.numeric(y))
    if (rho < 0) rho^round(lag) else rho^lag
  }

  # Cov(eps_a - eps_b, eps_c - eps_d) / sigma2 under the chosen AR(1) target.
  # At rho = 0 this is exactly the old i.i.d. indicator kernel (0^0 == 1).
  kern <- function(a, b, cc, d) {
    corr(a, cc) - corr(a, d) - corr(b, cc) + corr(b, d)
  }

  S <- matrix(0, nrow = H, ncol = H)
  for (j in seq_len(H)) {
    gp_j <- pairs$gp[j]; tp_j <- pairs$tpre[j]
    for (k in j:H) {
      gp_k <- pairs$gp[k]; tp_k <- pairs$tpre[k]
      # Term A: treated-group difference (target_t, period_1) x itself
      val <- kern(target_t, t1, target_t, t1) / n_g
      # Term B: never-treated cross-covariance (t, tpre_j) x (t, tpre_k)
      val <- val + kern(target_t, tp_j, target_t, tp_k) / n_inf
      # Terms C: treated x own-cohort comparison (only when gp == target_g)
      if (is.finite(gp_j) && gp_j == target_g) {
        val <- val - kern(target_t, t1, tp_j, t1) / n_g
      }
      if (is.finite(gp_k) && gp_k == target_g) {
        val <- val - kern(target_t, t1, tp_k, t1) / n_g
      }
      # Term D: shared comparison cohort (tpre_j, 1) x (tpre_k, 1)
      if (gp_j == gp_k) {
        val <- val + kern(tp_j, t1, tp_k, t1) / n_gp[j]
      }
      S[j, k] <- val
      S[k, j] <- val
    }
  }
  S
}

#' Per-unit moment influence matrix for a no-covariate cell (PT-All)
#'
#' Returns the \eqn{n \times H} matrix \eqn{\psi} with
#' \deqn{\psi_{ij} = \frac{1\{G_i=g\}}{\pi_g}(\Delta^g_i - \bar\Delta^g)
#'   - \frac{1\{G_i=\infty\}}{\pi_\infty}(\Delta^{\infty,j}_i - \bar\Delta^{\infty,j})
#'   - \frac{1\{G_i=g'_j\}}{\pi_{g'_j}}(\Delta^{g'_j}_i - \bar\Delta^{g'_j}),}
#' the influence vector of moment \eqn{j}'s group means, mirroring
#' \code{compute_eif_nocov_edid()} pair by pair without weights. Because
#' \code{cov_nn_edid()} divides by the group size, the exact finite-sample
#' identity \code{compute_omega_star_nocov_edid() == crossprod(psi) / n^2}
#' holds (regression-tested); \eqn{\psi} therefore supplies the per-unit
#' entry-variance estimate the Ledoit-Wolf intensity needs.
#'
#' @inheritParams compute_pole_structure_nocov_edid
#' @return numeric matrix n x H
#' @keywords internal
compute_psi_moments_nocov_edid <- function(target_g, target_t, pairs, panel_obj) {
  n  <- panel_obj$n
  H  <- nrow(pairs)
  ow <- panel_obj$outcome_wide

  mask_g   <- panel_obj$cohort_masks[[as.character(target_g)]]
  mask_inf <- panel_obj$never_treated_mask
  pi_g     <- panel_obj$cohort_fractions[[as.character(target_g)]]
  pi_inf   <- sum(mask_inf) / n

  col_t <- .col(panel_obj, target_t)
  col_1 <- .col(panel_obj, panel_obj$period_1)

  delta_g <- ow[mask_g, col_t] - ow[mask_g, col_1]
  ctr_g   <- (delta_g - mean(delta_g)) / pi_g

  psi <- matrix(0, nrow = n, ncol = H)
  for (j in seq_len(H)) {
    gp_j    <- pairs$gp[j]
    col_pre <- .col(panel_obj, pairs$tpre[j])

    psi[mask_g, j] <- ctr_g

    delta_inf <- ow[mask_inf, col_t] - ow[mask_inf, col_pre]
    psi[mask_inf, j] <- psi[mask_inf, j] -
      (delta_inf - mean(delta_inf)) / pi_inf

    mask_gp  <- panel_obj$cohort_masks[[as.character(gp_j)]]
    pi_gp    <- panel_obj$cohort_fractions[[as.character(gp_j)]]
    delta_gp <- ow[mask_gp, col_pre] - ow[mask_gp, col_1]
    psi[mask_gp, j] <- psi[mask_gp, j] -
      (delta_gp - mean(delta_gp)) / pi_gp
  }
  psi
}

#' Ledoit-Wolf shrinkage of Omega* toward a serial-correlation target
#'
#' Implements the \code{nocov_shrink} option of \code{\link{edid}} for one
#' (g, t) cell on the no-covariate PT-All path. The default target is the
#' closed-form i.i.d. pole covariance at the sample shares,
#' \eqn{T = \hat\sigma^2 S} with \eqn{S} from
#' \code{compute_pole_structure_nocov_edid(rho = 0)} and the method-of-moments
#' scale \eqn{\hat\sigma^2 = \langle\hat\Omega, S\rangle_F / \langle S, S\rangle_F}
#' (the Frobenius least-squares projection, i.e. the \eqn{\sigma^2} minimizing
#' \eqn{\|\hat\Omega - \sigma^2 S\|_F}). The intensity is the standard
#' Ledoit-Wolf ratio (variance-of-entries over distance-to-target, clamped to
#' \eqn{[0, 1]}):
#' \deqn{\hat\lambda = \min\!\Big(1, \frac{\bar b^2}{d^2}\Big), \qquad
#'   \bar b^2 = \frac{1}{n^2}\sum_i \|B_i - \hat\Omega\|_F^2, \quad
#'   d^2 = \|\hat\Omega - T\|_F^2,}
#' where \eqn{B_i = \psi_i\psi_i'/n} is unit \eqn{i}'s contribution
#' (\eqn{\hat\Omega = n^{-1}\sum_i B_i} exactly). Off the pole
#' \eqn{d^2 \to \|\Omega - T\|^2_F > 0} while \eqn{\bar b^2 = O_p(1/n)} times
#' the entry scale, so \eqn{\hat\lambda \to 0} and the asymptotic weights (and
#' gains) are unchanged; at the pole the target is consistent for the truth, so
#' a large \eqn{\hat\lambda} costs nothing asymptotically and removes the
#' finite-sample weight-estimation noise. For diagnostics,
#' \code{options(edid_nocov_shrink_target = "ar1", edid_nocov_ar1_rho = r)}
#' replaces the i.i.d. target with \eqn{S_r}; this is experimental and not an
#' automatic selector.
#'
#' Returns the input unchanged (with \code{lambda = NA}) for degenerate inputs
#' (H < 2, non-finite or all-zero \code{omega}, non-positive projection scale,
#' or an exactly-zero structure matrix).
#'
#' @param omega numeric H x H matrix from \code{compute_omega_star_nocov_edid()}
#' @inheritParams compute_pole_structure_nocov_edid
#' @return list with \code{omega} (the shrunk matrix), \code{lambda} (the
#'   intensity in \eqn{[0,1]}, or \code{NA} when shrinkage did not apply),
#'   \code{sigma2} (the method-of-moments scale), \code{target}, and \code{rho}
#' @keywords internal
shrink_omega_nocov_edid <- function(omega, target_g, target_t, pairs, panel_obj) {
  cfg <- nocov_shrink_target_edid()
  no_op <- list(omega = omega, lambda = NA_real_, sigma2 = NA_real_,
                target = cfg$target, rho = cfg$rho)
  H <- nrow(omega)
  if (is.null(H) || H < 2L || any(!is.finite(omega)) || all(omega == 0)) return(no_op)

  S <- compute_pole_structure_nocov_edid(target_g, target_t, pairs, panel_obj, rho = cfg$rho)
  ss <- sum(S * S)
  if (!is.finite(ss) || ss <= 0) return(no_op)
  sigma2 <- sum(omega * S) / ss
  if (!is.finite(sigma2) || sigma2 <= 0) return(no_op)

  target <- sigma2 * S
  d2 <- sum((omega - target)^2)
  # Omega already (numerically) equals its pole projection: shrinking is a no-op.
  if (d2 <= .Machine$double.eps * max(sum(omega * omega), .Machine$double.xmin)) {
    return(list(omega = omega, lambda = 0, sigma2 = sigma2,
                target = cfg$target, rho = cfg$rho))
  }

  n   <- panel_obj$n
  psi <- compute_psi_moments_nocov_edid(target_g, target_t, pairs, panel_obj)
  # b_bar^2 = (1/n^2) sum_i ||psi_i psi_i'/n - Omega||_F^2
  #         = (1/n^2) [ sum_i ||psi_i||_2^4 / n^2 - n ||Omega||_F^2 ]
  # (cross term collapses because Omega = (1/n) sum_i psi_i psi_i' / n exactly).
  q4 <- sum(rowSums(psi * psi)^2)
  b2 <- (q4 / n^2 - n * sum(omega * omega)) / n^2
  if (!is.finite(b2)) return(no_op)
  lambda <- min(1, max(0, b2) / d2)

  list(omega = (1 - lambda) * omega + lambda * target,
       lambda = lambda, sigma2 = sigma2,
       target = cfg$target, rho = cfg$rho)
}

# ---------------------------------------------------------------------------
# Second-order weight-estimation variance correction (no-covariate path)
# Engaged by estimation_effect = TRUE on a no-covariate fit; see edid().
# ---------------------------------------------------------------------------

#' Closed-form second-order weight-estimation variance correction (no-covariate cell)
#'
#' On the no-covariate PT-All path the cell estimator is
#' \eqn{\hat\theta = w(\hat\Omega)'\hat m}, where \eqn{\hat m} is the H-vector of
#' generated-outcome means and \eqn{w(\Omega) = \Omega^{-1}\mathbf 1 /
#' (\mathbf 1'\Omega^{-1}\mathbf 1)} the efficient-weight map (optionally through
#' the pole-target shrinkage \code{\link{shrink_omega_nocov_edid}}). The plug-in
#' variance is \eqn{\widehat V_{plug} = \hat w'\hat\Omega\hat w} (the empirical
#' variance of the realized weighted IF; \eqn{\hat\Omega = \Psi'\Psi/n^2} exactly,
#' with \eqn{\Psi} the per-unit moment influence matrix of
#' \code{\link{compute_psi_moments_nocov_edid}}). It accounts for nothing about the
#' estimation of \eqn{\hat\Omega}, which both \emph{generates the weights} and is
#' \emph{evaluated by the same minimized quadratic} -- the audited source of the
#' small-n no-covariate SE shortfall.
#'
#' \strong{Decomposition (what is actually missing).} Under parallel trends every
#' moment has mean \eqn{ATT}, so \eqn{m = ATT\cdot\mathbf 1} and
#' \eqn{E[\hat\theta \mid \hat\Omega] = \hat w' m = ATT} whenever
#' \eqn{\hat\Omega \perp \hat m}; with (approximately) Gaussian shocks the group
#' means and group-demeaned covariances are independent, so by conditioning on
#' \eqn{\hat\Omega}:
#' \deqn{\mathrm{Var}(\hat\theta) \;=\; E\big[\hat w'\,\Omega\,\hat w\big]
#'   \qquad(\text{cross term } \mathrm{Cov}(T_1, T_2) = 0 \text{ and the
#'   weight-noise variance } \mathrm{Var}(T_2) \text{ both fold in}).}
#' The plug-in replaces \eqn{\Omega} by \eqn{\hat\Omega} \emph{evaluated at the
#' weights chosen to minimize it}, so its bias is the in-sample optimism
#' \deqn{E[\widehat V_{plug}] - \mathrm{Var}(\hat\theta)
#'   = E\big[\hat w'(\hat\Omega - \Omega)\hat w\big]
#'   = -\Delta_{DF} \;-\; 2\,Q \;+\; O(n^{-3/2}\,\mathrm{rel.}),}
#' with two closed-form pieces this function adds back
#' (\eqn{\mathrm{Var}_{add} = \Delta_{DF} + 2\hat Q}):
#' \describe{
#'   \item{Bessel piece \eqn{\Delta_{DF}}}{\eqn{\hat\Omega}'s group covariances
#'     divide by the group size \eqn{m_\gamma}, so
#'     \eqn{E[\hat\Omega] = \Omega - \sum_\gamma \Omega_\gamma/m_\gamma}. Because each
#'     unit's \eqn{\psi_i} loads on exactly one cohort, the unit-level closed form is
#'     \eqn{\Delta_{DF} = n^{-2}\sum_i a_i^2/(m_{\gamma(i)} - 1)}, \eqn{a_i = \hat w'\psi_i}.}
#'   \item{Optimization optimism \eqn{2\hat Q}}{second-order in
#'     \eqn{dE = \hat\Omega - \Omega}: \eqn{Q = -E[(J[dE])' dE\, \hat w]} with
#'     \eqn{J} the weight-map Jacobian below. With the per-unit directions
#'     \eqn{d_i = J[v_i]}, \eqn{v_i = \psi_i\psi_i'/n - \hat\Omega} (exactly
#'     mean-zero, and \eqn{\sum_i d_i = 0} exactly), the estimate collapses to
#'     \eqn{\hat Q = -n^{-3}\sum_i a_i (d_i'\psi_i) = -\widehat{\mathrm{cov}}_{lead}}.
#'     On the unshrunk path \eqn{d_i'\psi_i = -(a_i/n)\,\psi_i'B\psi_i} with
#'     \eqn{B = A - (\mathbf 1'A\mathbf 1) w w' \succeq 0} (\eqn{A = \hat\Omega^{-1}};
#'     \eqn{B\mathbf 1 = 0}, \eqn{B\hat\Omega\hat w = 0}), so \eqn{\hat Q \ge 0}: the
#'     minimized quadratic is always optimistic.}
#' }
#' The J-linear third-moment cross term \eqn{\widehat{\mathrm{cov}}_{lead}
#' = n^{-3}\sum_i a_i(d_i'\psi_i)} and the degenerate-U weight-noise variance
#' \eqn{\widehat{\mathrm{Var}}(T_2) = n^{-2}[\sum_i d_i'\hat\Omega d_i +
#' \mathrm{tr}(\hat G^2)]}, \eqn{\hat G = n^{-1}\sum_i d_i\psi_i'}, are returned as
#' \emph{diagnostics} (\code{cov_lead}, \code{var_second}) but are NOT added:
#' a naive \eqn{\widehat V_{plug} + 2\widehat{\mathrm{cov}}_{lead} +
#' \widehat{\mathrm{Var}}(T_2)} assembly double-counts \eqn{\mathrm{Var}(T_2)}
#' (already inside \eqn{E[\widehat V_{plug}]} through the realized-weight wobble)
#' and keeps a truncated cross term whose higher-order (\eqn{J_2}) parts cancel it
#' under the Gaussian independence above -- a Monte Carlo channel decomposition at
#' the n = 50 i.i.d. pole confirms both (true \eqn{\mathrm{Cov}(T_1,T_2) \approx 0};
#' the naive assembly moves calibration the wrong way, while
#' \eqn{\Delta_{DF} + 2\hat Q} restores mean SE / MC SD to ~0.95). Under
#' non-Gaussian shocks the exact-zero cross term is approximate; the omitted
#' remainder is \eqn{O(n^{-3/2})} relative. The \eqn{O(1/n)} two-step \emph{bias}
#' of \eqn{\hat\theta} is exactly zero under the same independence (the estimator
#' is conditionally unbiased given \eqn{\hat\Omega} under PT).
#'
#' \strong{The Jacobian.} Matrix calculus on the normalized-inverse map gives
#' \deqn{dw[dE] = -B\, dE\, w, \qquad B = A - (\mathbf 1'A\mathbf 1)\, w w', \quad
#'   A = \Omega^{-1},}
#' so the perturbed weights keep summing to one (\eqn{B\mathbf 1 = 0}).
#'
#' \strong{Differentiating through the shrinkage} (\code{nocov_shrink}; engaged when
#' the cell's \code{shrink_lambda} is in \eqn{(0, 1]}): with the same target
#' \eqn{S} and \code{shrink_rho} the weight path used,
#' \eqn{\Omega_{sh}(\Omega) = (1-\lambda)\Omega + \lambda\sigma^2 S},
#' \eqn{\sigma^2 = \langle\Omega, S\rangle_F/\langle S,S\rangle_F}, and the
#' Ledoit-Wolf \eqn{\lambda = \min(1, \max(0, b^2)/d^2)} (holding \eqn{S}, \eqn{n},
#' and the fourth-moment statistic \eqn{q_4} fixed), the chain rule gives
#' \deqn{d\Omega_{sh}[dE] = (1-\lambda)\,dE
#'   + \lambda\,\frac{\langle dE, S\rangle_F}{\langle S, S\rangle_F}\,S
#'   + d\lambda[dE]\,(\sigma^2 S - \Omega),}
#' \deqn{d\lambda[dE] = \frac{-\tfrac{2}{n}\langle\Omega, dE\rangle_F
#'   - 2\lambda\,\langle\Omega - \sigma^2 S, dE\rangle_F}{d^2}
#'   \quad (\text{interior } \lambda; \; d\lambda = 0 \text{ at the clamps } 0, 1),}
#' using \eqn{d(b^2)[dE] = -(2/n)\langle\Omega,dE\rangle_F} (from
#' \eqn{b^2 = (q_4/n^2 - n\|\Omega\|_F^2)/n^2} with \eqn{q_4} fixed) and
#' \eqn{d(d^2)[dE] = 2\langle\Omega - \sigma^2 S, dE\rangle_F} (the \eqn{\sigma^2}
#' channel of \eqn{d^2} vanishes by the projection orthogonality
#' \eqn{\langle\Omega - \sigma^2 S, S\rangle_F = 0}). The data-dependence of
#' \eqn{\lambda} through \eqn{q_4} is omitted (it multiplies
#' \eqn{\sigma^2 S - \Omega}, which vanishes at the pole, while off the pole
#' \eqn{\lambda \to 0}; a genuinely higher-order channel). The composed Jacobian is
#' then \eqn{J[dE] = -B_{sh}\, d\Omega_{sh}[dE]\, w} with \eqn{B_{sh}} built from
#' \eqn{\Omega_{sh}^{-1}}. The whole map is finite-difference verified in
#' \code{test-edid-nocov-estimation-effect.R}.
#'
#' Returns \code{applied = FALSE} (with a reason; the cell then keeps the plug-in
#' SE) when the weights are not the smooth inverse-map weights -- the
#' pseudoinverse / uniform fallback of \code{\link{compute_efficient_weights_edid}}
#' breaks the premise of the derivative -- or for degenerate inputs. Derived under
#' unit-level sampling; with clustered fits the leading term is cluster-robust
#' while this additive term is the unit-level estimate.
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows (PT-All)
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param omega_raw H x H \emph{unshrunk} moment covariance from
#'   \code{compute_omega_star_nocov_edid()} (the exact \eqn{\Psi'\Psi/n^2})
#' @param omega_used H x H matrix the weights actually inverted (the shrunk matrix
#'   when \code{nocov_shrink} applied; \code{== omega_raw} otherwise)
#' @param weights numeric vector length H: the realized efficient weights
#' @param shrink_lambda the cell's Ledoit-Wolf intensity (\code{NA} or 0 when the
#'   shrinkage did not bind; the chain rule through the shrinkage is applied for
#'   \code{shrink_lambda > 0})
#' @param shrink_rho AR(1) target rho used by the shrinkage path. The default
#'   \code{0} is the i.i.d. pole.
#' @param return_D logical: include the n x H matrix of per-unit directions
#'   \eqn{d_i} in the result (tests / diagnostics only)
#'
#' @return list with \code{applied} (logical), \code{warn} (logical: numeric
#'   failure vs structural skip), \code{reason} (string or NA),
#'   \code{var_add} (\eqn{\Delta_{DF} + 2\hat Q}: the additive variance applied to
#'   the cell), \code{delta_df} (\eqn{\Delta_{DF}}), \code{q_opt} (\eqn{\hat Q}),
#'   the diagnostics \code{cov_lead} (\eqn{= -\hat Q}) and \code{var_second}
#'   (J-linear \eqn{\widehat{\mathrm{Var}}(T_2)}), and optionally \code{D}
#' @keywords internal
compute_nocov_ee_correction_edid <- function(
  target_g, target_t, pairs, panel_obj, omega_raw, omega_used, weights,
  shrink_lambda = NA_real_, shrink_rho = 0, return_D = FALSE
) {
  # warn = FALSE marks a STRUCTURAL skip: the cell's weights are not the smooth inverse-map
  # estimator (pseudoinverse/uniform fallback on an exactly-singular Omega, e.g. duplicated
  # moments in degenerate pre-period pair sets), so there is no smooth weight-estimation
  # channel to correct -- skipping IS the correct treatment, recorded per cell but not warned.
  # warn = TRUE marks a numeric failure on a cell that should have supported the correction.
  skip <- function(reason, warn = FALSE) list(applied = FALSE, reason = reason, warn = warn,
                                              var_add = NA_real_, delta_df = NA_real_,
                                              q_opt = NA_real_, cov_lead = NA_real_,
                                              var_second = NA_real_)
  H <- nrow(pairs)
  n <- panel_obj$n
  if (is.null(H) || H < 2L) return(skip("just-identified cell (H < 2): weights are not estimated"))
  if (any(!is.finite(omega_raw)) || any(!is.finite(omega_used)) || any(!is.finite(weights)))
    return(skip("non-finite omega/weights", warn = TRUE))

  # The Jacobian differentiates the SMOOTH normalized-inverse map. Reject cells whose stored
  # weights came from a different construction (pseudoinverse fallback, uniform fallback;
  # the no-covariate path has no eigenvalue floor, so these are the only non-smooth cases) by
  # recomputing the map from omega_used and requiring an exact match -- the same w_chk gate the
  # covariate psi_Omega channel uses.
  A <- tryCatch(solve(omega_used), error = function(e) NULL)
  if (is.null(A)) return(skip("omega not invertible (fallback weights)"))
  u <- drop(A %*% rep(1, H))
  s <- sum(u)
  if (!is.finite(s) || abs(s) < EDID_DENOM_EPS) return(skip("degenerate 1'Omega^{-1}1"))
  if (max(abs(u / s - weights)) > 1e-8 * (1 + max(abs(weights))))
    return(skip("weights are not the smooth inverse-map weights (fallback path)"))

  B <- A - s * tcrossprod(weights)                       # dw = -B dOmega w;  B 1 = 0 exactly
  psi  <- compute_psi_moments_nocov_edid(target_g, target_t, pairs, panel_obj)   # n x H
  a    <- drop(psi %*% weights)                          # leading per-unit IF (= the cell EIF)
  Bpsi <- psi %*% B                                      # rows (B psi_i)' (B symmetric)
  BOw  <- drop(B %*% (omega_raw %*% weights))            # exactly 0 when omega_used == omega_raw

  # Per-unit directions d_i = J vec(v_i), v_i = psi_i psi_i'/n - omega_raw (exactly mean-zero).
  # Plain map: B v_i w = (a_i/n) B psi_i - B Omega w  =>  d_i = -(a_i/n) (B psi)_i + B Omega w.
  use_chain <- is.finite(shrink_lambda) && shrink_lambda > 0
  if (!use_chain) {
    D <- -((a / n) * Bpsi - matrix(BOw, n, H, byrow = TRUE))
  } else {
    # Chain rule through Omega_sh = (1 - lambda) Omega + lambda sigma2 S (see roxygen above).
    S  <- compute_pole_structure_nocov_edid(target_g, target_t, pairs, panel_obj, rho = shrink_rho)
    ss <- sum(S * S)
    if (!is.finite(ss) || ss <= 0) return(skip("degenerate pole structure"))
    sigma2  <- sum(omega_raw * S) / ss
    lam     <- min(1, max(0, shrink_lambda))
    BSw     <- drop(B %*% (S %*% weights))
    alpha_i <- rowSums((psi %*% S) * psi) / n - sum(omega_raw * S)    # <v_i, S>_F
    D <- -((1 - lam) * ((a / n) * Bpsi - matrix(BOw, n, H, byrow = TRUE)) +
           (lam / ss) * outer(alpha_i, BSw))
    if (lam < 1) {
      # interior lambda: the d-lambda channel, direction (sigma2 S - Omega_raw)
      d2 <- sum((omega_raw - sigma2 * S)^2)
      if (!is.finite(d2) || d2 <= 0) return(skip("degenerate shrinkage distance"))
      beta_i  <- rowSums((psi %*% omega_raw) * psi) / n - sum(omega_raw * omega_raw)  # <v_i, Omega>_F
      gamma_i <- beta_i - sigma2 * alpha_i                                            # <v_i, Omega - sigma2 S>_F
      kappa_i <- -(2 / (n * d2)) * beta_i - (2 * lam / d2) * gamma_i                  # dlambda[v_i]
      Btw     <- drop(B %*% ((sigma2 * S - omega_raw) %*% weights))                   # B (sigma2 S - Omega) w
      D <- D - outer(kappa_i, Btw)
    }
  }

  # Variance assembly (see roxygen): Bessel piece + optimization optimism; the J-linear
  # cross term and degenerate-U variance are kept as diagnostics only.
  s_i        <- rowSums(D * psi)                          # d_i' psi_i
  cov_lead   <- sum(a * s_i) / n^3                        # J-linear Cov(T1, T2) (diagnostic)
  q_opt      <- -cov_lead                                 # optimism Q-hat = -cov_lead (sum_i d_i = 0 exactly)
  G          <- crossprod(D, psi) / n                     # (1/n) sum_i d_i psi_i'
  var_second <- (sum((D %*% omega_raw) * D) + sum(G * t(G))) / n^2   # J-linear Var(T2) (diagnostic)

  # Bessel piece: each unit's psi_i loads on exactly ONE cohort, so Omega-hat decomposes by
  # cohort and the unbiased-covariance gap is the unit-level closed form below. Cohort size
  # from the panel's unit_cohorts (Inf = never treated); m = 1 cohorts contribute through the
  # guard max(m - 1, 1) (their variance is not estimable; the thin-cohort guard upstream
  # normally prevents this).
  coh   <- panel_obj$unit_cohorts
  sizes <- table(coh)                                      # named by as.character (Inf included)
  m_i   <- as.numeric(sizes[as.character(coh)])
  delta_df <- sum(a^2 / pmax(m_i - 1, 1)) / n^2

  var_add <- delta_df + 2 * q_opt
  if (!is.finite(var_add)) return(skip("non-finite correction", warn = TRUE))

  out <- list(applied = TRUE, reason = NA_character_, warn = FALSE, var_add = var_add,
              delta_df = delta_df, q_opt = q_opt,
              cov_lead = cov_lead, var_second = var_second,
              # per-unit pieces for the CROSS-CELL assembly (nocov_ee_sigma_full_edid):
              # s_vec_i = d_i' psi_i; together with the cell EIF a_i these give the
              # cross-cell optimism in closed form without storing D or psi.
              s_vec = s_i)
  if (isTRUE(return_D)) out$D <- D
  out
}

#' Full K x K no-covariate weight-estimation covariance increment (cross-cell)
#'
#' The aggregate plug-in covariance \eqn{\widehat{\mathrm{Cov}}(\hat\theta_c,
#' \hat\theta_{c'}) = n^{-2}\sum_i a_i^c a_i^{c'}} (the EIF cross-products) is
#' optimism-biased for the same reason as the cell variances: under the Gaussian
#' independence of group means and group-demeaned covariances,
#' \eqn{\mathrm{Cov}(\hat\theta_c, \hat\theta_{c'}) = E[\hat w_c'\,\Omega_{cc'}\,
#' \hat w_{c'}]} while the plug-in evaluates \eqn{\widehat\Omega_{cc'}} at the
#' optimized weights. The closed-form second-order correction for entry
#' \eqn{(c, c')} is
#' \deqn{\Sigma_{cc'} = \Delta_{DF,cc'} - n^{-3}\textstyle\sum_i
#'   \big(s_i^c a_i^{c'} + a_i^c s_i^{c'}\big), \qquad
#'   \Delta_{DF,cc'} = n^{-2}\textstyle\sum_i \frac{a_i^c a_i^{c'}}{m_{\gamma(i)}-1},}
#' with \eqn{a_i^c} the cell EIFs, \eqn{s_i^c = d_i^{c\prime}\psi_i^c} the per-unit
#' Jacobian projections returned by \code{compute_nocov_ee_correction_edid()}
#' (\eqn{\sum_i d_i^c = 0} exactly kills the centering terms), and
#' \eqn{m_{\gamma(i)}} unit i's cohort size (the Bessel factor of the
#' group-demeaned cross covariances). The diagonal reproduces the cell formula
#' \eqn{\Delta_{DF} + 2\hat Q} exactly, so cell SEs, \code{sqrt(diag(Sig))} at the
#' band site, and the aggregate increments stay mutually consistent. Rows/columns
#' of cells without an applied correction are zero (those cells keep the plug-in
#' convention everywhere).
#'
#' @param eif_matrix n x K matrix of cell EIFs (\code{fit_edid_cells()} output)
#' @param s_matrix n x K matrix of per-unit projections \eqn{s_i^c} (NA columns
#'   for cells without an applied correction)
#' @param unit_cohorts length-n vector of unit cohort labels (Inf = never treated)
#' @return K x K matrix, or NULL when no cell carries an applied correction
#' @keywords internal
nocov_ee_sigma_full_edid <- function(eif_matrix, s_matrix, unit_cohorts) {
  if (is.null(eif_matrix) || is.null(s_matrix)) return(NULL)
  K <- ncol(s_matrix)
  ok <- which(vapply(seq_len(K), function(k) all(is.finite(s_matrix[, k])), logical(1L)))
  if (length(ok) == 0L) return(NULL)
  n <- nrow(s_matrix)
  sizes <- table(unit_cohorts)
  m_i   <- as.numeric(sizes[as.character(unit_cohorts)])
  A  <- eif_matrix[, ok, drop = FALSE]
  Sv <- s_matrix[, ok, drop = FALSE]
  out <- matrix(0, K, K)
  df_blk  <- crossprod(A / pmax(m_i - 1, 1), A) / n^2          # Delta_DF block (A' diag(1/(m-1)) A)
  sa      <- crossprod(Sv, A) / n^3                            # (1/n^3) sum_i s^c a^{c'}
  out[ok, ok] <- df_blk - (sa + t(sa))
  out
}

#' Assemble the per-cell no-covariate weight-estimation variance increments
#'
#' Returns the K x K \emph{diagonal} covariance increment whose entry k is cell
#' k's applied \code{nocov_ee$var_add} (0 where the correction did not apply), or
#' \code{NULL} when no cell carries an applied correction. Mirrors the
#' \code{sigma_quad_edid()} convention so the same consumers (the analytic cell
#' band in \code{edid()} and the aggregations in \code{aggte_edid()}) can add it
#' to the first-order covariance. Cross-cell second-order covariances are not
#' estimated (the omitted off-diagonal entries are higher-order for the
#' aggregate SEs in the same sense the own-cell term is for the cell SEs).
#'
#' @param cells list of \code{edid_cell_result} objects (cell order = ATT(g,t) order)
#' @return K x K diagonal matrix, or NULL
#' @keywords internal
nocov_ee_sigma_edid <- function(cells) {
  K <- length(cells)
  if (K == 0L) return(NULL)
  v <- vapply(cells, function(cc) {
    ee <- cc$nocov_ee
    if (is.null(ee) || !isTRUE(ee$applied) || !is.finite(ee$var_add)) 0 else ee$var_add
  }, numeric(1L))
  if (all(v == 0)) return(NULL)
  diag(v, nrow = K)
}

# ---------------------------------------------------------------------------
# Efficient weights
# ---------------------------------------------------------------------------

#' Compute efficient inverse-covariance weights
#'
#' Implements \eqn{w = (\Omega^{*-1} \mathbf{1}) / (\mathbf{1}' \Omega^{*-1} \mathbf{1})}
#' with fallback to uniform weights when the matrix is degenerate.
#'
#' @param omega_star numeric matrix H x H
#'
#' @return numeric vector length H, summing to 1
#' @keywords internal
compute_efficient_weights_edid <- function(omega_star) {
  H       <- nrow(omega_star)
  ones_H  <- rep(1, H)
  unif    <- ones_H / H

  # Degenerate: all zeros, or any non-finite entry (e.g. gmm cov() with NA moments
  # or a near-singular cell) -> uniform fallback instead of erroring in svd/eigen.
  if (all(omega_star == 0) || any(!is.finite(omega_star))) return(unif)

  # Check condition number
  kappa <- check_condition_edid(omega_star)
  if (!is.finite(kappa) || kappa > EDID_COND_THRESH) {
    inv_omega <- compute_pseudoinverse_edid(omega_star)
  } else {
    inv_omega <- tryCatch(
      solve(omega_star),
      error = function(e) compute_pseudoinverse_edid(omega_star)
    )
  }

  num   <- drop(inv_omega %*% ones_H)
  denom <- sum(num)
  if (!is.finite(denom) || abs(denom) < EDID_DENOM_EPS) return(unif)
  num / denom
}

# ---------------------------------------------------------------------------
# Generated outcomes (scalar moments)
# ---------------------------------------------------------------------------

#' Compute generated-outcome scalars for each valid pair
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param pt_assumption \code{"all"} or \code{"post"}
#'
#' @return numeric vector length H
#' @keywords internal
compute_generated_outcomes_nocov_edid <- function(
  target_g, target_t, pairs, panel_obj, pt_assumption
) {
  H   <- nrow(pairs)
  ow  <- panel_obj$outcome_wide

  mask_g   <- panel_obj$cohort_masks[[as.character(target_g)]]
  mask_inf <- panel_obj$never_treated_mask

  col_t  <- .col(panel_obj, target_t)

  if (pt_assumption == "post") {
    # PT-Post: one pair (Inf, base)
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)
    y_hat <- mean(ow[mask_g, col_t] - ow[mask_g, col_base]) -
             mean(ow[mask_inf, col_t] - ow[mask_inf, col_base])
    return(y_hat)  # scalar; will be treated as length-1 vector
  }

  # PT-All
  col_1 <- .col(panel_obj, panel_obj$period_1)
  term_g <- mean(ow[mask_g, col_t] - ow[mask_g, col_1])

  y_hat <- numeric(H)
  for (j in seq_len(H)) {
    gp_j    <- pairs$gp[j]
    tpre_j  <- pairs$tpre[j]
    col_pre <- .col(panel_obj, tpre_j)

    term_inf <- mean(ow[mask_inf, col_t] - ow[mask_inf, col_pre])

    mask_gp  <- panel_obj$cohort_masks[[as.character(gp_j)]]
    term_gp  <- mean(ow[mask_gp, col_pre] - ow[mask_gp, col_1])

    y_hat[j] <- term_g - term_inf - term_gp
  }
  y_hat
}

# ---------------------------------------------------------------------------
# Efficient Influence Function (per unit, length n)
# ---------------------------------------------------------------------------

#' Compute the no-covariate efficient influence function for a (g, t) cell
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#' @param weights numeric vector length H (efficient weights)
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param att_gt scalar ATT estimate for this cell
#' @param pt_assumption \code{"all"} or \code{"post"}
#'
#' @return numeric vector length n (zero-mean by construction)
#' @keywords internal
compute_eif_nocov_edid <- function(
  target_g, target_t, pairs, weights, panel_obj, att_gt, pt_assumption
) {
  n    <- panel_obj$n
  H    <- nrow(pairs)
  ow   <- panel_obj$outcome_wide

  mask_g   <- panel_obj$cohort_masks[[as.character(target_g)]]
  mask_inf <- panel_obj$never_treated_mask
  pi_g     <- panel_obj$cohort_fractions[[as.character(target_g)]]
  pi_inf   <- sum(mask_inf) / n

  col_t <- .col(panel_obj, target_t)

  eif <- numeric(n)

  if (pt_assumption == "post") {
    # PT-Post: one pair; base = g - 1 - anticipation
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)
    w_1      <- weights[1L]

    delta_g   <- ow[mask_g,   col_t] - ow[mask_g,   col_base]
    mean_g    <- mean(delta_g)
    eif[mask_g] <- eif[mask_g] + w_1 * (delta_g - mean_g) / pi_g

    delta_inf  <- ow[mask_inf, col_t] - ow[mask_inf, col_base]
    mean_inf   <- mean(delta_inf)
    eif[mask_inf] <- eif[mask_inf] - w_1 * (delta_inf - mean_inf) / pi_inf

    # gp_j == Inf: no comparison cohort term
  } else {
    # PT-All
    col_1    <- .col(panel_obj, panel_obj$period_1)
    col_base <- col_1   # base for treated group

    delta_g_t_base <- ow[mask_g, col_t] - ow[mask_g, col_base]
    mean_g_t_base  <- mean(delta_g_t_base)

    for (j in seq_len(H)) {
      w_j     <- weights[j]
      gp_j    <- pairs$gp[j]
      tpre_j  <- pairs$tpre[j]
      col_pre <- .col(panel_obj, tpre_j)

      # Treated group contribution (always present)
      # Y_hat_j enters via centering: phi_{ij} = I(G=g)/pi_g * (delta - Y_hat_j)
      # But the EIF is: score - att_gt, and score = sum w_j phi_{ij}
      # phi_{ij} for treated group = I(G=g)/pi_g * (delta_{g,t,base} - Y_hat_j)
      # We compute: sum_j w_j * I(G=g)/pi_g * (delta - Y_hat_j)
      # = I(G=g)/pi_g * [ (delta - mean_g) + (mean_g - sum_j w_j Y_hat_j) ]
      # But it's simpler to accumulate directly:
      eif[mask_g] <- eif[mask_g] +
        w_j * (delta_g_t_base - mean_g_t_base) / pi_g

      # Never-treated contribution (subtract)
      delta_inf_t_pre <- ow[mask_inf, col_t] - ow[mask_inf, col_pre]
      mean_inf_t_pre  <- mean(delta_inf_t_pre)
      eif[mask_inf] <- eif[mask_inf] -
        w_j * (delta_inf_t_pre - mean_inf_t_pre) / pi_inf

      # Comparison cohort contribution (subtract)
      mask_gp          <- panel_obj$cohort_masks[[as.character(gp_j)]]
      pi_gp            <- panel_obj$cohort_fractions[[as.character(gp_j)]]
      delta_gp_pre_base <- ow[mask_gp, col_pre] - ow[mask_gp, col_base]
      mean_gp_pre_base  <- mean(delta_gp_pre_base)
      eif[mask_gp] <- eif[mask_gp] -
        w_j * (delta_gp_pre_base - mean_gp_pre_base) / pi_gp
    }
  }

  # The score accumulated above already has mean 0 by construction:
  # each group's contribution is demeaned (delta_i - mean(delta)).
  # The EIF is the score itself; do NOT subtract att_gt again.
  eif
}
