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

  # Observation weights (NULL => unweighted; wvar_term_edid then reduces to cov_nn/n_g exactly).
  uw     <- panel_obj$unit_weights
  w_g    <- if (is.null(uw)) NULL else uw[mask_g]
  w_inf  <- if (is.null(uw)) NULL else uw[mask_inf]

  col_t  <- .col(panel_obj, target_t)
  col_1  <- .col(panel_obj, panel_obj$period_1)

  if (pt_assumption == "post") {
    # PT-Post: 1x1 matrix = var of standard (weighted) DiD moment
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)

    delta_g   <- ow[mask_g,   col_t] - ow[mask_g,   col_base]
    delta_inf <- ow[mask_inf, col_t] - ow[mask_inf, col_base]

    omega <- matrix(
      wvar_term_edid(delta_g, delta_g, w_g) +
        wvar_term_edid(delta_inf, delta_inf, w_inf),
      nrow = 1L, ncol = 1L
    )
    return(omega)
  }

  # ---------------------------------------------------------------------------
  # PT-All: H x H matrix, entry-by-entry
  # ---------------------------------------------------------------------------
  # Each term is a within-group (cross-)covariance of difference vectors divided by the
  # group size; wvar_term_edid(.) carries that division (n_g unweighted; W_g^2 / sum(w^2)
  # weighted) so the whole builder is weight-aware via the per-group weight vectors below.
  # Pre-compute treated-group change (same for all j, k)
  delta_g_t_1 <- ow[mask_g, col_t] - ow[mask_g, col_1]

  # Per-group weight vectors aligned to each cached difference vector (NULL when unweighted).
  w_gp_cache <- list()
  w_for <- function(gp_val) {
    if (is.null(uw)) return(NULL)
    key <- as.character(gp_val)
    if (is.null(w_gp_cache[[key]])) w_gp_cache[[key]] <<- uw[panel_obj$cohort_masks[[key]]]
    w_gp_cache[[key]]
  }

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
      delta_inf_k <- delta_inf_cache[[as.character(tpre_k)]]
      delta_gp_k  <- delta_gp_cache[[key_k]]

      # Term A: treated group variance (always present; same for all j, k)
      term_a <- wvar_term_edid(delta_g_t_1, delta_g_t_1, w_g)

      # Term B: never-treated cross-covariance
      term_b <- wvar_term_edid(delta_inf_j, delta_inf_k, w_inf)

      # Term C_j: non-zero only if gp_j == target_g
      term_cj <- 0
      if (is.finite(gp_j) && gp_j == target_g) {
        term_cj <- wvar_term_edid(delta_g_t_1, delta_gp_j, w_g)
      }

      # Term C_k: non-zero only if gp_k == target_g
      term_ck <- 0
      if (is.finite(gp_k) && gp_k == target_g) {
        term_ck <- wvar_term_edid(delta_g_t_1, delta_gp_k, w_g)
      }

      # Term D: non-zero only if gp_j == gp_k
      term_d <- 0
      if (gp_j == gp_k) {  # works for both finite and Inf
        term_d <- wvar_term_edid(delta_gp_j, delta_gp_k, w_for(gp_j))
      }

      omega[j, k] <- term_a + term_b - term_cj - term_ck + term_d
    }
  }

  omega
}

# ---------------------------------------------------------------------------
# Pole-target Ledoit-Wolf shrinkage of Omega* (no-covariate path; nocov_shrink)
# ---------------------------------------------------------------------------

#' i.i.d.-pole structure matrix for a no-covariate cell's moment covariance
#'
#' Builds the \eqn{H \times H} matrix \eqn{S} such that under i.i.d. shocks
#' \eqn{\varepsilon_{i,t}} with variance \eqn{\sigma^2} (plus arbitrary unit
#' effects and deterministic period effects, which difference out), the
#' population covariance of the cell's identifying moments is exactly
#' \eqn{\sigma^2 S} at the sample cohort sizes. It is the term-by-term mirror
#' of \code{compute_omega_star_nocov_edid()} (PT-All branch) with every
#' empirical covariance \code{cov_nn_edid(delta_a_b, delta_c_d)} replaced by
#' the i.i.d.-shock kernel
#' \deqn{Cov(\varepsilon_a-\varepsilon_b, \varepsilon_c-\varepsilon_d)/\sigma^2
#'   = 1\{a=c\} - 1\{a=d\} - 1\{b=c\} + 1\{b=d\},}
#' so entries depend only on the pair set and the group sizes (shares), per
#' the paper's closed-form pole covariance (the imputation/network algebra).
#' All edge cases (\code{tpre == period_1} degenerate self pairs, shared base
#' periods) are handled by the kernel mechanically, exactly as the empirical
#' builder handles them through zero/overlapping difference vectors.
#'
#' @param target_g scalar cohort value
#' @param target_t scalar time period
#' @param pairs data.frame with columns \code{gp} and \code{tpre}; H rows
#'   (PT-All enumeration: \code{gp} finite)
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#'
#' @return numeric matrix H x H (unit-\eqn{\sigma^2} pole covariance)
#' @keywords internal
compute_pole_structure_nocov_edid <- function(target_g, target_t, pairs, panel_obj) {
  H  <- nrow(pairs)
  t1 <- panel_obj$period_1
  # Group "effective inverse size": 1/m unweighted; the weighted Hajek design factor
  # sum(w^2)/(sum w)^2 weighted (the same wvar_term normalization the empirical Omega uses).
  # The pole structure then has the SAME 1/size prefactors as the empirical builder, so the
  # i.i.d.-pole target is the structured covariance at the actual (weighted) shares.
  uw <- panel_obj$unit_weights
  .inv_size <- function(mask) {
    if (is.null(uw)) return(1 / sum(mask))
    ww <- uw[mask]; sw <- sum(ww); sum(ww * ww) / (sw * sw)
  }
  inv_g   <- .inv_size(panel_obj$cohort_masks[[as.character(target_g)]])
  inv_inf <- .inv_size(panel_obj$never_treated_mask)
  inv_gp  <- vapply(pairs$gp, function(gp)
    .inv_size(panel_obj$cohort_masks[[as.character(gp)]]), numeric(1L))

  # Cov(eps_a - eps_b, eps_c - eps_d) / sigma2 under i.i.d. shocks
  kern <- function(a, b, cc, d) {
    (a == cc) - (a == d) - (b == cc) + (b == d)
  }

  S <- matrix(0, nrow = H, ncol = H)
  for (j in seq_len(H)) {
    gp_j <- pairs$gp[j]; tp_j <- pairs$tpre[j]
    for (k in j:H) {
      gp_k <- pairs$gp[k]; tp_k <- pairs$tpre[k]
      # Term A: treated-group difference (target_t, period_1) x itself
      val <- kern(target_t, t1, target_t, t1) * inv_g
      # Term B: never-treated cross-covariance (t, tpre_j) x (t, tpre_k)
      val <- val + kern(target_t, tp_j, target_t, tp_k) * inv_inf
      # Terms C: treated x own-cohort comparison (only when gp == target_g)
      if (is.finite(gp_j) && gp_j == target_g) {
        val <- val - kern(target_t, t1, tp_j, t1) * inv_g
      }
      if (is.finite(gp_k) && gp_k == target_g) {
        val <- val - kern(target_t, t1, tp_k, t1) * inv_g
      }
      # Term D: shared comparison cohort (tpre_j, 1) x (tpre_k, 1)
      if (gp_j == gp_k) {
        val <- val + kern(tp_j, t1, tp_k, t1) * inv_gp[j]
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

  # Observation weights. With weights the per-unit moment influence carries an explicit
  # w_i factor and a weighted (Hajek) centering, so crossprod(psi)/n^2 == the weighted
  # Omega* (= wvar_term_edid sums). With NULL weights, w_g/w_inf/w_gp are all-ones via the
  # `* 1` and the centering is the plain mean, so psi is byte-identical to the legacy form.
  uw     <- panel_obj$unit_weights
  w_g    <- if (is.null(uw)) rep(1, sum(mask_g))   else uw[mask_g]
  w_inf  <- if (is.null(uw)) rep(1, sum(mask_inf)) else uw[mask_inf]
  pi_inf <- if (is.null(uw)) sum(mask_inf) / n else sum(w_inf) / n

  col_t <- .col(panel_obj, target_t)
  col_1 <- .col(panel_obj, panel_obj$period_1)

  delta_g <- ow[mask_g, col_t] - ow[mask_g, col_1]
  # ctr_g = w_i * (delta_i - wmean) / pi_g  (NULL weights: w_i = 1, wmean = mean => legacy)
  ctr_g   <- w_g * (delta_g - wmean_edid(delta_g, if (is.null(uw)) NULL else w_g)) / pi_g

  psi <- matrix(0, nrow = n, ncol = H)
  for (j in seq_len(H)) {
    gp_j    <- pairs$gp[j]
    col_pre <- .col(panel_obj, pairs$tpre[j])

    psi[mask_g, j] <- ctr_g

    delta_inf <- ow[mask_inf, col_t] - ow[mask_inf, col_pre]
    psi[mask_inf, j] <- psi[mask_inf, j] -
      w_inf * (delta_inf - wmean_edid(delta_inf, if (is.null(uw)) NULL else w_inf)) / pi_inf

    mask_gp  <- panel_obj$cohort_masks[[as.character(gp_j)]]
    pi_gp    <- panel_obj$cohort_fractions[[as.character(gp_j)]]
    w_gp     <- if (is.null(uw)) rep(1, sum(mask_gp)) else uw[mask_gp]
    delta_gp <- ow[mask_gp, col_pre] - ow[mask_gp, col_1]
    psi[mask_gp, j] <- psi[mask_gp, j] -
      w_gp * (delta_gp - wmean_edid(delta_gp, if (is.null(uw)) NULL else w_gp)) / pi_gp
  }
  psi
}

#' Ledoit-Wolf shrinkage of Omega* toward its i.i.d.-pole structure
#'
#' Implements the \code{nocov_shrink} option of \code{\link{edid}} for one
#' (g, t) cell on the no-covariate PT-All path. The target is the closed-form
#' pole covariance at the sample shares,
#' \eqn{T = \hat\sigma^2 S} with \eqn{S} from
#' \code{compute_pole_structure_nocov_edid()} and the method-of-moments scale
#' \eqn{\hat\sigma^2 = \langle\hat\Omega, S\rangle_F / \langle S, S\rangle_F}
#' (the Frobenius least-squares projection, i.e. the \eqn{\sigma^2} minimizing
#' \eqn{\|\hat\Omega - \sigma^2 S\|_F}). The intensity is the standard
#' Ledoit-Wolf ratio (variance-of-entries over distance-to-target, clamped to
#' \eqn{[0, 1]}):
#' \deqn{\hat\lambda = \min\!\Big(1, \frac{\bar b^2}{d^2}\Big), \qquad
#'   \bar b^2 = \frac{\hat\pi}{n_{\mathrm{eff}}}, \quad
#'   \hat\pi = \frac{1}{n}\sum_i \|B_i - \hat\Omega\|_F^2, \quad
#'   d^2 = \|\hat\Omega - T\|_F^2,}
#' where \eqn{B_i = \psi_i\psi_i'/n} is unit \eqn{i}'s contribution
#' (\eqn{\hat\Omega = n^{-1}\sum_i B_i} exactly) and \eqn{n_{\mathrm{eff}}} is the
#' Kish effective sample size of the units active in this cell's weighted
#' \eqn{\hat\Omega} (\code{\link{n_eff_edid}}). Unweighted
#' \eqn{n_{\mathrm{eff}} = n} exactly, so
#' \eqn{\bar b^2 = (q_4/n^2 - n\|\hat\Omega\|_F^2)/n^2} bit-for-bit (the legacy
#' form); under dispersed observation weights the heavily-weighted units dominate
#' \eqn{\hat\Omega}, so the raw \eqn{n} would under-shrink by
#' \eqn{n/n_{\mathrm{eff}}} -- only the OUTER averaging factor \eqn{1/n_{\mathrm{eff}}}
#' (the variance-of-the-average) changes; the internal \eqn{1/n} of \eqn{B_i}
#' carries the fixed \eqn{1/\pi_g} scale of \eqn{\hat\Omega} and stays \eqn{n}.
#' Off the pole \eqn{d^2 \to \|\Omega - T\|^2_F > 0} while
#' \eqn{\bar b^2 = O_p(1/n_{\mathrm{eff}})} times the entry scale, so
#' \eqn{\hat\lambda \to 0} and the asymptotic weights (and gains) are unchanged;
#' at the pole the target is consistent for the truth, so a large
#' \eqn{\hat\lambda} costs nothing asymptotically and removes the finite-sample
#' weight-estimation noise.
#'
#' Returns the input unchanged (with \code{lambda = NA}) for degenerate inputs
#' (H < 2, non-finite or all-zero \code{omega}, non-positive projection scale,
#' or an exactly-zero structure matrix).
#'
#' @param omega numeric H x H matrix from \code{compute_omega_star_nocov_edid()}
#' @param cl_metric_on logical; \code{TRUE} when \code{omega} is the CLUSTER moment
#'   covariance \eqn{\Sigma_{cl}} (its i.i.d. sampling units are the \eqn{G}
#'   clusters, not the \eqn{n} units), in which case the LW averaging factor uses
#'   the cluster ESS \code{cl_n_eff} rather than the unit Kish ESS. Default
#'   \code{FALSE} (unit metric), byte-identical to the legacy call.
#' @param cl_n_eff numeric; effective number of clusters (Kish ESS of the active
#'   clusters' total weights) entering \eqn{\Sigma_{cl}}, used only when
#'   \code{cl_metric_on}. Mirrors the cluster-ESS switch of the no-cov ridge.
#' @inheritParams compute_pole_structure_nocov_edid
#' @return list with \code{omega} (the shrunk matrix), \code{lambda} (the
#'   intensity in \eqn{[0,1]}, or \code{NA} when shrinkage did not apply), and
#'   \code{sigma2} (the method-of-moments scale)
#' @keywords internal
shrink_omega_nocov_edid <- function(omega, target_g, target_t, pairs, panel_obj,
                                    cl_metric_on = FALSE, cl_n_eff = NA_real_) {
  no_op <- list(omega = omega, lambda = NA_real_, sigma2 = NA_real_)
  H <- nrow(omega)
  if (is.null(H) || H < 2L || any(!is.finite(omega)) || all(omega == 0)) return(no_op)

  S <- compute_pole_structure_nocov_edid(target_g, target_t, pairs, panel_obj)
  ss <- sum(S * S)
  if (!is.finite(ss) || ss <= 0) return(no_op)
  sigma2 <- sum(omega * S) / ss
  if (!is.finite(sigma2) || sigma2 <= 0) return(no_op)

  target <- sigma2 * S
  d2 <- sum((omega - target)^2)
  # Omega already (numerically) equals its pole projection: shrinking is a no-op.
  if (d2 <= .Machine$double.eps * max(sum(omega * omega), .Machine$double.xmin)) {
    return(list(omega = omega, lambda = 0, sigma2 = sigma2))
  }

  n   <- panel_obj$n
  psi <- compute_psi_moments_nocov_edid(target_g, target_t, pairs, panel_obj)
  # Ledoit-Wolf intensity b_bar^2 / d^2. The standard LW estimator is b_bar^2 = pi_hat / n_eff, where
  # pi_hat = (1/n) sum_{all i} ||psi_i psi_i'/n - Omega||_F^2 is the empirical mean squared entry-
  # deviation and the OUTER 1/n_eff is the variance-of-the-average factor (it shrinks like one over
  # the number of independent contributions). The legacy form b2 = (q4/n^2 - n ||Omega||^2)/n^2 IS
  # pi_hat/n with n_eff = the FULL n (the averaging is over all n units; inactive units contribute
  # ||Omega||_F^2 each, the source of the -n||Omega||^2 term). The internal psi_i psi_i'/n carries the
  # FIXED 1/pi_g = n/W_g scale of Omega-hat (NOT a count to replace), so only the OUTER averaging
  # factor reflects the effective sample size. Under dispersed observation weights the heavily-
  # weighted units dominate Omega-hat, so the raw n under-shrinks by n/n_eff; n_eff is the Kish ESS of
  # this cell's active units (the SAME denominator as the no-cov ridge). BYTE-IDENTITY: keep the
  # legacy expression verbatim, apply the correction factor n / n_eff, which is EXACTLY 1.0 unweighted
  # (n_eff_edid returns the same full n there, identical doubles), so b2 is bit-for-bit the legacy
  # value.
  q4        <- sum(rowSums(psi * psi)^2)
  b2_legacy <- (q4 / n^2 - n * sum(omega * omega)) / n^2          # = pi_hat_full / n (verbatim legacy)
  # The OUTER averaging factor reflects the number of INDEPENDENT contributions to Omega-hat. For the unit
  # metric that is the unit Kish ESS; for the CLUSTER metric Sig_cl the i.i.d. sampling units are the G
  # clusters, so the averaging is over the cluster ESS cl_n_eff (mirrors the cluster-ESS ridge switch at
  # edid-fit.R cl_metric_on). The internal psi_i psi_i'/n still carries the FIXED 1/pi_g scale of Omega-hat
  # (unchanged). For the unit metric (cl_metric_on = FALSE) this is byte-identical: n/n_eff == 1.0 unweighted.
  n_eff     <- if (isTRUE(cl_metric_on) && is.finite(cl_n_eff) && cl_n_eff > 0)
                 cl_n_eff
               else
                 n_eff_edid(panel_obj$unit_weights,
                            active_mask_nocov_edid(target_g, pairs, panel_obj), n)
  b2        <- b2_legacy * (n / n_eff)                            # n/n_eff == 1.0 exactly unit-unweighted
  if (!is.finite(b2)) return(no_op)
  lambda <- min(1, max(0, b2) / d2)

  list(omega = (1 - lambda) * omega + lambda * target,
       lambda = lambda, sigma2 = sigma2)
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
#' the cell's \code{shrink_lambda} is in \eqn{(0, 1]}): with
#' \eqn{\Omega_{sh}(\Omega) = (1-\lambda)\Omega + \lambda\sigma^2 S},
#' \eqn{\sigma^2 = \langle\Omega, S\rangle_F/\langle S,S\rangle_F}, and the
#' Ledoit-Wolf \eqn{\lambda = \min(1, \max(0, b^2)/d^2)} (holding \eqn{S}, \eqn{n},
#' and the fourth-moment statistic \eqn{q_4} fixed), the chain rule gives
#' \deqn{d\Omega_{sh}[dE] = (1-\lambda)\,dE
#'   + \lambda\,\frac{\langle dE, S\rangle_F}{\langle S, S\rangle_F}\,S
#'   + d\lambda[dE]\,(\sigma^2 S - \Omega),}
#' \deqn{d\lambda[dE] = \frac{-\tfrac{2}{n_{\mathrm{eff}}}\langle\Omega, dE\rangle_F
#'   - 2\lambda\,\langle\Omega - \sigma^2 S, dE\rangle_F}{d^2}
#'   \quad (\text{interior } \lambda; \; d\lambda = 0 \text{ at the clamps } 0, 1),}
#' using \eqn{d(b^2)[dE] = -(2/n_{\mathrm{eff}})\langle\Omega,dE\rangle_F} (from
#' \eqn{b^2 = q_4/(n^3 n_{\mathrm{eff}}) - \|\Omega\|_F^2/n_{\mathrm{eff}}} with
#' \eqn{q_4} fixed; \eqn{n_{\mathrm{eff}}} the cell's Kish ESS, \eqn{= n}
#' unweighted) and
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
#' @param return_D logical: include the n x H matrix of per-unit directions
#'   \eqn{d_i} in the result (tests / diagnostics only)
#' @param cluster_indices length-n cluster id vector, or \code{NULL} (i.i.d.). When
#'   supplied, the weights invert the CLUSTER moment covariance
#'   \eqn{\widehat\Sigma_{cl} = \mathrm{crossprod}(\mathrm{rowsum}(\psi, cl))/n^2}
#'   (== \code{omega_raw}/\code{omega_used} here), so the weight-estimation channel
#'   is driven by the per-CLUSTER moment IF \eqn{\Psi_g = \sum_{i\in g}\psi_i}. The
#'   function then returns: (1) a per-UNIT first-order \code{psi_omega} (the cluster
#'   misspecification IF, built from the cluster-broadcast EIF \eqn{a_{g(i)}}); and
#'   (2) a per-CLUSTER second-order \code{var_add} \eqn{= (G/(G-1))\,2\hat Q},
#'   \eqn{\hat Q = -(Gn^2)^{-1}\sum_g a_g (d_g'\Psi_g)}, \eqn{d_g = -B v_g w},
#'   \eqn{v_g = (G/n^2)\Psi_g\Psi_g' - \widehat\Sigma_{cl}}. There is no separate
#'   \eqn{\Delta_{DF}}: the cohort demeaning leaves the single constraint
#'   \eqn{\sum_g\Psi_g = 0} at the cluster-sum level, which the leading SE's CR1
#'   factor \eqn{G/(G-1)} already restores. \code{s_vec} is then per-CLUSTER (length
#'   G) for the cross-cell increment. All terms reduce to the i.i.d. forms below at
#'   clusters==units (\eqn{G=n}, \eqn{\Psi_g=\psi_i}); validated by FD oracle (the
#'   \eqn{d_g} map) and MC calibration on clustered DGPs. For
#'   \code{omega_cov_shrink = "ledoit_wolf"} the plain map retains the leading
#'   optimism (B from the shrunk \eqn{\widehat\Sigma_{cl}}); only the LW-intensity
#'   data-dependence (the \eqn{d\lambda} chain) is omitted (a higher-order term with
#'   no production consumer).
#' @param mbar numeric vector length H, or \code{NULL}: the cell's moment vector
#'   (the long-difference contrasts \eqn{\bar m}, \code{== compute_generated_outcomes_nocov_edid()}).
#'   When supplied, the result also carries \code{psi_omega}, the FIRST-ORDER
#'   misspecification weight-estimation influence function
#'   \eqn{\psi_{\Omega,i} = (D\,\bar m)_i} (the no-covariate sibling of the
#'   covariate \code{psi_Omega}), for the caller to fold into the cell EIF.
#'
#' @return list with \code{applied} (logical), \code{warn} (logical: numeric
#'   failure vs structural skip), \code{reason} (string or NA),
#'   \code{var_add} (\eqn{\Delta_{DF} + 2\hat Q}: the additive variance applied to
#'   the cell), \code{delta_df} (\eqn{\Delta_{DF}}), \code{q_opt} (\eqn{\hat Q}),
#'   the diagnostics \code{cov_lead} (\eqn{= -\hat Q}) and \code{var_second}
#'   (J-linear \eqn{\widehat{\mathrm{Var}}(T_2)}), optionally \code{D}, and --
#'   when \code{mbar} is supplied -- the per-unit first-order misspecification IF
#'   \code{psi_omega} (\eqn{= D\,\bar m}; mean-zero, exactly zero under correct
#'   specification)
#' @keywords internal
compute_nocov_ee_correction_edid <- function(
  target_g, target_t, pairs, panel_obj, omega_raw, omega_used, weights,
  shrink_lambda = NA_real_, return_D = FALSE, mbar = NULL, cluster_indices = NULL
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

  # ---- CLUSTER METRIC BRANCH ------------------------------------------------------------------------
  # When cluster_indices is supplied the weights invert the CLUSTER moment covariance
  # Sig_cl = crossprod(rowsum(psi, cluster))/n^2 (== omega_raw/omega_used here), so the weight-estimation
  # channel is driven by the per-CLUSTER moment IF Psi_g = sum_{i in g} psi_i (rows of psi_cl), not psi_i.
  # Two outputs, both reducing to the IID forms below at clusters==units (G=n, Psi_g=psi_i):
  #  (1) FIRST-ORDER misspecification IF psi_omega (per-UNIT; folds into the cell EIF). The per-unit IF of
  #      Sig_cl is t_i = (1/n) psi_i Psi_{g(i)}' (so (1/n) sum_i t_i = Sig_cl exactly), giving
  #      psi_omega_i = -mbar'B(t_i - Sig_cl)w = -(1/n)(mbar'B psi_i) a_{g(i)} + mbar'B Sig_cl w, with
  #      a_{g(i)} = Psi_{g(i)}'w the unit's cluster-sum EIF. Mean-zero (sum_i psi_omega_i = 0) and exactly
  #      0 under correct specification (mbar in span(1) => D 1 = 0). This is the deck/paper channel
  #      (misspec_robust = TRUE; reached on the ridge/none path where shrink_lambda is NA).
  #  (2) SECOND-ORDER var_add (per-CLUSTER; additive to the cell SE). Clusters are the iid sampling units:
  #      per-cluster IF v_g = (G/n^2) Psi_g Psi_g' - Sig_cl, direction d_g = -B v_g w =
  #      -(G a_g/n^2) B Psi_g + B Sig_cl w; optimism Q = -(G n^2)^{-1} sum_g a_g (d_g' Psi_g). The reported
  #      leading SE already carries the cluster Bessel G/(G-1) (safe_inference_edid), and the cohort
  #      demeaning leaves exactly the single constraint sum_g Psi_g = 0 at the cluster-sum level, so that
  #      CR1 factor restores the lost df: there is NO separate Delta_DF, and the additive correction is
  #      var_add = (G/(G-1)) * 2 Q. (Validated by MC calibration on a clustered DGP.)
  if (!is.null(cluster_indices)) {
    # The PLAIN map (B from omega_used) is exact for omega_cov_shrink in {none, ridge} (ridge is constant in
    # Omega-hat, d/dOmega = I). For ledoit_wolf the plain map retains the LEADING optimism (B is built from the
    # actually-inverted shrunk Sig_cl); only the LW intensity's data-dependence (the dlambda chain term) is
    # omitted -- a higher-order refinement with no deck/paper consumer (documented, NEWS). The first-order
    # psi_omega is ALWAYS computed below, so the misspec_robust channel is never dropped.
    ci_idx <- as.integer(factor(cluster_indices))        # compact 1..G codes (row order of rowsum)
    psi_cl <- rowsum(psi, ci_idx)                        # G x H cluster sums Psi_g
    Gn     <- nrow(psi_cl)
    if (Gn < 2L) return(skip("fewer than 2 clusters: cluster weight-estimation channel undefined"))
    a_cl   <- drop(psi_cl %*% weights)                   # length G: cluster-sum EIFs a_g
    a_brd  <- a_cl[ci_idx]                               # length n: each unit's own-cluster EIF a_{g(i)}
    Bpsi   <- psi %*% B                                  # n x H  (B symmetric: rows (B psi_i)')
    BOw    <- drop(B %*% (omega_raw %*% weights))        # = B Sig_cl w  (omega_raw == Sig_cl here)
    cr1    <- Gn / (Gn - 1)
    # (1) first-order per-unit direction -> psi_omega
    D_fo   <- -((a_brd / n) * Bpsi - matrix(BOw, n, H, byrow = TRUE))            # n x H
    # (2) second-order per-cluster direction -> var_add
    Bpsi_cl <- psi_cl %*% B                              # G x H
    d_cl    <- -((Gn * a_cl / n^2) * Bpsi_cl - matrix(BOw, Gn, H, byrow = TRUE)) # G x H  (d_g)
    s_g      <- rowSums(d_cl * psi_cl)                   # d_g' Psi_g
    cov_lead <- sum(a_cl * s_g) / (Gn * n^2)             # diagnostic (J-linear Cov(T1,T2))
    q_opt    <- cr1 * (-cov_lead)                        # CR1-scaled optimism (matches the reported leading)
    delta_df <- 0                                        # cluster Bessel absorbed by the leading G/(G-1)
    var_add  <- delta_df + 2 * q_opt
    if (!is.finite(var_add)) return(skip("non-finite correction", warn = TRUE))
    out <- list(applied = TRUE, reason = NA_character_, warn = FALSE, var_add = var_add,
                delta_df = delta_df, q_opt = q_opt, cov_lead = cov_lead, var_second = NA_real_,
                s_vec = s_g)                             # per-CLUSTER s_g for the cross-cell increment
    if (!is.null(mbar)) out$psi_omega <- drop(D_fo %*% mbar)
    if (isTRUE(return_D)) out$D <- D_fo
    return(out)
  }

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
    S  <- compute_pole_structure_nocov_edid(target_g, target_t, pairs, panel_obj)
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
      # d(b^2)[dE] = -(2/n_eff)<Omega, dE>_F now (b^2 = q4/(n^3 n_eff) - ||Omega||^2/n_eff with q4
      # fixed). n_eff is the SAME Kish ESS the LW intensity uses (full n unweighted, so the n_full
      # argument makes this EXACTLY 2/(n*d2) at w-equal), keeping the chain rule consistent with the
      # lambda actually applied. The d-lambda-through-Omega channel (the 2*lam/d2 term) is unchanged
      # (it differentiates d^2, not b^2).
      n_eff   <- n_eff_edid(panel_obj$unit_weights,
                            active_mask_nocov_edid(target_g, pairs, panel_obj), n)
      kappa_i <- -(2 / (n_eff * d2)) * beta_i - (2 * lam / d2) * gamma_i              # dlambda[v_i]
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
  # cohort and the unbiased-covariance gap is the unit-level closed form below. The biased
  # group covariance E[Omega-hat_g] = (1 - sum_i p_i^2) Omega_g (p_i = w_i / W_g the within-
  # cohort weight share), so the gap fraction the plug-in misses is f_g = s2_g / (1 - s2_g),
  # s2_g = sum_{i in g} w_i^2 / W_g^2 (the weighted Bessel factor). Unweighted (w_i = 1):
  # s2_g = m_g / m_g^2 = 1/m_g and f_g = 1/(m_g - 1), reproducing the legacy closed form
  # sum_i a_i^2/(m_{gamma(i)} - 1) bit-for-bit. The fraction is constant within a cohort, so it
  # rides per unit exactly as 1/(m-1) did. m = 1 cohorts (s2 = 1) are guarded to f = 1 (their
  # variance is not estimable; the thin-cohort guard upstream normally prevents this).
  coh <- panel_obj$unit_cohorts
  uw  <- panel_obj$unit_weights
  if (is.null(uw)) {
    sizes <- table(coh)                                    # named by as.character (Inf included)
    m_i   <- as.numeric(sizes[as.character(coh)])
    f_i   <- 1 / pmax(m_i - 1, 1)
  } else {
    # per-cohort weighted Bessel fraction f_g = s2_g / (1 - s2_g), s2_g = sum w^2 / (sum w)^2
    W_g   <- tapply(uw, coh, sum)
    W2_g  <- tapply(uw * uw, coh, sum)
    s2_g  <- as.numeric(W2_g / (W_g * W_g))
    names(s2_g) <- names(W_g)
    f_g   <- ifelse(s2_g >= 1, 1, s2_g / (1 - s2_g))       # s2 = 1 => single-unit cohort guard
    f_i   <- as.numeric(f_g[as.character(coh)])
  }
  delta_df <- sum(a^2 * f_i) / n^2

  var_add <- delta_df + 2 * q_opt
  if (!is.finite(var_add)) return(skip("non-finite correction", warn = TRUE))

  out <- list(applied = TRUE, reason = NA_character_, warn = FALSE, var_add = var_add,
              delta_df = delta_df, q_opt = q_opt,
              cov_lead = cov_lead, var_second = var_second,
              # per-unit pieces for the CROSS-CELL assembly (nocov_ee_sigma_full_edid):
              # s_vec_i = d_i' psi_i; together with the cell EIF a_i these give the
              # cross-cell optimism in closed form without storing D or psi.
              s_vec = s_i)
  # First-order MISSPECIFICATION weight-estimation influence function (the no-covariate analogue of the
  # covariate psi_Omega channel). theta_w = w'mbar is the weighted pseudo-estimand; estimating Omega -> w
  # contributes psi_omega_i = mbar' J[phi_i] = (D %*% mbar)_i, with D the per-unit Jacobian directions above
  # and phi_i the per-unit IF of Omega-hat. It is mean-zero (sum_i d_i = 0) and EXACTLY zero under correct
  # specification (then mbar lies in span(1) and D %*% 1 = 0 by the sum-to-one weight constraint -- the
  # optimal-weight FOC), so it only bites under misspecification, where it restores coverage of theta_w.
  # Unlike var_add it is a genuine per-unit IF, so the caller folds it into the cell EIF and it propagates
  # to the clustered covariance, the aggregations, the sup-t bands, and the multiplier bootstrap.
  if (!is.null(mbar)) out$psi_omega <- drop(D %*% mbar)
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
#' @param unit_weights length-n vector of per-unit observation weights, or
#'   \code{NULL} (unweighted). Sets the cohort Bessel factor \eqn{f_i}: the
#'   unweighted \eqn{1/(m_{\gamma(i)}-1)} when \code{NULL}, else the weighted
#'   fraction \eqn{s^2_g/(1-s^2_g)} with \eqn{s^2_g = \sum w^2/(\sum w)^2}, the
#'   same weighted Bessel fraction as the per-cell \code{delta_df}.
#' @param cluster_indices length-n cluster id vector, or \code{NULL} (i.i.d.).
#'   When supplied, the increment switches to the CLUSTER metric: cell EIFs are
#'   cluster-summed to \eqn{a_g^c}, there is no \eqn{\Delta_{DF}} block, and the
#'   increment is the CR1-scaled cross-cell optimism
#'   \eqn{\Sigma_{cc'} = -\frac{1}{(G-1)n^2}\sum_g (s_g^c a_g^{c'} + a_g^c s_g^{c'})}.
#' @return K x K matrix, or NULL when no cell carries an applied correction
#' @keywords internal
nocov_ee_sigma_full_edid <- function(eif_matrix, s_matrix, unit_cohorts, unit_weights = NULL,
                                     cluster_indices = NULL) {
  if (is.null(eif_matrix) || is.null(s_matrix)) return(NULL)
  K <- ncol(s_matrix)
  ok <- which(vapply(seq_len(K), function(k) all(is.finite(s_matrix[, k])), logical(1L)))
  if (length(ok) == 0L) return(NULL)
  out <- matrix(0, K, K)
  if (!is.null(cluster_indices)) {
    # CLUSTER metric. s_matrix rows are CLUSTERS (the per-cluster s_g = d_g' Psi_g from
    # compute_nocov_ee_correction_edid); cluster-sum the cell EIFs to a_g^c. The reported leading
    # covariance already carries the cluster Bessel G/(G-1), and the cohort demeaning leaves the single
    # constraint sum_g Psi_g = 0 at the cluster-sum level, so there is NO Delta_DF block: the increment is
    # purely the CR1-scaled cross-cell optimism Sigma_cc' = -(1/((G-1)n^2)) sum_g (s_g^c a_g^{c'} +
    # a_g^c s_g^{c'}). Its diagonal reproduces the per-cell var_add = (G/(G-1)) 2 Q exactly.
    ci   <- as.integer(factor(cluster_indices))
    G    <- length(unique(ci))
    if (G < 2L) return(NULL)
    n    <- length(ci)
    A_cl <- rowsum(eif_matrix[, ok, drop = FALSE], ci)     # G x |ok|: cluster-sum EIFs a_g^c
    Sv   <- s_matrix[, ok, drop = FALSE]                    # G x |ok|: per-cluster s_g^c
    if (nrow(A_cl) != nrow(Sv)) return(NULL)
    sa   <- crossprod(Sv, A_cl) / ((G - 1) * n^2)
    out[ok, ok] <- -(sa + t(sa))
    return(out)
  }
  n <- nrow(s_matrix)
  # Per-unit Bessel factor f_i (constant within cohort): 1/(m-1) unweighted, s2/(1-s2) weighted
  # (s2_g = sum w^2 / (sum w)^2) -- the same weighted Bessel fraction as the per-cell delta_df.
  if (is.null(unit_weights)) {
    sizes <- table(unit_cohorts)
    f_i   <- 1 / pmax(as.numeric(sizes[as.character(unit_cohorts)]) - 1, 1)
  } else {
    W_g  <- tapply(unit_weights, unit_cohorts, sum)
    W2_g <- tapply(unit_weights * unit_weights, unit_cohorts, sum)
    s2_g <- as.numeric(W2_g / (W_g * W_g)); names(s2_g) <- names(W_g)
    f_g  <- ifelse(s2_g >= 1, 1, s2_g / (1 - s2_g))
    f_i  <- as.numeric(f_g[as.character(unit_cohorts)])
  }
  A  <- eif_matrix[, ok, drop = FALSE]
  Sv <- s_matrix[, ok, drop = FALSE]
  out <- matrix(0, K, K)
  df_blk  <- crossprod(A * f_i, A) / n^2                       # Delta_DF block (A' diag(f) A)
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

  # Observation weights (NULL => wmean_edid reduces to mean() bit-for-bit).
  uw     <- panel_obj$unit_weights
  w_g    <- if (is.null(uw)) NULL else uw[mask_g]
  w_inf  <- if (is.null(uw)) NULL else uw[mask_inf]

  col_t  <- .col(panel_obj, target_t)

  if (pt_assumption == "post") {
    # PT-Post: one pair (Inf, base); weighted (Hajek) group means
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)
    y_hat <- wmean_edid(ow[mask_g, col_t] - ow[mask_g, col_base], w_g) -
             wmean_edid(ow[mask_inf, col_t] - ow[mask_inf, col_base], w_inf)
    return(y_hat)  # scalar; will be treated as length-1 vector
  }

  # PT-All
  col_1 <- .col(panel_obj, panel_obj$period_1)
  term_g <- wmean_edid(ow[mask_g, col_t] - ow[mask_g, col_1], w_g)

  y_hat <- numeric(H)
  for (j in seq_len(H)) {
    gp_j    <- pairs$gp[j]
    tpre_j  <- pairs$tpre[j]
    col_pre <- .col(panel_obj, tpre_j)

    term_inf <- wmean_edid(ow[mask_inf, col_t] - ow[mask_inf, col_pre], w_inf)

    mask_gp  <- panel_obj$cohort_masks[[as.character(gp_j)]]
    w_gp     <- if (is.null(uw)) NULL else uw[mask_gp]
    term_gp  <- wmean_edid(ow[mask_gp, col_pre] - ow[mask_gp, col_1], w_gp)

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

  # Observation weights. Each demeaned group contribution carries a w_i factor and a weighted
  # (Hajek) centering; pi = W/n. With NULL weights the multiplicand w is all-ones and the
  # centering is the plain mean, so the EIF is byte-identical to the legacy form. By construction
  # this EIF equals compute_psi_moments_nocov_edid() %*% weights (mirrored term by term).
  uw       <- panel_obj$unit_weights
  w_g_obs  <- if (is.null(uw)) rep(1, sum(mask_g))   else uw[mask_g]
  w_inf_obs<- if (is.null(uw)) rep(1, sum(mask_inf)) else uw[mask_inf]
  pi_inf   <- if (is.null(uw)) sum(mask_inf) / n else sum(w_inf_obs) / n

  col_t <- .col(panel_obj, target_t)

  eif <- numeric(n)

  if (pt_assumption == "post") {
    # PT-Post: one pair; base = g - 1 - anticipation
    tpre_val <- pairs$tpre[1L]
    col_base <- .col(panel_obj, tpre_val)
    w_1      <- weights[1L]

    delta_g   <- ow[mask_g,   col_t] - ow[mask_g,   col_base]
    mean_g    <- wmean_edid(delta_g, if (is.null(uw)) NULL else w_g_obs)
    eif[mask_g] <- eif[mask_g] + w_1 * w_g_obs * (delta_g - mean_g) / pi_g

    delta_inf  <- ow[mask_inf, col_t] - ow[mask_inf, col_base]
    mean_inf   <- wmean_edid(delta_inf, if (is.null(uw)) NULL else w_inf_obs)
    eif[mask_inf] <- eif[mask_inf] - w_1 * w_inf_obs * (delta_inf - mean_inf) / pi_inf

    # gp_j == Inf: no comparison cohort term
  } else {
    # PT-All
    col_1    <- .col(panel_obj, panel_obj$period_1)
    col_base <- col_1   # base for treated group

    delta_g_t_base <- ow[mask_g, col_t] - ow[mask_g, col_base]
    mean_g_t_base  <- wmean_edid(delta_g_t_base, if (is.null(uw)) NULL else w_g_obs)

    for (j in seq_len(H)) {
      w_j     <- weights[j]
      gp_j    <- pairs$gp[j]
      tpre_j  <- pairs$tpre[j]
      col_pre <- .col(panel_obj, tpre_j)

      # Treated group contribution (always present)
      # Y_hat_j enters via centering: phi_{ij} = I(G=g) w_i/pi_g * (delta - Y_hat_j)
      # We accumulate sum_j w_j * I(G=g) w_i/pi_g * (delta - wmean) directly.
      eif[mask_g] <- eif[mask_g] +
        w_j * w_g_obs * (delta_g_t_base - mean_g_t_base) / pi_g

      # Never-treated contribution (subtract)
      delta_inf_t_pre <- ow[mask_inf, col_t] - ow[mask_inf, col_pre]
      mean_inf_t_pre  <- wmean_edid(delta_inf_t_pre, if (is.null(uw)) NULL else w_inf_obs)
      eif[mask_inf] <- eif[mask_inf] -
        w_j * w_inf_obs * (delta_inf_t_pre - mean_inf_t_pre) / pi_inf

      # Comparison cohort contribution (subtract)
      mask_gp          <- panel_obj$cohort_masks[[as.character(gp_j)]]
      pi_gp            <- panel_obj$cohort_fractions[[as.character(gp_j)]]
      w_gp_obs         <- if (is.null(uw)) rep(1, sum(mask_gp)) else uw[mask_gp]
      delta_gp_pre_base <- ow[mask_gp, col_pre] - ow[mask_gp, col_base]
      mean_gp_pre_base  <- wmean_edid(delta_gp_pre_base, if (is.null(uw)) NULL else w_gp_obs)
      eif[mask_gp] <- eif[mask_gp] -
        w_j * w_gp_obs * (delta_gp_pre_base - mean_gp_pre_base) / pi_gp
    }
  }

  # The score accumulated above already has mean 0 by construction:
  # each group's contribution is demeaned (delta_i - mean(delta)).
  # The EIF is the score itself; do NOT subtract att_gt again.
  eif
}
