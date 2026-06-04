# edid-fit.R
# Outer (g, t) cell loop for the EDiD estimator.

#' Fit all (g, t) cells for the EDiD estimator
#'
#' Iterates over all treatment cohorts and all time periods (excluding
#' \code{period_1}), computes point estimates, EIFs, and analytical SEs for
#' each cell.
#'
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param pt_assumption character: \code{"all"} or \code{"post"}
#' @param alpha significance level in (0, 1)
#' @param store_eif logical: if TRUE, include EIF vectors in returned cells
#' @param xformla one-sided formula or NULL: covariate formula (routed to
#'   covariate path when non-trivial and \code{panel_obj$covariate_matrix}
#'   is non-NULL)
#' @param need_eif logical: if TRUE, always store EIF regardless of store_eif
#'   (used internally when \code{n_bootstrap > 0})
#'
#' @return list with elements:
#'   \describe{
#'     \item{\code{cells}}{list of \code{edid_cell_result} objects}
#'     \item{\code{eif_matrix}}{n x n_valid_cells numeric matrix, or NULL}
#'     \item{\code{cell_index}}{data.frame: group, time, cell_id, is_pre}
#'   }
#' @keywords internal
fit_edid_cells <- function(
  panel_obj, pt_assumption, alpha, store_eif, xformla = NULL, seed = NULL,
  need_eif = FALSE, weight_method = c("efficient", "averaged", "gmm", "uniform"),
  estimation_effect = FALSE, higher_order = FALSE, misspec_robust = FALSE
) {
  weight_method <- match.arg(weight_method)
  higher_order  <- isTRUE(higher_order)
  misspec_robust <- isTRUE(misspec_robust)
  # Determine if covariate path is active
  is_trivial_xformla <- is.null(xformla) ||
    identical(deparse(xformla, width.cutoff = 500L), "~1")
  use_cov_path <- !is_trivial_xformla && !is.null(panel_obj$covariate_matrix)

  # Curse of dimensionality guard for the pointwise efficient weights. The kernel
  # Omega*(X) is consistent only while its local effective sample size n * prod(h_k) ~
  # n^{(5-d)/5} grows, i.e. d < 5. At d >= 5 the admissible floor band (0, (5-d)/10) is
  # empty, Omega_hat*(X) is not consistent, and the regularized weights collapse toward
  # uniform (the estimator stays CONSISTENT via the valid DR moments, but is no longer
  # semiparametrically efficient). Recommend the constant-weight schemes, whose pooled
  # Omega-bar is sqrt(n)-consistent with no curse of dimensionality.
  if (weight_method == "efficient" && use_cov_path) {
    d_cov <- ncol(panel_obj$covariate_matrix)
    if (d_cov >= 5L) {
      warning(sprintf(paste0(
        "weights='efficient' with %d continuous covariates: the pointwise kernel ",
        "Omega*(X) is not consistently estimable (curse of dimensionality; admissible ",
        "floor band is empty for d>=5), so the efficient weights collapse toward uniform ",
        "and the efficiency guarantee is lost (ATT stays consistent). Use ",
        "weights='averaged' (sqrt(n)-consistent, no curse) or condition on a low-",
        "dimensional index (e.g. the propensity score)."), d_cov), call. = FALSE)
    }
  }

  # The "gmm" scheme inverts the unconditional covariance of the generated outcomes, which is
  # estimated from the same data as the moments. This induces a finite-sample two-step bias of
  # order O(1/n) that can be sizeable with strong treatment effects, many moments, or small
  # cohorts. Asymptotically it coincides with "averaged" and is never more efficient, so the
  # constant-weight default "averaged" (or the asymptotically efficient "efficient") is preferred.
  if (weight_method == "gmm") {
    warning(paste0(
      "weights = 'gmm' inverts the unconditional moment covariance and carries a finite-sample ",
      "O(1/n) bias that can be large under strong effects, many periods, or small cohorts; it is ",
      "asymptotically equal to 'averaged' but never more efficient. Prefer 'averaged' or 'efficient'."),
      call. = FALSE)
  }

  # Nuisance functions are estimated by plug-in (train = test = full sample), as in the paper's
  # remark that the efficient-influence-function moments are Neyman orthogonal, so the first-step
  # nuisance estimates do not affect the first-order asymptotic variance. Cross-fitting (K > 1) is
  # not used: the estimated weights and nuisances are already first-order negligible, while
  # sample-splitting can inflate the finite-sample variance of the just-identified DR estimator.
  K_use <- 1L
  if (use_cov_path) {
    fold_id <- rep(1L, panel_obj$n)
  } else {
    fold_id <- NULL
  }

  # The ACH (Ackerberg, Chen & Hahn 2012) first-step correction requires the plug-in
  # (train = test = full sample) M-estimator pieces; it is derived for K = 1. Guard
  # defensively in case cross-fitting is ever enabled upstream.
  if (isTRUE(estimation_effect)) {
    if (!use_cov_path) {
      warning("estimation_effect has no effect without covariates (no first-step nuisances).", call. = FALSE)
      estimation_effect <- FALSE
    } else if (K_use > 1L) {
      stop("estimation_effect = TRUE is only supported with plug-in nuisances (K = 1).", call. = FALSE)
    }
  }

  # The higher-order ("Wick") variance refinement reuses the SAME plug-in first-step M-estimator pieces
  # (the per-nuisance B / score / H_inv that feed the joint coefficient covariance V and the per-cell
  # Hessian). It is meaningful only on the covariate path -- with no covariates the nuisances are
  # unconditional means with no sieve coefficients, so the higher-order term is exactly zero -- and, like
  # the ACH correction, is derived for K = 1. edid() enforces the covariate-path requirement (stops on
  # xformla = NULL); guard defensively here too.
  if (higher_order) {
    if (!use_cov_path) {
      warning("higher_order has no effect without covariates (no first-step sieve nuisances).", call. = FALSE)
      higher_order <- FALSE
    } else if (K_use > 1L) {
      stop("higher_order = TRUE is only supported with plug-in nuisances (K = 1).", call. = FALSE)
    }
  }
  # The misspecification-robust SE augments the EIF with the weight-estimation channel psi_Omega (the sibling
  # of the ACH nuisance correction that estimation_effect explicitly leaves out). It needs estimated weights
  # (no covariates => weights not estimated from X; uniform => fixed weights, channel is exactly zero) and the
  # plug-in M-estimator aux (K = 1), like the ACH correction.
  if (misspec_robust) {
    if (!use_cov_path) {
      warning("misspec_robust has no effect without covariates (weights are not estimated from X).", call. = FALSE)
      misspec_robust <- FALSE
    } else if (weight_method == "uniform") {
      warning("misspec_robust has no effect for weights = 'uniform' (fixed weights have no estimation channel).",
              call. = FALSE)
      misspec_robust <- FALSE
    } else if (K_use > 1L) {
      stop("misspec_robust = TRUE is only supported with plug-in nuisances (K = 1).", call. = FALSE)
    }
  }
  # The first-step aux pieces (B / score / H_inv) are needed whenever EITHER the ACH correction OR the
  # higher-order refinement is requested -- and (experimental) for the gmm weight-channel correction, whose
  # quadratic moment u'Cw inherits the (r, m) nuisance estimation.
  want_aux <- isTRUE(estimation_effect) || higher_order ||
    ((isTRUE(getOption("edid_store_psiomega")) || misspec_robust) && weight_method == "gmm")

  tgroups   <- panel_obj$treatment_groups
  tperiods  <- panel_obj$time_periods
  period_1  <- panel_obj$period_1
  n         <- panel_obj$n

  # Periods to iterate over: all except period_1
  iter_periods <- tperiods[tperiods != period_1]

  # Pre-allocate cell list
  n_cells <- length(tgroups) * length(iter_periods)
  cells   <- vector("list", n_cells)

  # cell_index tracking
  ci_group  <- numeric(n_cells)
  ci_time   <- numeric(n_cells)
  ci_cell_id <- integer(n_cells)
  ci_is_pre  <- logical(n_cells)

  keep_eif <- store_eif || need_eif
  eif_list <- if (keep_eif) vector("list", n_cells) else NULL

  cell_id <- 0L
  n_extreme_ratio_instances <- 0L  # accumulate extreme-ratio warnings; emit once at end

  # Hoist the CELL-INVARIANT Nadaraya-Watson kernel: the n x n weight matrix K and the bandwidths depend only on
  # the full covariate matrix, so they are identical for every (g,t) cell. Build ONCE and reuse, instead of
  # rebuilding (an O(d*n^2) outer/dnorm + an n x n allocation) inside compute_omega_star_cov_edid on every cell --
  # the dominant cost of the covariate path. Numerically identical; only the covariate path needs it.
  kern_bw <- NULL; kern_K <- NULL
  if (use_cov_path) {
    kk_full <- build_kernel_weights_edid(panel_obj$covariate_matrix)
    kern_bw <- kk_full$bw; kern_K <- kk_full$K
  }

  for (g in tgroups) {
    for (t in iter_periods) {
      cell_id <- cell_id + 1L
      is_pre  <- (t < g)

      # Step 1: enumerate valid pairs
      pairs <- enumerate_valid_pairs_edid(
        target_g          = g,
        treatment_groups  = tgroups,
        time_periods      = tperiods,
        period_1          = period_1,
        pt_assumption     = pt_assumption,
        anticipation      = panel_obj$anticipation
      )

      # NA cell if no valid pairs
      if (nrow(pairs) == 0L) {
        cells[[cell_id]] <- list(
          group           = g,
          time            = t,
          att             = NA_real_,
          se              = NA_real_,
          ci_lower        = NA_real_,
          ci_upper        = NA_real_,
          t_stat          = NA_real_,
          p_value         = NA_real_,
          n_pairs         = 0L,
          weights         = NULL,
          condition_num   = NA_real_,
          is_pre          = is_pre,
          inference_valid = FALSE,
          eif             = NULL,
          ho              = NULL
        )
        ci_group[cell_id]   <- g
        ci_time[cell_id]    <- t
        ci_cell_id[cell_id] <- cell_id
        ci_is_pre[cell_id]  <- is_pre
        if (keep_eif) eif_list[[cell_id]] <- rep(NA_real_, n)
        next
      }

      # Covariate nuisance estimates (per cell)
      prop_ratios  <- NULL
      cond_means   <- NULL
      inv_propensities <- NULL
      ho_cell      <- NULL   # higher-order per-cell pieces (nuisance blocks + Hessian), set on cov path

      if (use_cov_path) {
        # Build nuisance estimation pairs.
        # For Eq. (4.4), we need nuisances for:
        #   - r[g, Inf] always (never-treated propensity ratio)
        #   - r[g, gp] for each cross-cohort gp
        #   - m_{Inf, t} and m_{Inf, tpre} for never-treated conditional means
        #   - m_{gp, tpre} for each cross-cohort gp's pretrend
        # Build pairs_for_nuisance that includes Inf + all cross-cohort gps
        pairs_for_nuisance <- pairs
        # Self-comparison pairs use Inf as comparison
        self_cmp <- is.finite(pairs_for_nuisance$gp) & (pairs_for_nuisance$gp == g)
        if (any(self_cmp)) {
          pairs_for_nuisance$gp[self_cmp] <- Inf
        }
        # Ensure Inf is always present for never-treated nuisances
        # (needed even for cross-cohort pairs)
        cross_pairs <- pairs[is.finite(pairs$gp) & pairs$gp != g, , drop = FALSE]
        if (nrow(cross_pairs) > 0L) {
          # Add Inf-based pairs for the same tpre values (for m_{Inf,tpre})
          inf_pairs <- data.frame(gp = Inf, tpre = unique(cross_pairs$tpre))
          pairs_for_nuisance <- rbind(pairs_for_nuisance, inf_pairs)
          pairs_for_nuisance <- unique(pairs_for_nuisance)
        }

        prop_ratios <- withCallingHandlers(
          estimate_all_propensity_ratios(
            panel_obj = panel_obj,
            g         = g,
            pairs     = pairs_for_nuisance,
            bs_df    = 4L,
            K_folds  = K_use,
            fold_id  = fold_id,
            return_aux = want_aux
          ),
          warning = function(w) {
            if (grepl("Extreme propensity ratios", conditionMessage(w), fixed = TRUE)) {
              n_extreme_ratio_instances <<- n_extreme_ratio_instances + 1L
              invokeRestart("muffleWarning")
            }
          }
        )
        cond_means <- estimate_all_conditional_means(
          panel_obj = panel_obj,
          pairs     = pairs_for_nuisance,
          t_val    = t,
          bs_df    = 4L,
          K_folds  = K_use,
          fold_id  = fold_id,
          return_aux = want_aux
        )
        # Split predictions from the first-step aux pieces (the *_aux return path wraps both). Requested by
        # the ACH correction (estimation_effect) and/or the higher-order refinement (higher_order).
        r_aux <- NULL; m_aux <- NULL
        if (want_aux) {
          r_aux <- prop_ratios$aux; prop_ratios <- prop_ratios$predictions
          m_aux <- cond_means$aux;  cond_means  <- cond_means$predictions
        }
        inv_propensities <- estimate_all_inverse_propensities(
          panel_obj = panel_obj,
          g         = g,
          pairs     = pairs,
          bs_df     = 4L,
          K_folds   = K_use,
          fold_id   = fold_id,
          return_aux = (isTRUE(getOption("edid_store_psiomega")) || misspec_robust) &&  # inv_p M-pieces for
                       weight_method %in% c("averaged", "efficient")  # the Sigma_Omega corr channel (psi consumers)
        )
      }

      # Steps 2-6: dispatch on covariate vs. no-covariate path
      if (use_cov_path) {
        # --- Covariate path ---
        gen_out_mat <- compute_generated_outcomes_cov_edid(
          panel_obj     = panel_obj,
          g             = g,
          t             = t,
          pairs         = pairs,
          prop_ratios   = prop_ratios,
          cond_means    = cond_means,
          pt_assumption = pt_assumption
        )
        if (isTRUE(getOption("edid_store_genout"))) {            # diagnostic: expose per-cell generated outcomes
          acc <- getOption("edid_genout_acc", list()); acc[[paste0(g, "_", t)]] <- gen_out_mat
          options(edid_genout_acc = acc)
        }
        if (weight_method == "efficient") {
          # Paper's pointwise efficient weights w(X_i)=Omega*(X_i)^{-1}1/(1'Omega*(X_i)^{-1}1),
          # with a dimension-aware eigenvalue-floor regularization (a=0.7*(5-d)/10) that is
          # asymptotically negligible yet dominates the NW estimation noise for stability.
          omega_arr <- compute_omega_star_cov_edid(panel_obj, g, t, pairs,
                                                   prop_ratios, cond_means,
                                                   inv_propensities, bw = kern_bw, K_mat = kern_K,
                                                   return_pointwise = TRUE)
          # Fuse the per-unit adjoint q_i into the weights' eigen pass when the weight-estimation channel is needed
          # (one eigendecomposition per unit instead of two). Q_pw is from the un-frozen weights; the .fwpw freeze
          # (research jackknife) is incompatible with the psi channel and is guarded with a stop below.
          want_q   <- isTRUE(getOption("edid_store_psiomega")) || misspec_robust
          pw_res   <- compute_pointwise_weights_edid(omega_arr, d = ncol(panel_obj$covariate_matrix),
                        gen_out_mat = if (want_q) gen_out_mat else NULL)        # n x H (+ Q_pw when want_q)
          if (want_q) { W_pw <- pw_res$W; Q_pw <- pw_res$Q } else W_pw <- pw_res
          # Research hook (default OFF, NOT on the PR): getOption("edid_fixed_wpw") = list("g_t" = list(ids=, W=))
          # FREEZES the per-unit pointwise weights at supplied (id-keyed) values instead of re-estimating them. Used by
          # the efficient weight-estimation-channel jackknife (refit nuisances but FREEZE W_pw) to isolate Sigma_Omega.
          .fwpw <- getOption("edid_fixed_wpw", NULL)
          .fwpw <- if (!is.null(.fwpw)) .fwpw[[paste0(g, "_", t)]] else NULL
          if (!is.null(.fwpw)) {
            mi <- match(panel_obj$all_units, .fwpw$ids)
            if (!anyNA(mi) && ncol(.fwpw$W) == ncol(W_pw)) W_pw <- .fwpw$W[mi, , drop = FALSE]
            else warning(sprintf("edid_fixed_wpw: freeze skipped for cell (%s,%s) (unmatched id or H change); W_pw re-estimated.",
                                 g, t), call. = FALSE)
          }
          cond_num <- tryCatch(check_condition_edid(apply(omega_arr, c(2, 3), mean)),
                               error = function(e) NA_real_)
          wY_i     <- rowSums(gen_out_mat * W_pw)
          if (anyNA(wY_i))                                                    # NA would split the support:
            stop(sprintf(paste0("edid cell (%s,%s): NA in weighted generated outcomes (typically a ",
              "non-invertible local covariance or a missing covariate prediction). The point estimate ",
              "would use complete cases while the EIF spans all units, breaking the empirical mean-zero ",
              "identity and invalidating the SE. edid requires complete cases: drop incomplete units ",
              "before calling, or inspect cell (%s,%s)."), g, t, g, t), call. = FALSE)
          att_gt   <- mean(wY_i, na.rm = TRUE)                               # E_n[w(X_i)' Ytilde_i]
          eif_gt   <- compute_eif_cov_edid(panel_obj, gen_out_mat, W_pw, att_gt, g)
          if (isTRUE(estimation_effect))                                    # ACH first-step correction (W frozen)
            eif_gt <- eif_gt - compute_ach_correction_cov_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, W_pw, m_aux, r_aux, pt_assumption)
          if (higher_order)                                                  # per-cell Hessian (W = W_pw frozen)
            ho_cell <- compute_cell_hessian_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, W_pw, m_aux, r_aux, pt_assumption)
          weights  <- colMeans(W_pw, na.rm = TRUE)                            # store mean weight per pair
          if (isTRUE(getOption("edid_store_wpw"))) {            # research diagnostic (OFF): full per-unit W_pw + ids
            acc <- getOption("edid_wpw_acc", list())            # (seeds the efficient frozen-W_pw jackknife)
            acc[[paste0(g, "_", t)]] <- list(ids = panel_obj$all_units, W = W_pw)
            options(edid_wpw_acc = acc)
          }
          # Weight-estimation channel psi_Omega. Computed when EITHER the research diagnostic (edid_store_psiomega)
          # OR the production misspec_robust SE is requested; STORED to the accumulator only for the diagnostic;
          # FOLDED into eif_gt (eif_gt + psi_Omega, mirroring the ACH subtraction above) only for misspec_robust.
          if (isTRUE(getOption("edid_store_psiomega")) || misspec_robust) {   # efficient pointwise Sigma_Omega
            # The pointwise weight-estimation IF is the per-unit five-term kernel IF: term_psi with per-unit coupling
            # coup_i = Q[i,j] W_pw[i,k] (q_i'1 = 0 per unit, so Term 1 cancels pointwise), plus the same analytic inv_p
            # correction (the per-unit coupled_C). wY_i is already NA-checked (stop above), so gen_out_mat is complete.
            if (K_use > 1L)
              stop("the misspec_robust / Sigma_Omega channel requires plug-in nuisances (K = 1).", call. = FALSE)
            if (!is.null(.fwpw))                                              # frozen W is not Minv-consistent =>
              stop("edid_store_psiomega is incompatible with edid_fixed_wpw (frozen W_pw breaks q_i'1 = 0).",
                   call. = FALSE)                                            # ...the Term-1 cancellation is invalid
            # Q_pw was computed in the fused weights pass above (q_i'1 = 0 per unit, same Minv as w_i).
            stopifnot(max(abs(rowSums(Q_pw))) < 1e-6 * (1 + max(abs(Q_pw))))  # defensive: pointwise Term-1 premise
            lam_cell <- attr(omega_arr, "shrink_lambda")                      # shrinkage intensity (psi omits its O(lam) IF)
            if (is.finite(lam_cell) && lam_cell > 0.05)                       # large lambda => omission non-negligible
              warning(sprintf("edid cell (%s,%s): efficient weight-channel SE omits the shrinkage IF and lambda=%.3f > 0.05; the misspec-robust SE is only an approximation to the weight-channel variance here.",
                              g, t, lam_cell), call. = FALSE)
            po    <- compute_omega_star_cov_edid(panel_obj, g, t, pairs, prop_ratios, cond_means,
                       inv_propensities, bw = kern_bw, K_mat = kern_K, return_pointwise = TRUE,
                       psi_qw = list(pointwise = TRUE, Q = Q_pw, W = W_pw))
            corr_an <- compute_invp_correction_analytic_cov_edid(panel_obj$n, attr(inv_propensities, "aux"),
                                                                 po$coupled_C)
            psi_i <- po$psi - corr_an                                         # weight-estimation IF (data - corr)
            if (isTRUE(getOption("edid_store_psiomega"))) {                   # research diagnostic accumulator
              acc <- getOption("edid_psiomega_acc", list())
              acc[[paste0(g, "_", t)]] <- list(data = po$psi, corr = corr_an, lambda = lam_cell)
              options(edid_psiomega_acc = acc)
            }
            if (misspec_robust) eif_gt <- eif_gt + psi_i                      # production: fold the channel into the EIF
          }
          # Diagnostic only (default off): per-component cross-unit SD of the pointwise weights
          # quantifies how much Omega*(X) shape-varies; ~0 => weights ~constant => efficient ~ averaged.
          if (isTRUE(getOption("edid_diag_wpw")))
            message(sprintf("WPWDIAG %d_%d wsd=[%s] meanw=[%s] maxw=%.3f",
              g, t, paste(round(apply(W_pw, 2, stats::sd, na.rm = TRUE), 4), collapse = ","),
              paste(round(colMeans(W_pw, na.rm = TRUE), 3), collapse = ","), max(abs(W_pw), na.rm = TRUE)))
        } else {
          # Constant-weight schemes (valid but not pointwise-efficient):
          #   averaged = invert kernel Omega-bar; gmm = invert unconditional S_hat; uniform = 1/H.
          psi_const <- NULL   # weight-estimation IF (averaged kernel channel / gmm sample-cov channel); NULL => fold +0
          omega    <- compute_omega_star_cov_edid(panel_obj, g, t, pairs,
                                                  prop_ratios, cond_means,
                                                  inv_propensities, bw = kern_bw, K_mat = kern_K)
          if (isTRUE(getOption("edid_store_omega"))) {       # research diagnostic (default OFF, NOT on the PR):
            arr <- compute_omega_star_cov_edid(panel_obj, g, t, pairs, prop_ratios, cond_means,  # per-unit Omega*(X_i)
                     inv_propensities, bw = kern_bw, K_mat = kern_K, return_pointwise = TRUE)     # array for the averaged
            acc <- getOption("edid_omega_acc", list()); acc[[paste0(g, "_", t)]] <- list(omega = omega, arr = arr)
            options(edid_omega_acc = acc)                    # Sigma_Omega derivation (IF_Omegabar(i) = Omega*(X_i) - Omegabar)
          }
          cond_num <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)
          H_local  <- nrow(omega)
          # Research hook (default OFF, NOT committed to the PR): getOption("edid_fixed_weights") =
          # list("g_t" = weight vector) FREEZES the constant weight at a supplied value instead of re-estimating it.
          # Used by the weight-estimation-channel bootstrap (resample + refit nuisances but FREEZE W) to isolate
          # Sigma_Omega vs a refit-all bootstrap. Inert unless the option is set; falls back to estimation if H changed.
          .fw <- getOption("edid_fixed_weights", NULL); .fw <- if (!is.null(.fw)) .fw[[paste0(g, "_", t)]] else NULL
          weights  <- if (!is.null(.fw) && length(.fw) == H_local) .fw else switch(weight_method,
            uniform = rep(1 / H_local, H_local),
            # pairwise.complete.obs so a single NA moment does not NA-poison the whole
            # covariance (compute_efficient_weights_edid then guards any residual non-finite).
            gmm     = compute_efficient_weights_edid(stats::cov(gen_out_mat, use = "pairwise.complete.obs")),
            compute_efficient_weights_edid(omega))   # "averaged"
          att_gt   <- sum(weights * colMeans(gen_out_mat, na.rm = TRUE))
          if ((isTRUE(getOption("edid_store_psiomega")) || misspec_robust) && weight_method == "averaged") {
            # Sigma_Omega = psi_data (kernel five-term IF of Omegabar) - corr (inv_p prefactor channel). The IF holds
            # only when q and w share omega's inverse (=> q'1 = 0, so Eq.(3.12) Term 1 cancels and is rightly skipped).
            # Skip cells where the stored weights came from a DIFFERENT construction -- frozen weights
            # (edid_fixed_weights), the pseudoinverse/uniform fallback, or an NA-incomplete gen_out_mat (complete-case
            # mbar would be inconsistent with the all-units kernel terms). One w_chk match covers all three; a skip
            # leaves psi_const = NULL so the misspec_robust fold adds +0 (that cell falls back to the plug-in SE).
            if (K_use > 1L)
              stop("the misspec_robust / Sigma_Omega channel requires plug-in nuisances (K = 1).", call. = FALSE)
            C_inv <- NULL
            if (is.null(.fw) && !anyNA(gen_out_mat)) {
              kappa <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)
              C_inv <- if (all(omega == 0) || any(!is.finite(omega))) NULL
                       else if (!is.finite(kappa) || kappa > EDID_COND_THRESH) compute_pseudoinverse_edid(omega)
                       else tryCatch(solve(omega), error = function(e) compute_pseudoinverse_edid(omega))
            }
            w_chk <- if (is.null(C_inv)) NULL else {                          # the efficient weights from THIS inverse
              num <- drop(C_inv %*% rep(1, ncol(gen_out_mat))); d <- sum(num)
              if (!is.finite(d) || abs(d) < EDID_DENOM_EPS) NULL else num / d }
            if (!is.null(w_chk) && max(abs(weights - w_chk)) < 1e-8) {        # weights ARE those efficient weights
              mbar  <- colMeans(gen_out_mat)
              q_vec <- drop(C_inv %*% (mbar - att_gt))                        # same inverse as the weights => q'1 = 0
              stopifnot(abs(sum(q_vec)) < 1e-6 * (1 + max(abs(q_vec))))       # Term-1 cancellation premise (defensive)
              po    <- compute_omega_star_cov_edid(panel_obj, g, t, pairs, prop_ratios, cond_means,
                         inv_propensities, bw = kern_bw, K_mat = kern_K, psi_qw = list(q = q_vec, w = weights))
              invp_aux <- attr(inv_propensities, "aux")
              corr_an  <- compute_invp_correction_analytic_cov_edid(panel_obj$n, invp_aux, po$coupled_C)  # optimized
              psi_const <- po$psi - corr_an                                  # captured for the misspec_robust fold
              if (isTRUE(getOption("edid_store_psiomega"))) {                # research diagnostic accumulator
                corr_fd  <- if (isTRUE(getOption("edid_psiomega_fd")))        # FD oracle (validation only)
                  compute_invp_correction_cov_edid(panel_obj, g, t, pairs, prop_ratios, cond_means,
                    inv_propensities, invp_aux, weights, mbar, bw = kern_bw, K_mat = kern_K) else NULL
                acc <- getOption("edid_psiomega_acc", list())
                acc[[paste0(g, "_", t)]] <- list(data = po$psi, corr = corr_an, corr_fd = corr_fd)  # = data - corr
                options(edid_psiomega_acc = acc)
              }
            }
          }
          if ((isTRUE(getOption("edid_store_psiomega")) || misspec_robust) && weight_method == "gmm") {
            # gmm weight-estimation IF: the gmm weight inverts C = cov(Ytilde), the unconditional sample covariance.
            # The channel is the sample-cov two-step IF psi_plug = -(d.u)(d.w) + u'Cw (d_i = Ytilde_i - mbar,
            # u = (C^-1 - w 1'C^-1)'mbar) PLUS the ACH correction for the quadratic moment u'Cw -- because a covariance
            # is NOT protected by the att moment's Neyman orthogonality, C inherits the first-step (r, m) nuisance
            # estimation. psi_Omega = psi_plug + nuis_corr (sign jackknife-locked, cor 1.000). w_chk gate as elsewhere.
            if (K_use > 1L)
              stop("the misspec_robust / Sigma_Omega channel requires plug-in nuisances (K = 1).", call. = FALSE)
            if (is.null(.fw) && !anyNA(gen_out_mat)) {
              Cmat  <- stats::cov(gen_out_mat, use = "pairwise.complete.obs")     # the SAME C the gmm weights invert
              C_inv <- if (any(!is.finite(Cmat))) NULL else {
                kappa <- tryCatch(check_condition_edid(Cmat), error = function(e) NA_real_)
                if (!is.finite(kappa) || kappa > EDID_COND_THRESH) compute_pseudoinverse_edid(Cmat)
                else tryCatch(solve(Cmat), error = function(e) compute_pseudoinverse_edid(Cmat)) }
              w_chk <- if (is.null(C_inv)) NULL else {
                num <- drop(C_inv %*% rep(1, ncol(gen_out_mat))); dd <- sum(num)
                if (!is.finite(dd) || abs(dd) < EDID_DENOM_EPS) NULL else num / dd }
              if (!is.null(w_chk) && max(abs(weights - w_chk)) < 1e-8) {
                mbar <- colMeans(gen_out_mat)
                B    <- C_inv - outer(weights, drop(crossprod(rep(1, length(weights)), C_inv)))  # C^-1 - w 1'C^-1
                u    <- drop(crossprod(B, mbar))
                dc   <- sweep(gen_out_mat, 2L, mbar, "-")
                psi_plug  <- -(as.numeric(dc %*% u) * as.numeric(dc %*% weights)) +
                             as.numeric(crossprod(u, Cmat %*% weights))           # exp07 plug-in sample-cov IF
                nuis_corr <- compute_gmm_weight_correction_cov_edid(panel_obj, g, t, pairs, prop_ratios,
                               cond_means, u, weights, m_aux, r_aux, pt_assumption)  # ACH correction for u'Cw
                psi_const <- psi_plug + nuis_corr                                 # captured for the misspec_robust fold
                if (isTRUE(getOption("edid_store_psiomega"))) {                   # research diagnostic accumulator
                  acc <- getOption("edid_psiomega_acc", list())                  # corr = -nuis_corr keeps the uniform
                  acc[[paste0(g, "_", t)]] <- list(data = psi_plug, corr = -nuis_corr)  # psi_Omega = data - corr convention
                  options(edid_psiomega_acc = acc)
                }
              }
            }
          }
          eif_gt   <- compute_eif_cov_edid(panel_obj, gen_out_mat, weights, att_gt, g)
          if (isTRUE(estimation_effect))                                    # ACH first-step correction (w frozen)
            eif_gt <- eif_gt - compute_ach_correction_cov_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, weights, m_aux, r_aux, pt_assumption)
          if (higher_order)                                                  # per-cell Hessian (w frozen)
            ho_cell <- compute_cell_hessian_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, weights, m_aux, r_aux, pt_assumption)
          if (misspec_robust && !is.null(psi_const))                         # fold the weight channel AFTER the ACH
            eif_gt <- eif_gt + psi_const                                     # subtraction (averaged/gmm); NULL => +0
        }
      } else {
        # --- No-covariate path ---
        y_hat    <- compute_generated_outcomes_nocov_edid(g, t, pairs, panel_obj, pt_assumption)
        omega    <- compute_omega_star_nocov_edid(g, t, pairs, panel_obj, pt_assumption)
        cond_num <- tryCatch(check_condition_edid(omega), error = function(e) NA_real_)
        # Honor weight_method: with no covariates there is no X-variation, so efficient,
        # averaged, and gmm all invert the same unconditional Omega* and coincide; only
        # uniform differs. (Previously this path silently ignored weight_method.)
        H_local  <- length(y_hat)
        weights  <- if (weight_method == "uniform") rep(1 / H_local, H_local) else compute_efficient_weights_edid(omega)
        att_gt   <- sum(weights * y_hat)
        eif_gt   <- compute_eif_nocov_edid(g, t, pairs, weights, panel_obj, att_gt, pt_assumption)
      }

      # Step 7: SE and inference
      inf_res <- safe_inference_edid(eif_gt, panel_obj$cluster_indices, alpha, att_gt)

      # Step 8: store. The per-cell EIF is intentionally NOT kept on the cell list: it is stored once in eif_list
      # -> eif_matrix -> fit$eif (the only place ever read, by as_MP_edid/aggregation), so a per-cell copy would
      # just hold the n x n_cells influence functions a second time. eif_gt still flows into eif_list below.
      cells[[cell_id]] <- list(
        group           = g,
        time            = t,
        att             = att_gt,
        se              = inf_res$se,
        ci_lower        = inf_res$ci_lower,
        ci_upper        = inf_res$ci_upper,
        t_stat          = inf_res$t_stat,
        p_value         = inf_res$p_value,
        n_pairs         = nrow(pairs),
        weights         = weights,
        condition_num   = cond_num,
        is_pre          = is_pre,
        inference_valid = inf_res$inference_valid,
        ho              = ho_cell   # higher-order pieces (NULL unless higher_order on the covariate path)
      )

      ci_group[cell_id]   <- g
      ci_time[cell_id]    <- t
      ci_cell_id[cell_id] <- cell_id
      ci_is_pre[cell_id]  <- is_pre

      if (keep_eif) eif_list[[cell_id]] <- eif_gt
    }
  }

  if (n_extreme_ratio_instances > 0L) {
    warning(sprintf(
      "Extreme propensity ratios detected (max > 100) in %d estimation step(s). Results may be unstable.",
      n_extreme_ratio_instances
    ))
  }

  # Build EIF matrix if needed
  eif_matrix <- NULL
  if (keep_eif) {
    # Stack EIF vectors as columns: n x n_cells matrix
    eif_matrix <- do.call(cbind, eif_list)
    if (is.null(dim(eif_matrix))) {
      eif_matrix <- matrix(eif_matrix, nrow = n)
    }
  }

  cell_index <- data.frame(
    group   = ci_group,
    time    = ci_time,
    cell_id = ci_cell_id,
    is_pre  = ci_is_pre,
    stringsAsFactors = FALSE
  )

  # Diagnostic: if every post-treatment cell is NA (no admissible pairs), warn instead
  # of silently returning NA. A common cause under pt_assumption='post' is a
  # non-consecutive/irregular time grid where the literal pre-period (g-1-anticipation)
  # is unobserved, so every (g,t) cell yields zero pairs.
  post_idx <- which(!ci_is_pre)
  if (length(post_idx)) {
    post_atts <- vapply(post_idx, function(k) {
      a <- cells[[ci_cell_id[k]]]$att; if (is.null(a)) NA_real_ else a
    }, numeric(1L))
    if (all(!is.finite(post_atts))) {
      warning(sprintf(paste0(
        "All post-treatment ATT(g,t) cells are NA (no admissible pairs). Under ",
        "pt_assumption='%s' the comparison base period is the period immediately before ",
        "treatment; on a non-consecutive/irregular time grid that period may be ",
        "unobserved. Check time spacing and anticipation."), pt_assumption), call. = FALSE)
    }
  }

  list(
    cells      = cells,
    eif_matrix = eif_matrix,
    cell_index = cell_index,
    misspec_robust = misspec_robust   # the EFFECTIVE flag (downgraded to FALSE by the guards above if applicable)
  )
}
