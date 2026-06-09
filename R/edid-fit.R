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
  estimation_effect = FALSE, higher_order = FALSE, misspec_robust = FALSE,
  trim_level = Inf
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
  n_shrink_approx <- 0L; max_shrink_lambda <- 0; shrink_eg_g <- NA; shrink_eg_t <- NA  # shrinkage-IF approx; emit once

  # Hoist the CELL-INVARIANT Nadaraya-Watson kernel: the n x n weight matrix K and the bandwidths depend only on
  # the full covariate matrix, so they are identical for every (g,t) cell. Build ONCE and reuse, instead of
  # rebuilding (an O(d*n^2) outer/dnorm + an n x n allocation) inside compute_omega_star_cov_edid on every cell --
  # the dominant cost of the covariate path. Numerically identical; only the covariate path needs it.
  kern_bw <- NULL; kern_K <- NULL
  # Conditional-covariance smoother for Omega*(X), two scenarios:
  #   "kernel" (default): O(n^2) NW kernel, mean-cached fast build (compute_omega_star_kernel_fast_edid);
  #                       "kernel_orig" forces the original per-(j,k) build (exact reference).
  #   "sieve":            O(n*p) series build, NO n x n matrix (scales past the kernel's memory wall).
  .omega_method <- getOption("edid_omega_method", "kernel")
  .omega_fun <- switch(.omega_method,
    sieve       = compute_omega_star_sieve_edid,
    kernel_orig = compute_omega_star_cov_edid,
    compute_omega_star_kernel_fast_edid)                  # "kernel" (default) -> mean-cached fast build
  if (use_cov_path && !identical(.omega_method, "sieve")) {  # any kernel variant needs the n x n weight matrix
    kk_full <- build_kernel_weights_edid(panel_obj$covariate_matrix)
    kern_bw <- kk_full$bw; kern_K <- kk_full$K
    # m_eff (Kish effective local sample size) is a function of K_mat ALONE => cell-invariant. Compute it ONCE
    # here and carry it on the matrix so the per-cell shrinkage step reuses it instead of re-summing the n x n
    # K_mat AND re-allocating the n x n K_mat^2 temporary on every (g,t). Read back via attr(K_mat,"edid_m_eff").
    .ks <- rowSums(kern_K); .ksq <- rowSums(kern_K^2)
    attr(kern_K, "edid_m_eff") <- stats::median(.ks^2 / pmax(.ksq, .Machine$double.eps))
  }

  # Per-cohort nuisance cache. The valid pairs, propensity ratios r_{g,.}, inverse propensities, and the
  # overlap-trim mask depend ONLY on the target cohort g, NOT the period t -- yet the per-(g,t) loop re-estimated
  # them for every t (a T-fold redundancy). Estimate them ONCE per g here (plug-in nuisances are deterministic =>
  # bit-identical), then each cell worker reads them and only fits the t-dependent conditional means. Built before
  # the parallel dispatch so the cache is shared copy-on-write across forks.
  # The per-cohort nuisance estimation (esp. the return_aux M-estimator pieces that the misspec_robust default
  # needs) is the dominant cost of the default path, so it is itself run in PARALLEL across cohorts (mclapply),
  # not just hoisted -- otherwise it would be a serial Amdahl bottleneck before the parallel cell loop.
  .gbuild <- function(g) {
    pairs_g <- enumerate_valid_pairs_edid(target_g = g, treatment_groups = tgroups, time_periods = tperiods,
                 period_1 = period_1, pt_assumption = pt_assumption, anticipation = panel_obj$anticipation)
    gb <- list(pairs = pairs_g, prop_ratios = NULL, r_aux = NULL, inv_propensities = NULL,
               trim_keep = NULL, pairs_for_nuisance = NULL, n_extreme = 0L)
    if (nrow(pairs_g) > 0L && use_cov_path) {
      pfn <- pairs_g
      self_cmp <- is.finite(pfn$gp) & (pfn$gp == g); if (any(self_cmp)) pfn$gp[self_cmp] <- Inf
      cross_pairs <- pairs_g[is.finite(pairs_g$gp) & pairs_g$gp != g, , drop = FALSE]
      if (nrow(cross_pairs) > 0L)
        pfn <- unique(rbind(pfn, data.frame(gp = Inf, tpre = unique(cross_pairs$tpre))))
      gb$pairs_for_nuisance <- pfn
      ne <- 0L
      pr <- withCallingHandlers(
        estimate_all_propensity_ratios(panel_obj = panel_obj, g = g, pairs = pfn, bs_df = 4L,
                                       K_folds = K_use, fold_id = fold_id, return_aux = want_aux),
        warning = function(w) {
          if (grepl("Extreme propensity ratios", conditionMessage(w), fixed = TRUE)) { ne <<- ne + 1L; invokeRestart("muffleWarning") }
        })
      if (want_aux) { gb$r_aux <- pr$aux; gb$prop_ratios <- pr$predictions } else gb$prop_ratios <- pr
      gb$n_extreme <- ne
      gb$inv_propensities <- estimate_all_inverse_propensities(panel_obj = panel_obj, g = g, pairs = pairs_g,
        bs_df = 4L, K_folds = K_use, fold_id = fold_id,
        return_aux = (isTRUE(getOption("edid_store_psiomega")) || misspec_robust) && weight_method %in% c("averaged", "efficient"))
      if (is.finite(trim_level) && (!is.null(gb$prop_ratios) || !is.null(gb$inv_propensities))) {
        ks <- union(names(gb$prop_ratios), names(gb$inv_propensities))
        gb$trim_keep <- stats::setNames(lapply(ks, function(k) {
          keep <- rep(TRUE, panel_obj$n)
          rr <- gb$prop_ratios[[k]];      if (!is.null(rr)) keep <- keep & (abs(rr) < trim_level)
          ip <- gb$inv_propensities[[k]]; if (!is.null(ip)) keep <- keep & (abs(ip) < trim_level)
          keep
        }), ks)
      }
    }
    gb
  }
  .mc_pre <- max(1L, suppressWarnings(as.integer(getOption("edid_mc_cores", 1L)))); if (is.na(.mc_pre)) .mc_pre <- 1L
  .glist <- if (.mc_pre > 1L && .Platform$OS.type != "windows")
              parallel::mclapply(tgroups, .gbuild, mc.cores = .mc_pre, mc.preschedule = FALSE)
            else lapply(tgroups, .gbuild)
  for (.gg in .glist) if (inherits(.gg, "try-error")) stop("edid: per-cohort nuisance precompute failed: ", as.character(.gg))
  gcache <- stats::setNames(.glist, as.character(tgroups))
  n_extreme_ratio_instances <- n_extreme_ratio_instances +
    sum(vapply(gcache, function(b) as.integer(b$n_extreme), integer(1)))  # counted once per cohort now

  # Global conditional-mean cache. m_{gp,period}(X) = E[Y_period - Y_1 | G=gp, X] depends ONLY on (gp, period),
  # not on the target cell -- but was re-fit per (g,t) cell (many cells share the same m). Fit each distinct
  # (gp, period) ONCE here (deterministic plug-in => bit-identical), in parallel over comparison cohorts; each
  # worker then subsets the cache to the m's its moment actually uses (preserving the original key order, so the
  # result is identical). Built before the cell dispatch => shared copy-on-write across forks; no duplicated work.
  mcache_pred <- NULL; mcache_aux <- NULL
  if (use_cov_path) {
    # EXACT set of (gp, period) any cell uses = U_g [ {(gp, t) : gp in pfn_g, t in iter_periods} U {(gp, tpre)} ].
    combo_set <- unique(do.call(rbind, lapply(tgroups, function(g) {
      pfn <- gcache[[as.character(g)]]$pairs_for_nuisance; if (is.null(pfn)) return(NULL)
      rbind(expand.grid(gp = unique(pfn$gp), period = iter_periods, KEEP.OUT.ATTRS = FALSE),
            data.frame(gp = pfn$gp, period = pfn$tpre))
    })))
    # One m per combo, parallelized over COMBOS (fine granularity, matches the cell loop) -- preschedule chunks them.
    .mfit1 <- function(i) estimate_all_conditional_means(
      panel_obj = panel_obj, pairs = data.frame(gp = combo_set$gp[i], tpre = combo_set$period[i]),
      t_val = combo_set$period[i], bs_df = 4L, K_folds = K_use, fold_id = fold_id, return_aux = want_aux)
    idx <- seq_len(nrow(combo_set))
    .mlist <- if (.mc_pre > 1L && .Platform$OS.type != "windows")
                parallel::mclapply(idx, .mfit1, mc.cores = .mc_pre, mc.preschedule = TRUE)
              else lapply(idx, .mfit1)
    for (.mm in .mlist) if (inherits(.mm, "try-error")) stop("edid: conditional-mean precompute failed: ", as.character(.mm))
    if (want_aux) {
      mcache_pred <- do.call(c, lapply(.mlist, `[[`, "predictions"))
      mcache_aux  <- do.call(c, lapply(.mlist, `[[`, "aux"))
    } else {
      mcache_pred <- do.call(c, .mlist)
    }
  }

  # Cell specs in the original (g outer, t inner) order so cell_id matches the serial loop exactly.
  .specs <- vector("list", n_cells); .k <- 0L
  for (g in tgroups) for (t in iter_periods) { .k <- .k + 1L; .specs[[.k]] <- list(g = g, t = t, cell_id = .k) }

  # Per-cell worker: the former double-loop body hoisted into a closure so cells can run in parallel
  # (parallel::mclapply) or serially (lapply). Cells are independent; the n x n kernel (kern_K) is shared
  # copy-on-write across forks. The shared warning accumulators become per-cell return values, reduced below.
  .fit_one_cell <- function(.sp) {
    g <- .sp$g; t <- .sp$t; cell_id <- .sp$cell_id
    is_pre  <- (t < g)
    n_extreme <- 0L; n_shrink <- 0L; max_lambda <- 0; shrink_g <- NA; shrink_t <- NA

      # Step 1: pairs from the per-cohort cache (g-only; computed once above)
      gb    <- gcache[[as.character(g)]]
      pairs <- gb$pairs

      # NA cell if no valid pairs
      if (nrow(pairs) == 0L) {
        return(list(cell = list(
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
        ),
        eif = rep(NA_real_, n),
        ci = c(g = g, time = t, cell_id = cell_id, is_pre = is_pre),
        n_extreme = 0L, n_shrink = 0L, max_lambda = 0, shrink_g = NA, shrink_t = NA))
      }

      # g-only nuisances from the per-cohort cache: propensity ratios, inverse propensities, trim mask, r-aux.
      prop_ratios      <- gb$prop_ratios
      inv_propensities <- gb$inv_propensities
      r_aux            <- gb$r_aux
      trim_keep        <- gb$trim_keep
      cond_means       <- NULL; m_aux <- NULL
      ho_cell          <- NULL   # higher-order per-cell pieces (nuisance blocks + Hessian), set on cov path

      if (use_cov_path) {
        # Conditional means: subset the global cache to exactly this cell's (gp, period) combos, in the SAME order
        # estimate_all_conditional_means would have produced them -- bit-identical to the per-cell estimate, but
        # computed once globally instead of once per cell.
        pfn     <- gb$pairs_for_nuisance
        .combos <- unique(rbind(data.frame(gp = pfn$gp, period = t),
                                data.frame(gp = pfn$gp, period = pfn$tpre)))
        .mk        <- paste0(.combos$gp, "_", .combos$period)
        cond_means <- mcache_pred[.mk]
        if (want_aux) m_aux <- mcache_aux[.mk]
      }

      # Steps 2-6: dispatch on covariate vs. no-covariate path
      if (use_cov_path) {
        # --- Covariate path ---
        go_res <- compute_generated_outcomes_cov_edid(
          panel_obj        = panel_obj,
          g                = g,
          t                = t,
          pairs            = pairs,
          prop_ratios      = prop_ratios,
          cond_means       = cond_means,
          pt_assumption    = pt_assumption,
          trim_keep        = trim_keep,
          return_trim_info = TRUE
        )
        gen_out_mat <- go_res$gen_out
        # Per-pair kept-treated masks + masses so the EIF centers on pi_g,kept (= m_kept_j), not pi_g, whenever
        # overlap trimming bit (NULL when trim_level = Inf => EIF takes its byte-identical no-trim path).
        eif_keep   <- go_res$keep
        eif_mkept  <- go_res$m_kept
        # One shared per-cell cache for the cell-invariant per-group kernel slices (K_mat[,idx] + row sums):
        # the array build and the psi pass slice the SAME groups from the SAME K_mat, so memoizing once here
        # halves get_kp's column-subset cost. Cleared each cell (reassigned next iteration) => memory-neutral.
        kp_cache <- new.env(parent = emptyenv())
        if (weight_method == "efficient") {
          # Paper's pointwise efficient weights w(X_i)=Omega*(X_i)^{-1}1/(1'Omega*(X_i)^{-1}1),
          # with a dimension-aware eigenvalue-floor regularization (a=0.7*(5-d)/10) that is
          # asymptotically negligible yet dominates the NW estimation noise for stability.
          omega_arr <- .omega_fun(panel_obj, g, t, pairs,
                                                   prop_ratios, cond_means,
                                                   inv_propensities, bw = kern_bw, K_mat = kern_K,
                                                   return_pointwise = TRUE, kp_cache = kp_cache)
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
          eif_gt   <- compute_eif_cov_edid(panel_obj, gen_out_mat, W_pw, att_gt, g, eif_keep, eif_mkept)
          if (isTRUE(estimation_effect))                                    # ACH first-step correction (W frozen)
            eif_gt <- eif_gt - compute_ach_correction_cov_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, W_pw, m_aux, r_aux, pt_assumption,
              trim_keep = trim_keep)
          if (higher_order)                                                  # per-cell Hessian (W = W_pw frozen)
            ho_cell <- compute_cell_hessian_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, W_pw, m_aux, r_aux, pt_assumption,
              trim_keep = trim_keep, keep_mat = eif_keep, m_kept = eif_mkept)
          weights  <- colMeans(W_pw, na.rm = TRUE)                            # store mean weight per pair
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
            if (is.finite(lam_cell) && lam_cell > 0.05) {                     # large lambda => omission; collect, warn once at end
              n_shrink <- n_shrink + 1L
              if (lam_cell > max_lambda) { max_lambda <- lam_cell; shrink_g <- g; shrink_t <- t }
            }
            po    <- compute_omega_star_cov_edid(panel_obj, g, t, pairs, prop_ratios, cond_means,
                       inv_propensities, bw = kern_bw, K_mat = kern_K, return_pointwise = TRUE,
                       psi_qw = list(pointwise = TRUE, Q = Q_pw, W = W_pw), kp_cache = kp_cache)
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
          omega    <- .omega_fun(panel_obj, g, t, pairs,
                                                  prop_ratios, cond_means,
                                                  inv_propensities, bw = kern_bw, K_mat = kern_K,
                                                  kp_cache = kp_cache)
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
                         inv_propensities, bw = kern_bw, K_mat = kern_K, psi_qw = list(q = q_vec, w = weights),
                         kp_cache = kp_cache)
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
                               cond_means, u, weights, m_aux, r_aux, pt_assumption,
                               trim_keep = trim_keep)                              # ACH correction for u'Cw
                psi_const <- psi_plug + nuis_corr                                 # captured for the misspec_robust fold
                if (isTRUE(getOption("edid_store_psiomega"))) {                   # research diagnostic accumulator
                  acc <- getOption("edid_psiomega_acc", list())                  # corr = -nuis_corr keeps the uniform
                  acc[[paste0(g, "_", t)]] <- list(data = psi_plug, corr = -nuis_corr)  # psi_Omega = data - corr convention
                  options(edid_psiomega_acc = acc)
                }
              }
            }
          }
          eif_gt   <- compute_eif_cov_edid(panel_obj, gen_out_mat, weights, att_gt, g, eif_keep, eif_mkept)
          if (isTRUE(estimation_effect))                                    # ACH first-step correction (w frozen)
            eif_gt <- eif_gt - compute_ach_correction_cov_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, weights, m_aux, r_aux, pt_assumption,
              trim_keep = trim_keep)
          if (higher_order)                                                  # per-cell Hessian (w frozen)
            ho_cell <- compute_cell_hessian_edid(
              panel_obj, g, t, pairs, prop_ratios, cond_means, weights, m_aux, r_aux, pt_assumption,
              trim_keep = trim_keep, keep_mat = eif_keep, m_kept = eif_mkept)
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
      the_cell <- list(
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

      list(cell = the_cell,
           eif = if (keep_eif) eif_gt else NULL,
           ci = c(g = g, time = t, cell_id = cell_id, is_pre = is_pre),
           n_extreme = n_extreme, n_shrink = n_shrink, max_lambda = max_lambda,
           shrink_g = shrink_g, shrink_t = shrink_t)
  }  # end .fit_one_cell

  # Dispatch: parallel over independent cells when edid_mc_cores > 1 (fork; not on Windows), else serial lapply
  # (byte-identical to the original loop). mc.preschedule = FALSE load-balances the very uneven per-cell costs.
  .mc <- suppressWarnings(as.integer(getOption("edid_mc_cores", 1L)))
  if (length(.mc) != 1L || is.na(.mc) || .mc < 1L) .mc <- 1L
  .results <- if (.mc > 1L && .Platform$OS.type != "windows")
                parallel::mclapply(.specs, .fit_one_cell, mc.cores = .mc, mc.preschedule = FALSE)
              else lapply(.specs, .fit_one_cell)

  # Reduce per-cell results back into the pre-allocated structures (identical layout to the serial loop).
  for (.r in .results) {
    if (inherits(.r, "try-error")) stop("edid: a parallel cell worker failed: ", as.character(.r))
    if (is.null(.r)) next
    .cid <- as.integer(.r$ci[["cell_id"]])
    cells[[.cid]]    <- .r$cell
    ci_group[.cid]   <- .r$ci[["g"]]
    ci_time[.cid]    <- .r$ci[["time"]]
    ci_cell_id[.cid] <- .cid
    ci_is_pre[.cid]  <- as.logical(.r$ci[["is_pre"]])
    if (keep_eif) eif_list[[.cid]] <- if (is.null(.r$eif)) rep(NA_real_, n) else .r$eif
    n_extreme_ratio_instances <- n_extreme_ratio_instances + .r$n_extreme
    if (.r$n_shrink > 0L) {
      n_shrink_approx <- n_shrink_approx + .r$n_shrink
      if (.r$max_lambda > max_shrink_lambda) {
        max_shrink_lambda <- .r$max_lambda; shrink_eg_g <- .r$shrink_g; shrink_eg_t <- .r$shrink_t
      }
    }
  }

  if (n_extreme_ratio_instances > 0L) {
    warning(sprintf(
      "Extreme propensity ratios detected (max > 100) in %d estimation step(s). Results may be unstable.",
      n_extreme_ratio_instances
    ))
  }

  if (n_shrink_approx > 0L) {
    warning(sprintf(
      paste0("misspec_robust: the efficient weight-channel SE is an approximation in %d cell(s) where the ",
             "pointwise shrinkage lambda > 0.05 (max %.3f, e.g. cell (%s,%s)); it omits the O(lambda) ",
             "shrinkage influence function there. Use misspec_robust = FALSE for the plug-in SE if exact ",
             "weight-channel SEs are needed."),
      n_shrink_approx, max_shrink_lambda, shrink_eg_g, shrink_eg_t), call. = FALSE)
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
