# edid-cov-sieve.R  (optimized SIEVE Omega build)
# Series-sieve estimator of the conditional covariance Omega*(X): a drop-in alternative to the Nadaraya-Watson
# kernel in compute_omega_star_cov_edid(). Same Eq.(3.12) 5-term structure + inv_p prefactors + eigenfloor /
# pointwise-shrinkage post-processing, so the downstream weight code is unchanged. The ONLY change is the
# conditional-covariance smoother:
#   kernel:  O(n^2 * H^2), materializes an n x n kernel matrix  (memory wall ~ n=8-16k)
#   sieve:   O(n * p * H^2), NO n x n matrix                    (scales to n >> 1e4)
# where p = additive B-spline basis dim (~ d*bs_df). As in the fast kernel build, each distinct difference
# vector's conditional mean is fit ONCE per (vector, group) and reused; only the per-pair cross-moment stays in
# the (j,k) loop. Estimates the SAME object as the kernel; values are close (not byte-identical) and both are
# consistent as n grows. Not wired for the psi/misspec channel. Selected via options(edid_omega_method="sieve").

#' Per-group sieve basis pieces (fit on the group, predicted for ALL units), memoized per group.
#' @keywords internal
#' @noRd
.sieve_group_pieces <- function(X_mat, grp_mask, bs_df) {
  idx <- which(grp_mask)
  if (length(idx) < 2L) return(list(ok = FALSE))
  Bobj   <- build_basis_matrix_edid(X_mat[idx, , drop = FALSE], bs_df)
  B_grp  <- unclass(Bobj); attr(B_grp, "bs_objects") <- NULL
  B_all  <- predict_basis_edid(attr(Bobj, "bs_objects"), X_mat)   # n x p (extrapolates outside group range)
  BtB    <- crossprod(B_grp)                                      # p x p
  BtB_inv <- tryCatch(chol2inv(chol(BtB)), error = function(e) compute_pseudoinverse_edid(BtB))
  list(ok = TRUE, idx = idx, B_grp = B_grp, B_all = B_all, BtB_inv = BtB_inv)
}

#' Series estimator of Omega*(X). Drop-in signature-compatible with compute_omega_star_cov_edid().
#' @keywords internal
#' @noRd
compute_omega_star_sieve_edid <- function(panel_obj, g, t, pairs,
                                          prop_ratios, cond_means,
                                          inv_propensities = NULL,
                                          bw = NULL, K_mat = NULL,
                                          return_pointwise = FALSE,
                                          psi_qw = NULL, bs_df = 4L,
                                          kp_cache = NULL) {
  # kp_cache is accepted for signature-compatibility with the kernel builders (fit_edid_cells passes it
  # uniformly via .omega_fun); the series scenario builds no n x n kernel matrix, so it is intentionally unused.
  # psi_qw triggers the weight-estimation channel (Sigma_Omega): handled below via the sieve OLS-projection IF.
  X_mat <- panel_obj$covariate_matrix
  n <- nrow(X_mat); H <- nrow(pairs); ow <- panel_obj$outcome_wide
  col_t <- panel_obj$period_to_col[[as.character(t)]]
  col_1 <- panel_obj$period_to_col[[as.character(panel_obj$period_1)]]
  mask_g <- panel_obj$cohort_masks[[as.character(g)]]; mask_inf <- panel_obj$never_treated_mask
  pi_g <- panel_obj$cohort_fractions[[as.character(g)]]; pi_inf <- sum(mask_inf) / n
  if (!is.null(inv_propensities)) {
    inv_pg_vec <- inv_propensities[[as.character(g)]]; if (is.null(inv_pg_vec)) inv_pg_vec <- rep(1/pi_g, n)
    inv_pinf_vec <- inv_propensities[["Inf"]];          if (is.null(inv_pinf_vec)) inv_pinf_vec <- rep(1/pi_inf, n)
  } else { inv_pg_vec <- rep(1/pi_g, n); inv_pinf_vec <- rep(1/pi_inf, n) }

  gp_cache <- new.env(parent = emptyenv())
  get_pieces <- function(mask, key) {
    if (exists(key, envir = gp_cache, inherits = FALSE)) return(get(key, envir = gp_cache))
    p <- .sieve_group_pieces(X_mat, mask, bs_df); assign(key, p, envir = gp_cache); p
  }
  # cached conditional mean of a difference vector over a group (one OLS fit, predicted for all n)
  mu_cache <- new.env(parent = emptyenv())
  cmean <- function(v, vkey, pc, gkey) {
    if (!isTRUE(pc$ok)) return(list(ok = FALSE))
    key <- paste0(vkey, "@", gkey)
    if (exists(key, envir = mu_cache, inherits = FALSE)) return(get(key, envir = mu_cache))
    vc <- v - mean(v[pc$idx])                                            # group-centered (shift-invariant)
    mu <- drop(pc$B_all %*% (pc$BtB_inv %*% crossprod(pc$B_grp, vc[pc$idx])))
    out <- list(ok = TRUE, mu = mu, vc = vc, pc = pc); assign(key, out, envir = mu_cache); out
  }
  ccov <- function(a, b) {
    if (!isTRUE(a$ok) || !isTRUE(b$ok)) return(rep(0, n))
    prod_idx <- (a$vc * b$vc)[a$pc$idx]
    muAB <- drop(a$pc$B_all %*% (a$pc$BtB_inv %*% crossprod(a$pc$B_grp, prod_idx)))
    cv <- muAB - a$mu * b$mu; cv[!is.finite(cv)] <- 0; cv
  }

  pc_g <- get_pieces(mask_g, as.character(g)); pc_inf <- get_pieces(mask_inf, "Inf")
  w_vec <- ow[, col_t] - ow[, col_1]
  cm_w_g <- cmean(w_vec, paste0("w", col_t), pc_g, as.character(g))
  term1_const <- inv_pg_vec * ccov(cm_w_g, cm_w_g)

  gp <- pairs$gp; tpre <- pairs$tpre; is_self <- is.finite(gp) & gp == g

  # ---------------------------------------------------------------------------
  # Weight-estimation channel (Sigma_Omega) for the SIEVE smoother
  # ---------------------------------------------------------------------------
  # Mirrors compute_omega_star_cov_edid()'s do_psi path term-for-term, replacing the kernel-local IF of
  # Cov(A,B|X) with the OLS-PROJECTION IF. Each Eq.(3.12) covariance term is Cov(A,B|X_i) = E[AB|X_i] -
  # E[A|X_i]E[B|X_i], all three group-OLS fits mu(X) = B_all beta, beta = BtB_inv B_grp' Z[idx]. The two-step
  # IF of att w.r.t. these coefficients is D' IF_beta with D = (1/n) sum_i s_i dCov_i/dbeta (s_i = pref_i*coup_i,
  # pref/coup/sign copied from the kernel term_psi) and IF_beta(l) = n BtB_inv B_l e_l (l in group); the 1/n and
  # n cancel, so per moment the contribution is e_l * B_l' BtB_inv (sum_i s_i dCov_i/dbeta) -- a length-p
  # accumulator summed over cells i FIRST (O(n p), no n x n), then one basis-dot per group unit. Product rule on
  # E[A|X]E[B|X] gives the -mu_B, -mu_A weights. Term 1 (cell-constant) is skipped: q_i'1 = 0 makes its coupling
  # sum to 0 per unit (the cancellation the kernel relies on). Centering of A,B is OLS-span-invariant => uses raw
  # fits. Returns list(psi, coupled_C); the inv-p correction and the eif fold (psi - corr) are smoother-agnostic.
  if (!is.null(psi_qw)) {
    pw_psi <- isTRUE(psi_qw$pointwise)
    if (pw_psi) { Q_mat <- psi_qw$Q; W_mat <- psi_qw$W } else { q_vec <- psi_qw$q; w_av <- psi_qw$w }
    # Eigen-floor-aware coupling gradient: preferred over the smooth Q/W or q/w adjoint when present. A 3D array
    # (n x H x H) is the per-unit (efficient) gradient; a 2D matrix (H x H) is the pooled (averaged) gradient,
    # broadcast as a constant coupling across cells. Built from the Daleckii-Krein derivative of the FLOORED
    # inverse so the floored directions are clamped (the smooth adjoint over-states psi_Omega where the floor binds).
    C_arr  <- psi_qw$C
    C_pooled <- !is.null(C_arr) && length(dim(C_arr)) == 2L
    # Leading-order shrinkage correction. The per-unit Omega is regularized to Omega^shrunk = (1-lam)Omega_i +
    # lam*Omega_bar before inversion, and the adjoint q (from W_mat/Q_mat) is computed on Omega^shrunk. The
    # covariance term enters Omega_i with coefficient (1-lam), so dtheta = -(1-lam) q'(dOmega_i)w. The sieve's
    # per-unit Omega is poorly conditioned => lam is large => omitting this factor over-states the channel
    # several-fold (the shrinkage is a big variance reduction, NOT asymptotically negligible here). Scale the
    # per-cell sensitivity by (1-lam). (The kernel path keeps lam~0, so this is a no-op there.) The data-driven
    # dlam and dOmega_bar terms are higher-order and omitted, as for the kernel.
    .lam_shr <- suppressWarnings(as.numeric(psi_qw$lambda))
    if (length(.lam_shr) != 1L || !is.finite(.lam_shr)) .lam_shr <- 0
    .shr <- min(1, max(0, 1 - .lam_shr))
    Wg     <- w_vec                                            # Y_t - Y_1 (the self-term left vector)
    psi_om <- numeric(n)
    cpl    <- new.env(parent = emptyenv())                     # coupled_C per inv_p group (smoother-agnostic)
    rf     <- new.env(parent = emptyenv())                     # raw (uncentered) group-OLS fits, cached per (vkey,group)
    rawfit <- function(v, vkey, gkey, pc) {
      key <- if (is.null(vkey)) NULL else paste0(vkey, "@@", gkey)
      if (!is.null(key) && exists(key, envir = rf, inherits = FALSE)) return(get(key, envir = rf))
      beta <- pc$BtB_inv %*% crossprod(pc$B_grp, v[pc$idx])
      mu   <- drop(pc$B_all %*% beta)
      out  <- list(mu = mu, e = v[pc$idx] - mu[pc$idx])        # mu length n; residual e length n_grp (on idx)
      if (!is.null(key)) assign(key, out, envir = rf)
      out
    }
    add_term <- function(A, B, pc, pref_vec, coup, grp_sign, gkey, akey, bkey) {
      if (!isTRUE(pc$ok)) return(invisible(NULL))
      fa <- rawfit(A, akey, gkey, pc); fb <- rawfit(B, bkey, gkey, pc); fab <- rawfit(A * B, NULL, gkey, pc)
      cov_vals <- fab$mu - fa$mu * fb$mu                        # conditional Cov(A,B|X), length n
      s   <- pref_vec * coup                                    # per-unit att-sensitivity (coup pointwise len-n or scalar)
      aAB <- drop(pc$BtB_inv %*% crossprod(pc$B_all, s))        # BtB_inv (sum_i s_i B_i); NO factor n (cancels)
      aA  <- drop(pc$BtB_inv %*% crossprod(pc$B_all, s * fb$mu))# weight = mu_B (product rule)
      aB  <- drop(pc$BtB_inv %*% crossprod(pc$B_all, s * fa$mu))# weight = mu_A
      idx <- pc$idx
      psi_om[idx] <<- psi_om[idx] +
        (-drop(pc$B_grp %*% aAB) * fab$e + drop(pc$B_grp %*% aA) * fa$e + drop(pc$B_grp %*% aB) * fb$e)
      cur <- if (exists(gkey, envir = cpl, inherits = FALSE)) get(gkey, envir = cpl) else numeric(n)
      assign(gkey, cur + (grp_sign * coup) * cov_vals, envir = cpl)
      invisible(NULL)
    }
    inv_pgp_psi <- function(j) {
      if (is.infinite(gp[j])) return(inv_pinf_vec)
      gpk <- as.character(gp[j])
      if (!is.null(inv_propensities) && !is.null(inv_propensities[[gpk]])) return(inv_propensities[[gpk]])
      pi_gp <- panel_obj$cohort_fractions[[gpk]]; if (!is.null(pi_gp) && pi_gp > 1e-15) rep(1/pi_gp, n) else rep(0, n)
    }
    for (j in seq_len(H)) {
      ctj <- panel_obj$period_to_col[[as.character(tpre[j])]]
      Uj  <- ow[, col_t] - ow[, ctj]; Vj <- ow[, ctj] - ow[, col_1]
      for (k in j:H) {
        ctk  <- panel_obj$period_to_col[[as.character(tpre[k])]]
        Uk   <- ow[, col_t] - ow[, ctk]; Vk <- ow[, ctk] - ow[, col_1]
        coup <- if (C_pooled)       { if (j == k) -C_arr[j, j]    else -2 * C_arr[j, k] }    # pooled (averaged) gradient: scalar
                else if (!is.null(C_arr)) { if (j == k) -C_arr[, j, j] else -2 * C_arr[, j, k] }  # per-unit (efficient) gradient: len-n
                else if (pw_psi)     { if (j == k) Q_mat[, j] * W_mat[, j] else Q_mat[, j] * W_mat[, k] + Q_mat[, k] * W_mat[, j] }
                else                 { if (j == k) q_vec[j] * w_av[j]      else q_vec[j] * w_av[k]      + q_vec[k] * w_av[j] }
        coup <- coup * .shr                                       # leading-order (1-lam) shrinkage-IF correction
        add_term(Uj, Uk, pc_inf, inv_pinf_vec, coup, 1, "Inf", paste0("u", ctj), paste0("u", ctk))             # T2
        if (is_self[j]) add_term(Wg, Vj, pc_g, -inv_pg_vec, coup, -1, as.character(g), "w", paste0("v", ctj))  # T3
        if (is_self[k]) add_term(Wg, Vk, pc_g, -inv_pg_vec, coup, -1, as.character(g), "w", paste0("v", ctk))  # T4
        if (identical(gp[j], gp[k])) {                                                                          # T5
          gpk <- as.character(gp[j])
          pc5 <- if (is.infinite(gp[j])) pc_inf else { mgp <- panel_obj$cohort_masks[[gpk]]; if (is.null(mgp)) list(ok = FALSE) else get_pieces(mgp, gpk) }
          add_term(Vj, Vk, pc5, inv_pgp_psi(j), coup, 1, gpk, paste0("v", ctj), paste0("v", ctk))
        }
      }
    }
    return(list(psi = psi_om, coupled_C = as.list(cpl)))
  }

  cm_u_inf <- vector("list", H); cm_v_g <- vector("list", H); cm_v_gp <- vector("list", H)
  for (j in seq_len(H)) {
    ctp <- panel_obj$period_to_col[[as.character(tpre[j])]]
    u_j <- ow[, col_t] - ow[, ctp]; v_j <- ow[, ctp] - ow[, col_1]
    cm_u_inf[[j]] <- cmean(u_j, paste0("u", ctp), pc_inf, "Inf")
    if (is_self[j]) cm_v_g[[j]] <- cmean(v_j, paste0("v", ctp), pc_g, as.character(g))
    gpk <- as.character(gp[j])
    pc_gp <- if (is.infinite(gp[j])) pc_inf else {
      mgp <- panel_obj$cohort_masks[[gpk]]; if (is.null(mgp)) list(ok = FALSE) else get_pieces(mgp, gpk) }
    cm_v_gp[[j]] <- cmean(v_j, paste0("v", ctp), pc_gp, gpk)
  }
  inv_pgp_of <- function(j) {
    if (is.infinite(gp[j])) return(inv_pinf_vec)
    gpk <- as.character(gp[j])
    if (!is.null(inv_propensities) && !is.null(inv_propensities[[gpk]])) return(inv_propensities[[gpk]])
    pi_gp <- panel_obj$cohort_fractions[[gpk]]; if (!is.null(pi_gp) && pi_gp > 1e-15) rep(1/pi_gp, n) else rep(0, n)
  }
  term34_g <- vector("list", H)
  for (j in seq_len(H)) term34_g[[j]] <- if (is_self[j]) ccov(cm_w_g, cm_v_g[[j]]) else NULL

  Omega_array <- NULL; Omega_hat <- matrix(0, H, H)
  if (return_pointwise) {
    Omega_array <- array(0, dim = c(n, H, H))
    for (j in seq_len(H)) for (k in j:H) {
      o <- term1_const + inv_pinf_vec * ccov(cm_u_inf[[j]], cm_u_inf[[k]])
      if (is_self[j]) o <- o - inv_pg_vec * term34_g[[j]]
      if (is_self[k]) o <- o - inv_pg_vec * term34_g[[k]]
      if (identical(gp[j], gp[k])) o <- o + inv_pgp_of(j) * ccov(cm_v_gp[[j]], cm_v_gp[[k]])
      ojk <- mean(o); Omega_hat[j, k] <- ojk; if (k != j) Omega_hat[k, j] <- ojk  # needed by the shrinkage step
      Omega_array[, j, k] <- o; if (k != j) Omega_array[, k, j] <- o
    }
  } else {
    # Averaged Omega-bar via symmetric crossprods (series analogue of the kernel batch). For the sieve smoother
    # E[V|X] = B_all (B_grp'B_grp)^{-1} B_grp' V, mean_i prefac_i E[V_jV_k|X_i] = (1/n) Vc' diag(h) Vc with
    # h = B_grp (B_grp'B_grp)^{-1} (B_all' prefac); one crossprod per term/group instead of H(H+1)/2 regressions.
    avg_block <- function(prefac, Vc, M, pc) {
      cB <- colSums(prefac * pc$B_all)
      h  <- drop(pc$B_grp %*% (pc$BtB_inv %*% cB))
      (crossprod(Vc, h * Vc) - crossprod(M, prefac * M)) / n
    }
    Omega_hat <- matrix(mean(term1_const), H, H)
    c3 <- vapply(seq_len(H), function(j) if (is_self[j]) mean(-inv_pg_vec * term34_g[[j]]) else 0, numeric(1))
    Omega_hat <- Omega_hat + outer(c3, rep(1, H)) + outer(rep(1, H), c3)
    if (isTRUE(pc_inf$ok)) {
      Uc <- vapply(cm_u_inf, function(z) z$vc[pc_inf$idx], numeric(length(pc_inf$idx)))
      Mu <- vapply(cm_u_inf, function(z) z$mu, numeric(n))
      Omega_hat <- Omega_hat + avg_block(inv_pinf_vec, Uc, Mu, pc_inf)
    }
    for (gpv in unique(gp)) {
      S <- which(gp == gpv); if (!length(S) || !isTRUE(cm_v_gp[[S[1]]]$ok)) next
      pc5 <- cm_v_gp[[S[1]]]$pc
      Vc <- vapply(S, function(j) cm_v_gp[[j]]$vc[pc5$idx], numeric(length(pc5$idx)))
      Mv <- vapply(S, function(j) cm_v_gp[[j]]$mu, numeric(n))
      Omega_hat[S, S] <- Omega_hat[S, S] + avg_block(inv_pgp_of(S[1]), Vc, Mv, pc5)
    }
  }

  if (return_pointwise) {
    Hh <- dim(Omega_array)[2]
    lam_opt <- suppressWarnings(as.numeric(getOption("edid_shrink_lambda", NA_real_)))
    if (length(lam_opt) == 1L && is.finite(lam_opt)) lam <- min(1, max(0, lam_opt)) else {
      shape_var <- mean(apply(Omega_array, c(2, 3), stats::var)); dg <- diag(Omega_hat)
      p_dim <- ncol(pc_g$B_all); m_eff <- max(2, stats::median(c(sum(mask_g), sum(mask_inf))) / max(p_dim, 1))
      samp_var <- mean(outer(dg, dg) + Omega_hat^2) / max(m_eff, 1)
      lam <- min(1, max(0, samp_var / max(shape_var, .Machine$double.eps)))
    }
    if (lam > 0) for (jj in seq_len(Hh)) for (kk in seq_len(Hh))
      Omega_array[, jj, kk] <- (1 - lam) * Omega_array[, jj, kk] + lam * Omega_hat[jj, kk]
    attr(Omega_array, "shrink_lambda") <- lam
    attr(Omega_array, "omega_bar") <- Omega_hat   # pooled (PSD after flooring): target for per-unit PD-blend
    return(Omega_array)
  }
  eig <- eigen(Omega_hat, symmetric = TRUE)
  d_cov <- ncol(X_mat); a_floor <- 0.7 * (5 - min(as.integer(d_cov), 4L)) / 10
  mx <- max(eig$values); floor_v <- if (is.finite(mx) && mx > 0) mx * n^(-a_floor) else 1e-12
  lam_raw <- eig$values                                       # raw (pre-floor) eigenvalues, for the coupling IF below
  eig$values <- pmax(eig$values, floor_v)
  out <- eig$vectors %*% diag(eig$values, nrow = H) %*% t(eig$vectors)
  # Attach the raw eigendecomposition so the AVERAGED+sieve weight channel can build the eigen-floor-aware
  # coupling (Daleckii-Krein derivative of the FLOORED inverse). The pooled Omega-bar floor binds in high-H
  # cells (floor ~ mx*n^-a, e.g. 18% of mx at n=600); there the smooth -sym(q w') adjoint over-states psi_Omega
  # (the floored directions do not respond to dOmega). Inert for kernel (different builder) and for the estimate
  # itself (attributes are stripped by solve/%*%); read only by compute_obar_coupling_edid.
  attr(out, "eig_floor") <- list(values = lam_raw, vectors = eig$vectors, floor = floor_v)
  out
}
