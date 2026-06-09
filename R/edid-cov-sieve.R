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
  if (!is.null(psi_qw))
    stop("compute_omega_star_sieve_edid: the psi/misspec channel is not implemented in the sieve scenario.")
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
  eig$values <- pmax(eig$values, floor_v)
  eig$vectors %*% diag(eig$values, nrow = H) %*% t(eig$vectors)
}
