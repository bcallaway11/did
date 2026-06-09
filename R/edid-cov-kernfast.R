# edid-cov-kernfast.R  (optimized KERNEL Omega build)
# Same Eq.(3.12) 5-term kernel estimator as compute_omega_star_cov_edid(), but the conditional MEANS of each
# distinct outcome-difference vector are computed ONCE per (vector, group) and reused, instead of being
# recomputed inside every (j,k) pair (the original calls kernel_cond_cov_kp per pair, so each mean is rebuilt
# ~H times). Only the per-pair cross-moment E_K[A_j A_k|X] (terms 2 and 5, the O(H^2) part) stays in the loop.
# Arithmetic is the SAME (same centering, same drop(Kg %*% .)/Ks, same accumulation order) => bit-identical to
# the kernel path for both the averaged H x H and the efficient per-unit array. Not wired for the psi/misspec
# channel (that stays on compute_omega_star_cov_edid). Selected via options(edid_omega_method = "kernel").

#' Optimized kernel Omega*(X): drop-in, signature-compatible with compute_omega_star_cov_edid().
#' @keywords internal
compute_omega_star_kernel_fast_edid <- function(panel_obj, g, t, pairs,
                                                prop_ratios, cond_means,
                                                inv_propensities = NULL,
                                                bw = NULL, K_mat = NULL,
                                                return_pointwise = FALSE,
                                                psi_qw = NULL,
                                                kp_cache = NULL) {
  if (!is.null(psi_qw))
    stop("compute_omega_star_kernel_fast_edid: the psi/misspec channel stays on compute_omega_star_cov_edid.")
  X_mat <- panel_obj$covariate_matrix
  n <- nrow(X_mat); H <- nrow(pairs); ow <- panel_obj$outcome_wide
  if (is.null(K_mat)) { kk <- build_kernel_weights_edid(X_mat, bw); bw <- kk$bw; K_mat <- kk$K }
  col_t <- panel_obj$period_to_col[[as.character(t)]]
  col_1 <- panel_obj$period_to_col[[as.character(panel_obj$period_1)]]
  mask_g <- panel_obj$cohort_masks[[as.character(g)]]; mask_inf <- panel_obj$never_treated_mask
  pi_g <- panel_obj$cohort_fractions[[as.character(g)]]; pi_inf <- sum(mask_inf) / n
  if (!is.null(inv_propensities)) {
    inv_pg_vec <- inv_propensities[[as.character(g)]]; if (is.null(inv_pg_vec)) inv_pg_vec <- rep(1/pi_g, n)
    inv_pinf_vec <- inv_propensities[["Inf"]];          if (is.null(inv_pinf_vec)) inv_pinf_vec <- rep(1/pi_inf, n)
  } else { inv_pg_vec <- rep(1/pi_g, n); inv_pinf_vec <- rep(1/pi_inf, n) }

  # per-group kernel pieces (idx, Kg, Ks), memoized exactly as the original get_kp. A shared kp_cache (from
  # fit_edid_cells) lets the psi pass reuse these byte-identical slices instead of re-cutting K_mat[,idx].
  kpiece <- if (is.null(kp_cache)) new.env(parent = emptyenv()) else kp_cache
  get_kp <- function(mask, key) {
    if (exists(key, envir = kpiece, inherits = FALSE)) return(get(key, envir = kpiece))
    idx <- which(mask)
    kp <- if (length(idx) < 2L) list(ok = FALSE) else {
      Kg <- K_mat[, idx, drop = FALSE]; Ks <- rowSums(Kg); Ks[Ks < 1e-15] <- NA_real_
      list(ok = TRUE, idx = idx, Kg = Kg, Ks = Ks)
    }
    assign(key, kp, envir = kpiece); kp
  }
  # cached conditional mean of a difference vector over a group: returns mu (length n) and centered vc (over idx).
  # SAME arithmetic as kernel_cond_cov_kp's mu_A: A_c = A[idx]-mean(A[idx]); mu = drop(Kg %*% A_c)/Ks.
  mu_cache <- new.env(parent = emptyenv())
  cmean <- function(v, vkey, kp, gkey) {
    if (!isTRUE(kp$ok)) return(list(ok = FALSE))
    key <- paste0(vkey, "@", gkey)
    if (exists(key, envir = mu_cache, inherits = FALSE)) return(get(key, envir = mu_cache))
    vc <- v[kp$idx] - mean(v[kp$idx]); mu <- drop(kp$Kg %*% vc) / kp$Ks
    out <- list(ok = TRUE, mu = mu, vc = vc, kp = kp); assign(key, out, envir = mu_cache); out
  }
  # conditional covariance from two cached means (cross-moment computed here; the O(H^2) per-pair part)
  ccov <- function(a, b) {
    if (!isTRUE(a$ok) || !isTRUE(b$ok)) return(rep(0, n))
    cv <- drop(a$kp$Kg %*% (a$vc * b$vc)) / a$kp$Ks - a$mu * b$mu
    cv[is.na(cv)] <- 0; cv
  }

  kp_g <- get_kp(mask_g, as.character(g)); kp_inf <- get_kp(mask_inf, "Inf")
  w_vec <- ow[, col_t] - ow[, col_1]                                   # Y_t - Y_1
  cm_w_g <- cmean(w_vec, paste0("w", col_t), kp_g, as.character(g))    # mu_w in group g
  term1_const <- inv_pg_vec * ccov(cm_w_g, cm_w_g)                     # cell-constant term 1

  # Precompute per-pair info + cached means: u_j (=Y_t-Y_tpre_j) in Inf; v_j (=Y_tpre_j-Y_1) in g (self) and in gp_j.
  gp <- pairs$gp; tpre <- pairs$tpre; is_self <- is.finite(gp) & gp == g
  cm_u_inf <- vector("list", H); cm_v_g <- vector("list", H); cm_v_gp <- vector("list", H)
  col_tp <- integer(H)
  for (j in seq_len(H)) {
    col_tp[j] <- panel_obj$period_to_col[[as.character(tpre[j])]]
    u_j <- ow[, col_t] - ow[, col_tp[j]]; v_j <- ow[, col_tp[j]] - ow[, col_1]
    cm_u_inf[[j]] <- cmean(u_j, paste0("u", col_tp[j]), kp_inf, "Inf")
    if (is_self[j]) cm_v_g[[j]] <- cmean(v_j, paste0("v", col_tp[j]), kp_g, as.character(g))
    # term5 group gp_j (only needed when some k shares gp; compute lazily but cache by (vkey, gpkey))
    gpk <- as.character(gp[j])
    if (is.infinite(gp[j])) { kp_gp <- kp_inf } else {
      mgp <- panel_obj$cohort_masks[[gpk]]; kp_gp <- if (is.null(mgp)) list(ok = FALSE) else get_kp(mgp, gpk)
    }
    cm_v_gp[[j]] <- cmean(v_j, paste0("v", col_tp[j]), kp_gp, gpk)
  }
  # inv_p prefactor per gp group for term 5 (per-unit vectors), matching the original
  inv_pgp_of <- function(j) {
    if (is.infinite(gp[j])) return(inv_pinf_vec)
    gpk <- as.character(gp[j])
    if (!is.null(inv_propensities) && !is.null(inv_propensities[[gpk]])) return(inv_propensities[[gpk]])
    pi_gp <- panel_obj$cohort_fractions[[gpk]]; if (!is.null(pi_gp) && pi_gp > 1e-15) rep(1/pi_gp, n) else rep(0, n)
  }
  # term 3/4 depend only on the single pair index: cov(w, v_j | g), precomputed per self pair
  term34_g <- vector("list", H)
  for (j in seq_len(H)) term34_g[[j]] <- if (is_self[j]) ccov(cm_w_g, cm_v_g[[j]]) else NULL

  Omega_array <- NULL; Omega_hat <- matrix(0, H, H)
  if (return_pointwise) {
    # Per-unit array: irreducible n x H x H. Keep the per-(j,k) accumulation (bit-identical arithmetic; the
    # eigen-inversion downstream is summation-order sensitive, so we do NOT reorder this path).
    Omega_array <- array(0, dim = c(n, H, H))
    for (j in seq_len(H)) for (k in j:H) {
      o <- term1_const + inv_pinf_vec * ccov(cm_u_inf[[j]], cm_u_inf[[k]])    # term1 + term2
      if (is_self[j]) o <- o - inv_pg_vec * term34_g[[j]]                     # term3
      if (is_self[k]) o <- o - inv_pg_vec * term34_g[[k]]                     # term4
      if (identical(gp[j], gp[k])) o <- o + inv_pgp_of(j) * ccov(cm_v_gp[[j]], cm_v_gp[[k]])  # term5
      ojk <- mean(o); Omega_hat[j, k] <- ojk; if (k != j) Omega_hat[k, j] <- ojk  # needed by the shrinkage step
      Omega_array[, j, k] <- o; if (k != j) Omega_array[, k, j] <- o
    }
  } else {
    # Averaged Omega-bar via SYMMETRIC crossprods: each (term, group) block is mean_i prefac_i (E_K[V_jV_k|X_i]
    # - muV_j muV_k), which collapses to one BLAS-3 crossprod over the n_grp x H_block matrix of centered diffs
    # (and one over the n x H_block means), instead of H(H+1)/2 separate matrix-vector smooths. Exploits Omega's
    # symmetry (crossprod returns a symmetric block) and that all pairs in a group share the same kernel weights.
    # With prefac_i = prefactor and cWp[l] = sum_i prefac_i Kg[i,l]/Ks_i:
    #   mean_i prefac_i E_K[V_jV_k|X_i] = (1/n) Vc' diag(cWp) Vc ;  mean_i prefac_i muV_j muV_k = (1/n) M' diag(prefac) M
    avg_block <- function(prefac, Vc, M, kp) {
      wrow <- prefac / kp$Ks; wrow[is.na(wrow)] <- 0
      cWp  <- colSums(wrow * kp$Kg); Ms <- M; Ms[is.na(Ms)] <- 0
      (crossprod(Vc, cWp * Vc) - crossprod(Ms, prefac * Ms)) / n
    }
    Omega_hat <- matrix(mean(term1_const), H, H)                              # term1 (constant block)
    c3 <- vapply(seq_len(H), function(j) if (is_self[j]) mean(-inv_pg_vec * term34_g[[j]]) else 0, numeric(1))
    Omega_hat <- Omega_hat + outer(c3, rep(1, H)) + outer(rep(1, H), c3)      # term3 + term4 (rank-2)
    if (isTRUE(kp_inf$ok)) {                                                  # term2 (group Inf, all H pairs)
      Uc <- vapply(cm_u_inf, function(z) z$vc, numeric(length(kp_inf$idx)))
      Mu <- vapply(cm_u_inf, function(z) z$mu, numeric(n))
      Omega_hat <- Omega_hat + avg_block(inv_pinf_vec, Uc, Mu, kp_inf)
    }
    for (gpv in unique(gp)) {                                                 # term5 (per gp group sub-block)
      S <- which(gp == gpv); if (!length(S) || !isTRUE(cm_v_gp[[S[1]]]$ok)) next
      kp5 <- cm_v_gp[[S[1]]]$kp
      Vc <- vapply(S, function(j) cm_v_gp[[j]]$vc, numeric(length(kp5$idx)))
      Mv <- vapply(S, function(j) cm_v_gp[[j]]$mu, numeric(n))
      Omega_hat[S, S] <- Omega_hat[S, S] + avg_block(inv_pgp_of(S[1]), Vc, Mv, kp5)
    }
  }

  # ---- identical post-processing to compute_omega_star_cov_edid ----
  if (return_pointwise) {
    Hh <- dim(Omega_array)[2]
    lam_opt <- suppressWarnings(as.numeric(getOption("edid_shrink_lambda", NA_real_)))
    if (length(lam_opt) == 1L && is.finite(lam_opt)) lam <- min(1, max(0, lam_opt)) else {
      m_eff <- attr(K_mat, "edid_m_eff")                          # cell-invariant; precomputed once by fit_edid_cells
      if (is.null(m_eff)) { ksum <- rowSums(K_mat); ksq <- rowSums(K_mat^2)  # standalone fallback (same value)
        m_eff <- stats::median(ksum^2 / pmax(ksq, .Machine$double.eps)) }
      shape_var <- mean(apply(Omega_array, c(2, 3), stats::var))
      dg <- diag(Omega_hat); samp_var <- mean(outer(dg, dg) + Omega_hat^2) / max(m_eff, 1)
      lam <- min(1, max(0, samp_var / max(shape_var, .Machine$double.eps)))
    }
    if (lam > 0) for (jj in seq_len(Hh)) for (kk in seq_len(Hh))
      Omega_array[, jj, kk] <- (1 - lam) * Omega_array[, jj, kk] + lam * Omega_hat[jj, kk]
    if (isTRUE(getOption("edid_diag_lambda"))) options(edid_lambda_acc = c(getOption("edid_lambda_acc", numeric(0)), lam))
    attr(Omega_array, "shrink_lambda") <- lam
    attr(Omega_array, "omega_bar") <- Omega_hat   # pooled (PSD after flooring): target for per-unit PD-blend
    return(Omega_array)
  }
  eig <- eigen(Omega_hat, symmetric = TRUE)
  d_cov <- ncol(X_mat); a_floor <- 0.7 * (5 - min(as.integer(d_cov), 4L)) / 10
  mx <- max(eig$values); floor_v <- if (is.finite(mx) && mx > 0) mx * panel_obj$n^(-a_floor) else 1e-12
  eig$values <- pmax(eig$values, floor_v)
  eig$vectors %*% diag(eig$values, nrow = H) %*% t(eig$vectors)
}
