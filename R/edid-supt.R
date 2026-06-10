# edid-supt.R
# Analytic simultaneous (sup-t) uniform confidence bands for edid.
#
# The uniform-band critical value is the equicoordinate (1 - alpha) quantile of max_k |Z_k|, with
# Z ~ N(0, corr(Sigma)) (Montiel Olea & Plagborg-Moller 2019). It is a pure function of the coefficient
# covariance matrix Sigma, which lets edid produce uniform bands WITHOUT a bootstrap and -- the reason
# this path exists -- lets the higher-order ("Wick") variance refinement enter through Sigma even though
# it is a degenerate second-order U-statistic and so cannot be carried by the IF-based multiplier
# bootstrap (`mboot`). The crit is computed by a fast base-R Monte Carlo (eigen square-root of the
# correlation matrix + a vectorized row-max); pure `stats`, no new package dependency.

#' Cluster-robust covariance of the columns of an influence-function matrix
#'
#' \eqn{\Sigma_{1,jk} = n^{-2}\sum_i \mathrm{IF}_{ij}\,\mathrm{IF}_{ik}} (i.i.d.), or the cluster-summed sandwich with the
#' G/(G-1) finite-cluster correction when \code{cluster_indices} is supplied. This is the analytic
#' first-order coefficient covariance; \code{sqrt(diag(.))} reproduces \code{safe_inference_edid()}'s SE.
#'
#' @param M n x K influence-function matrix.
#' @param cluster_indices length-n cluster id vector (1..G), or NULL for i.i.d.
#' @param n number of units (sample size).
#' @return K x K covariance matrix.
#' @keywords internal
cluster_cov_edid <- function(M, cluster_indices, n) {
  M <- as.matrix(M)
  if (is.null(cluster_indices)) return(crossprod(M) / n^2)
  G <- length(unique(cluster_indices))
  if (G <= 1L) return(matrix(NA_real_, ncol(M), ncol(M)))
  CS <- rowsum(M, cluster_indices)
  (G / (G - 1)) * crossprod(CS) / n^2
}

#' Analytic sup-t critical value from a coefficient covariance matrix
#'
#' Returns `c` such that the simultaneous band `theta_hat_k +/- c * se_k` (se_k = sqrt(diag(Sigma))) has
#' joint coverage `1 - alp`, i.e. the `(1 - alp)` quantile of `max_k |Z_k|`, `Z ~ N(0, corr(Sigma))`. Never
#' returns below the pointwise `qnorm(1 - alp/2)`. With < 2 non-degenerate coordinates it returns the
#' pointwise value.
#'
#' @param Sigma K x K coefficient covariance matrix.
#' @param alp significance level (two-sided simultaneous coverage 1 - alp). Default 0.05.
#' @param B number of Monte Carlo draws. Default 1e5.
#' @param seed optional integer for reproducibility (restores the RNG state on exit).
#' @return scalar critical value (>= qnorm(1 - alp/2)).
#' @keywords internal
supt_crit_edid <- function(Sigma, alp = 0.05, B = 1e5L, seed = NULL) {
  pointwise <- stats::qnorm(1 - alp / 2)
  Sigma <- as.matrix(Sigma)
  d  <- sqrt(diag(Sigma)); ok <- is.finite(d) & d > 0; p <- sum(ok)
  # Exclude coordinates with any non-finite covariance row/column entry so a single bad off-diagonal
  # entry does not break the whole familywise critical-value simulation.
  row_ok <- vapply(seq_len(nrow(Sigma)), function(i) all(is.finite(Sigma[i, ]) & is.finite(Sigma[, i])), logical(1L))
  ok <- ok & row_ok
  p <- sum(ok)
  if (p < 2L) return(pointwise)
  R <- Sigma[ok, ok, drop = FALSE] / tcrossprod(d[ok])
  R <- (R + t(R)) / 2                                        # symmetrize away roundoff
  # CANONICAL square root via Cholesky (UNIQUE for PD R), not the eigen-decomposition. The eigenvectors of a
  # near-degenerate R are not uniquely determined, so an eps-level change in R (an equivalent but FP-reordered
  # upstream computation -- e.g. a BLAS vs per-dimension kernel build) rotates them and, for a FIXED rng seed,
  # produces different draws and a crit that wobbles at ~1e-3 even though R itself moved only ~1e-9. The Cholesky
  # factor is a continuous, canonical function of R, so the seeded crit is reproducible and build-invariant. A
  # tiny relative ridge guarantees PD (R is only PSD; any rank-deficient coordinate then draws at sqrt(ridge)
  # scale => a negligible contribution to max_k |Z_k|). Falls back to the PSD eigen root if Cholesky still fails.
  rg <- 1e-10 * max(1, mean(diag(R)))
  U  <- tryCatch(chol(R + diag(rg, nrow(R))),               # upper-triangular: (R + ridge) = U'U
                 error = function(e2) { e <- eigen(R, symmetric = TRUE); sqrt(pmax(e$values, 0)) * t(e$vectors) })
  # ALWAYS restore the caller's RNG state: the B x p rnorm draws below must not perturb the user's stream.
  # The analytic cband is the DEFAULT, so a bare edid() call would otherwise silently advance .Random.seed.
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))             # reproducible when a seed is supplied
  Z <- matrix(stats::rnorm(B * p), B, p) %*% U               # B x p ~ N(0, R)
  m <- abs(Z[, 1L]); for (j in seq_len(p)[-1]) m <- pmax(m, abs(Z[, j]))   # vectorized row-max |Z|
  crit <- as.numeric(stats::quantile(m, 1 - alp, names = FALSE))
  max(crit, pointwise)
}

#' Higher-order ("Wick") covariance Sigma_quad of the cell ATT(g,t) vector
#'
#' Returns the K x K matrix \eqn{\Sigma_{quad}} whose \eqn{(k,j)} entry is the degenerate second-order
#' U-statistic ("Isserlis/Wick") covariance contributed by first-step sieve-nuisance estimation,
#' \deqn{\Sigma_{quad,kj} = \tfrac12\,\mathrm{tr}(H_k V H_j V),}
#' where \eqn{V} is the JOINT stacked-coefficient covariance across all cells' nuisance blocks and
#' \eqn{H_k} is cell \eqn{k}'s Hessian of \eqn{att} in those coefficients, embedded block-sparse in the
#' joint coefficient space (cell \eqn{k}'s \eqn{att} depends only on its own block, so off-diagonal cross-cell
#' entries come for free from \eqn{V}'s off-diagonal blocks -- the covariance of the two cells' scores over
#' their common units). \eqn{V} is the HC2-leverage-corrected, cluster-robust sandwich
#' \eqn{H^{-1}_{blk}\,(\sum_c S_c'S_c)\,H^{-1}_{blk}/n^2} with the \eqn{G/(G-1)} finite-cluster factor, the
#' stacked scores \eqn{S} corrected by \eqn{1/\sqrt{1-h}} (leverage \eqn{h} capped at 0.5). Adding
#' \eqn{\Sigma_{quad}} to the first-order \code{cluster_cov_edid()} covariance gives the higher-order-aware
#' Sigma the sup-t crit and SEs are read from. Mirrors the validated prototype
#' \code{exp10_vroute_supt.R::make_Sigma} exactly. Cells without an estimated Hessian (no covariates /
#' fallback nuisances; \code{ho$H = NULL} or 0 x 0) contribute zero rows and columns.
#'
#' @param cells list of \code{edid_cell_result} objects; each higher-order cell carries
#'   \code{$ho$blocks} (ordered nuisance blocks with \code{B}, \code{score_mat}, \code{H_inv}, \code{p}) and
#'   \code{$ho$H} (its P_k x P_k Hessian). Order must match the cell order of the ATT(g,t) vector.
#' @param cluster_indices length-n cluster id vector (1..G), or NULL for i.i.d.
#' @param n number of units (sample size).
#' @return K x K \eqn{\Sigma_{quad}} matrix (PSD up to roundoff).
#' @keywords internal
sigma_quad_edid <- function(cells, cluster_indices, n) {
  K <- length(cells)
  Sigma_quad <- matrix(0, K, K)

  # Per-cell stacked-coefficient dimension; cells with no estimated blocks get P_k = 0.
  Pk <- vapply(cells, function(cc) {
    h <- cc$ho
    if (is.null(h) || is.null(h$blocks) || length(h$blocks) == 0L) return(0L)
    sum(vapply(h$blocks, function(b) b$p, 1L))
  }, integer(1L))
  P_tot <- sum(Pk)
  if (P_tot == 0L) return(Sigma_quad)                       # no covariate cell -> Sigma_quad = 0

  cstart <- cumsum(c(0L, Pk[-K]))                            # 0-based joint-coef offset of each cell

  # HC2 leverage-corrected stacked scores S_all (n x P_tot). H^{-1} is block-diagonal: keep it as a LIST of
  # (index-range, H_inv) blocks rather than a dense P_tot x P_tot matrix -- the dense form makes the V build below
  # an O(P_tot^3) matmul of a mostly-zero matrix (the measured hotspot at large cell counts).
  S_all   <- matrix(0, n, P_tot)
  blk_idx <- list(); blk_Hinv <- list(); bi <- 0L
  for (k in seq_len(K)) {
    if (Pk[k] == 0L) next
    blocks <- cells[[k]]$ho$blocks
    o2 <- 0L
    for (b in blocks) {
      # leverage h = diag(B (B'B)^{-1} B') = diag(B (H_inv/n) B'); in-sample (score nonzero) units only.
      h  <- rowSums((b$B %*% (b$H_inv / n)) * b$B)
      nz <- rowSums(b$score_mat^2) > 0
      h  <- ifelse(nz, pmin(pmax(h, 0), 0.5), 0)             # cap at 0.5 (HC blow-up guard)
      jj <- cstart[k] + o2 + seq_len(b$p)
      S_all[, jj] <- b$score_mat / sqrt(1 - h)               # HC2 correction
      bi <- bi + 1L; blk_idx[[bi]] <- jj; blk_Hinv[[bi]] <- b$H_inv
      o2 <- o2 + b$p
    }
  }

  # Cluster-robust joint coefficient covariance V (G/(G-1) finite-cluster correction).
  if (is.null(cluster_indices)) {
    Ssum  <- S_all
    cfac  <- 1
  } else {
    Ssum  <- rowsum(S_all, cluster_indices)
    G     <- length(unique(cluster_indices))
    cfac  <- if (G > 1L) G / (G - 1) else 1
  }
  # V = cfac * D C D / n^2 with D = blockdiag(H_inv). Since D is block-diagonal, D %*% C %*% D scales C's block
  # rows then block cols by the per-block H_inv -- no dense P_tot x P_tot Hinv matmul. (Hinv_blk %*% C)[jj, ] =
  # H_inv_jj %*% C[jj, ]; identical arithmetic, O(P_tot^2 * p) instead of O(P_tot^3).
  V <- crossprod(Ssum)                                       # C (P_tot x P_tot score covariance)
  for (bi in seq_along(blk_idx)) { jj <- blk_idx[[bi]]; V[jj, ] <- blk_Hinv[[bi]] %*% V[jj, , drop = FALSE] }
  for (bi in seq_along(blk_idx)) { jj <- blk_idx[[bi]]; V[, jj] <- V[, jj, drop = FALSE] %*% blk_Hinv[[bi]] }
  V <- V * (cfac / (n^2))

  # 0.5 tr(H_k V H_j V), exploiting that H_k is nonzero ONLY in its own cell block idx_k. The original embedded
  # each H_k in a dense P_tot x P_tot Hbig and formed Hbig %*% V (a P_tot^3 matmul of a mostly-zero matrix) per
  # cell, then a dense t(HV) transpose per (k,j) -- O(K * P_tot^3) + O(K^2 * P_tot^2), quartic in the cell count.
  # Equivalent block-sparse form: HVb_k = H_k %*% V[idx_k, ] is Pk x P_tot (the only nonzero rows of H_k V), and
  #   tr(H_k V H_j V) = sum( HVb_k[, idx_j] * t(HVb_j[, idx_k]) )   (Pk x Pj small blocks).
  # Bit-identical arithmetic; no P_tot x P_tot allocations, matmuls, or transposes.
  idxs <- lapply(seq_len(K), function(k) if (Pk[k] == 0L) integer(0) else cstart[k] + seq_len(Pk[k]))
  HVb  <- vector("list", K)
  for (k in seq_len(K)) {
    if (Pk[k] == 0L) next
    HVb[[k]] <- cells[[k]]$ho$H %*% V[idxs[[k]], , drop = FALSE]   # Pk x P_tot
  }
  for (k in seq_len(K)) {
    if (is.null(HVb[[k]])) next
    for (j in k:K) {
      if (is.null(HVb[[j]])) next
      val <- 0.5 * sum(HVb[[k]][, idxs[[j]], drop = FALSE] * t(HVb[[j]][, idxs[[k]], drop = FALSE]))
      Sigma_quad[k, j] <- val
      Sigma_quad[j, k] <- val
    }
  }
  Sigma_quad
}

#' Analytic simultaneous bands for a vector of estimates from its covariance
#'
#' Helper that turns a covariance matrix into (se, crit, lower, upper). When \code{cband = FALSE} the crit
#' is the pointwise \code{qnorm(1 - alp/2)} (no simulation).
#'
#' @param att numeric vector of estimates.
#' @param Sigma covariance matrix of \code{att} (same order); \code{sqrt(diag)} gives the SEs.
#' @param alp significance level. @param cband logical: simultaneous (TRUE) vs pointwise (FALSE).
#' @param seed optional integer for the sup-t simulation.
#' @return list(se, crit, ci_lower, ci_upper).
#' @keywords internal
analytic_bands_edid <- function(att, Sigma, alp = 0.05, cband = TRUE, seed = NULL) {
  Sigma <- as.matrix(Sigma)
  se    <- sqrt(diag(Sigma))
  # Coordinates with degenerate (zero / non-finite) variance carry no band and -- crucially -- are
  # EXCLUDED from supt_crit_edid()'s simultaneous family (its internal ok = is.finite(d) & d > 0).
  # Emit the band over EXACTLY that family so the uniform guarantee is not silently claimed over more
  # coordinates than were simulated; degenerate coordinates get NA bounds rather than a spurious
  # zero-width / NaN interval. On non-degenerate input (the normal path) `good` is all-TRUE, so this is
  # a no-op: same se, same crit, same bounds.
  row_ok <- vapply(seq_len(nrow(Sigma)), function(i) all(is.finite(Sigma[i, ]) & is.finite(Sigma[, i])), logical(1L))
  good <- is.finite(se) & se > 0 & row_ok
  if (isTRUE(cband)) {
    crit <- supt_crit_edid(Sigma, alp = alp, seed = seed)
    if (any(!good)) {
      warning(sprintf(
        paste0("sup-t band: %d of %d coordinate(s) have degenerate variance; they are excluded from ",
               "the simultaneous critical value and their bands are returned as NA."),
        sum(!good), length(se)), call. = FALSE)
    }
  } else {
    crit <- stats::qnorm(1 - alp / 2)
  }
  ci_lower <- att - crit * se
  ci_upper <- att + crit * se
  ci_lower[!good] <- NA_real_
  ci_upper[!good] <- NA_real_
  list(se = se, crit = crit, ci_lower = ci_lower, ci_upper = ci_upper)
}
