# edid-cov.R
# Nuisance estimation functions for the EDiD covariate path.
# Implements sieve (B-spline) estimation of propensity ratios and conditional
# means. The package currently uses plug-in (K=1, no sample splitting)
# nuisance estimation, matching the paper's main-text proposal.

# ---------------------------------------------------------------------------
# Fold assignment
# ---------------------------------------------------------------------------

#' Generate cross-fitting fold assignments
#'
#' Assigns each of \code{n} units to one of \code{K} folds via simple
#' round-robin ordering (after optional random shuffling).
#'
#' @param n positive integer: number of units
#' @param K positive integer: number of folds (default 5)
#' @param seed integer or NULL: if not NULL, set.seed() is called and restored
#'
#' @return integer vector length \code{n}, values in \code{1:K}
#' @keywords internal
build_crossfit_folds_edid <- function(n, K = 5L, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  # Balanced folds: a shuffled round-robin so fold sizes differ by at most one (sampling with
  # replacement could leave a fold empty). Reduces to all-ones at K=1.
  sample(rep_len(seq_len(K), n))
}

# ---------------------------------------------------------------------------
# B-spline basis construction and prediction
# ---------------------------------------------------------------------------

#' Build B-spline basis matrix for a covariate matrix
#'
#' For the first covariate column, fits a B-spline basis with intercept
#' (\code{bs_df} columns). For each additional column, fits without intercept
#' (\code{bs_df - 1} columns, to avoid collinearity). Falls back to a linear
#' basis (intercept + raw column) if \code{splines::bs()} fails.
#'
#' @param X_mat numeric matrix, n x d. May also be a numeric vector (treated
#'   as n x 1).
#' @param bs_df positive integer: degrees of freedom for the B-spline basis
#'   (default 4)
#'
#' @return numeric matrix n x p, with attribute \code{"bs_objects"}: a list of
#'   length d, each element the fitted \code{bs} object for that column (used
#'   by \code{predict_basis_edid()} to evaluate on new data).
#' @keywords internal
build_basis_matrix_edid <- function(X_mat, bs_df = 4L) {
  if (is.vector(X_mat)) X_mat <- matrix(X_mat, ncol = 1L)
  n <- nrow(X_mat)
  d <- ncol(X_mat)

  blocks    <- vector("list", d)
  bs_objects <- vector("list", d)

  for (k in seq_len(d)) {
    xk <- X_mat[, k]
    use_intercept <- (k == 1L)
    # All warnings from splines::bs() are benign in cross-fitting contexts:
    # "boundary knots" fires when test covariates are outside the training range;
    # "interior knots match" fires for binary/factor-derived dummy columns.
    # Errors are caught and fall back to a linear basis.
    bs_result <- tryCatch(
      suppressWarnings(splines::bs(xk, df = bs_df, intercept = use_intercept)),
      error = function(e) NULL
    )
    if (!is.null(bs_result)) {
      bs_objects[[k]] <- bs_result
      blocks[[k]]     <- as.matrix(bs_result)
    } else {
      warning(sprintf("B-spline basis failed for covariate column %d; using linear basis.", k))
      bs_objects[[k]] <- list(fallback = TRUE, use_intercept = use_intercept)
      blocks[[k]]     <- if (use_intercept) cbind(1, xk) else matrix(xk, ncol = 1L)
    }
  }

  B <- do.call(cbind, blocks)
  attr(B, "bs_objects") <- bs_objects
  B
}

#' Predict B-spline basis at new data using stored knot information
#'
#' Evaluates the basis used during training (stored as \code{bs} objects) at
#' new covariate values. When the training basis fell back to a linear basis,
#' returns the linear approximation.
#'
#' @param bs_obj_list list of length d: the \code{"bs_objects"} attribute from
#'   \code{build_basis_matrix_edid()}
#' @param X_new_mat numeric matrix, n_test x d
#'
#' @return numeric matrix n_test x p (same column count as training basis)
#' @keywords internal
predict_basis_edid <- function(bs_obj_list, X_new_mat) {
  if (is.vector(X_new_mat)) X_new_mat <- matrix(X_new_mat, ncol = 1L)
  d      <- length(bs_obj_list)
  blocks <- vector("list", d)

  for (k in seq_len(d)) {
    xk  <- X_new_mat[, k]
    bsk <- bs_obj_list[[k]]
    # NULL or fallback sentinel -> linear basis
    is_fallback <- is.null(bsk) || (is.list(bsk) && isTRUE(bsk$fallback))
    if (is_fallback) {
      use_intercept <- if (is.list(bsk) && !is.null(bsk$use_intercept)) bsk$use_intercept else (k == 1L)
      blocks[[k]] <- if (use_intercept) cbind(1, xk) else matrix(xk, ncol = 1L)
    } else {
      # Suppress the splines::bs() "beyond boundary knots" warning that fires
      # whenever a cross-fitting test fold contains covariate values outside the
      # training fold's knot range.  This is a normal artifact of random splits
      # and does not affect prediction correctness.  All other warnings propagate.
      blocks[[k]] <- withCallingHandlers(
        predict(bsk, newx = xk),
        warning = function(w) {
          if (grepl("boundary knots", conditionMessage(w), fixed = TRUE))
            invokeRestart("muffleWarning")
        }
      )
    }
  }

  do.call(cbind, blocks)
}

# ---------------------------------------------------------------------------
# IC-based sieve-dimension selection (opt-in via edid(bs_df = "ic"))
# ---------------------------------------------------------------------------

#' Select the B-spline df by the paper's information criterion
#'
#' Implements the sieve-index selection rule of Chen, Sant'Anna & Xie (2025)
#' (the display after Eq. (4.2)):
#' \deqn{\widehat{K} = \arg\min_K \ 2\,\mathbb{E}_n[\,\ell_K(\widehat\beta_K)\,]
#'   + C_n K / n, \qquad C_n = \log(n)\ \text{(BIC flavor)},}
#' where \eqn{\ell_K} is the estimator's own convex loss evaluated at the fitted
#' sieve coefficients (for the propensity ratio \eqn{r}:
#' \eqn{\mathbb{E}_n[r^2 G_{g'} - 2 r G_g]}; for the inverse propensity \eqn{s}:
#' \eqn{\mathbb{E}_n[s^2 G_{g'} - 2 s]}; for the conditional mean \eqn{m}: the
#' least-squares loss \eqn{\mathbb{E}_n[G_{g'} (Y_\Delta - m)^2]}), \eqn{K} is the
#' TOTAL basis dimension \code{ncol(B)} implied by the candidate df (the paper's
#' \eqn{\psi^K} dimension, not the per-covariate df), and \eqn{\mathbb{E}_n}
#' averages over the full training sample (\eqn{n = n_{train}}; plug-in regime,
#' the only one \code{edid()} uses). Candidate dfs are \code{grid}; infeasible
#' candidates (\code{fit_loss} returns \code{NULL}, errors, or gives a non-finite
#' loss) are skipped; ties keep the smaller df (parsimony); if every candidate is
#' infeasible the package default \code{4L} is returned.
#'
#' @param fit_loss function(df) returning \code{c(loss = , K = )} (empirical
#'   loss and total basis dimension), or \code{NULL} when the candidate df is
#'   infeasible for this fit
#' @param n training-sample size (the \eqn{\mathbb{E}_n} and penalty denominator)
#' @param grid integer vector of candidate B-spline dfs (default \code{3:8})
#' @return the selected df (integer scalar)
#' @keywords internal
select_bs_df_ic_edid <- function(fit_loss, n, grid = 3:8) {
  best_df <- NA_integer_
  best_ic <- Inf
  Cn      <- log(n)
  for (df_k in grid) {
    fl <- tryCatch(fit_loss(df_k), error = function(e) NULL)
    if (is.null(fl) || !is.finite(fl[["loss"]]) || !is.finite(fl[["K"]])) next
    ic <- 2 * fl[["loss"]] + Cn * fl[["K"]] / n
    if (is.finite(ic) && ic < best_ic) {
      best_ic <- ic
      best_df <- df_k
    }
  }
  if (is.na(best_df)) 4L else as.integer(best_df)
}

# ---------------------------------------------------------------------------
# Propensity ratio estimation
# ---------------------------------------------------------------------------

#' Estimate the propensity ratio r(X) = P(G=g|X) / P(G=g'|X)
#'
#' Implements the sieve (B-spline) estimator for the propensity ratio from
#' Chen, Sant'Anna & Xie (2025) Eq. (4.1)-(4.2). The ratio is estimated via
#' OLS minimising \eqn{E[r(X)^2 G_{g'} - 2 r(X) G_g]}.
#'
#' Closed form:
#' \deqn{\hat\beta = [B_{g'}' B_{g'}]^{-1} \sum_{i: G_i = g} B(X_i)}
#' Then \eqn{\hat r(X_i) = B(X_i)' \hat\beta}.
#'
#' @param X_train numeric matrix n_train x d
#' @param G_train numeric vector n_train: cohort values (Inf for never-treated)
#' @param X_test  numeric matrix n_test x d
#' @param g  scalar: target treatment cohort
#' @param gp scalar: comparison cohort (may be Inf for never-treated)
#' @param bs_df integer B-spline degrees of freedom (default 4), or \code{"ic"}
#'   for the per-fit information-criterion selection of
#'   \code{\link{select_bs_df_ic_edid}} (the paper's BIC-flavored rule)
#'
#' @return numeric vector length n_test: estimated r(X) values. Under
#'   \code{bs_df = "ic"} the selected df is attached as attribute
#'   \code{"edid_bs_df"} (or list element \code{bs_df} when
#'   \code{return_aux = TRUE}).
#' @keywords internal
estimate_propensity_ratio_edid <- function(X_train, G_train, X_test, g, gp,
                                           bs_df = 4L, return_aux = FALSE) {
  n_test <- nrow(X_test)

  # Masks for g and g' units in training data
  mask_gp <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  mask_g  <- (G_train == g)

  n_gp <- sum(mask_gp)
  n_g  <- sum(mask_g)

  if (n_gp < 2L) {
    warning(sprintf(
      "estimate_propensity_ratio_edid: fewer than 2 units in g'=%g training fold; returning 0.", gp
    ))
    fb <- rep(0, n_test)
    return(if (return_aux) list(pred = fb, is_fallback = TRUE) else fb)
  }
  if (n_g < 1L) {
    warning(sprintf(
      "estimate_propensity_ratio_edid: 0 units in g=%g training fold; returning 0.", g
    ))
    fb <- rep(0, n_test)
    return(if (return_aux) list(pred = fb, is_fallback = TRUE) else fb)
  }

  # Opt-in IC sieve-dimension selection (edid(bs_df = "ic")): pick the df in 3:8
  # minimizing the paper's criterion 2*E_n[r^2 G_gp - 2 r G_g] + log(n)*K/n
  # (K = total basis dimension); the winner then takes the standard fitting path
  # below unchanged. The integer default never enters this branch (byte-identical).
  ic_pick <- NULL
  if (identical(bs_df, "ic")) {
    bs_df <- select_bs_df_ic_edid(function(df_k) {
      Bo  <- build_basis_matrix_edid(X_train, df_k)
      B   <- unclass(Bo); attr(B, "bs_objects") <- NULL
      bg  <- B[mask_gp, , drop = FALSE]
      bet <- as.vector(compute_pseudoinverse_edid(t(bg) %*% bg) %*%
                         colSums(B[mask_g, , drop = FALSE]))
      r_in <- drop(B %*% bet)                              # in-sample (train = test under plug-in)
      c(loss = mean(mask_gp * r_in^2 - 2 * mask_g * r_in), K = ncol(B))
    }, n = nrow(X_train))
    ic_pick <- bs_df
  }

  # Build basis on training data
  B_train_obj <- build_basis_matrix_edid(X_train, bs_df)
  B_train     <- unclass(B_train_obj)
  attr(B_train, "bs_objects") <- NULL

  B_gp         <- B_train[mask_gp, , drop = FALSE]
  col_sums_g   <- colSums(B_train[mask_g, , drop = FALSE])

  # beta_hat = [B_gp' B_gp]^{-1} * col_sums_g
  BtB_gp   <- t(B_gp) %*% B_gp
  beta_hat <- as.vector(compute_pseudoinverse_edid(BtB_gp) %*% col_sums_g)

  # Evaluate the fitted sieve basis on the evaluation sample.
  B_test <- predict_basis_edid(attr(B_train_obj, "bs_objects"), X_test)

  # Predict the propensity ratio. The sieve estimate is left unconstrained (no truncation at
  # zero and no link function): truncating would break the linear first-order condition that
  # delivers Neyman orthogonality, so we keep the raw projection. Under good overlap the estimate
  # is positive; large or negative values arise only under weak overlap and are flagged by
  # check_condition_edid().
  r_hat <- drop(B_test %*% beta_hat)

  if (!return_aux) {
    if (!is.null(ic_pick)) attr(r_hat, "edid_bs_df") <- ic_pick
    return(r_hat)
  }

  # ACH (Ackerberg, Chen & Hahn 2012) first-step pieces. The estimating equation is
  #   sum_i [ G_{gp,i} B_i (B_i'beta) - G_{g,i} B_i ] = 0,
  # so the per-unit M-estimator score is s_i = B_i (G_{gp,i} r_i - G_{g,i}) and the
  # Hessian H = (B_gp'B_gp)/n, i.e. H^{-1} = n * pinv(B_gp'B_gp) (same pseudoinverse used
  # for beta_hat, so the score columns sum to ~0 and the EIF stays mean-zero). Valid only
  # for the plug-in (K=1, train=test=full) regime fit_edid_cells enforces with this flag.
  mask_gp_t <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  mask_g_t  <- (G_train == g)
  score_mat <- B_test * (mask_gp_t * r_hat - mask_g_t)   # n x p (row i scaled by G_gp,i r_i - G_g,i)
  H_inv     <- n_test * compute_pseudoinverse_edid(BtB_gp)
  out <- list(pred = r_hat, B_test = B_test, score_mat = score_mat, H_inv = H_inv, is_fallback = FALSE)
  if (!is.null(ic_pick)) out$bs_df <- ic_pick
  out
}

#' Estimate the inverse propensity s(X) = 1 / P(G=g'|X)
#'
#' Implements the sieve estimator for the inverse propensity from
#' Chen, Sant'Anna & Xie (2025) Eq. after (4.2). Estimated via
#' minimising \eqn{E[s(X)^2 G_{g'} - 2 s(X)]}.
#'
#' Closed form:
#' \deqn{\hat\beta = [B_{g'}' B_{g'}]^{-1} \sum_{i=1}^{n} B(X_i)}
#' Then \eqn{\hat s(X_i) = B(X_i)' \hat\beta}, clipped to [0, Inf).
#'
#' @param X_train numeric matrix n_train x d
#' @param G_train numeric vector n_train: cohort values (Inf for never-treated)
#' @param X_test  numeric matrix n_test x d
#' @param gp scalar: cohort whose inverse propensity to estimate
#' @param bs_df integer B-spline degrees of freedom (default 4), or \code{"ic"}
#'   for the per-fit information-criterion selection (see
#'   \code{\link{select_bs_df_ic_edid}}; loss \eqn{\mathbb{E}_n[s^2 G_{g'} - 2 s]})
#'
#' @return numeric vector length n_test: estimated inverse propensities 1/p_g'(X), >= 0.
#'   Under \code{bs_df = "ic"} the selected df is attached as attribute
#'   \code{"edid_bs_df"} (or list element \code{bs_df} when \code{return_aux = TRUE}).
#' @keywords internal
estimate_inverse_propensity_edid <- function(X_train, G_train, X_test, gp,
                                             bs_df = 4L, return_aux = FALSE) {
  n_test <- nrow(X_test)

  mask_gp <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  n_gp <- sum(mask_gp)

  if (n_gp == 0L) {
    stop(sprintf(paste0(
      "estimate_inverse_propensity_edid: 0 units in comparison cohort g'=%g; cannot estimate ",
      "1/p_{g'}(X) (the fallback would divide by zero). Check the comparison group."), gp))
  }
  if (n_gp < 2L) {
    warning(sprintf(
      "estimate_inverse_propensity_edid: fewer than 2 units in g'=%g; returning 1/pi_g'.", gp
    ))
    sv <- rep(length(G_train) / n_gp, n_test)
    return(if (return_aux) list(s_hat = sv, is_fallback = TRUE) else sv)
  }

  # Opt-in IC sieve-dimension selection (edid(bs_df = "ic")): the s-loss is the
  # analogue of the ratio loss with the SECOND term unindicated (min E[s^2 G_gp - 2 s],
  # the paper's display below Eq. (4.2)). Selection uses the raw (unclamped) sieve s,
  # consistent with the linear first-order condition the loss encodes.
  ic_pick <- NULL
  if (identical(bs_df, "ic")) {
    bs_df <- select_bs_df_ic_edid(function(df_k) {
      Bo  <- build_basis_matrix_edid(X_train, df_k)
      B   <- unclass(Bo); attr(B, "bs_objects") <- NULL
      bg  <- B[mask_gp, , drop = FALSE]
      bet <- as.vector(compute_pseudoinverse_edid(t(bg) %*% bg) %*% colSums(B))
      s_in <- drop(B %*% bet)
      c(loss = mean(mask_gp * s_in^2 - 2 * s_in), K = ncol(B))
    }, n = nrow(X_train))
    ic_pick <- bs_df
  }

  B_train_obj <- build_basis_matrix_edid(X_train, bs_df)
  B_train     <- unclass(B_train_obj)
  attr(B_train, "bs_objects") <- NULL

  B_gp       <- B_train[mask_gp, , drop = FALSE]
  col_sums_all <- colSums(B_train)

  BtB_pinv <- compute_pseudoinverse_edid(t(B_gp) %*% B_gp)
  beta_hat <- as.vector(BtB_pinv %*% col_sums_all)

  B_test <- predict_basis_edid(attr(B_train_obj, "bs_objects"), X_test)
  s_raw  <- drop(B_test %*% beta_hat)
  s_hat  <- pmax(s_raw, 0)

  if (return_aux) {
    # M-estimator pieces for the inv_p Sigma_Omega channel (SAME convention as estimate_propensity_ratio_edid).
    # Moment E[B (s 1{gp} - 1)] = 0 (FOC of min E[s^2 G_gp - 2s]); per-unit score s_i = B_i (1{i in gp} s_raw_i - 1)
    # (s_raw = B beta, before the >=0 clamp). Hessian H = E[B B' 1{gp}] = (B_gp'B_gp)/n, so H^{-1} = n * pinv(B_gp'B_gp)
    # (the n factor matches the ratio aux). IF of beta at unit l = -H_inv score_l.
    score_mat <- (as.numeric(mask_gp) * s_raw - 1) * B_test
    out <- list(s_hat = s_hat, B_test = B_test, score_mat = score_mat, H_inv = n_test * BtB_pinv,
                s_pos = s_raw > 0, is_fallback = FALSE)
    if (!is.null(ic_pick)) out$bs_df <- ic_pick
    return(out)
  }

  if (!is.null(ic_pick)) attr(s_hat, "edid_bs_df") <- ic_pick
  s_hat
}


#' Check for extreme propensity ratios and warn once
#' @keywords internal
.check_extreme_ratios_edid <- function(r_vec, g, gp) {
  if (any(is.finite(r_vec)) && max(r_vec, na.rm = TRUE) > 100) {
    warning(sprintf(
      "Extreme propensity ratios detected for g=%g, g'=%g (max > 100). Results may be unstable.",
      g, gp
    ))
  }
}

#' Build the per-comparison overlap-trim masks (ratio-targeted)
#'
#' One \{TRUE, FALSE\} n-vector per comparison key, consumed by
#' \code{edid_cell_trim_structure} (a pair's own mask is \code{keep_inf} for
#' self/two-period pairs and \code{keep_inf * keep_gp} for cross pairs).
#'
#' \strong{Ratio-targeted semantics.} A finite comparison cohort's mask keys on the
#' propensity RATIO \eqn{r_{g,g'}(X)} only -- the actual reweighting factor of that
#' pair's moment (Eq. 4.4 term 3). The inverse propensity \eqn{1/p_{g'}(X)} is a
#' variance-channel object (an Omega* prefactor, never a moment weight) whose
#' absolute scale is \eqn{\approx 1/\pi_{g'}}: thresholding it at a fixed
#' \code{trim_level} mechanically removes (nearly) every unit from any pair whose
#' comparison cohort is small (\eqn{1/\pi_{g'} \ge} \code{trim_level}), killing
#' small-cohort pairs regardless of actual covariate overlap -- the audited
#' dead-pair pathology. The never-treated mask (key \code{"Inf"}) retains the
#' legacy definition (ratio AND inverse propensity), which keeps the PT-Post and
#' self-pair behavior byte-identical.
#'
#' @param prop_ratios named list of n-vectors (keys \code{"Inf"} and finite gp's)
#' @param inv_propensities named list of n-vectors (same key space, superset)
#' @param trim_level positive scalar; \code{Inf} disables trimming
#' @param n number of units
#' @return named list of logical n-vectors, or \code{NULL} when trimming is off
#' @keywords internal
#' @noRd
build_trim_keep_edid <- function(prop_ratios, inv_propensities, trim_level, n) {
  if (!is.finite(trim_level)) return(NULL)
  if (is.null(prop_ratios) && is.null(inv_propensities)) return(NULL)
  ks <- union(names(prop_ratios), names(inv_propensities))
  if (!length(ks)) return(NULL)
  stats::setNames(lapply(ks, function(k) {
    keep <- rep(TRUE, n)
    rr <- prop_ratios[[k]]
    if (!is.null(rr)) keep <- keep & (abs(rr) < trim_level)
    if (identical(k, "Inf")) {                       # legacy never-treated mask: ratio AND 1/p_NT
      ip <- inv_propensities[[k]]
      if (!is.null(ip)) keep <- keep & (abs(ip) < trim_level)
    }
    keep
  }), ks)
}

# ---------------------------------------------------------------------------
# (Removed 2026-06-12) The "coherent" multinomial-logit sieve engine -- the
# ridge-logistic per-cohort system fitter, its joint stacked-coefficient aux
# assembler, and the cross-cohort ratio / inverse-propensity aux packers -- was
# evaluated against the default "exp" engine and removed: a thin-cohort-share
# Monte Carlo (quality_reports/drafts/gate_runs/ratio_method_thinshares_mc.md)
# found it hard-fails ~6.7% of draws (fitted p = 0 -> 1/p = Inf) and is
# anti-conservative in the body, with uniformly worse CI coverage than "exp";
# the comfortable-share comparison (ratio_method_comparison.md) showed no
# offsetting advantage. The general first-step infrastructure the coherent audit
# produced (bootstrap coef_id dedup, trim/keep threading, correlation-scale eigen
# floors, the $args refit snapshot, link-aware FD steps) is retained.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Exponential-link Riesz regressions (ratio_method = "exp")
# ---------------------------------------------------------------------------

#' Newton solver for the tailored exponential-link Riesz loss
#'
#' Minimizes the globally convex tailored loss
#' \deqn{L(\beta) = \mathbb{E}_n[\,\mathrm{comp}_i\, e^{\psi_i'\beta}\,] - t_{col}'\beta
#'       \;(+\tfrac{\lambda}{2}\,\beta' \mathrm{diag}(pen)\,\beta\ \text{on rescue}),}
#' whose first-order condition is EXACT basis-mean balancing:
#' \eqn{\mathbb{E}_n[\psi_i\, e^{\psi_i'\beta}\, \mathrm{comp}_i] = t_{col}}. For the
#' propensity ratio \eqn{r_{g,g'}} take \eqn{\mathrm{comp} = G_{g'}} and
#' \eqn{t_{col} = \mathbb{E}_n[\psi\, G_g]} (population minimizer
#' \eqn{\psi'\beta^* = \log(p_g/p_{g'})}); for the inverse propensity \eqn{s_{g'}} take
#' \eqn{\mathrm{comp} = G_{g'}} and \eqn{t_{col} = \mathbb{E}_n[\psi]} (population minimizer
#' \eqn{\log(1/p_{g'})}).
#'
#' Safeguards: (i) basis columns with no comparison-side support (\code{colSums(|B| comp) ~ 0})
#' cannot be balanced -- the loss is unbounded below along them -- so they are pinned at 0
#' ("dead" columns; the optimization runs on the live block); (ii) Newton steps are
#' step-halved on the loss (up to 30 halvings); (iii) exp overflow is capped inside the
#' objective and treated as non-descent; (iv) if the unpenalized problem does not converge
#' (balancing infeasible / quasi-separation), a scale-normalized ridge
#' \eqn{pen = \lambda\,\mathrm{colMeans}(B^2)/n} is escalated over \code{ridge_grid} until
#' the (then strictly convex, coercive) problem converges -- the \eqn{1/n} scaling keeps an
#' \eqn{O(1)} penalty against the \eqn{O(n)}-scale criterion,
#' so a fixed \eqn{\lambda} stays asymptotically negligible. The returned \code{pen} vector
#' is 0 except on rescue, and is folded into the aux score/Hessian by the callers so the
#' M-estimator pieces stay mean-zero at the fitted coefficients.
#'
#' @param B numeric matrix n x p (sieve basis at the training sample)
#' @param comp 0/1 comparison-group indicator length n
#' @param tcol length-p target basis means (see above)
#' @param beta_init optional warm start (length p)
#' @param maxit,tol Newton controls (tol is relative on the gradient sup-norm)
#' @param ridge_grid increasing ridge scales tried after the unpenalized fit fails
#' @return list(beta, pen, converged, lambda, n_iter, live)
#' @keywords internal
#' @noRd
fit_exp_riesz_edid <- function(B, comp, tcol, beta_init = NULL, maxit = 200L, tol = 1e-8,
                               ridge_grid = c(0, 1, 100, 1e4)) {
  n <- nrow(B); p <- ncol(B)
  comp <- as.numeric(comp)
  live <- colSums(abs(B) * comp) > 1e-12          # balanceable directions (comparison-side support)
  scale_t <- 1 + max(abs(tcol))
  best <- NULL
  for (lam in ridge_grid) {
    pen <- lam * colMeans(B * B) / n              # O(1) penalty vs the O(n)-scale criterion
    pen[!live] <- 0
    beta <- if (!is.null(beta_init) && length(beta_init) == p && all(is.finite(beta_init))) beta_init else numeric(p)
    beta[!live] <- 0
    loss_fn <- function(b) {
      eta <- drop(B %*% b)
      if (max(eta) > 350) return(Inf)             # exp overflow guard (squares of exp(350) stay finite)
      mean(comp * exp(eta)) - sum(tcol * b) + 0.5 * sum(pen * b * b)
    }
    l0 <- loss_fn(beta)
    if (!is.finite(l0)) { beta <- numeric(p); l0 <- loss_fn(beta) }
    converged <- FALSE; it_used <- 0L
    for (it in seq_len(maxit)) {
      it_used <- it
      eta  <- drop(B %*% beta)
      w    <- comp * exp(pmin(eta, 350))
      grad <- drop(crossprod(B, w)) / n - tcol + pen * beta
      if (max(abs(grad[live])) < tol * scale_t) { converged <- TRUE; break }
      Hm   <- crossprod(B, w * B) / n
      if (lam > 0) diag(Hm) <- diag(Hm) + pen
      Hl   <- Hm[live, live, drop = FALSE]
      step <- tryCatch(-solve(Hl, grad[live]),
                       error = function(e) -drop(compute_pseudoinverse_edid(Hl) %*% grad[live]))
      if (!all(is.finite(step))) break
      ok <- FALSE; fac <- 1; lnew <- l0
      for (hh in 1:30) {                          # step-halving on the (penalized) loss
        bnew <- beta; bnew[live] <- beta[live] + fac * step
        lnew <- loss_fn(bnew)
        if (is.finite(lnew) && lnew <= l0 + 1e-12 * (1 + abs(l0))) { ok <- TRUE; break }
        fac <- fac / 2
      }
      if (!ok) break                              # no descent direction left: stop (converged stays FALSE)
      beta[live] <- beta[live] + fac * step
      moved <- max(abs(fac * step))
      l0 <- lnew
      if (moved < 1e-12 * (1 + max(abs(beta)))) { # negligible step: accept if the gradient is small-ish
        eta  <- drop(B %*% beta)
        w    <- comp * exp(pmin(eta, 350))
        grad <- drop(crossprod(B, w)) / n - tcol + pen * beta
        converged <- max(abs(grad[live])) < sqrt(tol) * scale_t
        break
      }
    }
    if (converged && all(is.finite(beta))) return(list(beta = beta, pen = pen, converged = TRUE,
                                                       lambda = lam, n_iter = it_used, live = live))
    if (is.null(best) && all(is.finite(beta))) best <- list(beta = beta, pen = pen, converged = FALSE,
                                                            lambda = lam, n_iter = it_used, live = live)
  }
  if (is.null(best)) best <- list(beta = numeric(p), pen = numeric(p), converged = FALSE,
                                  lambda = NA_real_, n_iter = 0L, live = live)
  best
}

#' Warm starts for the exponential-link Riesz fit
#'
#' Returns the better (lower tailored loss) of two candidates: (a) the constant fit at the
#' closed-form level \code{const_level} (= \eqn{n_g/n_{g'}} for the ratio, \eqn{n/n_{g'}} for
#' the inverse propensity), represented as \code{log(const_level) * q} with \code{q} the LS
#' projection of the constant 1 onto the basis (skipped when constants are not in the span);
#' (b) the log of the CLIPPED per-target LS sieve fit (the paper's closed-form linear fit,
#' floored at a small positive value) projected back onto the basis. Either may be the zero
#' vector when degenerate; the zero start is always a valid fallback (the loss is globally
#' convex, so the warm start affects speed and overflow risk, not the optimum).
#' @keywords internal
#' @noRd
exp_riesz_warmstart_edid <- function(B, comp, tcol, const_level) {
  n <- nrow(B); p <- ncol(B)
  BtB_pinv <- compute_pseudoinverse_edid(crossprod(B))
  cands <- list(numeric(p))
  # (a) constant log level through the basis (B-spline first block is a partition of unity)
  if (is.finite(const_level) && const_level > 0) {
    q <- drop(BtB_pinv %*% colSums(B))
    if (max(abs(drop(B %*% q) - 1)) < 0.01) cands[[length(cands) + 1L]] <- log(const_level) * q
  }
  # (b) log of the clipped linear (paper LS) fit
  beta_ls <- tryCatch(
    drop(compute_pseudoinverse_edid(crossprod(B, comp * B)) %*% (n * tcol)),
    error = function(e) NULL)
  if (!is.null(beta_ls) && all(is.finite(beta_ls))) {
    r_ls <- drop(B %*% beta_ls)
    clip <- max(1e-6, 1e-3 * stats::median(abs(r_ls)))
    bl   <- drop(BtB_pinv %*% crossprod(B, log(pmax(r_ls, clip))))
    if (all(is.finite(bl))) cands[[length(cands) + 1L]] <- bl
  }
  loss0 <- function(b) {
    eta <- drop(B %*% b)
    if (max(eta) > 350) return(Inf)
    mean(comp * exp(eta)) - sum(tcol * b)
  }
  ls <- vapply(cands, loss0, numeric(1))
  cands[[which.min(replace(ls, !is.finite(ls), Inf))]]
}

#' Literal paper-loss refinement of an exponential-link fit (internal cross-check)
#'
#' Damped (Levenberg) Newton minimization of the PAPER's loss with the exp parametrization,
#' \eqn{L_2(\beta) = \mathbb{E}_n[\,\mathrm{comp}\, e^{2\psi'\beta} - 2\,\mathrm{tgt}\,
#' e^{\psi'\beta}\,]} (Eq. (4.1)'s \eqn{r^2 G_{g'} - 2 r G_g} with \eqn{r = e^{\psi'\beta}};
#' for \eqn{s}, \code{tgt = 1}), warm-started from the tailored solution. \eqn{L_2} is not
#' globally convex in \eqn{\beta} (difference of convex), and in finite samples it can be
#' UNBOUNDED below along directions that raise \eqn{\eta} where the target has basis mass
#' but the comparison has (essentially) none -- the multiplicative analogue of the LS sieve's
#' thin-support pathology, which an unconstrained quasi-Newton line search will find and
#' exploit. The refinement therefore stays in a TRUST REGION around the tailored solution
#' (steps that push \eqn{\max\eta} more than 5 above the warm start's are rejected), uses
#' loss-based step-halving with Levenberg damping when the Hessian is indefinite, and is
#' accepted only at an interior stationary point (small gradient); otherwise the caller keeps
#' the tailored fit and its aux. Under correct specification both losses share the population
#' minimizer, so the refinement converges and the two fits agree. Reached via
#' \code{options(edid_exp_loss = "paper")}; the default \code{"tailored"} never calls this.
#' @keywords internal
#' @noRd
exp_riesz_paper_refine_edid <- function(B, comp, tgt, beta0, maxit = 50L) {
  n <- nrow(B); p <- ncol(B)
  eta_max0 <- max(drop(B %*% beta0))
  eta_cap  <- eta_max0 + 5                         # trust bound against the empirical-divergence escape
  fn <- function(eta) mean(comp * exp(2 * eta) - 2 * tgt * exp(eta))
  beta <- beta0
  eta  <- drop(B %*% beta)
  f0   <- fn(eta)
  grad <- drop(crossprod(B, 2 * (comp * exp(2 * eta) - tgt * exp(eta)))) / n
  tol_g <- 1e-8 * (1 + max(abs(grad)))
  converged <- max(abs(grad)) < tol_g
  it <- 0L
  while (!converged && it < maxit) {
    it <- it + 1L
    Hm <- crossprod(B, (2 * (2 * comp * exp(2 * eta) - tgt * exp(eta))) * B) / n
    mu <- 0; accepted <- FALSE
    for (damp in 1:8) {                            # Levenberg escalation if the step is not descent
      Hd <- Hm; if (mu > 0) diag(Hd) <- diag(Hd) + mu * (1 + abs(diag(Hm)))
      step <- tryCatch(-solve(Hd, grad),
                       error = function(e) -drop(compute_pseudoinverse_edid(Hd) %*% grad))
      if (all(is.finite(step))) {
        fac <- 1
        for (hh in 1:25) {
          bnew <- beta + fac * step
          eta_new <- drop(B %*% bnew)
          if (max(eta_new) <= eta_cap) {           # stay in the trust region
            fnew <- fn(eta_new)
            if (is.finite(fnew) && fnew <= f0 - 1e-14 * (1 + abs(f0))) {
              beta <- bnew; eta <- eta_new; f0 <- fnew; accepted <- TRUE
              break
            }
          }
          fac <- fac / 2
        }
      }
      if (accepted) break
      mu <- if (mu == 0) 1e-4 else mu * 100
    }
    if (!accepted) break                           # no acceptable descent step: stop
    grad <- drop(crossprod(B, 2 * (comp * exp(2 * eta) - tgt * exp(eta)))) / n
    converged <- max(abs(grad)) < tol_g
  }
  if (!converged || !all(is.finite(beta))) return(list(beta = beta0, converged = FALSE))
  list(beta = beta, converged = TRUE)
}

#' Estimate the propensity ratio r(X) = p_g(X)/p_g'(X) by the exponential-link Riesz regression
#'
#' The \code{ratio_method = "exp"} engine (see \code{\link{edid}}): the paper-compatible
#' per-target alternative to the LS sieve of \code{estimate_propensity_ratio_edid}, with
#' \eqn{\hat r_{g,g'}(X) = \exp(\psi^K(X)'\hat\beta)} -- positive by construction -- fit
#' INDEPENDENTLY per (g, g') on the same B-spline machinery. The primary criterion is the
#' tailored convex loss \eqn{\mathbb{E}_n[e^{\psi'\beta} G_{g'} - (\psi'\beta) G_g]}
#' (\code{fit_exp_riesz_edid}): globally convex, FOC = exact basis-mean balancing
#' \eqn{\mathbb{E}_n[\psi\, \hat r\, G_{g'}] = \mathbb{E}_n[\psi\, G_g]}, population target
#' \eqn{\log(p_g/p_{g'})}. Unlike the per-pair LS sieve, no thin-denominator Gram is inverted
#' on the raw scale: the exp link rules out negative fitted "ratios" entirely.
#' \code{options(edid_exp_loss = "paper")} refines by the literal paper loss
#' (\code{exp_riesz_paper_refine_edid}) for cross-checking.
#'
#' \strong{Estimation-effect aux (full integration, no fallback-marking).} Under
#' \code{return_aux = TRUE} (plug-in regime) the M-estimator pieces are returned in the SAME
#' contract every correction consumes, with the exp-link chain rule baked in:
#' \itemize{
#'   \item \code{B_test} = \eqn{\partial \hat r/\partial\beta = \hat r\,\psi} (n x p Jacobian;
#'     consumers use \code{B_test} as the prediction-perturbation direction and as the Gamma
#'     basis \eqn{\Gamma = n^{-1} B_{test}'s}, both of which are exactly the coefficient
#'     derivative under this packing);
#'   \item \code{score_mat} = \eqn{\psi_i\,(G_{g',i}\hat r_i - G_{g,i})} (+ the ridge term on
#'     rescue), the tailored-loss score -- mean-zero at \eqn{\hat\beta};
#'   \item \code{H_inv} = \eqn{n\,[\Psi'\mathrm{diag}(G_{g'}\hat r)\Psi + n\,\mathrm{diag}(pen)]^{-}},
#'     the tailored-loss Hessian (positive semi-definite, same pseudoinverse convention as the
#'     LS aux). Under \code{options(edid_exp_loss = "paper")} the score/Hessian are the paper
#'     loss's: score \eqn{\psi\,\hat r\,(G_{g'}\hat r - G_g)}, Hessian
#'     \eqn{n^{-1}\Psi'\mathrm{diag}(2 G_{g'}\hat r^2 - G_g \hat r)\Psi}.
#' }
#' \code{link = "exp"} and the raw basis \code{B_raw} ride along for the FD oracles.
#'
#' @inheritParams estimate_propensity_ratio_edid
#' @return as \code{estimate_propensity_ratio_edid} (vector, or aux list under
#'   \code{return_aux = TRUE}); predictions are strictly positive
#' @keywords internal
estimate_propensity_ratio_exp_edid <- function(X_train, G_train, X_test, g, gp,
                                               bs_df = 4L, return_aux = FALSE) {
  n_test <- nrow(X_test)
  mask_gp <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  mask_g  <- (G_train == g)
  n_gp <- sum(mask_gp); n_g <- sum(mask_g)
  n_train <- length(G_train)

  if (n_gp < 2L) {
    warning(sprintf(
      "estimate_propensity_ratio_exp_edid: fewer than 2 units in g'=%g training fold; returning 0.", gp))
    fb <- rep(0, n_test)
    return(if (return_aux) list(pred = fb, is_fallback = TRUE) else fb)
  }
  if (n_g < 1L) {
    warning(sprintf(
      "estimate_propensity_ratio_exp_edid: 0 units in g=%g training fold; returning 0.", g))
    fb <- rep(0, n_test)
    return(if (return_aux) list(pred = fb, is_fallback = TRUE) else fb)
  }

  comp <- as.numeric(mask_gp)
  fit_at_df <- function(df_k) {
    Bo <- build_basis_matrix_edid(X_train, df_k)
    B  <- unclass(Bo); attr(B, "bs_objects") <- NULL
    tcol <- colSums(B * mask_g) / n_train
    w0 <- exp_riesz_warmstart_edid(B, comp, tcol, const_level = n_g / n_gp)
    list(Bo = Bo, B = B, tcol = tcol, ft = fit_exp_riesz_edid(B, comp, tcol, beta_init = w0))
  }

  # Opt-in IC sieve-dimension selection (edid(bs_df = "ic")): the estimator's OWN convex loss
  # (the tailored loss, unpenalized) in the paper's criterion 2*loss + log(n)*K/n.
  ic_pick <- NULL
  if (identical(bs_df, "ic")) {
    bs_df <- select_bs_df_ic_edid(function(df_k) {
      fk <- fit_at_df(df_k)
      if (!fk$ft$converged || !all(is.finite(fk$ft$beta))) return(NULL)
      eta <- drop(fk$B %*% fk$ft$beta)
      c(loss = mean(comp * exp(pmin(eta, 350)) - mask_g * eta), K = ncol(fk$B))
    }, n = n_train)
    ic_pick <- bs_df
  }

  fk <- fit_at_df(bs_df)
  ft <- fk$ft
  if (!all(is.finite(ft$beta)) || (!ft$converged && is.na(ft$lambda))) {
    warning(sprintf(
      "estimate_propensity_ratio_exp_edid: exp-link fit failed for g=%g vs g'=%g; using the constant share ratio.",
      g, gp))
    fb <- rep(n_g / n_gp, n_test)
    return(if (return_aux) list(pred = fb, is_fallback = TRUE) else fb)
  }
  beta <- ft$beta
  pen  <- ft$pen

  # Literal paper-loss refinement (internal cross-check flag; default "tailored" skips this).
  # A rejected refinement keeps the tailored fit AND its aux (the score/Hessian must encode
  # the loss whose stationary point beta-hat actually is).
  loss_used <- getOption("edid_exp_loss", "tailored")
  if (identical(loss_used, "paper")) {
    pf <- exp_riesz_paper_refine_edid(fk$B, comp, as.numeric(mask_g), beta)
    if (pf$converged) {
      beta <- pf$beta
      pen  <- numeric(length(beta))                # the refinement is unpenalized
    } else {
      loss_used <- "tailored"
    }
  }

  B_test <- predict_basis_edid(attr(fk$Bo, "bs_objects"), X_test)
  r_hat  <- exp(pmin(drop(B_test %*% beta), 350))

  if (!return_aux) {
    if (!is.null(ic_pick)) attr(r_hat, "edid_bs_df") <- ic_pick
    return(r_hat)
  }

  # M-estimator pieces (plug-in regime: train = test = full sample, as the LS aux).
  mask_gp_t <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  mask_g_t  <- (G_train == g)
  if (identical(loss_used, "paper")) {
    score_mat <- B_test * (r_hat * (mask_gp_t * r_hat - mask_g_t))
    Hmat      <- crossprod(B_test, (2 * mask_gp_t * r_hat^2 - mask_g_t * r_hat) * B_test)
  } else {
    score_mat <- B_test * (mask_gp_t * r_hat - mask_g_t)
    Hmat      <- crossprod(B_test, (mask_gp_t * r_hat) * B_test)
    if (any(pen > 0)) {                            # ridge rescue: keep the aux mean-zero at beta-hat
      score_mat <- score_mat + matrix(pen * beta, n_test, ncol(B_test), byrow = TRUE)
      Hmat      <- Hmat + n_test * diag(pen, ncol(B_test))
    }
  }
  H_inv <- n_test * compute_pseudoinverse_edid(Hmat)
  out <- list(pred = r_hat, B_test = B_test * r_hat, score_mat = score_mat, H_inv = H_inv,
              is_fallback = FALSE, link = "exp", B_raw = B_test, beta = beta,
              exp_converged = ft$converged, exp_lambda = ft$lambda, exp_loss = loss_used)
  if (!is.null(ic_pick)) out$bs_df <- ic_pick
  out
}

#' Estimate the inverse propensity 1/p_g'(X) by the exponential-link Riesz regression
#'
#' The \code{ratio_method = "exp"} engine for the FINITE-cohort inverse propensities (the
#' \eqn{\Omega^*} variance prefactors): \eqn{\hat s_{g'}(X) = \exp(\psi^K(X)'\hat\beta) > 0},
#' fit by the tailored convex loss \eqn{\mathbb{E}_n[e^{\psi'\beta} G_{g'} - \psi'\beta]}
#' (FOC: \eqn{\mathbb{E}_n[\psi\,\hat s\,G_{g'}] = \mathbb{E}_n[\psi]}; population target
#' \eqn{\log(1/p_{g'})}). Same solver, warm starts, ridge rescue, paper-loss flag, and
#' full-aux contract as \code{estimate_propensity_ratio_exp_edid} (chain rule
#' \eqn{\partial\hat s/\partial\beta = \hat s\,\psi} packed as \code{B_test}; tailored score
#' \eqn{\psi(G_{g'}\hat s - 1)}; Hessian \eqn{n^{-1}\Psi'\mathrm{diag}(G_{g'}\hat s)\Psi}).
#' \code{s_pos} is all-TRUE (the exp fit is never clamped), so the analytic inv-p
#' weight-channel correction covers this fit with no masked rows.
#'
#' @inheritParams estimate_inverse_propensity_edid
#' @return as \code{estimate_inverse_propensity_edid}; predictions strictly positive
#' @keywords internal
estimate_inverse_propensity_exp_edid <- function(X_train, G_train, X_test, gp,
                                                 bs_df = 4L, return_aux = FALSE) {
  n_test <- nrow(X_test)
  mask_gp <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  n_gp <- sum(mask_gp)
  n_train <- length(G_train)

  if (n_gp == 0L) {
    stop(sprintf(paste0(
      "estimate_inverse_propensity_exp_edid: 0 units in comparison cohort g'=%g; cannot estimate ",
      "1/p_{g'}(X). Check the comparison group."), gp))
  }
  if (n_gp < 2L) {
    warning(sprintf(
      "estimate_inverse_propensity_exp_edid: fewer than 2 units in g'=%g; returning 1/pi_g'.", gp))
    sv <- rep(n_train / n_gp, n_test)
    return(if (return_aux) list(s_hat = sv, is_fallback = TRUE) else sv)
  }

  comp <- as.numeric(mask_gp)
  fit_at_df <- function(df_k) {
    Bo <- build_basis_matrix_edid(X_train, df_k)
    B  <- unclass(Bo); attr(B, "bs_objects") <- NULL
    tcol <- colMeans(B)
    w0 <- exp_riesz_warmstart_edid(B, comp, tcol, const_level = n_train / n_gp)
    list(Bo = Bo, B = B, tcol = tcol, ft = fit_exp_riesz_edid(B, comp, tcol, beta_init = w0))
  }

  ic_pick <- NULL
  if (identical(bs_df, "ic")) {
    bs_df <- select_bs_df_ic_edid(function(df_k) {
      fk <- fit_at_df(df_k)
      if (!fk$ft$converged || !all(is.finite(fk$ft$beta))) return(NULL)
      eta <- drop(fk$B %*% fk$ft$beta)
      c(loss = mean(comp * exp(pmin(eta, 350)) - eta), K = ncol(fk$B))
    }, n = n_train)
    ic_pick <- bs_df
  }

  fk <- fit_at_df(bs_df)
  ft <- fk$ft
  if (!all(is.finite(ft$beta)) || (!ft$converged && is.na(ft$lambda))) {
    warning(sprintf(
      "estimate_inverse_propensity_exp_edid: exp-link fit failed for g'=%g; using 1/pi_g'.", gp))
    sv <- rep(n_train / n_gp, n_test)
    return(if (return_aux) list(s_hat = sv, is_fallback = TRUE) else sv)
  }
  beta <- ft$beta
  pen  <- ft$pen

  loss_used <- getOption("edid_exp_loss", "tailored")
  if (identical(loss_used, "paper")) {
    pf <- exp_riesz_paper_refine_edid(fk$B, comp, rep(1, n_train), beta)
    if (pf$converged) {
      beta <- pf$beta
      pen  <- numeric(length(beta))
    } else {
      loss_used <- "tailored"                      # rejected refinement: keep the tailored fit + aux
    }
  }

  B_test <- predict_basis_edid(attr(fk$Bo, "bs_objects"), X_test)
  s_hat  <- exp(pmin(drop(B_test %*% beta), 350))

  if (return_aux) {
    mask_gp_t <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
    if (identical(loss_used, "paper")) {
      score_mat <- B_test * (s_hat * (mask_gp_t * s_hat - 1))
      Hmat      <- crossprod(B_test, (2 * mask_gp_t * s_hat^2 - s_hat) * B_test)
    } else {
      score_mat <- B_test * (as.numeric(mask_gp_t) * s_hat - 1)
      Hmat      <- crossprod(B_test, (mask_gp_t * s_hat) * B_test)
      if (any(pen > 0)) {
        score_mat <- score_mat + matrix(pen * beta, n_test, ncol(B_test), byrow = TRUE)
        Hmat      <- Hmat + n_test * diag(pen, ncol(B_test))
      }
    }
    out <- list(s_hat = s_hat, B_test = B_test * s_hat, score_mat = score_mat,
                H_inv = n_test * compute_pseudoinverse_edid(Hmat),
                s_pos = rep(TRUE, n_test), is_fallback = FALSE, link = "exp",
                B_raw = B_test, beta = beta,
                exp_converged = ft$converged, exp_lambda = ft$lambda, exp_loss = loss_used)
    if (!is.null(ic_pick)) out$bs_df <- ic_pick
    return(out)
  }

  if (!is.null(ic_pick)) attr(s_hat, "edid_bs_df") <- ic_pick
  s_hat
}

#' Estimate the conditional mean \eqn{E[Y_s - Y_1 | G=g', X]}
#'
#' Fits an OLS B-spline regression of \code{Y_delta} on \code{B(X)} using only
#' units with \code{G_train == gp}, then predicts for all test units.
#'
#' @param X_train numeric matrix n_train x d
#' @param Y_delta_train numeric vector n_train: Y_s - Y_1 for all training units
#' @param G_train numeric vector n_train: cohort values (Inf for never-treated)
#' @param X_test  numeric matrix n_test x d
#' @param gp scalar: cohort to regress on (may be Inf)
#' @param bs_df integer B-spline degrees of freedom (default 4), or \code{"ic"}
#'   for the per-fit information-criterion selection (see
#'   \code{\link{select_bs_df_ic_edid}}; least-squares loss
#'   \eqn{\mathbb{E}_n[G_{g'} (Y_\Delta - m)^2]})
#'
#' @return numeric vector length n_test. Under \code{bs_df = "ic"} the selected
#'   df is attached as attribute \code{"edid_bs_df"} (or list element
#'   \code{bs_df} when \code{return_aux = TRUE}).
#' @keywords internal
estimate_conditional_mean_edid <- function(X_train, Y_delta_train, G_train,
                                           X_test, gp, bs_df = 4L, return_aux = FALSE) {
  n_test  <- nrow(X_test)

  mask_gp <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  n_gp    <- sum(mask_gp)

  if (n_gp < 2L) {
    fallback_val <- if (n_gp == 1L) Y_delta_train[mask_gp] else 0
    warning(sprintf(
      "estimate_conditional_mean_edid: fewer than 2 units in g'=%g training fold; using constant.", gp
    ))
    fb <- rep(fallback_val, n_test)
    return(if (return_aux) list(pred = fb, is_fallback = TRUE) else fb)
  }

  X_gp    <- X_train[mask_gp, , drop = FALSE]
  y_gp    <- Y_delta_train[mask_gp]

  # Opt-in IC sieve-dimension selection (edid(bs_df = "ic")): least-squares loss
  # E_n[G_gp (Y_delta - m)^2] over the FULL training sample (off-cohort terms are
  # zero), penalty log(n)*K/n. Candidates with more basis columns than cohort
  # observations are infeasible (mirrors the n_gp <= p_basis fallback below).
  ic_pick <- NULL
  if (identical(bs_df, "ic")) {
    n_train <- length(Y_delta_train)
    bs_df <- select_bs_df_ic_edid(function(df_k) {
      Bo <- build_basis_matrix_edid(X_gp, df_k)
      B  <- unclass(Bo); attr(B, "bs_objects") <- NULL
      if (n_gp <= ncol(B)) return(NULL)
      ft <- solve_ols_edid(B, y_gp)
      c(loss = sum(ft$residuals^2) / n_train, K = ncol(B))
    }, n = n_train)
    ic_pick <- bs_df
  }

  B_gp_train_obj <- build_basis_matrix_edid(X_gp, bs_df)
  B_gp_train     <- unclass(B_gp_train_obj)
  attr(B_gp_train, "bs_objects") <- NULL

  # Check minimum sample for basis dimension
  p_basis <- ncol(B_gp_train)
  if (n_gp <= p_basis) {
    # Fewer obs than basis columns: use simpler 1-column basis (intercept only)
    fit <- list(coef = mean(y_gp))
    m_hat <- rep(mean(y_gp), n_test)
    return(if (return_aux) list(pred = m_hat, is_fallback = TRUE) else m_hat)
  }

  fit   <- solve_ols_edid(B_gp_train, y_gp)
  B_test <- predict_basis_edid(attr(B_gp_train_obj, "bs_objects"), X_test)
  m_hat <- drop(B_test %*% fit$coef)

  if (!return_aux) {
    if (!is.null(ic_pick)) attr(m_hat, "edid_bs_df") <- ic_pick
    return(m_hat)
  }

  # ACH (Ackerberg, Chen & Hahn 2012) first-step pieces for the within-cohort OLS
  #   sum_{i: G=gp} B_i (Y_delta_i - B_i'beta) = 0.
  # SIGN: the OLS estimating moment B*resid has Jacobian E[d s/d beta'] = -E[G_gp BB'] (NEGATIVE),
  # so the first-step influence function of beta_hat is +H^{-1} s and the two-step correction must be
  # ADDED. The correction helper uses a uniform "psi - score %*% (H_inv %*% Gamma)" (subtract)
  # convention -- correct for the propensity RATIO, whose moment B*(G_gp r - G_g) has the OPPOSITE
  # (positive) Jacobian +E[G_gp BB']. To make the shared subtract convention correct for this OLS
  # channel too, the score carries a leading MINUS: score = -B*resid. (Validated: the corrected EIF
  # then matches the numerical two-step IF, cor = +1; with +B*resid it is exactly negated, cor = -1.)
  # Hessian H = (B_gp'B_gp)/n => H^{-1} = n * pinv(B_gp'B_gp), FULL-sample n (the score is zero
  # off-cohort, so the M-estimator average is over n, not n_gp); same pseudoinverse as beta. Plug-in only.
  resid_full           <- numeric(n_test)
  resid_full[mask_gp]  <- fit$residuals
  score_mat <- -B_test * resid_full                      # n x p (leading minus: see SIGN note above)
  H_inv     <- n_test * compute_pseudoinverse_edid(crossprod(B_gp_train))
  out <- list(pred = m_hat, B_test = B_test, score_mat = score_mat, H_inv = H_inv, is_fallback = FALSE)
  if (!is.null(ic_pick)) out$bs_df <- ic_pick
  out
}

# ---------------------------------------------------------------------------
# Cross-fitted nuisance estimation (aggregate over folds)
# ---------------------------------------------------------------------------

#' Estimate propensity ratios for all comparison cohorts via cross-fitting
#'
#' For each unique \code{gp} in \code{pairs}, produces a full-sample n-vector of
#' \eqn{\hat r_{g, g'}(X_i)}.
#'
#' \strong{Ratio construction (\code{ratio_method}).} The never-treated ratio
#' \eqn{r_{g,\infty}} is ALWAYS estimated by the paper's per-pair LS sieve
#' (\code{estimate_propensity_ratio_edid}; its denominator group is the large
#' never-treated pool, the well-conditioned case). For \emph{finite} comparison
#' cohorts \eqn{g' \ne g} (the cross-cohort pairs of the PT-All moment set):
#' \describe{
#'   \item{\code{"exp"}}{the exponential-link Riesz regression of
#'     \code{estimate_propensity_ratio_exp_edid} for EVERY \code{gp} -- including the
#'     never-treated pool: each ratio is an independent per-target fit
#'     \eqn{\hat r = \exp(\psi'\hat\beta)} of the tailored convex balancing loss
#'     (positive by construction, paper-loss-compatible). Under
#'     \code{return_aux = TRUE} every entry carries FULL M-estimator aux (chain-rule
#'     Jacobian \code{B_test}, tailored score, Hessian inverse) with
#'     \code{is_fallback = FALSE}: the ACH / higher-order / gmm / bootstrap first-step
#'     corrections COVER the cross-cohort channels (no fallback-skipping). This is
#'     \code{edid()}'s default.}
#'   \item{\code{"direct"}}{the paper's literal per-pair LS sieve for every
#'     \code{gp} (byte-identical legacy behavior; retained for forensics).}
#' }
#'
#' @param panel_obj panel object with \code{covariate_matrix} and
#'   \code{unit_cohorts}
#' @param g scalar: target treatment cohort
#' @param pairs data.frame with column \code{gp}
#' @param bs_df integer: B-spline df, or \code{"ic"}
#' @param K_folds integer: number of cross-fitting folds
#' @param fold_id integer vector length n: pre-generated fold assignments
#' @param ratio_method \code{"direct"} (default at the function level; legacy
#'   per-pair LS sieve for every comparison) or \code{"exp"} (per-target
#'   exponential-link Riesz regressions for every comparison, full
#'   estimation-effect aux; \code{edid()}'s default). The function-level default
#'   stays \code{"direct"} so existing direct callers and validation harnesses are
#'   unchanged; \code{fit_edid_cells} passes the user's choice explicitly.
#'
#' @return named list of n-vectors, keyed by \code{as.character(gp)}
#' @keywords internal
estimate_all_propensity_ratios <- function(panel_obj, g, pairs, bs_df,
                                           K_folds, fold_id, return_aux = FALSE,
                                           ratio_method = c("direct", "exp")) {
  ratio_method <- match.arg(ratio_method)
  n       <- panel_obj$n
  X_mat   <- panel_obj$covariate_matrix
  G_vec   <- panel_obj$unit_cohorts
  result  <- list()
  aux     <- list()
  ic_mode <- identical(bs_df, "ic")   # per-fit IC selection: record the selected dfs
  sel_key <- character(0L); sel_df <- integer(0L)

  unique_gps <- unique(pairs$gp)

  # Per-target fitter: the paper's LS sieve ("direct"), or (ratio_method = "exp") the
  # exponential-link Riesz regression -- identical signature/contract, so the fold loop and
  # the aux/ic bookkeeping are shared verbatim. Both fit EVERY gp (including the
  # never-treated pool) as an independent per-target regression.
  ratio_fitter <- if (identical(ratio_method, "exp")) estimate_propensity_ratio_exp_edid
                  else estimate_propensity_ratio_edid

  for (gp in unique_gps) {
    r_full <- numeric(n)

    for (ell in seq_len(K_folds)) {
      if (K_folds == 1L) {
        # Plug-in: train = test = full sample (paper's main-text proposal)
        test_idx  <- seq_len(n)
        train_idx <- seq_len(n)
      } else {
        test_idx  <- which(fold_id == ell)
        train_idx <- which(fold_id != ell)
      }
      if (length(test_idx) == 0L) next

      out_gp <- ratio_fitter(
        X_train = X_mat[train_idx, , drop = FALSE],
        G_train = G_vec[train_idx],
        X_test  = X_mat[test_idx,  , drop = FALSE],
        g       = g,
        gp      = gp,
        bs_df   = bs_df,
        return_aux = return_aux
      )
      if (return_aux) {
        r_full[test_idx] <- out_gp$pred              # aux only requested with K=1 (single fold)
        aux[[as.character(gp)]] <- out_gp
      } else {
        r_full[test_idx] <- out_gp
      }
      if (ic_mode) {                                 # selected df (per fit; K = 1 in edid())
        dfk <- if (return_aux) out_gp$bs_df else attr(out_gp, "edid_bs_df")
        if (!is.null(dfk)) { sel_key <- c(sel_key, as.character(gp)); sel_df <- c(sel_df, dfk) }
      }
    }

    .check_extreme_ratios_edid(r_full, g, gp)
    result[[as.character(gp)]] <- r_full
  }

  out <- if (return_aux) list(predictions = result, aux = aux) else result
  if (ic_mode && length(sel_key)) {
    attr(out, "bs_df_selected") <- data.frame(key = sel_key, bs_df = sel_df,
                                              stringsAsFactors = FALSE)
  }
  out
}

#' Estimate inverse propensities 1/p_g'(X) for all groups via cross-fitting
#'
#' For each group g needed in Omega* (the target group g, the never-treated,
#' and each comparison cohort g'), performs K-fold cross-fitting to produce
#' a full-sample n-vector of \eqn{\hat s_{g'}(X_i) = 1/\hat p_{g'}(X_i)}.
#'
#' \strong{Inverse-propensity construction (\code{ratio_method}).} The never-treated
#' inverse propensity \eqn{1/p_{NT}} is ALWAYS the LS sieve of
#' \code{estimate_inverse_propensity_edid} (its Gram uses the large never-treated
#' pool). The paper's per-cohort LS sieve (\code{"direct"}) minimizes
#' \eqn{\mathbb{E}_n[s^2 G_{g'} - 2 s]}, whose Gram uses only the \eqn{n_{g'}}
#' cohort observations: for small cohorts it is the s-channel instance of the same
#' thin-denominator explosion as the direct ratio fits (audited fitted \eqn{1/p}
#' of order \eqn{10^8} against a true scale of \eqn{10^2}, ~half the sample clamped
#' at 0), which poisons every Omega* prefactor it enters. Under \code{"exp"} (the
#' default), each FINITE cohort's inverse propensity is instead the per-target
#' exponential-link Riesz regression of \code{estimate_inverse_propensity_exp_edid}
#' (tailored loss \eqn{\mathbb{E}_n[e^{\psi'\beta} G_{g'} - \psi'\beta]}, strictly
#' positive fit), and the aux entries carry FULL M-estimator pieces with
#' \code{is_fallback = FALSE}, so the analytic inv-p weight-channel correction covers
#' every channel (no skipping). PT-Post fits are byte-invariant to this choice:
#' their single moment uses no \eqn{s}, their weight is identically 1
#' (\eqn{H = 1}), and the \eqn{H = 1} weight-channel coupling is exactly zero.
#'
#' @param panel_obj panel object
#' @param g scalar: target treatment cohort
#' @param pairs data.frame with column \code{gp}
#' @param bs_df integer: B-spline df
#' @param K_folds integer: number of cross-fitting folds
#' @param fold_id integer vector length n: pre-generated fold assignments
#' @param ratio_method \code{"direct"} (default at the function level; legacy LS sieve
#'   for every group) or \code{"exp"} (exponential-link Riesz regression for finite
#'   cohorts, full estimation-effect aux; \code{edid()}'s default). The never-treated
#'   \eqn{1/p_{NT}} is ALWAYS the LS sieve.
#'
#' @return named list of n-vectors, keyed by \code{as.character(group)}
#' @keywords internal
estimate_all_inverse_propensities <- function(panel_obj, g, pairs, bs_df,
                                              K_folds, fold_id, return_aux = FALSE,
                                              ratio_method = c("direct", "exp")) {
  ratio_method <- match.arg(ratio_method)
  n       <- panel_obj$n
  X_mat   <- panel_obj$covariate_matrix
  G_vec   <- panel_obj$unit_cohorts
  result  <- list()
  aux     <- if (return_aux) list() else NULL       # per-group M-estimator pieces for the inv_p Sigma_Omega channel
  ic_mode <- identical(bs_df, "ic")                 # per-fit IC selection: record the selected dfs
  sel_key <- character(0L); sel_df <- integer(0L)

  groups_needed <- unique(c(g, Inf, pairs$gp[is.finite(pairs$gp)]))

  for (gp in groups_needed) {
    s_full <- numeric(n)
    want_aux_gp <- return_aux && K_folds == 1L       # aux only in the plug-in (train = test) regime, like ACH
    # Per-target fitter: "exp" routes FINITE cohorts to the exponential-link Riesz regression
    # (full aux); the never-treated 1/p_NT always keeps the paper's LS sieve.
    invp_fitter <- if (identical(ratio_method, "exp") && is.finite(gp)) estimate_inverse_propensity_exp_edid
                   else estimate_inverse_propensity_edid

    for (ell in seq_len(K_folds)) {
      if (K_folds == 1L) {
        # Plug-in: train = test = full sample (paper's main-text proposal)
        test_idx  <- seq_len(n)
        train_idx <- seq_len(n)
      } else {
        test_idx  <- which(fold_id == ell)
        train_idx <- which(fold_id != ell)
      }
      if (length(test_idx) == 0L) next

      res_ell <- invp_fitter(
        X_train = X_mat[train_idx, , drop = FALSE],
        G_train = G_vec[train_idx],
        X_test  = X_mat[test_idx,  , drop = FALSE],
        gp      = gp,
        bs_df   = bs_df,
        return_aux = want_aux_gp
      )
      if (want_aux_gp) { s_full[test_idx] <- res_ell$s_hat; aux[[as.character(gp)]] <- res_ell }
      else             s_full[test_idx] <- res_ell
      if (ic_mode) {                                 # selected df (per fit; K = 1 in edid())
        dfk <- if (want_aux_gp) res_ell$bs_df else attr(res_ell, "edid_bs_df")
        if (!is.null(dfk)) { sel_key <- c(sel_key, as.character(gp)); sel_df <- c(sel_df, dfk) }
      }
    }

    result[[as.character(gp)]] <- s_full
  }

  if (return_aux) attr(result, "aux") <- aux
  if (ic_mode && length(sel_key)) {
    attr(result, "bs_df_selected") <- data.frame(key = sel_key, bs_df = sel_df,
                                                 stringsAsFactors = FALSE)
  }
  result
}


#' Estimate conditional means for all (g', period) combinations via cross-fitting
#'
#' For each unique (gp, period) pair needed by the cell, performs K-fold
#' cross-fitting to produce a full-sample n-vector of
#' \eqn{\hat m_{g', \text{period}, 1}(X_i)}.
#'
#' @param panel_obj panel object with \code{covariate_matrix}, \code{unit_cohorts},
#'   \code{outcome_wide}, and \code{period_to_col}
#' @param pairs data.frame with columns \code{gp} and \code{tpre}
#' @param t_val scalar: target time period for this cell
#' @param bs_df integer: B-spline df
#' @param K_folds integer: number of cross-fitting folds
#' @param fold_id integer vector length n: pre-generated fold assignments
#'
#' @return named list of n-vectors, keyed by \code{paste0(gp, "_", period)}
#' @keywords internal
estimate_all_conditional_means <- function(panel_obj, pairs, t_val, bs_df,
                                           K_folds, fold_id, return_aux = FALSE) {
  n           <- panel_obj$n
  X_mat       <- panel_obj$covariate_matrix
  G_vec       <- panel_obj$unit_cohorts
  ow          <- panel_obj$outcome_wide
  t1_col      <- panel_obj$period_to_col[[as.character(panel_obj$period_1)]]
  result      <- list()
  aux         <- list()
  ic_mode     <- identical(bs_df, "ic")   # per-fit IC selection: record the selected dfs
  sel_key     <- character(0L); sel_df <- integer(0L)

  # Collect unique (gp, period) combinations needed
  # We need m_{gp, t_val, 1}(X) and m_{gp, tpre, 1}(X) for each pair
  unique_gps  <- unique(pairs$gp)
  unique_tpre <- unique(pairs$tpre)

  # Build full list of (gp, period) to estimate
  combos <- unique(rbind(
    data.frame(gp = pairs$gp, period = t_val),
    data.frame(gp = pairs$gp, period = pairs$tpre)
  ))

  for (ii in seq_len(nrow(combos))) {
    gp     <- combos$gp[ii]
    period <- combos$period[ii]
    key    <- paste0(gp, "_", period)

    if (!is.null(result[[key]])) next  # already computed

    period_col <- panel_obj$period_to_col[[as.character(period)]]
    if (is.null(period_col)) {
      warning(sprintf("estimate_all_conditional_means: period %g not found in panel.", period))
      result[[key]] <- rep(NA_real_, n)
      next
    }

    Y_delta <- ow[, period_col] - ow[, t1_col]  # Y_{period} - Y_1
    m_full  <- numeric(n)

    for (ell in seq_len(K_folds)) {
      if (K_folds == 1L) {
        # Plug-in: train = test = full sample (paper's main-text proposal)
        test_idx  <- seq_len(n)
        train_idx <- seq_len(n)
      } else {
        test_idx  <- which(fold_id == ell)
        train_idx <- which(fold_id != ell)
      }
      if (length(test_idx) == 0L) next

      out_m <- estimate_conditional_mean_edid(
        X_train       = X_mat[train_idx, , drop = FALSE],
        Y_delta_train = Y_delta[train_idx],
        G_train       = G_vec[train_idx],
        X_test        = X_mat[test_idx, , drop = FALSE],
        gp            = gp,
        bs_df         = bs_df,
        return_aux    = return_aux
      )
      if (return_aux) {
        m_full[test_idx] <- out_m$pred               # aux only requested with K=1 (single fold)
        aux[[key]] <- out_m
      } else {
        m_full[test_idx] <- out_m
      }
      if (ic_mode) {                                 # selected df (per fit; K = 1 in edid())
        dfk <- if (return_aux) out_m$bs_df else attr(out_m, "edid_bs_df")
        if (!is.null(dfk)) { sel_key <- c(sel_key, key); sel_df <- c(sel_df, dfk) }
      }
    }

    result[[key]] <- m_full
  }

  out <- if (return_aux) list(predictions = result, aux = aux) else result
  if (ic_mode && length(sel_key)) {
    attr(out, "bs_df_selected") <- data.frame(key = sel_key, bs_df = sel_df,
                                              stringsAsFactors = FALSE)
  }
  out
}
