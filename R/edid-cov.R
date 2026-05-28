# edid-cov.R
# Nuisance estimation functions for the EDiD covariate path.
# Implements sieve (B-spline) estimation of propensity ratios and conditional
# means with K-fold cross-fitting.

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
  sample(seq_len(K), n, replace = TRUE)
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
#' Then \eqn{\hat r(X_i) = B(X_i)' \hat\beta}, clipped to [0, Inf).
#'
#' @param X_train numeric matrix n_train x d
#' @param G_train numeric vector n_train: cohort values (Inf for never-treated)
#' @param X_test  numeric matrix n_test x d
#' @param g  scalar: target treatment cohort
#' @param gp scalar: comparison cohort (may be Inf for never-treated)
#' @param bs_df integer: B-spline degrees of freedom (default 4)
#'
#' @return numeric vector length n_test: estimated r(X) values, >= 0
#' @keywords internal
estimate_propensity_ratio_edid <- function(X_train, G_train, X_test, g, gp,
                                           bs_df = 4L) {
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
    return(rep(0, n_test))
  }
  if (n_g < 1L) {
    warning(sprintf(
      "estimate_propensity_ratio_edid: 0 units in g=%g training fold; returning 0.", g
    ))
    return(rep(0, n_test))
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

  # Predict on test data
  B_test <- predict_basis_edid(attr(B_train_obj, "bs_objects"), X_test)
  r_hat  <- pmax(drop(B_test %*% beta_hat), 0)

  r_hat
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
#' @param bs_df integer: B-spline degrees of freedom (default 4)
#'
#' @return numeric vector length n_test: estimated 1/p_{g'}(X) values, >= 0
#' @keywords internal
estimate_inverse_propensity_edid <- function(X_train, G_train, X_test, gp,
                                             bs_df = 4L) {
  n_test <- nrow(X_test)

  mask_gp <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  n_gp <- sum(mask_gp)

  if (n_gp < 2L) {
    warning(sprintf(
      "estimate_inverse_propensity_edid: fewer than 2 units in g'=%g; returning 1/pi_g'.", gp
    ))
    return(rep(length(G_train) / n_gp, n_test))
  }

  B_train_obj <- build_basis_matrix_edid(X_train, bs_df)
  B_train     <- unclass(B_train_obj)
  attr(B_train, "bs_objects") <- NULL

  B_gp       <- B_train[mask_gp, , drop = FALSE]
  col_sums_all <- colSums(B_train)

  BtB_gp   <- t(B_gp) %*% B_gp
  beta_hat <- as.vector(compute_pseudoinverse_edid(BtB_gp) %*% col_sums_all)

  B_test <- predict_basis_edid(attr(B_train_obj, "bs_objects"), X_test)
  s_hat  <- pmax(drop(B_test %*% beta_hat), 0)

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
#' @param bs_df integer: B-spline degrees of freedom (default 4)
#'
#' @return numeric vector length n_test
#' @keywords internal
estimate_conditional_mean_edid <- function(X_train, Y_delta_train, G_train,
                                           X_test, gp, bs_df = 4L) {
  n_test  <- nrow(X_test)

  mask_gp <- if (is.infinite(gp)) is.infinite(G_train) else (G_train == gp)
  n_gp    <- sum(mask_gp)

  if (n_gp < 2L) {
    fallback_val <- if (n_gp == 1L) Y_delta_train[mask_gp] else 0
    warning(sprintf(
      "estimate_conditional_mean_edid: fewer than 2 units in g'=%g training fold; using constant.", gp
    ))
    return(rep(fallback_val, n_test))
  }

  X_gp    <- X_train[mask_gp, , drop = FALSE]
  y_gp    <- Y_delta_train[mask_gp]

  B_gp_train_obj <- build_basis_matrix_edid(X_gp, bs_df)
  B_gp_train     <- unclass(B_gp_train_obj)
  attr(B_gp_train, "bs_objects") <- NULL

  # Check minimum sample for basis dimension
  p_basis <- ncol(B_gp_train)
  if (n_gp <= p_basis) {
    # Fewer obs than basis columns: use simpler 1-column basis (intercept only)
    fit <- list(coef = mean(y_gp))
    m_hat <- rep(mean(y_gp), n_test)
    return(m_hat)
  }

  fit   <- solve_ols_edid(B_gp_train, y_gp)
  B_test <- predict_basis_edid(attr(B_gp_train_obj, "bs_objects"), X_test)
  m_hat <- drop(B_test %*% fit$coef)

  m_hat
}

# ---------------------------------------------------------------------------
# Cross-fitted nuisance estimation (aggregate over folds)
# ---------------------------------------------------------------------------

#' Estimate propensity ratios for all comparison cohorts via cross-fitting
#'
#' For each unique \code{gp} in \code{pairs}, performs K-fold cross-fitting to
#' produce a full-sample n-vector of \eqn{\hat r_{g, g'}(X_i)}.
#'
#' @param panel_obj panel object with \code{covariate_matrix} and
#'   \code{unit_cohorts}
#' @param g scalar: target treatment cohort
#' @param pairs data.frame with column \code{gp}
#' @param bs_df integer: B-spline df
#' @param K_folds integer: number of cross-fitting folds
#' @param fold_id integer vector length n: pre-generated fold assignments
#'
#' @return named list of n-vectors, keyed by \code{as.character(gp)}
#' @keywords internal
estimate_all_propensity_ratios <- function(panel_obj, g, pairs, bs_df,
                                           K_folds, fold_id) {
  n       <- panel_obj$n
  X_mat   <- panel_obj$covariate_matrix
  G_vec   <- panel_obj$unit_cohorts
  result  <- list()

  unique_gps <- unique(pairs$gp)

  for (gp in unique_gps) {
    r_full <- numeric(n)

    for (ell in seq_len(K_folds)) {
      test_idx  <- which(fold_id == ell)
      train_idx <- which(fold_id != ell)
      if (length(test_idx) == 0L) next

      r_full[test_idx] <- estimate_propensity_ratio_edid(
        X_train = X_mat[train_idx, , drop = FALSE],
        G_train = G_vec[train_idx],
        X_test  = X_mat[test_idx,  , drop = FALSE],
        g       = g,
        gp      = gp,
        bs_df   = bs_df
      )
    }

    .check_extreme_ratios_edid(r_full, g, gp)
    result[[as.character(gp)]] <- r_full
  }

  result
}

#' Estimate inverse propensities 1/p_{g'}(X) for all groups via cross-fitting
#'
#' For each group g needed in Omega* (the target group g, the never-treated,
#' and each comparison cohort g'), performs K-fold cross-fitting to produce
#' a full-sample n-vector of \eqn{\hat s_{g'}(X_i) = 1/\hat p_{g'}(X_i)}.
#'
#' @param panel_obj panel object
#' @param g scalar: target treatment cohort
#' @param pairs data.frame with column \code{gp}
#' @param bs_df integer: B-spline df
#' @param K_folds integer: number of cross-fitting folds
#' @param fold_id integer vector length n: pre-generated fold assignments
#'
#' @return named list of n-vectors, keyed by \code{as.character(group)}
#' @keywords internal
estimate_all_inverse_propensities <- function(panel_obj, g, pairs, bs_df,
                                              K_folds, fold_id) {
  n       <- panel_obj$n
  X_mat   <- panel_obj$covariate_matrix
  G_vec   <- panel_obj$unit_cohorts
  result  <- list()

  groups_needed <- unique(c(g, Inf, pairs$gp[is.finite(pairs$gp)]))

  for (gp in groups_needed) {
    s_full <- numeric(n)

    for (ell in seq_len(K_folds)) {
      test_idx  <- which(fold_id == ell)
      train_idx <- which(fold_id != ell)
      if (length(test_idx) == 0L) next

      s_full[test_idx] <- estimate_inverse_propensity_edid(
        X_train = X_mat[train_idx, , drop = FALSE],
        G_train = G_vec[train_idx],
        X_test  = X_mat[test_idx,  , drop = FALSE],
        gp      = gp,
        bs_df   = bs_df
      )
    }

    result[[as.character(gp)]] <- s_full
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
                                           K_folds, fold_id) {
  n           <- panel_obj$n
  X_mat       <- panel_obj$covariate_matrix
  G_vec       <- panel_obj$unit_cohorts
  ow          <- panel_obj$outcome_wide
  t1_col      <- panel_obj$period_to_col[[as.character(panel_obj$period_1)]]
  result      <- list()

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
      test_idx  <- which(fold_id == ell)
      train_idx <- which(fold_id != ell)
      if (length(test_idx) == 0L) next

      m_full[test_idx] <- estimate_conditional_mean_edid(
        X_train       = X_mat[train_idx, , drop = FALSE],
        Y_delta_train = Y_delta[train_idx],
        G_train       = G_vec[train_idx],
        X_test        = X_mat[test_idx, , drop = FALSE],
        gp            = gp,
        bs_df         = bs_df
      )
    }

    result[[key]] <- m_full
  }

  result
}
