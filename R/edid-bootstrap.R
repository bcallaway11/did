# edid-bootstrap.R
# Multiplier bootstrap for the EDiD estimator.

#' Run the multiplier bootstrap for EDiD estimates
#'
#' Generates \code{n_bootstrap} perturbed versions of all cell-level ATTs
#' by multiplying stored EIF vectors with multiplier weights, then re-aggregates
#' using the same fixed cohort-share weights.
#'
#' @param cells list of \code{edid_cell_result} objects
#' @param eif_matrix n x n_cells numeric matrix of stored EIFs
#' @param cell_index data.frame with columns \code{group}, \code{time},
#'   \code{cell_id}, \code{is_pre}
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param n_bootstrap positive integer number of bootstrap draws
#' @param bootstrap_weights character: \code{"rademacher"}, \code{"mammen"},
#'   or \code{"webb"}
#' @param seed integer seed or NULL
#' @param aggregate character: which aggregations to return
#' @param balance_e integer or NULL
#' @param alpha significance level
#'
#' @return list with elements \code{overall_b}, \code{event_study_b},
#'   \code{group_b}, \code{n_bootstrap}, \code{weight_type}, \code{seed}
#' @keywords internal
run_multiplier_bootstrap_edid <- function(
  cells, eif_matrix, cell_index, panel_obj,
  n_bootstrap, bootstrap_weights = "rademacher", seed = NULL,
  aggregate = "all", balance_e = NULL, alpha = 0.05
) {
  n <- panel_obj$n

  # Generate multiplier weights: n x n_bootstrap matrix
  xi_mat <- generate_multiplier_weights_edid(
    n         = n,
    n_bootstrap = n_bootstrap,
    type        = bootstrap_weights,
    cluster_indices = panel_obj$cluster_indices,
    seed        = seed
  )

  # Post-treatment cell ids and their ATTs
  post_ci  <- cell_index[!cell_index$is_pre, , drop = FALSE]
  post_ids <- post_ci$cell_id
  valid_post <- vapply(post_ids, function(cid) {
    !is.null(cells[[cid]]$att) && is.finite(cells[[cid]]$att)
  }, logical(1L))
  post_ci_v  <- post_ci[valid_post, , drop = FALSE]
  post_ids_v <- post_ci_v$cell_id

  # Perturbed cell ATTs: ATT_b(g,t) = ATT_hat(g,t) + (1/n) * xi' * EIF
  # Build perturbed ATT matrix: n_valid_post x n_bootstrap
  att_hat_post <- vapply(post_ids_v, function(cid) cells[[cid]]$att, numeric(1L))
  eif_post     <- eif_matrix[, post_ids_v, drop = FALSE]  # n x n_valid_post
  # perturbation: (1/n) * t(xi_mat) %*% eif_post -> n_bootstrap x n_valid_post
  perturb_mat  <- (1 / n) * t(xi_mat) %*% eif_post
  att_boot_mat <- sweep(perturb_mat, 2, att_hat_post, `+`)
  # att_boot_mat: n_bootstrap x n_valid_post

  # Cohort-share weights for overall
  pi_g_k <- vapply(seq_along(post_ids_v), function(jj) {
    panel_obj$cohort_fractions[[as.character(post_ci_v$group[jj])]]
  }, numeric(1L))
  S_u    <- sum(pi_g_k)
  q_norm <- if (S_u > EDID_DENOM_EPS) pi_g_k / S_u else rep(1 / length(pi_g_k), length(pi_g_k))

  # -----------------------------------------------------------------------
  # Overall bootstrap draws
  # -----------------------------------------------------------------------
  overall_b <- NULL
  if (aggregate %in% c("all", "overall")) {
    overall_b <- drop(att_boot_mat %*% q_norm)  # n_bootstrap vector
  }

  # -----------------------------------------------------------------------
  # Event-study bootstrap draws
  # -----------------------------------------------------------------------
  event_study_b <- NULL
  if (aggregate %in% c("all", "event_study")) {
    cell_index_v <- cell_index[cell_index$cell_id %in% post_ids_v, , drop = FALSE]
    cell_index_v$e <- cell_index_v$time - cell_index_v$group
    if (!is.null(balance_e)) {
      cell_index_v <- cell_index_v[abs(cell_index_v$e) <= balance_e, , drop = FALSE]
    }
    unique_e <- sort(unique(cell_index_v$e))
    event_study_b <- vector("list", length(unique_e))
    names(event_study_b) <- as.character(unique_e)

    for (ii in seq_along(unique_e)) {
      e_val   <- unique_e[ii]
      e_mask  <- cell_index_v$e == e_val
      e_ids   <- cell_index_v$cell_id[e_mask]
      e_g     <- cell_index_v$group[e_mask]
      # map e_ids to column indices in att_boot_mat
      col_idx <- match(e_ids, post_ids_v)
      col_idx <- col_idx[!is.na(col_idx)]
      if (length(col_idx) == 0L) next
      e_pi    <- pi_g_k[col_idx]
      e_Su    <- sum(e_pi)
      if (e_Su < EDID_DENOM_EPS) next
      e_q     <- e_pi / e_Su
      event_study_b[[ii]] <- drop(att_boot_mat[, col_idx, drop = FALSE] %*% e_q)
    }
  }

  # -----------------------------------------------------------------------
  # Group bootstrap draws
  # -----------------------------------------------------------------------
  group_b <- NULL
  if (aggregate %in% c("all", "group")) {
    tgroups <- panel_obj$treatment_groups
    group_b <- vector("list", length(tgroups))
    names(group_b) <- as.character(tgroups)

    for (ii in seq_along(tgroups)) {
      g_val   <- tgroups[ii]
      g_mask  <- post_ci_v$group == g_val
      g_ids   <- post_ci_v$cell_id[g_mask]
      col_idx <- match(g_ids, post_ids_v)
      col_idx <- col_idx[!is.na(col_idx)]
      if (length(col_idx) == 0L) next
      m_g     <- length(col_idx)
      g_q     <- rep(1 / m_g, m_g)
      group_b[[ii]] <- drop(att_boot_mat[, col_idx, drop = FALSE] %*% g_q)
    }
  }

  list(
    overall_b       = overall_b,
    event_study_b   = event_study_b,
    group_b         = group_b,
    n_bootstrap     = n_bootstrap,
    weight_type     = bootstrap_weights,
    seed            = seed
  )
}

#' Generate multiplier bootstrap weights
#'
#' Returns an \code{n x n_bootstrap} matrix of multiplier weights drawn from
#' the specified distribution.  When \code{cluster_indices} is supplied,
#' weights are drawn at the cluster level (G x n_bootstrap) and then
#' expanded to unit level by repeating within cluster.
#'
#' @param n integer: number of units
#' @param n_bootstrap positive integer: number of bootstrap draws
#' @param type character: \code{"rademacher"} (default), \code{"mammen"},
#'   or \code{"webb"}
#' @param cluster_indices integer vector length n (values 1..G) or NULL
#' @param seed integer seed or NULL
#'
#' @return numeric matrix n x n_bootstrap
#' @keywords internal
generate_multiplier_weights_edid <- function(
  n, n_bootstrap, type = "rademacher", cluster_indices = NULL, seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Determine draw size
  if (!is.null(cluster_indices)) {
    G      <- length(unique(cluster_indices))
    draw_n <- G
  } else {
    G      <- NULL
    draw_n <- n
  }

  xi_raw <- switch(type,
    rademacher = {
      matrix(sample(c(-1, 1), draw_n * n_bootstrap, replace = TRUE),
             nrow = draw_n, ncol = n_bootstrap)
    },
    mammen = {
      p_pos <- (sqrt(5) + 1) / (2 * sqrt(5))
      p_neg <- 1 - p_pos
      v_pos <- (sqrt(5) + 1) / 2
      v_neg <- -(sqrt(5) - 1) / 2
      raw   <- matrix(
        sample(c(v_neg, v_pos), draw_n * n_bootstrap, replace = TRUE,
               prob = c(p_neg, p_pos)),
        nrow = draw_n, ncol = n_bootstrap
      )
      raw
    },
    webb = {
      webb_vals <- c(-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2))
      matrix(
        sample(webb_vals, draw_n * n_bootstrap, replace = TRUE),
        nrow = draw_n, ncol = n_bootstrap
      )
    },
    stop(sprintf("Unknown bootstrap_weights type: \"%s\". ",
                 type),
         "Choose one of \"rademacher\", \"mammen\", or \"webb\".")
  )

  # Expand cluster-level draws to unit level
  if (!is.null(cluster_indices)) {
    xi_unit <- xi_raw[cluster_indices, , drop = FALSE]
    return(xi_unit)
  }
  xi_raw
}

#' Compute bootstrap SE, CI, and p-value from a vector of bootstrap draws
#'
#' @param boot_draws numeric vector of length \code{n_bootstrap}
#' @param att_hat scalar point estimate
#' @param alpha significance level in (0, 1)
#'
#' @return named list: \code{se_boot}, \code{ci_lower}, \code{ci_upper},
#'   \code{p_value_boot}
#' @keywords internal
compute_bootstrap_stats_edid <- function(boot_draws, att_hat, alpha = 0.05) {
  se_boot   <- stats::sd(boot_draws - att_hat)
  ci_lower  <- unname(stats::quantile(boot_draws, alpha / 2))
  ci_upper  <- unname(stats::quantile(boot_draws, 1 - alpha / 2))
  p_value   <- mean(abs(boot_draws - att_hat) >= abs(att_hat))
  list(
    se_boot    = se_boot,
    ci_lower   = ci_lower,
    ci_upper   = ci_upper,
    p_value_boot = p_value
  )
}
