# edid-aggregate.R
# Aggregation functions for the EDiD estimator:
#   aggregate_overall_edid()
#   aggregate_event_study_edid()
#   aggregate_group_edid()
#   compute_wif_contribution_edid()

# ---------------------------------------------------------------------------
# Internal helper: pull EIF column from eif_matrix for a given cell_id
# ---------------------------------------------------------------------------
.eif_col <- function(eif_matrix, cell_id) {
  if (is.null(eif_matrix)) return(NULL)
  eif_matrix[, cell_id, drop = TRUE]
}

# ---------------------------------------------------------------------------
# WIF correction
# ---------------------------------------------------------------------------

#' Compute WIF (weight influence function) correction for aggregated EIF
#'
#' Corrects the aggregated EIF for estimation error in the cohort-share weights
#' \eqn{\pi_g = n_g / n}. The correction accounts for the fact that \eqn{\pi_g}
#' is estimated from data and therefore contributes to the variance of the
#' aggregated estimator.
#'
#' @param weight_fn function(cells, cell_index, panel_obj) that returns a named
#'   numeric vector of normalized weights \eqn{q_k} for all post-treatment cells;
#'   names are cell_ids (character)
#' @param cells list of \code{edid_cell_result} objects
#' @param eif_matrix n x n_cells numeric matrix (or NULL)
#' @param cell_index data.frame with columns \code{group}, \code{time},
#'   \code{cell_id}, \code{is_pre}
#' @param panel_obj panel object from \code{prepare_edid_panel()}
#' @param agg_att scalar: the already-computed aggregated ATT (used for WIF)
#'
#' @return numeric vector length n: WIF correction to add to aggregated EIF
#' @keywords internal
compute_wif_contribution_edid <- function(
  weight_fn, cells, eif_matrix, cell_index, panel_obj, agg_att
) {
  n   <- panel_obj$n
  wif <- numeric(n)

  if (is.null(eif_matrix)) return(wif)

  # Obtain normalized weights for all post-treatment cells
  q_vec <- weight_fn(cells, cell_index, panel_obj)  # named numeric

  # Get post-treatment cell ids and their ATTs
  post_ci   <- cell_index[!cell_index$is_pre, , drop = FALSE]
  post_ids  <- as.character(post_ci$cell_id)
  post_atts <- vapply(post_ids, function(cid) {
    idx <- which(post_ci$cell_id == as.integer(cid))
    cells[[post_ci$cell_id[idx]]]$att
  }, numeric(1L))

  # For each cohort g, compute the WIF contribution from pi_g uncertainty.
  # WIF_i = sum_k [ (att_k - agg_att) * d(q_k)/d(pi_g_k) * d(pi_g_k)/dI(G_i=g_k) ]
  # With q_k = pi_{g_k} / sum_{k'} pi_{g_{k'}} :
  # d(q_k)/d(pi_g) = (I(g_k==g) * S - pi_{g_k}) / S^2
  # where S = sum_{k' in post} pi_{g_{k'}}
  # d(pi_g)/dI(G_i=g) = (1/n) * (1 - I(G_i=g) * n_g_hat / n_g) ... simplified below

  S <- sum(q_vec[post_ids] * (q_vec[post_ids] > 0))  # sum of unnormalized weights
  # Actually q_vec already gives normalized weights. Rebuild unnorm (pi_g_k):
  # unnorm_k = pi_{g_k}; norm_k = pi_{g_k} / S_unnorm
  # We need: d(q_k_norm)/d(pi_{g}) = (I(g_k==g)*S_u - pi_{g_k}) / S_u^2
  # where S_u = sum pi_{g_{k'}}.

  # Collect pi_{g_k} for each post cell
  pi_by_cell <- vapply(post_ids, function(cid) {
    g_val <- post_ci$group[post_ci$cell_id == as.integer(cid)]
    panel_obj$cohort_fractions[[as.character(g_val)]]
  }, numeric(1L))
  S_u <- sum(pi_by_cell)  # sum of unnormalized (cohort share) weights

  if (S_u < EDID_DENOM_EPS) return(wif)

  for (g in panel_obj$treatment_groups) {
    mask_g <- panel_obj$cohort_masks[[as.character(g)]]
    n_g    <- sum(mask_g)
    if (n_g == 0L) next
    pi_g   <- panel_obj$cohort_fractions[[as.character(g)]]

    # Cells belonging to cohort g (among post cells)
    g_post_mask <- post_ci$group == g
    if (!any(g_post_mask)) next

    # sum_{k: g_k=g} (att_k - agg_att) * (S_u - pi_g * K_g) / S_u^2
    # where K_g = number of post cells for cohort g
    # More precisely, the partial derivative for cohort g:
    # dq_k/d(pi_g) = (I(g_k==g)*S_u - pi_g * sum_all 1) / S_u^2
    # But S_u = sum_k pi_{g_k}, so d(S_u)/d(pi_g) = K_g (# cells with g_k=g)
    K_g <- sum(g_post_mask)
    g_ids <- post_ids[g_post_mask]
    g_atts <- post_atts[g_post_mask]

    # For each post cell k:
    #   if g_k == g: dq_k/d(pi_g) = (S_u - pi_g * K_g) / S_u^2 + 0 ... simplified:
    #   dq_k/d(pi_g) = (I(g_k==g) * S_u - pi_g * K_g) / S_u^2
    #   if g_k != g: dq_k/d(pi_g) = -pi_{g_k} * K_g / S_u^2
    # Then d(agg_att)/d(pi_g) = sum_k att_k * dq_k/d(pi_g)
    # And WIF contribution for unit i in cohort g:
    #   d(agg_att)/d(pi_g) * d(pi_g)/d(I(G_i=g)) = d(agg_att)/d(pi_g) * (1 - pi_g * K_g / K_g?) ...
    # d(pi_g)/d(I(G_i=g)) = 1/n  (adding one unit to cohort g increases pi_g by 1/n)

    # Contribution from cells in cohort g (g_k == g):
    d_agg_d_pig_from_g  <- sum(g_atts) * (S_u - pi_g * K_g) / (S_u^2)
    # Contribution from cells NOT in cohort g (g_k != g):
    not_g_mask <- !g_post_mask
    if (any(not_g_mask)) {
      not_g_atts    <- post_atts[not_g_mask]
      not_g_pi      <- pi_by_cell[not_g_mask]
      d_agg_d_pig_from_notg <- sum(not_g_atts * (-not_g_pi * K_g) / (S_u^2))
    } else {
      d_agg_d_pig_from_notg <- 0
    }
    d_agg_d_pig <- d_agg_d_pig_from_g + d_agg_d_pig_from_notg

    # Units in cohort g: d(pi_g)/dI(G_i=g) = 1/n
    wif[mask_g] <- wif[mask_g] + d_agg_d_pig * (1 / n)

    # Units outside cohort g also have a contribution via d(pi_g)/dI(G_i=g):
    # d(pi_g)/dI(G_i=g) for i not in g = -pi_g / (n * (1 - pi_g))... but
    # using the simpler sample-level: d(n_g/n)/dI(G_i not g) = 0 exactly.
    # (Adding unit i outside cohort g does not change n_g.) So no contribution
    # for units outside cohort g from this cohort's weight derivative. Correct.
  }

  wif
}

# ---------------------------------------------------------------------------
# Overall ATT aggregation
# ---------------------------------------------------------------------------

#' Aggregate cell-level ATTs into an overall ATT
#'
#' Uses cohort-share weights \eqn{q_k = \pi_{g_k}} over post-treatment cells,
#' normalized to sum to 1.  Includes WIF correction for estimated weights.
#'
#' @param cells list of \code{edid_cell_result} objects (ordered by cell_id)
#' @param eif_matrix n x n_cells numeric matrix (or NULL if EIF not stored)
#' @param cell_index data.frame with columns \code{group}, \code{time},
#'   \code{cell_id}, \code{is_pre}
#' @param panel_obj panel object
#' @param alpha significance level
#'
#' @return named list: \code{att}, \code{se}, \code{ci_lower}, \code{ci_upper},
#'   \code{t_stat}, \code{p_value}, \code{eif_agg}
#' @keywords internal
aggregate_overall_edid <- function(cells, eif_matrix, cell_index, panel_obj, alpha) {

  # Post-treatment cells only
  post_ci  <- cell_index[!cell_index$is_pre, , drop = FALSE]
  post_ids <- post_ci$cell_id  # integer indices into cells list

  # Filter to cells with valid (non-NA) ATT
  valid_mask <- vapply(post_ids, function(cid) {
    !is.null(cells[[cid]]$att) && is.finite(cells[[cid]]$att)
  }, logical(1L))

  if (!any(valid_mask)) {
    return(list(att = NA_real_, se = NA_real_, ci_lower = NA_real_,
                ci_upper = NA_real_, t_stat = NA_real_, p_value = NA_real_,
                eif_agg = NULL))
  }

  post_ci_v  <- post_ci[valid_mask, , drop = FALSE]
  post_ids_v <- post_ci_v$cell_id

  # Cohort-share weights (unnormalized)
  pi_g_k <- vapply(post_ids_v, function(cid) {
    panel_obj$cohort_fractions[[as.character(post_ci_v$group[post_ci_v$cell_id == cid])]]
  }, numeric(1L))
  S_u    <- sum(pi_g_k)
  if (S_u < EDID_DENOM_EPS) {
    return(list(att = NA_real_, se = NA_real_, ci_lower = NA_real_,
                ci_upper = NA_real_, t_stat = NA_real_, p_value = NA_real_,
                eif_agg = NULL))
  }
  q_norm <- pi_g_k / S_u

  # Point estimate
  att_k  <- vapply(post_ids_v, function(cid) cells[[cid]]$att, numeric(1L))
  overall_att <- sum(q_norm * att_k)

  # Aggregated EIF (direct contribution)
  eif_agg <- NULL
  if (!is.null(eif_matrix)) {
    n       <- panel_obj$n
    eif_agg <- numeric(n)
    for (ii in seq_along(post_ids_v)) {
      cid    <- post_ids_v[ii]
      eif_c  <- .eif_col(eif_matrix, cid)
      if (!is.null(eif_c)) {
        eif_agg <- eif_agg + q_norm[ii] * eif_c
      }
    }

    # WIF correction: weight function for overall ATT
    weight_fn_overall <- function(cells_arg, cell_index_arg, panel_obj_arg) {
      ci_post <- cell_index_arg[!cell_index_arg$is_pre, , drop = FALSE]
      ids     <- ci_post$cell_id
      valid   <- vapply(ids, function(cid) {
        !is.null(cells_arg[[cid]]$att) && is.finite(cells_arg[[cid]]$att)
      }, logical(1L))
      ids_v   <- ids[valid]
      ci_v    <- ci_post[valid, , drop = FALSE]
      pg      <- vapply(ids_v, function(cid) {
        panel_obj_arg$cohort_fractions[[as.character(
          ci_v$group[ci_v$cell_id == cid])]]
      }, numeric(1L))
      Su <- sum(pg)
      if (Su < EDID_DENOM_EPS) {
        return(stats::setNames(rep(0, length(ids_v)), as.character(ids_v)))
      }
      stats::setNames(pg / Su, as.character(ids_v))
    }

    wif <- compute_wif_contribution_edid(
      weight_fn_overall, cells, eif_matrix, cell_index, panel_obj, overall_att
    )
    eif_agg <- eif_agg + wif
  }

  # SE and inference
  inf_res <- safe_inference_edid(
    eif_agg, panel_obj$cluster_indices, alpha, overall_att
  )

  list(
    att      = overall_att,
    se       = inf_res$se,
    ci_lower = inf_res$ci_lower,
    ci_upper = inf_res$ci_upper,
    t_stat   = inf_res$t_stat,
    p_value  = inf_res$p_value,
    eif_agg  = eif_agg
  )
}

# ---------------------------------------------------------------------------
# Event-study aggregation
# ---------------------------------------------------------------------------

#' Aggregate cell-level ATTs by relative time (event study)
#'
#' For each unique relative time \eqn{e = t - g}, computes the cohort-share-
#' weighted average ATT over all \code{(g, t)} cells with \eqn{t - g = e}.
#' Includes WIF correction.
#'
#' @param cells list of \code{edid_cell_result} objects
#' @param eif_matrix n x n_cells numeric matrix (or NULL)
#' @param cell_index data.frame with columns \code{group}, \code{time},
#'   \code{cell_id}, \code{is_pre}
#' @param panel_obj panel object
#' @param alpha significance level
#' @param balance_e integer or NULL: if not NULL, restrict to \code{[-balance_e, balance_e]}
#'
#' @return named list, one entry per unique relative time; each entry is a list
#'   with \code{e}, \code{att}, \code{se}, \code{ci_lower}, \code{ci_upper},
#'   \code{t_stat}, \code{p_value}, \code{eif_agg}
#' @keywords internal
aggregate_event_study_edid <- function(
  cells, eif_matrix, cell_index, panel_obj, alpha, balance_e = NULL
) {
  # Compute relative times
  cell_index$e <- cell_index$time - cell_index$group

  # Apply balance_e restriction
  if (!is.null(balance_e)) {
    cell_index <- cell_index[abs(cell_index$e) <= balance_e, , drop = FALSE]
  }

  unique_e <- sort(unique(cell_index$e))
  result   <- vector("list", length(unique_e))
  names(result) <- as.character(unique_e)

  n <- panel_obj$n

  for (ii in seq_along(unique_e)) {
    e_val   <- unique_e[ii]
    e_mask  <- cell_index$e == e_val
    e_ci    <- cell_index[e_mask, , drop = FALSE]
    e_ids   <- e_ci$cell_id

    # Valid cells
    valid   <- vapply(e_ids, function(cid) {
      !is.null(cells[[cid]]$att) && is.finite(cells[[cid]]$att)
    }, logical(1L))

    na_entry <- list(e = e_val, att = NA_real_, se = NA_real_,
                     ci_lower = NA_real_, ci_upper = NA_real_,
                     t_stat = NA_real_, p_value = NA_real_, eif_agg = NULL)
    if (!any(valid)) {
      result[[ii]] <- na_entry
      next
    }

    e_ci_v  <- e_ci[valid, , drop = FALSE]
    e_ids_v <- e_ci_v$cell_id

    pi_g_k  <- vapply(seq_along(e_ids_v), function(jj) {
      panel_obj$cohort_fractions[[as.character(e_ci_v$group[jj])]]
    }, numeric(1L))
    S_u    <- sum(pi_g_k)
    if (S_u < EDID_DENOM_EPS) { result[[ii]] <- na_entry; next }
    q_norm <- pi_g_k / S_u

    att_k  <- vapply(e_ids_v, function(cid) cells[[cid]]$att, numeric(1L))
    es_att <- sum(q_norm * att_k)

    eif_agg <- NULL
    if (!is.null(eif_matrix)) {
      eif_agg <- numeric(n)
      for (jj in seq_along(e_ids_v)) {
        cid   <- e_ids_v[jj]
        eif_c <- .eif_col(eif_matrix, cid)
        if (!is.null(eif_c)) eif_agg <- eif_agg + q_norm[jj] * eif_c
      }

      # Build a cell_index restricted to this e-group for WIF
      # WIF weight function: cohort-share weights restricted to this e
      e_cell_index_restricted <- cell_index
      # Mark non-e cells as pre so weight_fn ignores them
      e_cell_index_restricted$is_pre <- TRUE
      e_cell_index_restricted$is_pre[e_mask] <- e_ci$is_pre

      weight_fn_e <- .make_weight_fn_es(e_val)
      wif <- compute_wif_contribution_edid(
        weight_fn_e, cells, eif_matrix, e_cell_index_restricted, panel_obj, es_att
      )
      eif_agg <- eif_agg + wif
    }

    inf_res <- safe_inference_edid(eif_agg, panel_obj$cluster_indices, alpha, es_att)

    result[[ii]] <- list(
      e        = e_val,
      att      = es_att,
      se       = inf_res$se,
      ci_lower = inf_res$ci_lower,
      ci_upper = inf_res$ci_upper,
      t_stat   = inf_res$t_stat,
      p_value  = inf_res$p_value,
      eif_agg  = eif_agg
    )
  }
  result
}

# Helper: build weight function for event-study relative time e
.make_weight_fn_es <- function(e_val) {
  force(e_val)
  function(cells_arg, cell_index_arg, panel_obj_arg) {
    ci <- cell_index_arg
    ci$e_local <- ci$time - ci$group
    e_mask_local <- ci$e_local == e_val & !ci$is_pre
    ids_v <- ci$cell_id[e_mask_local]
    valid <- vapply(ids_v, function(cid) {
      !is.null(cells_arg[[cid]]$att) && is.finite(cells_arg[[cid]]$att)
    }, logical(1L))
    ids_v  <- ids_v[valid]
    ci_v   <- ci[e_mask_local, , drop = FALSE][valid, , drop = FALSE]
    pg     <- vapply(seq_along(ids_v), function(jj) {
      panel_obj_arg$cohort_fractions[[as.character(ci_v$group[jj])]]
    }, numeric(1L))
    Su     <- sum(pg)
    if (Su < EDID_DENOM_EPS) {
      return(stats::setNames(rep(0, length(ids_v)), as.character(ids_v)))
    }
    stats::setNames(pg / Su, as.character(ids_v))
  }
}

# ---------------------------------------------------------------------------
# Group aggregation
# ---------------------------------------------------------------------------

#' Aggregate cell-level ATTs by treatment cohort
#'
#' For each cohort \code{g}, computes the equal-time-weighted average ATT over
#' all post-treatment cells \code{(g, t)} with \code{t >= g}.
#' Group aggregation uses equal weights so there is no WIF correction.
#'
#' @param cells list of \code{edid_cell_result} objects
#' @param eif_matrix n x n_cells numeric matrix (or NULL)
#' @param cell_index data.frame with columns \code{group}, \code{time},
#'   \code{cell_id}, \code{is_pre}
#' @param panel_obj panel object
#' @param alpha significance level
#'
#' @return named list, one entry per cohort; each entry is a list with
#'   \code{group}, \code{att}, \code{se}, \code{ci_lower}, \code{ci_upper},
#'   \code{t_stat}, \code{p_value}, \code{eif_agg}
#' @keywords internal
aggregate_group_edid <- function(cells, eif_matrix, cell_index, panel_obj, alpha) {

  tgroups <- panel_obj$treatment_groups
  result  <- vector("list", length(tgroups))
  names(result) <- as.character(tgroups)
  n <- panel_obj$n

  for (ii in seq_along(tgroups)) {
    g_val  <- tgroups[ii]
    g_mask <- cell_index$group == g_val & !cell_index$is_pre
    g_ci   <- cell_index[g_mask, , drop = FALSE]
    g_ids  <- g_ci$cell_id

    na_entry <- list(group = g_val, att = NA_real_, se = NA_real_,
                     ci_lower = NA_real_, ci_upper = NA_real_,
                     t_stat = NA_real_, p_value = NA_real_, eif_agg = NULL)

    valid <- vapply(g_ids, function(cid) {
      !is.null(cells[[cid]]$att) && is.finite(cells[[cid]]$att)
    }, logical(1L))
    if (!any(valid)) { result[[ii]] <- na_entry; next }

    g_ids_v  <- g_ids[valid]
    m_g      <- length(g_ids_v)
    att_k    <- vapply(g_ids_v, function(cid) cells[[cid]]$att, numeric(1L))
    group_att <- mean(att_k)  # equal weights

    eif_agg <- NULL
    if (!is.null(eif_matrix)) {
      eif_agg <- numeric(n)
      for (cid in g_ids_v) {
        eif_c <- .eif_col(eif_matrix, cid)
        if (!is.null(eif_c)) eif_agg <- eif_agg + eif_c / m_g
      }
      # Equal weights: no WIF correction needed (weights don't depend on pi_g)
    }

    inf_res <- safe_inference_edid(eif_agg, panel_obj$cluster_indices, alpha, group_att)

    result[[ii]] <- list(
      group    = g_val,
      att      = group_att,
      se       = inf_res$se,
      ci_lower = inf_res$ci_lower,
      ci_upper = inf_res$ci_upper,
      t_stat   = inf_res$t_stat,
      p_value  = inf_res$p_value,
      eif_agg  = eif_agg
    )
  }
  result
}
