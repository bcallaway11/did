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

  # Participating post-treatment cells. (The event-study caller marks non-e cells
  # as is_pre, so this selects exactly the cells entering THIS aggregation.)
  # Restrict to finite-ATT cells -- an NA-att (empty) cell would NaN-poison the
  # aggregated EIF/SE and silently return an NA SE on a finite ATT.
  post_ci <- cell_index[!cell_index$is_pre, , drop = FALSE]
  if (nrow(post_ci) == 0L) return(wif)
  valid <- vapply(post_ci$cell_id, function(cid) {
    !is.null(cells[[cid]]$att) && is.finite(cells[[cid]]$att)
  }, logical(1L))
  post_ci <- post_ci[valid, , drop = FALSE]
  K <- nrow(post_ci)
  if (K == 0L) return(wif)

  att_k <- vapply(post_ci$cell_id, function(cid) cells[[cid]]$att, numeric(1L))
  # Unnormalized cohort-share weight pi_{g_k} per participating cell; S_u = sum.
  pgg <- vapply(seq_len(K), function(k) {
    panel_obj$cohort_fractions[[as.character(post_ci$group[k])]]
  }, numeric(1L))
  S_u <- sum(pgg)
  if (S_u < EDID_DENOM_EPS) return(wif)

  # Weight-influence-function (WIF) correction for the estimated cohort-share weights, following
  # the canonical construction of the did package (did::wif). The influence function of the
  # estimated share pi_hat_g is the O(1) quantity 1(G_i = g) - pi_g, so the cohort-share weights
  # contribute this term to the aggregated influence function; omitting it (or shrinking it by an
  # n^{-1} factor) makes the aggregated standard error anti-conservative under cohort-ATT
  # heterogeneity. Here ind[i, k] = 1(G_i = g_k).
  ind <- vapply(seq_len(K), function(k) {
    as.numeric(panel_obj$cohort_masks[[as.character(post_ci$group[k])]])
  }, numeric(n))
  if (is.null(dim(ind))) ind <- matrix(ind, nrow = n)
  centered <- sweep(ind, 2L, pgg, "-")                       # 1(G = g_k) - pi_{g_k}
  wif_mat  <- centered / S_u - outer(rowSums(centered), pgg / S_u^2)
  drop(wif_mat %*% att_k)                                    # WIF_i = sum_k att_k * wif_mat[i,k]
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
      # Mark non-e cells as pre so the WIF ignores them; cells AT this event time participate
      # regardless of pre/post. This includes LEADS (e < 0): their cohort shares pi_g are also
      # estimated, so the weight-influence term must enter their SE. (The previous code set this
      # to e_ci$is_pre, which is TRUE for leads and therefore zeroed the lead WIF.)
      e_cell_index_restricted$is_pre <- TRUE
      e_cell_index_restricted$is_pre[e_mask] <- FALSE

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

# ---------------------------------------------------------------------------
# Calendar-time aggregation
# ---------------------------------------------------------------------------

#' Aggregate cell-level ATTs by calendar time
#'
#' For each post-treatment calendar period \code{t}, computes the cohort-share-
#' weighted average ATT over all cells \code{(g, t)} with \code{g <= t} (cohorts
#' already treated by \code{t}), with weights \eqn{\pi_g / \sum_{g' \le t} \pi_{g'}}.
#' Includes the WIF correction for the estimated cohort-share weights. Mirrors
#' \code{did::aggte(type = "calendar")}.
#'
#' @inheritParams aggregate_event_study_edid
#' @return named list, one entry per post calendar time; each entry is a list with
#'   \code{time}, \code{att}, \code{se}, \code{ci_lower}, \code{ci_upper},
#'   \code{t_stat}, \code{p_value}, \code{eif_agg}
#' @keywords internal
aggregate_calendar_edid <- function(cells, eif_matrix, cell_index, panel_obj, alpha) {
  post_times <- sort(unique(cell_index$time[!cell_index$is_pre]))
  result <- vector("list", length(post_times))
  names(result) <- as.character(post_times)
  n <- panel_obj$n

  for (ii in seq_along(post_times)) {
    t_val  <- post_times[ii]
    t_mask <- cell_index$time == t_val & !cell_index$is_pre
    t_ci   <- cell_index[t_mask, , drop = FALSE]
    t_ids  <- t_ci$cell_id

    na_entry <- list(time = t_val, att = NA_real_, se = NA_real_,
                     ci_lower = NA_real_, ci_upper = NA_real_,
                     t_stat = NA_real_, p_value = NA_real_, eif_agg = NULL)

    valid <- vapply(t_ids, function(cid) {
      !is.null(cells[[cid]]$att) && is.finite(cells[[cid]]$att)
    }, logical(1L))
    if (!any(valid)) { result[[ii]] <- na_entry; next }

    t_ci_v  <- t_ci[valid, , drop = FALSE]
    t_ids_v <- t_ci_v$cell_id

    pi_g_k  <- vapply(seq_along(t_ids_v), function(jj) {
      panel_obj$cohort_fractions[[as.character(t_ci_v$group[jj])]]
    }, numeric(1L))
    S_u <- sum(pi_g_k)
    if (S_u < EDID_DENOM_EPS) { result[[ii]] <- na_entry; next }
    q_norm <- pi_g_k / S_u

    att_k   <- vapply(t_ids_v, function(cid) cells[[cid]]$att, numeric(1L))
    cal_att <- sum(q_norm * att_k)

    eif_agg <- NULL
    if (!is.null(eif_matrix)) {
      eif_agg <- numeric(n)
      for (jj in seq_along(t_ids_v)) {
        cid   <- t_ids_v[jj]
        eif_c <- .eif_col(eif_matrix, cid)
        if (!is.null(eif_c)) eif_agg <- eif_agg + q_norm[jj] * eif_c
      }
      # WIF: restrict cell_index to this calendar time (mark non-t cells as pre so the
      # cohort-share WIF is computed over exactly the (g, t_val) cells with g <= t_val).
      t_cell_index_restricted <- cell_index
      t_cell_index_restricted$is_pre <- TRUE
      t_cell_index_restricted$is_pre[t_mask] <- t_ci$is_pre
      wif <- compute_wif_contribution_edid(
        function(...) NULL, cells, eif_matrix, t_cell_index_restricted, panel_obj, cal_att
      )
      eif_agg <- eif_agg + wif
    }

    inf_res <- safe_inference_edid(eif_agg, panel_obj$cluster_indices, alpha, cal_att)

    result[[ii]] <- list(
      time     = t_val,
      att      = cal_att,
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
