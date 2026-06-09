#' @title Compute Group-Time Average Treatment Effects
#'
#' @description `compute.att_gt` does the main work for computing
#'  multiperiod group-time average treatment effects
#'
#' @param dp A DIDparams object
#'
#' @return a list with length equal to the number of groups times the
#'  number of time periods; each element of the list contains an
#'  object that contains group-time average treatment effect as well
#'  as which group it is for and which time period it is for. It also exports
#'  the influence function which is used externally to compute
#'  standard errors.
#'
#' @keywords internal
#'
#' @export
compute.att_gt <- function(dp) {
  #-----------------------------------------------------------------------------
  # unpack DIDparams
  #-----------------------------------------------------------------------------
  data <- data.table::as.data.table(dp$data)
  yname <- dp$yname
  tname <- dp$tname
  idname <- dp$idname
  xformla <- dp$xformla
  weightsname <- dp$weightsname
  est_method <- dp$est_method
  extra_args <- if (is.null(dp$extra_args)) list() else dp$extra_args
  base_period <- dp$base_period
  panel <- dp$panel
  true_repeated_cross_sections <- dp$true_repeated_cross_sections
  print_details <- dp$print_details
  control_group <- dp$control_group
  anticipation <- dp$anticipation
  gname <- dp$gname
  n <- dp$n
  nT <- dp$nT
  nG <- dp$nG
  tlist <- dp$tlist
  glist <- dp$glist

  #-----------------------------------------------------------------------------
  # main computations
  #-----------------------------------------------------------------------------

  # will populate with all att(g,t)
  attgt.list <- list()

  # place holder in lists
  counter <- 1

  # number of time periods
  tlist.length <- if (base_period != "universal") length(tlist) - 1L else length(tlist)
  tfac <- if (base_period != "universal") 1L else 0L

  # influence function dimensions
  inffunc_nrow <- n

  # list to collect sparse matrix updates (built into sparse matrix after loops)
  inffunc_updates <- list()
  # counter for keeping track of updates
  update_counter <- 1

  # never treated option
  nevertreated <- (control_group[1] == "nevertreated")

  fix_weights <- dp$fix_weights

  # Whether to compute influence functions (default TRUE; FALSE = point estimates
  # only). When FALSE we skip requesting influence functions from the 2x2 estimators,
  # skip the per-cell influence-function scaling, and skip storing inffunc_updates --
  # so a point-estimates-only run does not pay the influence-function compute/memory
  # cost. The att estimates are identical either way.
  do_inf <- is.null(dp$compute_inffunc) || isTRUE(dp$compute_inffunc)

  # Pre-extract columns to avoid repeated get() inside data.table (which is slow)
  g_col <- data[[gname]]
  t_col <- data[[tname]]

  # Build weight lookup by period for fix_weights options (balanced panel only)
  if (!is.null(fix_weights) && panel) {
    weights_by_period <- list()
    for (tp in seq_along(tlist)) {
      weights_by_period[[tp]] <- data[t_col == tlist[tp], .w]
    }
  }

  if (nevertreated) {
    set(data, j = ".C", value = as.integer(g_col == 0))
  }
  set(data, j = ".y", value = data[[yname]])

  # -- design matrix: built ONCE over the full data, then row-subset per (g,t) cell --
  # This (a) avoids rebuilding model.matrix() nG*nT times, and (b) gives every cell the
  # GLOBAL factor levels / GLOBAL transform basis -- exactly what the repeated-cross-
  # section path, the fast path, and manually-constructed dummy columns use. As a
  # result a factor covariate behaves identically to adding its dummies by hand: a
  # level absent from a particular 2x2 comparison becomes an all-zero column (handled
  # by the estimator's overlap/rank check, giving the same NA + warning as the manual
  # dummies) rather than a dropped contrast that changes the design across cells.
  #
  # For purely row-wise numeric formulae (e.g. ~X, ~X1+X2, ~X1*X2, ~I(X^2)) the slice
  # is bit-identical to a per-cell model.matrix(xformla, disdat); for factors it equals
  # the manual-dummy design; for data-dependent bases (poly/scale/ns/...) it uses the
  # global basis, which yields the same ATT/SE as a per-cell fit (the estimator is
  # invariant to a full-rank reparameterization of the covariates).
  #
  # `data` is already complete-cases and finite (see pre_process_did), so na.pass adds
  # no NA rows and full_mm aligns 1:1 (by row position) with `data`.
  full_mm <- model.matrix(xformla, data = data, na.action = na.pass)
  if (panel) {
    # per-period row slices keyed by idname, for the earlier-period covariates that
    # get_wide_data() retains. (idname is used rather than .rowid because .rowid is
    # only created on the RC/unbalanced path; in a balanced panel idname is unique
    # within each period.)
    period_mm <- vector("list", length(tlist))
    period_id <- vector("list", length(tlist))
    period_y  <- vector("list", length(tlist))
    period_w  <- vector("list", length(tlist))
    y_col <- data[[".y"]]
    w_col <- data[[".w"]]
    for (tp in seq_along(tlist)) {
      rows_tp <- which(t_col == tlist[tp])
      period_mm[[tp]] <- full_mm[rows_tp, , drop = FALSE]
      period_id[[tp]] <- data[[idname]][rows_tp]
      period_y[[tp]]  <- y_col[rows_tp]
      period_w[[tp]]  <- w_col[rows_tp]
    }
    # Precompute fast path: when the panel is balanced with every period sharing the
    # same units in the same (id-sorted) row order and default weights are used, each
    # 2x2 cell can be assembled by position from these per-period blocks -- skipping
    # the per-cell data[time_mask] subset + get_wide_data() reshape (the dominant
    # cost on the slow path). Falls through to the original get_wide_data() code,
    # bit-identically, for fix_weights / unaligned panels.
    panel_units_aligned <- all(vapply(period_id, function(x)
      length(x) == length(period_id[[1L]]) && all(x == period_id[[1L]]), logical(1)))
    g_unit_panel <- g_col[which(t_col == tlist[1L])]
    # options(did.disable_precompute = TRUE) forces the original get_wide_data() path
    # (escape hatch / used to verify bit-identical equivalence).
    use_precompute_panel <- is.null(fix_weights) && panel_units_aligned &&
      !isTRUE(getOption("did.disable_precompute"))
  } else {
    # repeated cross-section / unbalanced: row-subset the full design by POSITION
    # (.rowid is not unique across periods here).
    global_mm <- full_mm
    use_precompute_panel <- FALSE
  }

  # loop over groups
  for (g in 1:nG) {
    # Set up .G once (use pre-extracted g_col to avoid get() inside data.table)
    current_g <- glist[g]
    set(data, j = ".G", value = as.numeric(g_col == current_g))

    # pre-compute the universal/post-treatment base period for this group (used multiple times)
    idx_g <- which((tlist + anticipation) < glist[g])
    pret_g <- if (length(idx_g) == 0L) NA_integer_ else idx_g[length(idx_g)]

    # loop over time periods
    for (t in 1:tlist.length) {
      #-----------------------------------------------------------------------------
      # Set pret

      # varying base period
      pret <- t

      # universal base period
      if (base_period == "universal") {
        # use same base period as for post-treatment periods
        pret <- pret_g
      }

      # check if in post-treatment period
      if ((glist[g] <= tlist[(t + tfac)])) {
        # update pre-period if in post-treatment period to
        # be  period (g-delta-1)
        pret <- pret_g
      }

      # check if there are no pre-treatment periods
      if (is.na(pret)) {
        warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g], "\nUnits from this group are dropped"))
        break
      }

      # use "not yet treated as control"
      # that is, never treated + units that are eventually treated,
      # but not treated by the current period (+ anticipation)
      if (!nevertreated) {
        # Use pre-extracted g_col to avoid get() inside data.table
        time_threshold <- tlist[max(t, pret) + tfac] + anticipation
        set(data, j = ".C", value = as.integer((g_col == 0) |
          ((g_col > time_threshold) &
            (g_col != current_g))))
      }


      #-----------------------------------------------------------------------------
      # if we are in period (g-1), normalize results to be equal to 0
      # and break without computing anything
      if (base_period == "universal") {
        if (tlist[pret] == tlist[(t + tfac)]) {
          attgt.list[[counter]] <- list(att = 0, group = glist[g], year = tlist[(t + tfac)], post = 0)
          # inffunc[,counter] <- rep(0,n)
          # counter <- counter+1
          if (do_inf) inffunc_updates[[update_counter]] <- list(
            indices = integer(0), # No non-zero entries
            values = numeric(0)
          )

          # Update the counters
          update_counter <- update_counter + 1
          counter <- counter + 1
          next
        }
      }

      # print the details of which iteration we are on
      if (print_details) {
        cat(paste("current period:", tlist[(t + tfac)]), "\n")
        cat(paste("current group:", glist[g]), "\n")
        cat(paste("set pretreatment period to be", tlist[pret]), "\n")
      }

      #-----------------------------------------------------------------------------
      # results for the case with panel data
      #-----------------------------------------------------------------------------

      # post treatment dummy variable
      post.treat <- 1 * (glist[g] <= tlist[t + tfac])

      # total number of units (not just included in G or C)
      # disdat <- data[data[,tname] == tlist[t+tfac] | data[,tname] == tlist[pret],]
      target_times <- c(tlist[t + tfac], tlist[pret])
      time_mask <- t_col %in% target_times
      # The precompute panel path assembles each cell from per-period blocks and does
      # not need this long 2-period subset (the dominant per-cell cost); RC,
      # fix_weights, and unaligned panels still build disdat as before.
      if (!(panel && use_precompute_panel)) {
        disdat <- data[time_mask]
      }


      if (panel) {
        if (use_precompute_panel) {
          # ---- precompute fast path (balanced, id-aligned panel; default weights) --
          # Build the 2x2 cell directly from the per-period blocks. g is time-invariant
          # so .G / .C are unit-level; Ypost is always period (t+tfac) and Ypre always
          # period pret because get_wide_data's .y1/.y0 (later/earlier) swap exactly
          # cancels the universal-base ordering used below. The earlier period (whose
          # weights + covariates get_wide_data retains) is min(t + tfac, pret).
          G_full <- as.numeric(g_unit_panel == current_g)
          if (nevertreated) {
            C_full <- as.numeric(g_unit_panel == 0)
          } else {
            time_threshold <- tlist[max(t, pret) + tfac] + anticipation
            C_full <- as.numeric((g_unit_panel == 0) |
              ((g_unit_panel > time_threshold) & (g_unit_panel != current_g)))
          }
          disidx <- G_full == 1 | C_full == 1
          n <- length(g_unit_panel)
          kept <- which(disidx)
          n1 <- length(kept)
          G <- G_full[kept]
          C <- C_full[kept]
          Ypost <- period_y[[t + tfac]][kept]
          Ypre  <- period_y[[pret]][kept]
          earlier_idx <- min(t + tfac, pret)
          w <- period_w[[earlier_idx]][kept]
          covariates <- period_mm[[earlier_idx]][kept, , drop = FALSE]
          dimnames(covariates) <- list(NULL, colnames(covariates))
        } else {
          # ---- original path: per-cell get_wide_data() reshape ----
          # transform  disdat it into "cross-sectional" data where one of the columns
          # contains the change in the outcome over time.
          disdat <- get_wide_data(disdat, yname, idname, tname)

          # still total number of units (not just included in G or C)
          n <- nrow(disdat)

          # pick up the indices for units that will be used to compute ATT(g,t)
          disidx <- disdat$.G == 1 | disdat$.C == 1

          # pick up the data that will be used to compute ATT(g,t)
          disdat <- disdat[disidx]

          n1 <- nrow(disdat) # num obs. for computing ATT(g,t)

          # drop missing factors
          disdat <- droplevels(disdat)

          # give short names for data in this iteration
          G <- disdat$.G
          C <- disdat$.C

          # handle pre-treatment universal base period differently
          # we need to do this because panel2cs2 just puts later period
          # in .y1, but if we are in a pre-treatment period with a universal
          # base period, then the "base period" is actually the later period
          Ypre <- if (tlist[(t + tfac)] > tlist[pret]) disdat$.y0 else disdat$.y1
          Ypost <- if (tlist[(t + tfac)] > tlist[pret]) disdat$.y1 else disdat$.y0

          # Select weights based on fix_weights
          if (is.null(fix_weights)) {
            # Default: .w from get_wide_data (earlier period)
            w <- disdat$.w
          } else if (fix_weights == "base_period") {
            w <- weights_by_period[[pret_g]][disidx]
          } else if (fix_weights == "first_period") {
            w <- weights_by_period[[1L]][disidx]
          } else if (fix_weights == "varying") {
            w <- disdat$.w  # will be overridden below when switching to RC estimator
          } else {
            w <- disdat$.w
          }

          # matrix of covariates: index the precomputed earlier-period design by this
          # cell's units (idname). get_wide_data retained the earlier period, whose
          # index is min(t + tfac, pret).
          earlier_idx <- min(t + tfac, pret)
          covariates <- period_mm[[earlier_idx]][match(disdat[[idname]], period_id[[earlier_idx]]), , drop = FALSE]
          dimnames(covariates) <- list(NULL, colnames(covariates))
        }

        #-----------------------------------------------------------------------------
        # code for actually computing att(g,t)
        #-----------------------------------------------------------------------------

        attgt <- tryCatch({
          if (!is.null(fix_weights) && fix_weights == "varying") {
            # fix_weights = "varying": use RC estimators with per-period weights
            # Go back to long-format data for this (g,t) cell
            disdat_long <- data[time_mask]
            disdat_long_idx <- disdat_long$.G == 1 | disdat_long$.C == 1
            disdat_long <- droplevels(disdat_long[disdat_long_idx])
            Y_rc <- disdat_long[[yname]]
            G_rc <- disdat_long$.G
            post_rc <- as.numeric(disdat_long[[tname]] == tlist[t + tfac])
            w_rc <- disdat_long$.w
            # Use earlier-period covariates for all observations — fix_weights
            # only changes weights, not the covariate conditioning set.
            # Use min(pret, t+tfac) to match the panel estimator's convention:
            # with base_period="universal", pret can be later than t for placebo cells.
            earlier_idx_v <- min(pret, t + tfac)
            earlier_period <- tlist[earlier_idx_v]
            early_mask <- disdat_long[[tname]] == earlier_period
            disdat_early <- disdat_long[early_mask]
            cov_early <- period_mm[[earlier_idx_v]][match(disdat_early[[idname]], period_id[[earlier_idx_v]]), , drop = FALSE]
            dimnames(cov_early) <- list(NULL, colnames(cov_early))
            # Map each row in disdat_long to its unit's earlier-period covariates
            early_ids <- disdat_early[[idname]]
            all_ids <- disdat_long[[idname]]
            id_map <- match(all_ids, early_ids)
            covariates_rc <- cov_early[id_map, , drop = FALSE]

            # Run overlap/rank checks on RC data (not wide panel data)
            if (!is.function(est_method)) {
              if (est_method %in% c("dr", "ipw")) {
                preliminary_logit <- overlap_logit_fit(covariates_rc, G_rc)
                if (max(preliminary_logit$fitted.values) >= 0.999) {
                  warning(paste0("overlap condition violated for ", glist[g], " in time period ", tlist[t + tfac]))
                  stop("overlap")
                }
              }
              if (est_method %in% c("dr", "reg")) {
                control_covs_rc <- covariates_rc[G_rc == 0, , drop = FALSE]
                if (rcond(t(control_covs_rc) %*% control_covs_rc) < .Machine$double.eps) {
                  warning(paste0("Not enough control units for group ", glist[g], " in time period ", tlist[t + tfac], " to run specified regression"))
                  stop("singular")
                }
              }
            }

            if (inherits(est_method, "function")) {
              res <- do.call(est_method, c(list(
                y = Y_rc, post = post_rc,
                D = G_rc, covariates = covariates_rc,
                i.weights = w_rc, inffunc = do_inf
              ), extra_args))
            } else if (est_method == "ipw") {
              res <- DRDID::std_ipw_did_rc(Y_rc, post_rc, G_rc,
                covariates = covariates_rc,
                i.weights = w_rc, boot = FALSE, inffunc = do_inf)
            } else if (est_method == "reg") {
              res <- DRDID::reg_did_rc(Y_rc, post_rc, G_rc,
                covariates = covariates_rc,
                i.weights = w_rc, boot = FALSE, inffunc = do_inf)
            } else {
              res <- DRDID::drdid_rc(Y_rc, post_rc, G_rc,
                covariates = covariates_rc,
                i.weights = w_rc, boot = FALSE, inffunc = do_inf)
            }
          } else {
            # Panel path: run overlap/rank checks on panel data
            if (!is.function(est_method)) {
              if (est_method %in% c("dr", "ipw")) {
                preliminary_logit <- overlap_logit_fit(covariates, G)
                if (max(preliminary_logit$fitted.values) >= 0.999) {
                  warning(paste0("overlap condition violated for ", glist[g], " in time period ", tlist[t + tfac]))
                  stop("overlap")
                }
              }
              if (est_method %in% c("dr", "reg")) {
                control_covs <- covariates[G == 0, , drop = FALSE]
                if (rcond(t(control_covs) %*% control_covs) < .Machine$double.eps) {
                  warning(paste0("Not enough control units for group ", glist[g], " in time period ", tlist[t + tfac], " to run specified regression"))
                  stop("singular")
                }
              }
            }

            if (inherits(est_method, "function")) {
              # user-specified function
              res <- do.call(est_method, c(list(
                y1 = Ypost, y0 = Ypre,
                D = G,
                covariates = covariates,
                i.weights = w,
                inffunc = do_inf
              ), extra_args))
            } else if (est_method == "ipw") {
              # inverse-probability weights
              res <- DRDID::std_ipw_did_panel(Ypost, Ypre, G,
                covariates = covariates,
                i.weights = w,
                boot = FALSE, inffunc = do_inf
              )
            } else if (est_method == "reg") {
              # regression
              res <- DRDID::reg_did_panel(Ypost, Ypre, G,
                covariates = covariates,
                i.weights = w,
                boot = FALSE, inffunc = do_inf
              )
            } else {
              # doubly robust, this is default
              res <- DRDID::drdid_panel(Ypost, Ypre, G,
                covariates = covariates,
                i.weights = w,
                boot = FALSE, inffunc = do_inf
              )
            }
          }

          # adjust influence function to account for only using
          # subgroup to estimate att(g,t) (skipped for point estimates only)
          if (do_inf) {
            if (!is.null(fix_weights) && fix_weights == "varying") {
              # RC influence function has one entry per obs in disdat_long
              # (2 per unit: pre + post). Aggregate to unit level by ID,
              # independent of row ordering.
              res$att.inf.func <- as.numeric(rowsum(res$att.inf.func,
                                                    disdat_long[[idname]],
                                                    reorder = FALSE))
              res$att.inf.func <- (n / n1) * res$att.inf.func
            } else {
              res$att.inf.func <- (n / n1) * res$att.inf.func
            }
          }
          res
        }, error = function(e) {
          warning("Error computing internal 2x2 DiD for (g, t) = (", glist[g], ", ", tlist[t + tfac], "): ", e$message, ". The ATT for this cell will be set to NA.")
          NULL
        })

        if (is.null(attgt)) {
          attgt.list[[counter]] <- list(att = NA, group = glist[g], year = tlist[(t + tfac)], post = post.treat)
          if (do_inf) inffunc_updates[[update_counter]] <- list(
            indices = seq_len(n),
            values = rep(NA_real_, n)
          )
          update_counter <- update_counter + 1
          counter <- counter + 1
          next
        }
      } else { # repeated cross sections / unbalanced panel

        # pick up the indices for units that will be used to compute ATT(g,t)
        # these conditions are (1) you are observed in the right period and
        # (2) you are in the right group (it is possible to be observed in
        # the right period but still not be part of the treated or control
        # group in that period here
        rightids <- disdat$.rowid[disdat$.G == 1 | disdat$.C == 1]

        # this is the fix for unbalanced panels; 2nd criteria shouldn't do anything
        # with true repeated cross sections, but should pick up the right time periods
        # only with unbalanced panel
        disidx <- (data$.rowid %in% rightids) & ((t_col == tlist[t + tfac]) | (t_col == tlist[pret]))

        # pick up the data that will be used to compute ATT(g,t)
        disdat <- data[disidx, ]
        # positional index into `data` for slicing the global design matrix, kept in
        # sync with any later row drops below (droplevels does not reorder/drop rows).
        # global_mm[disdat_rows, ] is the design for the rows used in this cell.
        disdat_rows <- which(disidx)

        # drop missing factors
        disdat <- droplevels(disdat)

        # give short names for data in this iteration
        G <- disdat$.G
        C <- disdat$.C
        Y <- disdat[[yname]]
        post <- 1 * (disdat[[tname]] == tlist[t + tfac])
        # num obs. for computing ATT(g,t), have to be careful here
        n1 <- sum(G + C)

        # Handle fix_weights for RC/unbalanced panel
        if (!is.null(fix_weights) && fix_weights %in% c("base_period", "first_period")) {
          # Determine which period's weight to use
          if (fix_weights == "base_period") {
            target_period <- tlist[pret_g]
          } else {
            target_period <- tlist[1]
          }
          # Build lookup: weight from target period per unit
          target_rows <- data[t_col == target_period, ]
          # Map each observation's weight from the target period by integer .rowid
          # matching (equivalent to the previous named-character lookup: first match
          # wins, unmatched -> NA), without coercing every .rowid to character.
          w <- target_rows$.w[match(disdat$.rowid, target_rows$.rowid)]
          # Drop units not observed in the target period
          missing_w <- is.na(w)
          if (any(missing_w)) {
            n_dropped <- length(unique(disdat$.rowid[missing_w]))
            warning(paste0("Dropped ", n_dropped, " units not observed in ",
                           fix_weights, " (period ", target_period, ") ",
                           "for group ", glist[g], " in time period ", tlist[t + tfac]))
            disdat <- disdat[!missing_w, ]
            disdat_rows <- disdat_rows[!missing_w]
            G <- disdat$.G
            C <- disdat$.C
            Y <- disdat[[yname]]
            post <- 1 * (disdat[[tname]] == tlist[t + tfac])
            n1 <- sum(G + C)
            w <- w[!missing_w]
          }
        } else {
          w <- disdat$.w
        }

        #-----------------------------------------------------------------------------
        # checks to make sure that we have enough observations
        skip_this_att_gt <- FALSE
        if (sum(G * post) == 0) {
          warning(paste0("No units in group ", glist[g], " in time period ", tlist[t + tfac]))
          skip_this_att_gt <- TRUE
        }
        if (sum(G * (1 - post)) == 0) {
          warning(paste0("No units in group ", glist[g], " in time period ", tlist[t]))
          skip_this_att_gt <- TRUE
        }
        if (sum(C * post) == 0) {
          warning(paste0("No available control units for group ", glist[g], " in time period ", tlist[t + tfac]))
          skip_this_att_gt <- TRUE
        }
        if (sum(C * (1 - post)) == 0) {
          warning(paste0("No available control units for group ", glist[g], " in time period ", tlist[t]))
          skip_this_att_gt <- TRUE
        }

        if (skip_this_att_gt) {
          attgt.list[[counter]] <- list(att = NA, group = glist[g], year = tlist[(t + tfac)], post = post.treat)
          # inffunc[,counter] <- NA
          # counter <- counter+1
          if (do_inf) inffunc_updates[[update_counter]] <- list(
            indices = seq_len(n),
            values = rep(NA_real_, n)
          )

          # Update the counters
          update_counter <- update_counter + 1
          counter <- counter + 1
          next
        }


        # matrix of covariates: row-subset the global design by position.
        covariates <- global_mm[disdat_rows, , drop = FALSE]
        dimnames(covariates) <- list(NULL, colnames(covariates))

        #-----------------------------------------------------------------------------
        # more checks for enough observations in each group

        # if using custom estimation method, skip this part
        custom_est_method <- is.function(est_method)

        if (!custom_est_method) {
          pscore_problems_likely <- FALSE
          reg_problems_likely <- FALSE

          # checks for pscore based methods
          if (est_method %in% c("dr", "ipw")) {
            preliminary_logit <- overlap_logit_fit(covariates, G)
            preliminary_pscores <- preliminary_logit$fitted.values
            if (max(preliminary_pscores) >= 0.999) {
              pscore_problems_likely <- TRUE
              warning(paste0("overlap condition violated for ", glist[g], " in time period ", tlist[t + tfac]))
            }
          }

          # check if can run regression using control units
          if (est_method %in% c("dr", "reg")) {
            control_covs <- covariates[G == 0, , drop = FALSE]
            if (rcond(t(control_covs) %*% control_covs) < .Machine$double.eps) {
              reg_problems_likely <- TRUE
              warning(paste0("Not enough control units for group ", glist[g], " in time period ", tlist[t + tfac], " to run specified regression"))
            }
          }

          if (reg_problems_likely | pscore_problems_likely) {
            attgt.list[[counter]] <- list(att = NA, group = glist[g], year = tlist[(t + tfac)], post = post.treat)
            if (do_inf) inffunc_updates[[update_counter]] <- list(
              indices = seq_len(n),
              values = rep(NA_real_, n)
            )

            # Update the counters
            update_counter <- update_counter + 1
            counter <- counter + 1
            next
          }
        }

        #-----------------------------------------------------------------------------
        # code for actually computing att(g,t)
        #-----------------------------------------------------------------------------

        attgt <- tryCatch({
          if (inherits(est_method, "function")) {
            # user-specified function
            res <- do.call(est_method, c(list(
              y = Y,
              post = post,
              D = G,
              covariates = covariates,
              i.weights = w,
              inffunc = do_inf
            ), extra_args))
          } else if (est_method == "ipw") {
            # inverse-probability weights
            res <- DRDID::std_ipw_did_rc(
              y = Y,
              post = post,
              D = G,
              covariates = covariates,
              i.weights = w,
              boot = FALSE, inffunc = do_inf
            )
          } else if (est_method == "reg") {
            # regression
            res <- DRDID::reg_did_rc(
              y = Y,
              post = post,
              D = G,
              covariates = covariates,
              i.weights = w,
              boot = FALSE, inffunc = do_inf
            )
          } else {
            # doubly robust, this is default
            res <- DRDID::drdid_rc(
              y = Y,
              post = post,
              D = G,
              covariates = covariates,
              i.weights = w,
              boot = FALSE, inffunc = do_inf
            )
          }

          # n/n1 adjusts for estimating the att_gt only using observations from
          # groups G and C (skipped for point estimates only)
          if (do_inf) res$att.inf.func <- (n / n1) * res$att.inf.func

          # If ATT is NaN, replace it with NA, and mark influence function as missing
          if (is.nan(res$ATT)) {
            res$ATT <- NA
            if (do_inf) res$att.inf.func <- rep(NA_real_, length(res$att.inf.func))
          }
          res
        }, error = function(e) {
          warning("Error computing internal 2x2 DiD for (g, t) = (", glist[g], ", ", tlist[t + tfac], "): ", e$message, ". The ATT for this cell will be set to NA.")
          NULL
        })

        if (is.null(attgt)) {
          attgt.list[[counter]] <- list(att = NA, group = glist[g], year = tlist[(t + tfac)], post = post.treat)
          if (do_inf) inffunc_updates[[update_counter]] <- list(
            indices = seq_len(n),
            values = rep(NA_real_, n)
          )
          update_counter <- update_counter + 1
          counter <- counter + 1
          next
        }
      } # end panel if

      # save results for this att(g,t)
      attgt.list[[counter]] <- list(
        att = attgt$ATT, group = glist[g], year = tlist[(t + tfac)], post = post.treat
      )


      # populate the influence function in the right places (skipped for point
      # estimates only -- compute_inffunc = FALSE -- so we neither aggregate nor
      # store the per-cell influence function)
      if (do_inf) {
      if (panel) {
        # Store integer indices directly (avoids which() later)
        idx <- which(disidx)
        inffunc_updates[[update_counter]] <- list(
          indices = idx,
          values = attgt$att.inf.func
        )
      } else {
        # aggregate inf functions by id (order by id)
        # Use current disdat$.rowid (may differ from rightids if fix_weights dropped obs)
        # rowsum(reorder = TRUE) sums per id and returns rows in ascending-id order,
        # matching stats::aggregate() but ~40x faster (see test-slowpath-rowsum).
        current_ids <- disdat$.rowid[disdat$.G == 1 | disdat$.C == 1]
        rs <- rowsum(attgt$att.inf.func, current_ids, reorder = TRUE)
        rs_ids <- as.numeric(rownames(rs))
        idx <- which(unique(data$.rowid) %in% rs_ids)
        inffunc_updates[[update_counter]] <- list(
          indices = idx,
          values = as.numeric(rs[, 1])
        )
      }
      }


      # save it in influence function matrix
      # inffunc[g,t,] <- inf.func
      # inffunc[,counter] <- inf.func
      update_counter <- update_counter + 1

      # update counter
      counter <- counter + 1
    } # end looping over t
  } # end looping over g

  # Build the influence function sparse matrix directly from triplets.
  # Gather per-column indices/values into lists and concatenate once, rather than
  # growing the triplet vectors with c() on every iteration (which reallocates
  # O(ncol^2) in total). The triplet order is unchanged (column-major, original
  # within-column order), so the resulting sparse matrix is identical.
  # Skipped for point estimates only (compute_inffunc = FALSE): no matrix is returned.
  if (do_inf) {
    nz_list <- lapply(inffunc_updates, `[[`, "indices")
    val_list <- lapply(inffunc_updates, `[[`, "values")
    trip_i <- unlist(nz_list, use.names = FALSE)
    trip_x <- unlist(val_list, use.names = FALSE)
    trip_j <- rep.int(seq_along(inffunc_updates), lengths(nz_list))
    inffunc <- Matrix::sparseMatrix(
      i = trip_i,
      j = trip_j,
      x = trip_x,
      dims = c(inffunc_nrow, length(inffunc_updates))
    )
  } else {
    inffunc <- NULL
  }

  return(list(attgt.list = attgt.list, inffunc = inffunc))
}
