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
          inffunc_updates[[update_counter]] <- list(
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
      disdat <- data[time_mask]


      if (panel) {
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

        # matrix of covariates
        covariates <- model.matrix(xformla, data = disdat)

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
            covariates_rc <- model.matrix(xformla, data = disdat_long)

            # Run overlap/rank checks on RC data (not wide panel data)
            if (!is.function(est_method)) {
              if (est_method %in% c("dr", "ipw")) {
                preliminary_logit <- fastglm::fastglm(covariates_rc, G_rc, family = binomial())
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
                i.weights = w_rc, inffunc = TRUE
              ), extra_args))
            } else if (est_method == "ipw") {
              res <- DRDID::std_ipw_did_rc(Y_rc, post_rc, G_rc,
                covariates = covariates_rc,
                i.weights = w_rc, boot = FALSE, inffunc = TRUE)
            } else if (est_method == "reg") {
              res <- DRDID::reg_did_rc(Y_rc, post_rc, G_rc,
                covariates = covariates_rc,
                i.weights = w_rc, boot = FALSE, inffunc = TRUE)
            } else {
              res <- DRDID::drdid_rc(Y_rc, post_rc, G_rc,
                covariates = covariates_rc,
                i.weights = w_rc, boot = FALSE, inffunc = TRUE)
            }
          } else {
            # Panel path: run overlap/rank checks on panel data
            if (!is.function(est_method)) {
              if (est_method %in% c("dr", "ipw")) {
                preliminary_logit <- fastglm::fastglm(covariates, G, family = binomial())
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
                inffunc = TRUE
              ), extra_args))
            } else if (est_method == "ipw") {
              # inverse-probability weights
              res <- DRDID::std_ipw_did_panel(Ypost, Ypre, G,
                covariates = covariates,
                i.weights = w,
                boot = FALSE, inffunc = TRUE
              )
            } else if (est_method == "reg") {
              # regression
              res <- DRDID::reg_did_panel(Ypost, Ypre, G,
                covariates = covariates,
                i.weights = w,
                boot = FALSE, inffunc = TRUE
              )
            } else {
              # doubly robust, this is default
              res <- DRDID::drdid_panel(Ypost, Ypre, G,
                covariates = covariates,
                i.weights = w,
                boot = FALSE, inffunc = TRUE
              )
            }
          }

          # adjust influence function to account for only using
          # subgroup to estimate att(g,t)
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
          res
        }, error = function(e) {
          warning("Error computing internal 2x2 DiD for (g, t) = (", glist[g], ", ", tlist[t + tfac], "): ", e$message, ". The ATT for this cell will be set to NA.")
          NULL
        })

        if (is.null(attgt)) {
          attgt.list[[counter]] <- list(att = NA, group = glist[g], year = tlist[(t + tfac)], post = post.treat)
          inffunc_updates[[update_counter]] <- list(
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
          target_w <- stats::setNames(target_rows$.w, target_rows$.rowid)
          # Look up weight for each observation's unit
          w <- as.numeric(target_w[as.character(disdat$.rowid)])
          # Drop units not observed in the target period
          missing_w <- is.na(w)
          if (any(missing_w)) {
            n_dropped <- length(unique(disdat$.rowid[missing_w]))
            warning(paste0("Dropped ", n_dropped, " units not observed in ",
                           fix_weights, " (period ", target_period, ") ",
                           "for group ", glist[g], " in time period ", tlist[t + tfac]))
            disdat <- disdat[!missing_w, ]
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
          inffunc_updates[[update_counter]] <- list(
            indices = seq_len(n),
            values = rep(NA_real_, n)
          )

          # Update the counters
          update_counter <- update_counter + 1
          counter <- counter + 1
          next
        }


        # matrix of covariates
        covariates <- model.matrix(xformla, data = disdat)

        #-----------------------------------------------------------------------------
        # more checks for enough observations in each group

        # if using custom estimation method, skip this part
        custom_est_method <- is.function(est_method)

        if (!custom_est_method) {
          pscore_problems_likely <- FALSE
          reg_problems_likely <- FALSE

          # checks for pscore based methods
          if (est_method %in% c("dr", "ipw")) {
            preliminary_logit <- fastglm::fastglm(covariates, G, family = binomial())
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
            inffunc_updates[[update_counter]] <- list(
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
              inffunc = TRUE
            ), extra_args))
          } else if (est_method == "ipw") {
            # inverse-probability weights
            res <- DRDID::std_ipw_did_rc(
              y = Y,
              post = post,
              D = G,
              covariates = covariates,
              i.weights = w,
              boot = FALSE, inffunc = TRUE
            )
          } else if (est_method == "reg") {
            # regression
            res <- DRDID::reg_did_rc(
              y = Y,
              post = post,
              D = G,
              covariates = covariates,
              i.weights = w,
              boot = FALSE, inffunc = TRUE
            )
          } else {
            # doubly robust, this is default
            res <- DRDID::drdid_rc(
              y = Y,
              post = post,
              D = G,
              covariates = covariates,
              i.weights = w,
              boot = FALSE, inffunc = TRUE
            )
          }

          # n/n1 adjusts for estimating the
          # att_gt only using observations from groups
          # G and C
          res$att.inf.func <- (n / n1) * res$att.inf.func

          # If ATT is NaN, replace it with NA, and mark influence function as missing
          if (is.nan(res$ATT)) {
            res$ATT <- NA
            res$att.inf.func <- rep(NA_real_, length(res$att.inf.func))
          }
          res
        }, error = function(e) {
          warning("Error computing internal 2x2 DiD for (g, t) = (", glist[g], ", ", tlist[t + tfac], "): ", e$message, ". The ATT for this cell will be set to NA.")
          NULL
        })

        if (is.null(attgt)) {
          attgt.list[[counter]] <- list(att = NA, group = glist[g], year = tlist[(t + tfac)], post = post.treat)
          inffunc_updates[[update_counter]] <- list(
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


      # populate the influence function in the right places
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
        current_ids <- disdat$.rowid[disdat$.G == 1 | disdat$.C == 1]
        aggte_inffunc <- suppressWarnings(stats::aggregate(attgt$att.inf.func, list(current_ids), sum))
        idx <- which(unique(data$.rowid) %in% aggte_inffunc[, 1])
        inffunc_updates[[update_counter]] <- list(
          indices = idx,
          values = aggte_inffunc[, 2]
        )
      }


      # save it in influence function matrix
      # inffunc[g,t,] <- inf.func
      # inffunc[,counter] <- inf.func
      update_counter <- update_counter + 1

      # update counter
      counter <- counter + 1
    } # end looping over t
  } # end looping over g

  # Build the influence function sparse matrix directly from triplets
  trip_i <- integer(0)
  trip_j <- integer(0)
  trip_x <- numeric(0)
  for (j in seq_along(inffunc_updates)) {
    update <- inffunc_updates[[j]]
    idx <- update$indices  # already integer indices
    if (length(idx) > 0L) {
      trip_i <- c(trip_i, idx)
      trip_j <- c(trip_j, rep.int(j, length(idx)))
      trip_x <- c(trip_x, update$values)
    }
  }
  inffunc <- Matrix::sparseMatrix(
    i = trip_i,
    j = trip_j,
    x = trip_x,
    dims = c(inffunc_nrow, length(inffunc_updates))
  )

  return(list(attgt.list = attgt.list, inffunc = inffunc))
}
