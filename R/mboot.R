#' @title Multiplier Bootstrap
#'
#' @description A function to take an influence function and use the
#'  multiplier bootstrap to compute standard errors and critical values for
#'  uniform confidence bands.
#'
#' @param inf.func an influence function
#' @param DIDparams DIDparams object
#' @param pl whether or not to use parallel processing in the multiplier
#'  bootstrap, default=FALSE
#' @param cores the number of cores to use with parallel processing,
#'  default=1
#' @param return_V whether to compute and return the bootstrap variance
#'  matrix `V`. Default is `TRUE`. Internal callers that only consume
#'  `bres`, `se`, or `crit.val` set this to `FALSE` to skip the computation
#'  (it is the only O(biters x k^2) step in the function).
#'
#' @return list with elements
#' \item{bres}{results from each bootstrap iteration}
#' \item{V}{variance matrix (`NULL` when `return_V = FALSE`)}
#' \item{se}{standard errors}
#' \item{crit.val}{a critical value for computing uniform confidence bands}
#'
#' @export
mboot <- function(inf.func, DIDparams, pl = FALSE, cores = 1, return_V = TRUE) {

  # setup needed variables according to faster_mode; This returns different type of objects
  # depending on whether we are in faster_mode or not that has to be handled in the code below
  idname <- DIDparams$idname
  clustervars <- DIDparams$clustervars
  biters <- DIDparams$biters
  tname <- DIDparams$tname
  alp <- DIDparams$alp
  panel <- DIDparams$panel
  true_repeated_cross_sections <- DIDparams$true_repeated_cross_sections
  unbalanced_panel <- DIDparams$allow_unbalanced_panel
  # Reuse the per-unit cluster vector that att_gt() stored in DIDparams when it
  # aligns with the influence-function rows: it was already derived and validated
  # (time invariance) there, so the data load, the aggregate() check, and the
  # unique() derivation below can all be skipped. Fall back to deriving the
  # clusters from the data for (older) DIDparams objects lacking an aligned
  # cluster_vector -- the fallback keeps the time-varying-cluster stop.
  cluster_vector <- DIDparams$cluster_vector
  use_cluster_vector <- !is.null(cluster_vector) &&
    length(cluster_vector) == NROW(inf.func)
  # Avoid full data copy when possible; only needed for clustering without an
  # aligned cluster_vector
  need_data <- length(clustervars) > 0 && !use_cluster_vector
  if (need_data) {
    if (isTRUE(DIDparams$faster_mode) && !is.null(DIDparams$time_invariant_data)) {
      dta <- as.data.frame(DIDparams$time_invariant_data)
    } else {
      data <- as.data.frame(DIDparams$data)
      tlist <- if (isTRUE(DIDparams$faster_mode)) DIDparams$time_periods else sort(unique(data[,tname]))
      if (panel) {
        dta <- data[data[,tname] == tlist[1], ]
      } else {
        dta <- data
      }
    }
  }

  # Convert sparse matrix to dense for bootstrap computation
  inf.func <- as.matrix(inf.func)

  # set correct number of units
  n <- nrow(inf.func)

  # if include id as variable to cluster on
  # drop it as we do this automatically
  if (idname %in% clustervars) {
    clustervars <- clustervars[-which(clustervars==idname)]
  }

  if(!is.null(clustervars)){
    if(is.numeric(clustervars)){
      stop("'clustervars' must be a character string specifying the name of the clustering variable, not a numeric value.")
    }
  }
  # we can only handle up to 2-way clustering
  # (in principle could do more, but not high priority now)
  if (length(clustervars) > 1) {
    stop("At most one cluster variable (beyond 'idname') is supported. Please reduce to one.")
  }

  if (length(clustervars) > 0 && !use_cluster_vector) {
    if(!DIDparams$faster_mode){
      # check that cluster variable does not vary over time within unit
      # reuse 'data' already loaded above (need_data is TRUE on this fallback path)
      clust_tv <- aggregate(data[,clustervars], by=list(data[,idname]), function(rr) length(unique(rr))==1)
      if (!all(clust_tv[,2])) {
        stop("Time-varying cluster variables are not supported. Please provide a time-invariant cluster variable.")
      }
    }
  }

  # multiplier bootstrap
  n_clusters <- n
  if (length(clustervars)==0) {
    bres <- sqrt(n) * run_multiplier_bootstrap(inf.func, biters, pl, cores)
  } else {
    # Cluster-robust multiplier bootstrap following Callaway & Sant'Anna (2021, Remark 10): draw
    # cluster-specific multipliers and apply them to the influence function aggregated to cluster
    # *sums*, consistent with the (1/n) empirical average that defines the estimator. For
    # equal-sized clusters this coincides with the previous aggregation; for unbalanced clusters and
    # repeated cross-sections it keeps the bootstrap aligned with the cluster-robust variance.
    if (use_cluster_vector) {
      # att_gt() already built this vector aligned with the inf.func rows (and
      # validated time invariance when doing so); reusing it avoids reloading
      # the data and re-deriving the clusters on every call.
      cluster <- cluster_vector
      n_clusters <- length(unique(cluster))
    } else {
      n_clusters <- length(unique(dta[,clustervars]))
      cluster <- unique(dta[,c(idname,clustervars)])[,2]
    }
    cluster_sum_if <- rowsum(inf.func, cluster, reorder=TRUE)
    bres <- sqrt(n_clusters) * run_multiplier_bootstrap(cluster_sum_if, biters, pl, cores)
  }


  # handle vector and matrix case differently, so you get nxk matrix
  # ifelse(class(bres)=="matrix", bres <- t(bres), bres <- as.matrix(bres))

  if (isTRUE(class(bres) == "numeric")) bres <- as.matrix(bres)

  # Non-degenerate dimensions
  # ndg.dim <- (base::colSums(bres) != 0)
  ndg.dim <- (!is.na(colSums(bres))) & (base::colSums(bres^2) > sqrt(.Machine$double.eps)*10)
  # If NA, set it to false
  #ndg.dim[is.na(ndg.dim)] <- FALSE
  bres <- as.matrix(bres[ , ndg.dim])

  # bootstrap variance matrix (this matrix can be defective because of degenerate cases).
  # Computed only on request: no internal caller consumes it, and it is the only
  # O(biters * k^2) term in mboot.
  V <- if (return_V) cov(bres) else NULL
  # bootstrap standard error: IQR-based scale normalized to the SD of a normal.
  # The qnorm difference is a constant (hoisted out of the per-column loop), and
  # the type = 1 quantiles are exactly the order statistics x_(ceil(B * p)), so a
  # single sort per column replaces the quantile() calls bit-identically (the
  # post-sort length is the per-column non-NA count, preserving na.rm = TRUE;
  # sort.int() drops NAs by default (na.last = NA)).
  iqr_norm <- qnorm(.75) - qnorm(.25)
  bSigma <- apply(bres, 2, function(b) {
    b <- sort.int(b)
    nb <- length(b)
    (b[ceiling(.75 * nb)] - b[ceiling(.25 * nb)]) / iqr_norm
  })

  # critical value for uniform confidence band: per bootstrap draw, the max over
  # columns of |b / bSigma|. Vectorized as a column-fold of pmax (equivalent to the
  # previous row-wise apply(max, na.rm = TRUE): NA entries are ignored via -Inf, and
  # an all-NA row yields -Inf, which is dropped by the is.finite filter below).
  if (ncol(bres) == 0L) {
    # All bootstrap dimensions degenerate: every row has no finite t-statistic.
    # This matches the previous apply(max, na.rm = TRUE), which returned -Inf for
    # every row and then dropped them all via the is.finite filter.
    bT <- numeric(0)
  } else {
    scaled <- abs(bres / rep(bSigma, each = nrow(bres)))
    scaled[is.na(scaled)] <- -Inf
    bT <- scaled[, 1L]
    for (j in seq_len(ncol(scaled))[-1L]) bT <- pmax(bT, scaled[, j])
    bT <- bT[is.finite(bT)]
  }
  crit.val <- quantile(bT, 1-alp, type=1, na.rm = TRUE)

  #se <- rep(0, length(ndg.dim))
  se <- rep(NA, length(ndg.dim))
  # Standard error from the bootstrap interquartile-range scale. With cluster sums (above) this is the
  # cluster-robust SE (1/n) * sqrt(sum_c S_c^2); without clustering n_clusters = n and it reduces to the
  # i.i.d. form bSigma / sqrt(n). (Equivalent to the previous bSigma / sqrt(n_clusters) when clusters are
  # equal-sized; the change matters only for unbalanced clusters / repeated cross-sections.)
  se[ndg.dim] <- as.numeric(bSigma) * sqrt(n_clusters) / n
  #se[se==0] <- NA

  list(bres = bres, V = V, se = se, crit.val = crit.val)
}

run_multiplier_bootstrap <- function(inf.func, biters, pl = FALSE, cores = 1) {
  ngroups = ceiling(biters/cores)
  chunks = rep(ngroups, cores)
  # Round down so you end up with the right number of biters
  chunks[1] = chunks[1] + biters - sum(chunks)

  n <- nrow(inf.func)
  parallel.function <- function(biters) {
    BMisc::multiplier_bootstrap(inf.func, biters)
  }
  # From tests, this is about where it becomes worth it to parallelize
  if (.Platform$OS.type == "windows" && pl == TRUE && cores > 1) {
    warning("Parallel processing (pl=TRUE) is not supported on Windows. Using sequential processing instead.")
    pl <- FALSE
  }
  if(n > 2500 & pl == TRUE & cores > 1) {
    results = parallel::mclapply(
      chunks,
      FUN = parallel.function,
      mc.cores = cores
    )
    results = do.call(rbind, results)
  } else {
    results = BMisc::multiplier_bootstrap(inf.func, biters)
  }
  return(results)
}
