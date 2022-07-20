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
#'
#' @return list with elements
#' \item{bres}{results from each bootstrap iteration}
#' \item{V}{variance matrix}
#' \item{se}{standard errors}
#' \item{crit.val}{a critical value for computing uniform confidence bands}
#'
#' @export
mboot <- function(inf.func, DIDparams, pl = FALSE, cores = 1) {

  # setup needed variables
  data <- as.data.frame(DIDparams$data)
  idname <- DIDparams$idname
  clustervars <- DIDparams$clustervars
  biters <- DIDparams$biters
  tname <- DIDparams$tname
  tlist <- unique(data[,tname])[order(unique(data[,tname]))]
  alp <- DIDparams$alp
  panel <- DIDparams$panel
  true_repeated_cross_sections <- DIDparams$true_repeated_cross_sections

  # just get n obsevations (for clustering below...)
  ifelse(panel,
         dta <- data[ data[,tname]==tlist[1], ],
         dta <- data)
  # Make sure inf.func is matrix because we need this for computing n below
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
      stop("clustervars need to be the name of the clustering variable.")
    }
  }
  # we can only handle up to 2-way clustering
  # (in principle could do more, but not high priority now)
  if (length(clustervars) > 1) {
    stop("can't handle that many cluster variables")
  }

  if (length(clustervars) > 0) {
    # check that cluster variable does not vary over time within unit
    clust_tv <- aggregate(data[,clustervars], by=list(data[,idname]), function(rr) length(unique(rr))==1)
    if (!all(clust_tv[,2])) {
      stop("can't handle time-varying cluster variables")
    }
    ## # CHECK iF CLUSTERVAR is TIME-VARYING
    ## clust_tv = base::suppressWarnings(stats::aggregate(data[,clustervars], list((data[,idname])), sd))
    ## clust_tv$x[is.na(clust_tv$x)] <- 0
    ## if(any(clust_tv[,2]>.Machine$double.eps)){
    ##   stop("can't handle time-varying cluster variables")
    ## } else if (!panel){
    ##   # IF NOT, SUBSET DTA TO ONE VALUE PER ID
    ##   # Here we do not care about tname and yname as we do not use these
    ##   dta <- base::suppressWarnings(stats::aggregate(dta, list((data[,idname])), mean)[,-1])
    ## }

  }

  # multiplier bootstrap
  n_clusters <- n
  if (length(clustervars)==0) {
    bres <- sqrt(n) * run_multiplier_bootstrap(inf.func, biters, pl, cores)
  } else {
    n_clusters <- length(unique(data[,clustervars]))
    cluster <- unique(dta[,c(idname,clustervars)])[,2]
    cluster_n <- aggregate(cluster, by=list(cluster), length)[,2]
    cluster_mean_if <- rowsum(inf.func, cluster,reorder=TRUE) / cluster_n
    bres <- sqrt(n_clusters) * run_multiplier_bootstrap(cluster_mean_if, biters, pl, cores)
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

  # bootstrap variance matrix (this matrix can be defective because of degenerate cases)
  V <- cov(bres)
  # bootstrap standard error
  bSigma <- apply(bres, 2,
                  function(b) (quantile(b, .75, type=1, na.rm = T) -
                                 quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))

  # critical value for uniform confidence band
  bT <- base::suppressWarnings(apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = T)))
  bT <- bT[is.finite(bT)]
  crit.val <- quantile(bT, 1-alp, type=1, na.rm = T)

  #se <- rep(0, length(ndg.dim))
  se <- rep(NA, length(ndg.dim))
  se[ndg.dim] <- as.numeric(bSigma) / sqrt(n_clusters)
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
