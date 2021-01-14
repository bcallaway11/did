#' @title Multiplier Bootstrap
#'
#' @description A function to take an influence function and use the
#'  multiplier bootstrap to compute standard errors and critical values for
#'  uniform confidence bands.
#'
#' @param inf.func an influence function
#' @param DIDparams DIDparams object
#'
#' @return list with elements
#' \item{bres}{results from each bootstrap iteration}
#' \item{V}{variance matrix}
#' \item{se}{standard errors}
#' \item{crit.val}{a critical value for computing uniform confidence bands}
#'
#' @export
mboot <- function(inf.func, DIDparams) {

  # setup needed variables
  data <- DIDparams$data
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
    # CHECK iF CLUSTERVAR is TIME-VARYING
    clust_tv = stats::aggregate(data[,clustervars], list((data[,idname])), sd)
    if(any(clust_tv[,2]>0)){
      stop("can't handle time-varying cluster variables")
    } else if (!panel){
      # IF NOT, SUBSET DTA TO ONE VALUE PER ID
      # Here we do not care about tname and yname as we do not use these
      dta <- stats::aggregate(dta, list((data[,idname])), mean)[,-1]
    }

  }


  # bootstrap
  bout <- lapply(1:biters, FUN=function(b) {
    if (length(clustervars) > 0) {
      # draw Rademachar weights
      # these are the same within clusters
      # see paper for details
      n1 <- length(unique(dta[,clustervars]))
      Vb <- matrix(sample(c(-1,1), n1, replace=T))
      Vb <- cbind.data.frame(unique(dta[,clustervars]), Vb)
      Ub <- data.frame(dta[,clustervars])
      Ub <- Vb[match(Ub[,1], Vb[,1]),]
      Ub <- Ub[,-1]
    } else {
      Ub <- sample(c(-1,1), n, replace=T)
    }
    # multiply weights onto influence function
    Rb <- sqrt(n)*(apply(Ub*(inf.func), 2, mean))
    # return bootstrap draw
    Rb
  })
  # bootstrap results
  bres <- simplify2array(bout)
  # handle vector and matrix case differently, so you get nxk matrix
  ifelse(class(bres)=="matrix", bres <- t(bres), bres <- as.matrix(bres))

  # Non-degenerate dimensions
  # ndg.dim <- (base::colSums(bres) != 0)
  ndg.dim <- (!is.na(colSums(bres))) & (base::colSums(bres^2) > .Machine$double.eps/2)
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
  bT <- apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = T))
  crit.val <- quantile(bT, 1-alp, type=1, na.rm = T)

  se <- rep(0, length(ndg.dim))
  se[ndg.dim] <- as.numeric(bSigma)/sqrt(n)
  #se[se==0] <- NA

  list(bres = bres, V = V, se = se, crit.val = crit.val)
}
