#' @title Multiplier Bootstrap
#'
#' @description A function to take an influence function and use the
#'  multiplier bootsrap to compute standard errors and critical values for
#'  uniform confidence bands.
#'
#' @param inf.func an influence function
#' @param DIDparams DIDparams object
#'
#' @return list with elements
#' * bres returns from each bootstrap iteration
#' * V variance matrix
#' * se standard errors
#' * crit.val a critical value for computing uniform confidence bands
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
  # Make sure inf.func is matrix because we need this for computing n below
  inf.func <- as.matrix(inf.func)
  # just get n obsevations (for clustering below...)
  ifelse(panel,
         dta <- data[ data[,tname]==tlist[1], ],
         dta <- data)

  n <- nrow(inf.func) # this adjusts automatically to panel vs. repeated cross sections

  # if include id as variable to cluster on
  # drop it as we do this automatically
  if (idname %in% clustervars) {
    clustervars <- clustervars[-which(clustervars==idname)]
  }

  # we can only handle up to 2-way clustering
  # (in principle could do more, but not high priority now)
  if (length(clustervars) > 1) {
    stop("can't handle that many cluster variables")
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
  ndg.dim <- (base::colSums(bres) != 0)
  # If NA, set it to false
  ndg.dim[is.na(ndg.dim)] <- FALSE
  bres <- as.matrix(bres[ , ndg.dim])

  # bootstrap variance matrix (this matrix can be defective because of degenerate cases)
  V <- cov(bres)
  # bootstrap standard error
  bSigma <- apply(bres, 2,
                  function(b) (quantile(b, .75, type=1, na.rm = T) -
                                 quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))

  # critical value for uniform confidence band
  bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
  crit.val <- quantile(bT, 1-alp, type=1, na.rm = T)

  se <- rep(0, length(ndg.dim))
  se[ndg.dim] <- as.numeric(bSigma)/sqrt(n)
  #se[se==0] <- NA

  list(bres = bres, V = V, se = se, crit.val = crit.val)
}
