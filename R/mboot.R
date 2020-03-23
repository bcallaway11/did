mboot <- function(inf.func, DIDparams) {

  # setup needed variables
  data <- DIDparams$data
  tlist <- DIDparams$tlist
  idname <- DIDparams$idname
  clustervars <- DIDparams$clustervars
  biters <- DIDparams$biters
  
  # just get n obsevations (for clustering below...)
  dta <- data[ data[,tname]==tlist[1], ] 
  
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
    Rb <- sqrt(n)*(apply(Ub*(inffunc1), 2, mean))
    # return bootstrap draw
    Rb
  })
  # bootstrap results
  bres <- t(simplify2array(bout))
  # bootstrap variance matrix 
  V <- cov(bres)

  list(bres=bres, V=V)
}
