#' @title compute.aggte
#'
#' @description does the heavy lifting on computing aggregated group-time
#'  average treatment effects
#'
#' @inheritParams att_gt
#'
#' @return \code{AGGTE} object
#'
#' @keywords internal
#'
#' @export
compute.aggte <- function(flist, tlist, group, t, att, first.treat.name, inffunc1, n,
                          clustervars, dta, idname, bstrap, biters, alp, cband, maxe, mine) {

  ## if ( (length(clustervars) > 0) & !bstrap) {
  ##   warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  ## }

  #-----------------------------------------------------------------------------
  # Internal functions for getteing standard errors
  #-----------------------------------------------------------------------------

  # This formula, the argument is the relevant influence function. It return the standard errors
  getSE_inf <- function(thisinffunc) {
    if (bstrap) {
      if (idname %in% clustervars) {
        clustervars <- clustervars[-which(clustervars==idname)]
      }
      if (length(clustervars) > 1) {
        stop("can't handle that many cluster variables")
      }
      bout <- lapply(1:biters, FUN=function(b) {
        if (length(clustervars) > 0) {
          n1 <- length(unique(dta[,clustervars]))
          Vb <- matrix(sample(c(-1,1), n1, replace=T))
          Vb <- cbind.data.frame(unique(dta[,clustervars]), Vb)
          Ub <- data.frame(dta[,clustervars])
          Ub <- Vb[match(Ub[,1], Vb[,1]),]
          Ub <- Ub[,-1]
        } else {
          Ub <- sample(c(-1,1), n, replace=T)
        }
        Rb <- base::mean(Ub*(thisinffunc), na.rm = T)
        Rb
      })
      bres <- as.vector(base::simplify2array(bout))

      bSigma <- (quantile(bres, .75, type=1, na.rm = T) - quantile(bres, .25, type=1, na.rm = T)) /
        (qnorm(.75) - qnorm(.25))
      return(as.numeric(bSigma))
    } else {
      return(sqrt( mean((thisinffunc)^2)/n ))
    }
  }

  # internal function for computing standard errors
  #  this method is used across different types of
  #  aggregate treatment effect parameters and is just
  #  based on using the right influence function and weights
  #  -- these are specific to which aggregate treatment
  #  effect parameter is being considered.
  # @param wif is the influence function for the weights

  getSE <- function(inffunc1, whichones, weights, wif=NULL) {
    # enforce weights are in matrix form
    weights <- as.matrix(weights)
    # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
    thisinffunc <- inffunc1[,whichones]%*%weights
    # Incorporate influence function of the weights
    if (!is.null(wif)) {
      thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
    }
    # Now, compute the standard errror
    getSE_inf(thisinffunc)
  }

  #-----------------------------------------------------------------------------
  # data organization and recoding
  #-----------------------------------------------------------------------------

  # do some recoding to make sure time periods are 1 unit apart
  # and then put these back together at the end
  originalt <- t
  originalgroup <- group
  originalflist <- flist
  originaltlist <- tlist
  uniquet <- seq(1,length(unique(t)))
  # function to switch from "new" t values to  original t values
  t2orig <- function(t) {
    unique(c(originalt,0))[which(c(uniquet,0)==t)]
  }
  # function to switch between "original"
  #  t values and new t values
  orig2t <- function(orig) {
    c(uniquet,0)[which(unique(c(originalt,0))==orig)]
  }
  t <- sapply(originalt, orig2t)
  group <- sapply(originalgroup, orig2t)
  flist <- sapply(originalflist, orig2t)

  # Set the weights
  weights.agg = dta$w

  # some variables used throughout
  # Ever treated only among the units we actually compute the ATT(g,t)
  ever.treated <- 1 * (dta[,first.treat.name]>0)
  mean.w.ever.treated <- mean(weights.agg * ever.treated)

  browser()
  # Probability of being in group g (among the relevant ever-treated groups!)
  pg <- sapply(originalflist,
               function(g) mean(weights.agg * ever.treated * (dta[,first.treat.name]==g))/
                 mean(weights.agg * ever.treated))
  pgg <- pg
  pg <- pg[match(group, flist)] ## make it have the same length as att
  attg <- split(att, group)
  tg <- split(t, group)
  keepers <- which(group <= t)
  G <-  unlist(lapply(dta[,first.treat.name], orig2t))

  #-----------------------------------------------------------------------------
  # Compute the simple ATT summary
  #-----------------------------------------------------------------------------

  # simple att
  simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))
  # Estimation effect coming from P(G=g| Ever treated)
  # Part 1: est effect from P(G=g, ever treated = 1) treating P(ever treated) as known
  simple.oif1 <- sapply(keepers,
                        function(k) {
                          ( (weights.agg *(G==group[k]) -
                               mean(weights.agg * (G==group[k]))) /
                              mean.w.ever.treated
                          )
                        })
  # Part 2: est effect from  P(ever treated) treating P(G=g, ever treated = 1) as known
  #simple.oif2 <- sapply(keepers,
  #                      function(j) mean(weights.agg * (G==group[j])) *
  #                        apply(sapply(keepers, function(k) (weights.agg*(G==group[k]) - mean(weights.agg*(G==group[k])))),1,sum))
  simple.oif2 <- sapply(keepers,
                        function(j) {
                          ((mean(weights.agg * (G==group[j])) /
                              (mean.w.ever.treated^2)) *
                             (weights.agg * ever.treated - mean.w.ever.treated)
                          )
                        })
  
  # Estimation effect from numerator
  simple.oif <- (simple.oif1 - simple.oif2)/(sum(pg[keepers]))
  #Estimation effect from denominator of the weights (normalization)
  simple.oif3 <- base::rowSums(simple.oif) %*%  t(matrix(pg[keepers]/sum(pg[keepers])))
  #Estimation effect from estimated weights (in total)
  simple.oif <- simple.oif - simple.oif3

  # Get standard error for the simple ATT average
  simple.se <- getSE(inffunc1, keepers, pg[keepers]/sum(pg[keepers]), simple.oif)



  #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  # Make sure maximum e is feasible
  if(is.null(maxe)) maxe <- max(t-group)+1
  if (maxe > (max(t-group)+1)) maxe <- max(t-group)+1

  # Sequence of post-treatment e
  eseq <- seq(1,maxe)
  # Get event-study estimators
  dynamic.att.e <- sapply(eseq, function(e) {
    whiche <- which(t - group + 1 == e)
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })

  # Amount of time each unit is treated
  time.treated <- ((max(originalt)) - dta[,first.treat.name] +1) *
                    (dta[,first.treat.name]>0)

  # Get the influence function of the event-study estimators
  dynamic.e.inf.f <- sapply(eseq, function(e) {
    whiche <- which(t - group + 1 == e)
    pge <- pg[whiche]/sum(pg[whiche])

    ## some variables used
    atleast.e.treated <- 1 * (time.treated >=e)
    mean.w.atleast.e.treated <- mean(weights.agg * atleast.e.treated)

    # Estimation effect coming from P(G=g| treated for at least e periods)
    # Part 1: est effect from P(G=g, treatted at least e periods) treating P(treated for at least e periods) as known
    dynamic.oif1 <- sapply(whiche,
                           function(k) {
                             ( (weights.agg *(G==group[k]) -
                                  mean(weights.agg * (G==group[k]))) /
                                 mean.w.atleast.e.treated
                             )
                           })

    # Part 2: est effect from  P(treated for at least e periods) treating P(G=g) as known
    dynamic.oif2 <- sapply(whiche,
                           function(j) ((mean(weights.agg * (G==group[j]))/(mean.w.atleast.e.treated^2)) *
                                          (weights.agg * atleast.e.treated - mean.w.atleast.e.treated)
                           )
    )
    dynamic.oif <- (dynamic.oif1 - dynamic.oif2)

    inffunc1[,whiche] %*% as.matrix(pge) + dynamic.oif %*% as.matrix(att[whiche])

    #getSE(whiche, pge, dynamic.oif )

  })

  # Pre-treatment periods
  pre.treat <- which(group > t)
  if(is.null(mine)) mine <- 0
  if (mine < (min(t-group))) mine <- min(t-group)+1

  eseq.pre <- seq(mine, 0)


  dynamic.att.pre.e <- sapply(eseq.pre, function(e) {
    whiche <- which(t - group + 1 == e)
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })


  n.pre.treat <- (time.treated - max(t)+1) *
    (dta[,first.treat.name]>0)



  dynamic.pre.e.inf.f <- sapply(eseq.pre, function(e) {
    whiche <- which(t - group + 1 == e)
    pge <- pg[whiche]/sum(pg[whiche])

    ## some variables used
    atleast.e.pre.treated <- 1 * (n.pre.treat <=e) * (dta[,first.treat.name]>0)
    mean.w.atleast.e.pre.treated <- mean(weights.agg * atleast.e.pre.treated)

    # Estimation effect coming from P(G=g| treated for at least e periods)
    # Part 1: est effect from P(G=g) treating P(treated for at least e periods) as known
    dynamic.oif1 <- sapply(whiche,
                           function(k) ( (weights.agg *(G==group[k]) - mean(weights.agg * (G==group[k]))) /
                                           mean.w.atleast.e.pre.treated
                           )
    )

    # Part 2: est effect from  P(treated for at least e periods) treating P(G=g) as known
    dynamic.oif2 <- sapply(whiche,
                           function(j) ((mean(weights.agg * (G==group[j]))/(mean.w.atleast.e.pre.treated^2)) *
                                          (weights.agg * atleast.e.pre.treated - mean.w.atleast.e.pre.treated)
                           )
    )
    dynamic.oif <- (dynamic.oif1 - dynamic.oif2)

    inffunc1[,whiche] %*% as.matrix(pge) + dynamic.oif %*% as.matrix(att[whiche])
  })

  dynamic.se.e <- sqrt(base::colMeans(cbind(dynamic.pre.e.inf.f, dynamic.e.inf.f)^2)/n)

  # Average over all dynamic att(e), e>=1
  dynamic.att <- mean(dynamic.att.e)
  # Influence function of the average over all dynamic att(e), e>=1
  dynamic.if <- rowMeans(dynamic.e.inf.f)

  # std. error
  #dynamic.se <- sqrt((mean(dynamic.if^2)/n))
  dynamic.se <- getSE_inf(as.matrix(dynamic.if))

  c.dynamic <- qnorm(1 - alp/2)

  # Bootstrap for simulatanerous Conf. Int for the event study
  if (bstrap) {
    if (idname %in% clustervars) {
      clustervars <- clustervars[-which(clustervars==idname)]
    }
    if (length(clustervars) > 1) {
      stop("can't handle that many cluster variables")
    }
    ## new version
    bout <- lapply(1:biters, FUN=function(b) {

      if (length(clustervars) > 0) {
        n1 <- length(unique(dta[,clustervars]))
        Vb <- matrix(sample(c(-1,1), n1, replace=T))
        Vb <- cbind.data.frame(unique(dta[,clustervars]), Vb)
        Ub <- data.frame(dta[,clustervars])
        Ub <- Vb[match(Ub[,1], Vb[,1]),]
        Ub <- Ub[,-1]
      } else {
        Ub <- sample(c(-1,1), n, replace=T)
      }
      ##Ub <- sample(c(-1,1), n, replace=T)
      Rb <- (base::colMeans(Ub * cbind(dynamic.pre.e.inf.f, dynamic.e.inf.f), na.rm = T))
      Rb
    })
    bres <- t(simplify2array(bout))
    # Non-degenerate dimensions
    ndg.dim <- base::colSums(bres)!=0
    #V.dynamic <- cov(bres)
    bres <- bres[,ndg.dim]

    bSigma <- apply(bres, 2,
                    function(b) (quantile(b, .75, type=1,na.rm = T) -
                                   quantile(b, .25, type=1,na.rm = T))/(qnorm(.75) - qnorm(.25)))

    bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
    c.dynamic <- quantile(bT, 1-alp, type=1,na.rm = T)
    dynamic.se.e <- rep(0,length(ndg.dim))
    dynamic.se.e[ndg.dim] <- as.numeric(bSigma)
  }


  AGGTE(simple.att=simple.att, simple.se=simple.se,
        dynamic.att=dynamic.att, dynamic.se=dynamic.se,
        dynamic.att.e = c(dynamic.att.pre.e, dynamic.att.e), dynamic.se.e=dynamic.se.e,  c.dynamic=c.dynamic,
        group=originalflist,times=originaltlist,
        e = c(eseq.pre, eseq))
}
