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
compute.aggte <- function(glist, tlist, group, t, att, first.treat.name, inffunc1, n,
                          clustervars, dta, idname, bstrap, biters, alp, cband, maxe, mine) {

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

  getSE <- function(att, inffunc1, whichones, weights, wif=NULL) {
    # enforce weights are in matrix form
    weights <- as.matrix(weights)
    # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
    thisinffunc <- inffunc1[,whichones]%*%weights
    # Incorporate influence function of the weights
    if (!is.null(wif)) {
      thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
    }
    # Now, compute the standard errror
    return(list(se=getSE_inf(thisinffunc),inf.func=thisinffunc))
  }


  
  #-----------------------------------------------------------------------------
  # data organization and recoding
  #-----------------------------------------------------------------------------

  # do some recoding to make sure time periods are 1 unit apart
  # and then put these back together at the end
  originalt <- t
  originalgroup <- group
  originalglist <- glist
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
  glist <- sapply(originalglist, orig2t)

  # Set the weights
  weights.agg = dta$w

  # some variables used throughout
  # Ever treated only among the units we actually compute the ATT(g,t)
  ever.treated <- 1 * (dta[,first.treat.name]>0)
  mean.w.ever.treated <- mean(weights.agg * ever.treated)

  
  # Probability of being in group g (among the relevant ever-treated groups!)
  ## pg <- sapply(originalflist,
  ##              function(g) {
  ##                mean(weights.agg * ever.treated * (dta[,first.treat.name]==g)) /
  ##                  mean(weights.agg * ever.treated)
  ##              })
  
  # we can work in overall probabilities because conditioning will cancel out
  # cause it shows up in numerator and denominator
  pg <- sapply(originalglist, function(g) mean(weights.agg*dta[,first.treat.name]==g))
  pgg <- pg
  pg <- pg[match(group, glist)] ## make it have the same length as att
  attg <- split(att, group)
  tg <- split(t, group)
  keepers <- which(group <= t)
  G <-  unlist(lapply(dta[,first.treat.name], orig2t))

  #-----------------------------------------------------------------------------
  # Compute the simple ATT summary
  #-----------------------------------------------------------------------------

  ## # simple att
  ## simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))
  ## # Estimation effect coming from P(G=g| Ever treated)
  ## # Part 1: est effect from P(G=g, ever treated = 1) treating P(ever treated) as known
  ## simple.oif1 <- sapply(keepers,
  ##                       function(k) {
  ##                         ( (weights.agg *(G==group[k]) -
  ##                              mean(weights.agg * (G==group[k]))) /
  ##                             mean.w.ever.treated
  ##                         )
  ##                       })
  ## # Part 2: est effect from  P(ever treated) treating P(G=g, ever treated = 1) as known
  ## #simple.oif2 <- sapply(keepers,
  ## #                      function(j) mean(weights.agg * (G==group[j])) *
  ## #                        apply(sapply(keepers, function(k) (weights.agg*(G==group[k]) - mean(weights.agg*(G==group[k])))),1,sum))
  ## simple.oif2 <- sapply(keepers,
  ##                       function(j) {
  ##                         ((mean(weights.agg * (G==group[j])) /
  ##                             (mean.w.ever.treated^2)) *
  ##                            (weights.agg * ever.treated - mean.w.ever.treated)
  ##                         )
  ##                       })
  
  ## # Estimation effect from numerator
  ## simple.oif <- (simple.oif1 - simple.oif2)/(sum(pg[keepers]))
  ## #Estimation effect from denominator of the weights (normalization)
  ## simple.oif3 <- base::rowSums(simple.oif) %*%  t(matrix(pg[keepers]/sum(pg[keepers])))
  ## #Estimation effect from estimated weights (in total)
  ## simple.oif <- simple.oif - simple.oif3


  ## # Get standard error for the simple ATT average
  ## simple.se <- getSE(inffunc1, keepers, pg[keepers]/sum(pg[keepers]), simple.oif)

  
  # simple att
  simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))

  # account for having to estimate the weights in converting att(g,t)
  # effect of estimating weights in the numerator
  simple.oif1 <- sapply(keepers, function(k) {
    (weights.agg * 1*(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  simple.oif2 <- rowSums( sapply( keepers, function(k) {
    weights.agg*1*(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2)) 
  ## # effect of estimating weights in the denominator
  ## simple.oif2 <- rowSums( sapply( keepers, function(k) {
  ##   weights.agg*1*(G==group[k]) - mean(weights.agg*1*(G==group[k]))
  ## })) %*%
  ##   t(pg[keepers]/(sum(pg[keepers])^2)) 
  # combine everything using getSE function
  simple.se <- getSE(att, inffunc1, keepers, pg[keepers]/sum(pg[keepers]), simple.oif1 - simple.oif2)$se

  #-----------------------------------------------------------------------------
  # Compute the selective treatment timing estimators
  #-----------------------------------------------------------------------------


  # function to compute extra term
  # in influence function due to estimating the weights
  wif <- function(keepers, pg, weights.agg, G) {
    if1 <- sapply(keepers, function(k) {
      (weights.agg * 1*(G==group[k]) - pg[k]) /
        sum(pg[keepers])
    })
    # effect of estimating weights in the denominator
    if2 <- rowSums( sapply( keepers, function(k) {
      weights.agg*1*(G==group[k]) - pg[k]
    })) %*%
      t(pg[keepers]/(sum(pg[keepers])^2))

    if1 - if2
  }
  
  ## Selective Treatment Timing
  ## Note: for selective.att.g, don't need to adjust standard
  ##  errors for estimating weights because they are known
  selective.att.g <- sapply(glist, function(g) {
    whichg <- which( (group == g) & (g <= t))
    attg <- att[whichg]
    mean(attg)
  })
  selective.se.inner <- lapply(glist, function(g) {
        whichg <- which( (group == g) & (g <= t))
        getSE(att,inffunc1,whichg, pg[whichg]/sum(pg[whichg]))
  })
  
  selective.se.g <- unlist(getListElement(selective.se.inner, "se"))
  selective.inf.func.g <- simplify2array(getListElement(selective.se.inner, "inf.func"))[,1,]
  selective.att <- sum(selective.att.g * pgg)/sum(pgg)
  selective.wif <- wif(1:length(glist), pgg, weights.agg, G)
  selective.se <- getSE(selective.att.g, selective.inf.func.g, 1:length(glist), pgg/sum(pgg), selective.wif)$se
    
  #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  eseq <- unique(t-group) #seq(0,max(t-group)) # fix this to cover pre-treatment periods
  eseq <- eseq[order(eseq)]
  
  dynamic.att.e <- sapply(eseq, function(e) {
    whiche <- which(t - group == e)
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })

  dynamic.se.inner <- lapply(eseq, function(e) {
    whiche <- which(t - group == e)
    pge <- pg[whiche]/(sum(pg[whiche]))
    wif.e <- wif(whiche, pg, weights.agg, G)
    getSE(att, inffunc1, whiche, pge, wif.e)
  })

  dynamic.se.e <- unlist(getListElement(dynamic.se.inner, "se"))
  dynamic.inf.func.e <- simplify2array(getListElement(dynamic.se.inner, "inf.func"))[,1,]

  # get overall average treatment effect
  # by averaging over positive dynamics
  epos <- eseq >= 0
  dynamic.att <- mean(dynamic.att.e[epos])
  dynamic.se  <- getSE(att=dynamic.att.e[epos],
                       inffunc1=dynamic.inf.func.e[,epos],
                       whichones=(1:sum(epos)),
                       weights=(rep(1/sum(epos), sum(epos))))$se


  #-----------------------------------------------------------------------------
  # calendar time effects
  #-----------------------------------------------------------------------------


  
  browser()

  dynamic.se.e <- sapply(eseq, function(e) {
        whiche <- which(t - group + 1 == e)
        pge <- pg[whiche]/(sum(pg[whiche]))
        dynamic.oif1 <- sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k]))) / sum(pg[whiche]))
        dynamic.oif2 <- sapply(whiche, function(j) mean(1*(G==group[j])) * apply(sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k])))),1,sum))
        getSE(whiche, pge, dynamic.oif1 - dynamic.oif2)
    })

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
        group=originalglist,times=originaltlist,
        e = c(eseq.pre, eseq))
}



