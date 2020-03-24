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
compute.aggte <- function(MP, type="simple", balance.e=NULL) {

  #-----------------------------------------------------------------------------
  # process MP object
  #-----------------------------------------------------------------------------
  # load parameters
  group <- MP$group
  t <- MP$t
  att <- MP$att
  dp <- MP$DIDparams
  first.treat.name <- dp$first.treat.name
  inffunc1 <- MP$inffunc
  n <- MP$n
  clustervars <- dp$clustervars
  data <- dp$data
  tname <- dp$tname
  idname <- dp$idname
  bstrap <- dp$bstrap
  biters <- dp$biters
  alp <- dp$alp
  cband <- dp$cband
  maxe <- dp$maxe
  mine <- dp$mine
  # figure out the dates
  # list of dates from smallest to largest
  tlist <- unique(data[,tname])[order(unique(data[,tname]))] 
  # list of treated groups (by time) from smallest to largest
  glist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
  # Only the treated groups
  glist <- glist[glist>0]
  
  # data from first period
  dta <- data[ data[,tname]==tlist[1], ]
  dta$w <- 1 # TODO: fix this

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
  maxT <- max(t)
  
  # Set the weights
  weights.ind  <-  dta$w

  # some variables used throughout
  # Ever treated only among the units we actually compute the ATT(g,t)
  ever.treated <- 1 * (dta[,first.treat.name]>0)
  mean.w.ever.treated <- mean(weights.ind * ever.treated)

  
  # we can work in overall probabilities because conditioning will cancel out
  # cause it shows up in numerator and denominator
  pg <- sapply(originalglist, function(g) mean(weights.ind*dta[,first.treat.name]==g))
  pgg <- pg
  pg <- pg[match(group, glist)] ## make it have the same length as att
  attg <- split(att, group)
  tg <- split(t, group)
  keepers <- which(group <= t)
  G <-  unlist(lapply(dta[,first.treat.name], orig2t))

  #-----------------------------------------------------------------------------
  # Compute the simple ATT summary
  #-----------------------------------------------------------------------------

  if (type == "simple") {
    
    # simple att
    simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))
    simple.wif <- wif(keepers, pg, weights.ind, G, group)
    simple.if <- get_agg_inf_func(att=att,
                                  inffunc1=inffunc1,
                                  whichones=keepers,
                                  weights.agg=pg[keepers]/sum(pg[keepers]),
                                  wif=simple.wif)
    simple.se <- getSE(simple.if, dp)

    return(AGGTEobj(overall.att=simple.att, overall.se=simple.se, type=type))
  }

  #-----------------------------------------------------------------------------
  # Compute the selective treatment timing estimators
  #-----------------------------------------------------------------------------

  if (type == "selective") {
  
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
      inf.func.g <- get_agg_inf_func(att=att,
                                     inffunc1=inffunc1,
                                     whichones=whichg,
                                     weights.agg=pg[whichg]/sum(pg[whichg]),
                                     wif=NULL)
      se.g <- getSE(inf.func.g, dp)
      list(inf.func=inf.func.g, se=se.g)
    })
    
    selective.se.g <- unlist(getListElement(selective.se.inner, "se"))
    selective.inf.func.g <- simplify2array(getListElement(selective.se.inner, "inf.func"))[,1,]
    selective.crit.val <- mboot(selective.inf.func.g, dp)$crit.val
    selective.att <- sum(selective.att.g * pgg)/sum(pgg)
    selective.wif <- wif(keepers=1:length(glist),
                         pg=pgg,
                         weights.ind=weights.ind,
                         G=G,
                         group=group)
    selective.inf.func <- get_agg_inf_func(att=selective.att.g,
                                           inffunc1=selective.inf.func.g,
                                           whichones=(1:length(glist)),
                                           weights.agg=pgg/sum(pgg),
                                           wif=selective.wif)
    selective.se <- getSE(selective.inf.func, dp)
    
    return(AGGTEobj(overall.att=selective.att,
                    overall.se=selective.se,
                    type=type,
                    egt=originalglist,
                    att.egt=selective.att.g,
                    se.egt=selective.se.g,
                    crit.val.egt=selective.crit.val))

  }


  
  #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  if (type == "dynamic") {
    
    eseq <- unique(t-group) 
    eseq <- eseq[order(eseq)]

    # if we balance the sample with resepect to event time
    if (!is.null(balance.e)) { 
      eseq <- eseq[ (eseq <= balance.e) & (eseq >= balance.e - maxT + 1)]
    }

    # these are not currently used, but if we want to trim
    # out some lengths of exposure, we could use this
    # eseq <- eseq[ (eseq >= mine) & (eseq <= maxe) ]    
    # note that they would still be included in estimating overall effects
    
    dynamic.att.e <- sapply(eseq, function(e) {
      whiche <- which( (t - group == e) & (maxT - group >= balance.e) ) 
      atte <- att[whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(atte*pge)
    })

    dynamic.se.inner <- lapply(eseq, function(e) {
      whiche <- which( (t - group == e) & (maxT - group >= balance.e) ) 
      pge <- pg[whiche]/(sum(pg[whiche]))
      wif.e <- wif(whiche, pg, weights.ind, G, group)
      inf.func.e <- get_agg_inf_func(att=att,
                                     inffunc1=inffunc1,
                                     whichones=whiche,
                                     weights.agg=pge,
                                     wif=wif.e)
      se.e <- getSE(inf.func.e, dp)
      list(inf.func=inf.func.e, se=se.e)
    })

    dynamic.se.e <- unlist(getListElement(dynamic.se.inner, "se"))
    dynamic.inf.func.e <- simplify2array(getListElement(dynamic.se.inner, "inf.func"))[,1,]
    dynamic.crit.val <- mboot(dynamic.inf.func.e, dp)$crit.val

    # get overall average treatment effect
    # by averaging over positive dynamics
    epos <- eseq >= 0
    dynamic.att <- mean(dynamic.att.e[epos])
    dynamic.inf.func <- get_agg_inf_func(att=dynamic.att.e[epos],
                                         inffunc1=dynamic.inf.func.e[,epos],
                                         whichones=(1:sum(epos)),
                                         weights.agg=(rep(1/sum(epos), sum(epos))),
                                         wif=NULL)
    dynamic.se <- getSE(dynamic.inf.func, dp)
    
    return(AGGTEobj(overall.att=dynamic.att,
                    overall.se=dynamic.se,
                    type=type,
                    egt=eseq,
                    att.egt=dynamic.att.e,
                    se.egt=dynamic.se.e,
                    crit.val.egt=dynamic.crit.val))
  }

  #-----------------------------------------------------------------------------
  # calendar time effects
  #-----------------------------------------------------------------------------


  # TODO...
  
  
}

#-----------------------------------------------------------------------------
# Internal functions for getteing standard errors
#-----------------------------------------------------------------------------


# internal function for computing standard errors
#  this method is used across different types of
#  aggregate treatment effect parameters and is just
#  based on using the right influence function and weights
#  -- these are specific to which aggregate treatment
#  effect parameter is being considered.
# @param wif is the influence function for the weights

get_agg_inf_func <- function(att, inffunc1, whichones, weights.agg, wif=NULL) {
  # enforce weights are in matrix form
  weights.agg <- as.matrix(weights.agg)
  # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
  thisinffunc <- inffunc1[,whichones]%*%weights.agg
  # Incorporate influence function of the weights
  if (!is.null(wif)) {
    thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
  }
  # Now, compute the standard errror
  return(thisinffunc)
}


# This formula, the argument is the relevant influence function. It return the standard errors
getSE <- function(thisinffunc, DIDparams=NULL) {
  alp <- .05
  bstrap <- FALSE
  if (!is.null(DIDparams)) {
    bstrap <- DIDparams$bstrap
    alp <- DIDparams$alp
    cband <- DIDparams$cband
    n <- length(thisinffunc)
  }
  
  if (bstrap) {
    bout <- mboot(thisinffunc, DIDparams)
    return(bout$se)
  } else {
    return(sqrt( mean((thisinffunc)^2)/n ))
  }
}

# function to compute extra term
# in influence function due to estimating the weights
wif <- function(keepers, pg, weights.ind, G, group) {
  if1 <- sapply(keepers, function(k) {
    (weights.ind * 1*(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  if2 <- rowSums( sapply( keepers, function(k) {
    weights.ind*1*(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2))

  if1 - if2
}
