#' @title Compute Aggregated Treatment Effect Parameters
#'
#' @description Does the heavy lifting on computing aggregated group-time
#'  average treatment effects
#'
#' @inheritParams att_gt
#' @inheritParams aggte
#' @param call The function call to aggte
#'
#' @return \code{\link{AGGTEobj}} object
#'
#' @keywords internal
#'
#' @export
compute.aggte <- function(MP,
                          type = "group",
                          balance_e = NULL,
                          min_e = -Inf,
                          max_e = Inf,
                          na.rm = FALSE,
                          bstrap = NULL,
                          biters = NULL,
                          cband = NULL,
                          alp = NULL,
                          clustervars = NULL,
                          call = NULL) {

  #-----------------------------------------------------------------------------
  # unpack MP object
  #-----------------------------------------------------------------------------
  # load parameters
  group <- MP$group
  t <- MP$t
  att <- MP$att
  dp <- MP$DIDparams
  inffunc1 <- MP$inffunc
  n <- MP$n


  gname <- dp$gname
  data <- dp$data
  tname <- dp$tname
  idname <- dp$idname
  if(is.null(clustervars)){
    clustervars <- dp$clustervars
  }
  if(is.null(bstrap)){
    bstrap <- dp$bstrap
  }
  if(is.null(biters)){
    biters <- dp$biters
  }
  if(is.null(alp)){
    alp <- dp$alp
  }
  if(is.null(cband)){
    cband <- dp$cband
  }

  tlist <- dp$tlist
  glist <- dp$glist
  panel <- dp$panel

  # overwrite MP objects (so we can actually compute bootstrap)
  MP$DIDparams$clustervars <- clustervars
  MP$DIDparams$bstrap <- bstrap
  MP$DIDparams$biters <- biters
  MP$DIDparams$alp <- alp
  MP$DIDparams$cband <- cband
  dp <- MP$DIDparams


  if(na.rm){
    notna <- !is.na(att)
    group <- group[notna]
    t <- t[notna]
    att <- att[notna]
    inffunc1 <- inffunc1[, notna]
    #tlist <- sort(unique(t))
    glist <- sort(unique(group))
  }

  if((na.rm == FALSE) && base::anyNA(att)) stop("Missing values at att_gt found. If you want to remove these, set `na.rm = TRUE'.")

  # data from first period
  ifelse(panel,
         dta <- data[ data[,tname]==tlist[1], ],
         dta <- data
         )
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
  tlist <- unique(t)
  maxT <- max(t)

  # Set the weights
  weights.ind  <-  dta$w

  # we can work in overall probabilities because conditioning will cancel out
  # cause it shows up in numerator and denominator
  pg <- sapply(originalglist, function(g) mean(weights.ind*(dta[,gname]==g)))

  # length of this is equal to number of groups
  pgg <- pg

  # same but length is equal to the number of ATT(g,t)
  pg <- pg[match(group, glist)]

  # which group time average treatment effects are post-treatment
  keepers <- which(group <= t)

  # n x 1 vector of group variable
  G <-  unlist(lapply(dta[,gname], orig2t))

  #-----------------------------------------------------------------------------
  # Compute the simple ATT summary
  #-----------------------------------------------------------------------------

  if (type == "simple") {

    # simple att
    # averages all post-treatment ATT(g,t) with weights
    # given by group size
    simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))

    # get the part of the influence function coming from estimated weights
    simple.wif <- wif(keepers, pg, weights.ind, G, group)

    # get the overall influence function
    simple.if <- get_agg_inf_func(att=att,
                                  inffunc1=inffunc1,
                                  whichones=keepers,
                                  weights.agg=pg[keepers]/sum(pg[keepers]),
                                  wif=simple.wif)

    # get standard errors from overall influence function
    simple.se <- getSE(simple.if, dp)

    return(AGGTEobj(overall.att = simple.att,
                    overall.se = simple.se,
                    type = type,
                    inf.function = list(simple.att = simple.if),
                    call=call,
                    DIDparams=dp))
  }

  #-----------------------------------------------------------------------------
  # Compute the group (i.e., selective) treatment timing estimators
  #-----------------------------------------------------------------------------

  if (type == "group") {

    # get group specific ATTs
    # note: there are no estimated weights here
    selective.att.g <- sapply(glist, function(g) {
      # look at post-treatment periods for group g
      whichg <- which( (group == g) & (g <= t))
      attg <- att[whichg]
      mean(attg)
    })

    # get standard errors for each group specific ATT
    selective.se.inner <- lapply(glist, function(g) {
      whichg <- which( (group == g) & (g <= t))
      inf.func.g <- as.numeric(get_agg_inf_func(att=att,
                                     inffunc1=inffunc1,
                                     whichones=whichg,
                                     weights.agg=pg[whichg]/sum(pg[whichg]),
                                     wif=NULL))
      se.g <- getSE(inf.func.g, dp)
      list(inf.func=inf.func.g, se=se.g)
    })

    # recover standard errors separately by group
    selective.se.g <- unlist(BMisc::getListElement(selective.se.inner, "se"))

    # recover influence function separately by group
    selective.inf.func.g <- simplify2array(BMisc::getListElement(selective.se.inner, "inf.func"))

    # use multiplier bootstrap (across groups) to get critical value
    # for constructing uniform confidence bands
    selective.crit.val <- NULL
    if(dp$cband==TRUE){
      selective.crit.val <- mboot(selective.inf.func.g, dp)$crit.val
    }

    # get overall att under selective treatment timing
    # (here use pgg instead of pg because we can just look at each group)
    selective.att <- sum(selective.att.g * pgg)/sum(pgg)

    # account for having to estimate pgg in the influence function
    selective.wif <- wif(keepers=1:length(glist),
                         pg=pgg,
                         weights.ind=weights.ind,
                         G=G,
                         group=group)

    # get overall influence function
    selective.inf.func <- get_agg_inf_func(att=selective.att.g,
                                           inffunc1=selective.inf.func.g,
                                           whichones=(1:length(glist)),
                                           weights.agg=pgg/sum(pgg),
                                           wif=selective.wif)

    # get overall standard error
    selective.se <- getSE(selective.inf.func, dp)

    return(AGGTEobj(overall.att=selective.att,
                    overall.se=selective.se,
                    type=type,
                    egt=originalglist,
                    att.egt=selective.att.g,
                    se.egt=selective.se.g,
                    crit.val.egt=selective.crit.val,
                    inf.function = list(selective.inf.func.g = selective.inf.func.g,
                                        selective.inf.func = selective.inf.func),
                    call=call,
                    DIDparams=dp))

  }



  #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  if (type == "dynamic") {

    # event times
    # this looks at all available event times
    # note: event times can be negative here.
    # note: event time = 0 corresponds to "on impact"
    #eseq <- unique(t-group)
    eseq <- unique(originalt - originalgroup)
    eseq <- eseq[order(eseq)]

    # if the user specifies balance_e, then we are going to
    # drop some event times and some groups; if not, we just
    # keep everything (that is what this variable is for)
    include.balanced.gt <- rep(TRUE, length(originalgroup))

    # if we balance the sample with resepect to event time
    if (!is.null(balance_e)) {
      eseq <- eseq[ (eseq <= balance_e) & (eseq >= balance_e - t2orig(maxT) + t2orig(1))]
      include.balanced.gt <- (t2orig(maxT) - originalgroup >= balance_e)
    }

    # only looks at some event times 
     eseq <- eseq[ (eseq >= min_e) & (eseq <= max_e) ]

    # compute atts that are specific to each event time
    dynamic.att.e <- sapply(eseq, function(e) {
      # keep att(g,t) for the right g&t as well as ones that
      # are not trimmed out from balancing the sample
      whiche <- which( (originalt - originalgroup == e) & (include.balanced.gt) )
      atte <- att[whiche]
      pge <- pg[whiche]/(sum(pg[whiche]))
      sum(atte*pge)
    })

    # compute standard errors for dynamic effects
    dynamic.se.inner <- lapply(eseq, function(e) {
      whiche <- which( (originalt - originalgroup == e) & (include.balanced.gt) )
      pge <- pg[whiche]/(sum(pg[whiche]))
      wif.e <- wif(whiche, pg, weights.ind, G, group)
      inf.func.e <- as.numeric(get_agg_inf_func(att=att,
                                     inffunc1=inffunc1,
                                     whichones=whiche,
                                     weights.agg=pge,
                                     wif=wif.e))
      se.e <- getSE(inf.func.e, dp)
      list(inf.func=inf.func.e, se=se.e)
    })

    dynamic.se.e <- unlist(BMisc::getListElement(dynamic.se.inner, "se"))
    dynamic.inf.func.e <- simplify2array(BMisc::getListElement(dynamic.se.inner, "inf.func"))

    dynamic.crit.val <- NULL
    if(dp$cband==TRUE){
      dynamic.crit.val <- mboot(dynamic.inf.func.e, dp)$crit.val
    }

    # get overall average treatment effect
    # by averaging over positive dynamics
    epos <- eseq >= 0
    dynamic.att <- mean(dynamic.att.e[epos])
    dynamic.inf.func <- get_agg_inf_func(att=dynamic.att.e[epos],
                                         inffunc1=as.matrix(dynamic.inf.func.e[,epos]),
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
                    crit.val.egt=dynamic.crit.val,
                    inf.function = list(dynamic.inf.func.e = dynamic.inf.func.e,
                                        dynamic.inf.func = dynamic.inf.func),
                    call=call,
                    min_e=min_e,
                    max_e=max_e,
                    balance_e=balance_e,
                    DIDparams=dp
                    ))
  }

  #-----------------------------------------------------------------------------
  # calendar time effects
  #-----------------------------------------------------------------------------

  if (type == "calendar") {

    # drop time periods where no one is treated yet
    # (can't get treatment effects in those periods)
    minG <- min(group)
    calendar.tlist <- tlist[tlist>=minG]

    # calendar time specific atts
    calendar.att.t <- sapply(calendar.tlist, function(t1) {
      # look at post-treatment periods for group g
      whicht <- which( (t == t1) & (group <= t))
      attt <- att[whicht]
      mean(attt)
    })

    # get standard errors and influence functions
    # for each time specific att
    calendar.se.inner <- lapply(calendar.tlist, function(t1) {
      whicht <- which( (t == t1) & (group <= t))
      wif.t <- wif(keepers=whicht,
                   pg=pg,
                   weights.ind=weights.ind,
                   G=G,
                   group=group)
      inf.func.t <- as.numeric(get_agg_inf_func(att=att,
                                     inffunc1=inffunc1,
                                     whichones=whicht,
                                     weights.agg=pg[whicht]/sum(pg[whicht]),
                                     wif=wif.t))
      se.t <- getSE(inf.func.t, dp)
      list(inf.func=inf.func.t, se=se.t)
    })

    # recover standard errors separately by time
    calendar.se.t <- unlist(BMisc::getListElement(calendar.se.inner, "se"))

    # recover influence function separately by time
    calendar.inf.func.t <- simplify2array(BMisc::getListElement(calendar.se.inner, "inf.func"))

    # use multiplier boostrap (across groups) to get critical value
    # for constructing uniform confidence bands
    calendar.crit.val <- NULL
    if(dp$cband==TRUE){
      calendar.crit.val <- mboot(calendar.inf.func.t, dp)$crit.val
    }

    # get overall att under calendar time effects
    # this is just average over all time periods
    calendar.att <- mean(calendar.att.t)

    # get overall influence function
    calendar.inf.func <- get_agg_inf_func(att=calendar.att.t,
                                           inffunc1=calendar.inf.func.t,
                                           whichones=(1:length(calendar.tlist)),
                                           weights.agg=rep(1/length(calendar.tlist), length(calendar.tlist)),
                                           wif=NULL)

    # get overall standard error
    calendar.se <- getSE(calendar.inf.func, dp)

    return(AGGTEobj(overall.att=calendar.att,
                    overall.se=calendar.se,
                    type=type,
                    egt=sapply(calendar.tlist,t2orig),
                    att.egt=calendar.att.t,
                    se.egt=calendar.se.t,
                    crit.val.egt=calendar.crit.val,
                    inf.function = list(calendar.inf.func.t = calendar.inf.func.t,
                                        calendar.inf.func = calendar.inf.func),
                    call=call,
                    DIDparams=dp
                    ))

  }


}

#-----------------------------------------------------------------------------
# Internal functions for getteing standard errors
#-----------------------------------------------------------------------------

#' @title Compute extra term in influence function due to estimating weights
#'
#' @description A function to compute the extra term that shows up in the
#'  influence function for aggregated treatment effect parameters
#'  due to estimating the weights
#'
#' @param keepers a vector of indices for which group-time average
#'  treatment effects are used to compute a particular aggregated parameter
#' @param pg a vector with same length as total number of group-time average
#'  treatment effects that contains the probability of being in particular group
#' @param weights.ind additional sampling weights (nx1)
#' @param G vector containing which group a unit belongs to (nx1)
#' @param group vector of groups
#'
#' @return nxk influence function matrix
#'
#' @keywords internal
wif <- function(keepers, pg, weights.ind, G, group) {
  # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
  # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights.ind * 1*(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  if2 <- rowSums( sapply( keepers, function(k) {
    weights.ind*1*(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2))

  # return the influence function for the weights
  if1 - if2
}


#' @title Get an influence function for particular aggregate parameters
#'
#' @title This is a generic internal function for combining influence
#'  functions across ATT(g,t)'s to return an influence function for
#'  various aggregated treatment effect parameters.
#'
#' @param att vector of group-time average treatment effects
#' @param inffunc1 influence function for all group-time average treatment effects
#'  (matrix)
#' @param whichones which elements of att will be used to compute the aggregated
#'  treatment effect parameter
#' @param weights.agg the weights to apply to each element of att[whichones];
#'  should have the same dimension as att[whichones]
#' @param wif extra influence function term coming from estimating the weights;
#'  should be n x k matrix where k is dimension of whichones
#'
#' @return nx1 influence function
#'
#' @keywords internal
get_agg_inf_func <- function(att, inffunc1, whichones, weights.agg, wif=NULL) {
  # enforce weights are in matrix form
  weights.agg <- as.matrix(weights.agg)

  # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
  thisinffunc <- inffunc1[,whichones]%*%weights.agg

  # Incorporate influence function of the weights
  if (!is.null(wif)) {
    thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
  }

  # return influence function
  return(thisinffunc)
}


#' @title Take influence function and return standard errors
#'
#' @description Function to take an nx1 influence function and return
#'  a standard error
#'
#' @param thisinffunc An influence function
#' @inheritParams compute.aggte
#'
#' @return scalar standard error
#'
#' @keywords internal
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

