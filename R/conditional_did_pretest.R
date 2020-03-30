## The idea here is to combine the weighting function with Y and run the previous
## code for computing group-time average treatment effects
#' @title mp.spatt.test
#'
#' @description integrated moments test for conditional common trends holding in all pre-treatment time
#'  periods across all groups
#'
#' @inheritParams mp.spatt
#' @param weightfun A function that takes in two arguments, X and u, to compute
#'  the weighting function for the test.  The default is \code{1*(X <= u)}
#' @param xformlalist A list of formulas for the X variables.  This allows to
#'  test using different specifications for X, if desired
#' @param clustervarlist A list of cluster variables.  This allows to conduct
#'  the test using different levels of clustering, if desired.
#'
#' @examples
#' \dontrun{
#' data(mpdta)
#' mptest <- mp.spatt.test(lemp ~ treat, xformlalist=list(~lpop), data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year", clustervarlist=list(NULL))
#' summary(mptest[[1]])
#' }
#'
#' data(mpdta)
#' mptest <- mp.spatt.test(lemp ~ treat, xformlalist=list(NULL), data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year", clustervarlist=list(NULL))
#' summary(mptest[[1]])
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#'
#' @return list containing test results
#' @export
conditional_did_pretest <- function(yname, 
                                    tname,
                                    idname=NULL,
                                    first.treat.name,
                                    xformla,
                                    data,
                                    control.group=c("nevertreated","notyettreated"),
                                    aggte=TRUE,
                                    maxe = NULL,
                                    mine = NULL,
                                    w=NULL,
                                    alp=0.05,
                                    bstrap=T, biters=1000, clustervars=NULL,
                                    cband=T,
                                    printdetails=TRUE,
                                    seedvec=NULL, pl=FALSE, cores=1,method="logit",
                                    estMethod="dr", panel=TRUE,
                                    weightfun=NULL) {


  #-----------------------------------------------------------------------------
  # Data pre-processing and error checking
  #-----------------------------------------------------------------------------

  # set control group
  control.group <- control.group[1]
  
  # store parameters for passing around later
  dp <- DIDparams(yname=yname,
                  tname=tname,
                  idname=idname,
                  first.treat.name=first.treat.name,
                  xformla=xformla,
                  data=data,
                  control.group=control.group,
                  maxe=maxe,
                  mine=mine,
                  w=w,
                  alp=alp,
                  bstrap=bstrap,
                  biters=biters,
                  clustervars=clustervars,
                  cband=cband,
                  printdetails=printdetails,
                  seedvec=seedvec,
                  pl=pl,
                  cores=cores,
                  method=method,
                  estMethod=estMethod,
                  panel=panel)
  
  # make sure dataset is a data.frame
  # this gets around RStudio's default of reading data as tibble
  if (!all( class(data) == "data.frame")) {
    #warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  # weights if null
  if(is.character(w))  w <- data[, as.character(w)]
  if(is.null(w)) {
    w <- as.vector(rep(1, nrow(data)))
  } else if(min(w) < 0) stop("'w' must be non-negative")
  data$w <- w

  # Outcome variable will be denoted by y
  data$y <- data[, yname]
  
  # figure out the dates
  # list of dates from smallest to largest
  tlist <- unique(data[,tname])[order(unique(data[,tname]))] 
  # list of treated groups (by time) from smallest to largest
  glist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]

  # Check if there is a never treated grup
  if ( length(glist[glist==0]) == 0) {
    if(control.group=="nevertreated"){
      stop("It seems you do not have a never-treated group in the data. If you do have a never-treated group in the data, make sure to set data[,first.treat.name] = 0 for the observation in this group. Otherwise, select control.group = \"notyettreated\" so you can use the not-yet treated units as a comparison group.")
    } else {
      warning("It seems like that there is not a never-treated group in the data. In this case, we cannot identity the ATT(g,t) for the group that is treated las, nor any ATT(g,t) for t higher than or equal to the largest g.\n \nIf you do have a never-treated group in the data, make sure to set data[,first.treat.name] = 0 for the observation in this group.")
      # Drop all time periods with time periods >= latest treated
      data <- base::subset(data,(data[,tname] < max(glist)))
      # Replace last treated time with zero
      lines.gmax = data[,first.treat.name]==max(glist)
      data[lines.gmax,first.treat.name] <- 0

      ##figure out the dates
      tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest
      # Figure out the groups
      glist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
    }
  }

  # Only the treated groups
  glist <- glist[glist>0]
  
  # check for groups treated in the first period and drop these
  mint <- tlist[1]
  nfirstperiod <- nrow( data[ data[,first.treat.name] == mint, ] )
  if ( nfirstperiod > 0 ) {
    warning(paste0("dropping ", nfirstperiod, " units that were already treated in the first period...this is normal"))
    data <- data[ data[,first.treat.name] != mint, ]
  }

  
  if (!pl) {
    cores <- 1
  }

  ##################################
  ## do some error checking
  if (!is.numeric(tlist)) {
    warning("not guaranteed to order time periods correclty if they are not numeric")
  }
  ## check that first.treat doesn't change across periods for particular individuals
  if (!all(sapply( split(data, data[,idname]), function(df) {
    length(unique(df[,first.treat.name]))==1
  }))) {
    stop("Error: the value of first.treat must be the same across all periods for each particular individual.")
  }
  ####################################
  # How many time periods
  nT <- length(tlist)
  # How many treated groups
  nG <- length(glist)

  
  #################################################################
  # Size of each group (divide by lenth(tlist because we are working with long data))
  #gsize <- aggregate(data[,"w"], by=list(data[,first.treat.name]),
  #                   function(x) sum(x)/length(tlist))
  #################################################################

  # setup data in panel case
  if (panel) {
    # make it a balanced data set
    data <- BMisc::makeBalancedPanel(data, idname, tname)

    # create an n-row data.frame to hold the influence function later
    #dta <- data[ data[,tname]==tlist[1], ]  

    n <- nrow(data[ data[,tname]==tlist[1], ]) # use this for influence function
    
    # check that first.treat doesn't change across periods for particular individuals
    if (!all(sapply( split(data, data[,idname]), function(df) {
      length(unique(df[,first.treat.name]))==1
    }))) {
      stop("Error: the value of first.treat must be the same across all periods for each particular individual.")
    }
  } else {

    # n-row data.frame to hold the influence function
    data$rowid <- seq(1:nrow(data))
    n <- nrow(data)
    #dta <- data
  }

  # put in blank xformla if no covariates
  if (is.null(xformla)) {
    xformla <- ~1
  }

  if (is.null(weightfun)) {
    weightfun <- indicator
  }

  #xformlalist <- lapply(xformlalist, function(ff) if (is.null(ff)) ~1 else ff)

  thecount <- 1
  innercount <- 1

  #outlist <- lapply(xformlalist, function(xformla) {

  ##
  ##X1 <- apply(model.matrix(xformla, data), 2, function(col) {
  ##    pnorm(col)
  ##})## for the entire dataset

  dta <- data[ data[,tname]==tlist[1], ]
  
  ## X <- apply(model.matrix(xformla, dta), 2, function(col) {
  ##   ##pnorm(col)
  ##   ecdf(col)(col)
  ## })

  ## #for debugging:
  ## #X <- as.matrix(X[1:100,])

  ## #X <- unique(X) # do we want to do this...
  
  ## X1 <- apply(model.matrix(xformla, dta), 2, function(col) {
  ##   ecdf(col)(col)
  ##   ##pnorm(col)
  ## })

  X <- model.matrix(xformla, dta)
  X1 <- X

  #for debugging:
  # X <- as.matrix(X[1:100,])
  
  ## thetlist <- list()
  ## pscorelist <- list()
  ## for (g in 1:nG) {

  ##   disdat <- data[(data[,tname]==tlist[1]),]

  ##   disdat$C <- 1*(disdat[,first.treat.name] == 0)

  ##   disdat$G <- 1*(disdat[,first.treat.name] == flist[g])

  ##   disdat <- droplevels(disdat)

  ##   pformla <- xformla
  ##   pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))##formula.tools::lhs(pformla) <- as.name("G")
  ##   pscore.reg <- glm(pformla, family=binomial(link="logit"),
  ##                     data=subset(disdat, C+G==1))
  ##   thetlist[[g]] <- coef(pscore.reg)
  ##   pscorelist[[g]] <- predict(pscore.reg, newdata=disdat, type="response")
  ## }   

  
  
  thecount <<- thecount+1
  n.preperiods <- sum(sapply(glist, function(g) 1*(g>tlist[-1])))
  inffunc1 <- array(data=0, dim=c(n.preperiods, n, n))


  cat("Step 1 of 2: Computing test statistic....\n")
  out <- pbapply::pblapply(1:nrow(X), function(i) {
    # these are the weights for the conditional moment test
    # all this is super hack and specific to indicator weighting
    www <- as.numeric(weightfun(X1, X[i,]))
    rightids <- dta[,idname][www==1]
    thisdata <- data
    thisdata[,yname] <- 0
    thisdata[ thisdata[,idname] %in% rightids,yname] <- data[ thisdata[,idname] %in% rightids,yname] # this is only going to work if use indicator weights, I think (otherwise you will multiply by www twice

    Jres <- compute.att_gt(nG=nG,
                           nT=nT,
                           glist=glist,
                           tlist=tlist,
                           data=thisdata,
                           n=n,
                           first.treat.name=first.treat.name,
                           yname=yname,
                           tname=tname,
                           w=w,
                           idname=idname,
                           xformla=xformla,
                           method=method,
                           seedvec=seedvec,
                           pl=FALSE,
                           cores=1,
                           printdetails=FALSE,
                           control.group="nevertreated",
                           estMethod=estMethod,
                           panel=panel)


    J.results <- process.attgt(Jres)
    group <- J.results$group
    att <- J.results$att
    tt <- J.results$tt
    inf.func <- J.results$inf.func

    list(J=att, group=group, t=tt, inf.func=inf.func)
  }, cl=cores)

  

  # outinffunc <- lapply(out, function(o) t(o$inf.func))
  Jinf.func <- simplify2array(getListElement(out, "inf.func"))
  
  J <- t(sapply(out, function(o) o$J))
  #KS <- sqrt(n) * sum(apply(J,2,function(j) max(abs(j))))
  #CvM <- n*sum(apply(J, 2, function(j) mean( j^2 )))
  # average J^2 across x, and then sum across g and t
  CvM <- n*sum(apply(J^2, 2, mean)) 
  

  
  #boot.res <- test.empboot(dp)
  cat("Step 2 of 2: Simulating limiting distribution of test statistic....\n")
  boot.res <- test.mboot(Jinf.func, dp, cores=cores)


  ## #some debugging code
  ## ddd <- 5
  ## bout <- sapply(1:100, function(b) {

  ##   Jstar <- apply(sample(c(-1,1), size=n, replace=TRUE)*Jinf.func[,ddd,], 2, mean)
  ##   mean((sqrt(n)*Jstar)^2)
  ## })

  ## ts <- mean((sqrt(n)*J[,ddd])^2)
  ## 1-ecdf(bout)(ts)


  ## ts2.inner <- apply( (sqrt(n)*J)^2, 2, mean)
  ## ts2 <- sum(ts2.inner)
  ## bout2 <- sapply(1:100, function(b) {

  ##   Jstar <- t(apply(sample(c(-1,1), size=n, replace=TRUE)*Jinf.func[,,], c(2,3), mean))
  ##   cvm.inner <- apply((sqrt(n)*Jstar)^2, 2, mean)
  ##   sum(cvm.inner)
  ## })
  ## 1-ecdf(bout2)(ts2)
 
  
  CvMb <- boot.res$bres
  ## CvMb <- unlist(lapply(bout, function(b1) b1$CvMb))
  CvM.crit.val <- boot.res$crit.val
  CvMpval <- 1-ecdf(CvMb)(CvM)

  ## KSb <- unlist(lapply(bout, function(b1) b1$KSb))
  ## KSocval <- quantile(KSb, probs=(1-alp), type=1)
  ## KSpval <- 1-ecdf(KSb)(KS)


  ## bres <- lapply(bout, simplify2array)
  ## CvMb <- sapply(bres, function(b) apply(b, 1, function(bb) mean(bb^2)))
  ## CvMb <- n*apply(CvMb, 2, sum)
  ## CvMocval <- quantile(CvMb, probs=(1-alp), type=1)
  ## CvMpval <- 1-ecdf(CvMb)(CvM)
  ## bres <- sapply(bres, function(b) apply(b, 1, function(bb) max(abs(bb))))
  ## KSb <- sqrt(n)*apply(bres, 2, sum)
  ## KSocval <- quantile(KSb, probs=(1-alp), type=1)
  ## KSpval=1-ecdf(KSb)(KS)
  out <- list(CvM=CvM, CvMb=CvMb, CvMcval=CvM.crit.val, CvMpval=CvMpval)#KS=, KSb=KSb, KScval=KSocval, KSpval=KSpval, clustervars=clustervars, xformla=xformla)
  out
}



#' @title expf
#'
#' @description exponential weighting function
#'
#' @param X matrix of X's from the data
#' @param u a particular value to multiply times the X's
#'
#' @return numeric vector
#' @examples
#' data(mpdta)
#' dta <- subset(mpdta, year==2007)
#' X <- model.matrix(~lpop, data=dta)
#' X <- expf(X, X[1,])
#'
#' @export
expf <- function(X, u) {
    exp(X%*%u)
}

#' @title indicator
#'
#' @description indicator weighting function
#'
#' @param X matrix of X's from the data
#' @param u a particular value to compare X's to
#'
#' @return numeric vector
#'
#' @examples
#' data(mpdta)
#' dta <- subset(mpdta, year==2007)
#' X <- model.matrix(~lpop, data=dta)
#' X <- indicator(X, X[1,])
#'
#' @export
indicator <- function(X, u) {
  # check if each element in each row of X <= corresponding element in u
  cond <- t(apply( X, 1, function(x) x <= u))
  # check if entire row of X <= entire row of u
  1*apply(cond, 1, all)
  #apply(X <= u, 1, function(b) 1*all(b))
}


test.mboot <- function(inf.func, DIDparams, cores=1) {

  # setup needed variables
  data <- DIDparams$data
  idname <- DIDparams$idname
  clustervars <- DIDparams$clustervars
  biters <- DIDparams$biters
  tname <- DIDparams$tname
  tlist <- unique(data[,tname])[order(unique(data[,tname]))]
  alp <- DIDparams$alp
  
  # just get n obsevations (for clustering below...)
  dta <- data[ data[,tname]==tlist[1], ]
  n <- nrow(dta)
  
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
  bout <- pbapply::pbsapply(1:biters, FUN=function(b) {
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
    Jb <- t(apply(Ub*inf.func, c(2,3), mean))
    CvMb <- n*sum(apply(Jb^2, 2, mean))
    # return bootstrap draw
    CvMb
  }, cl=cores)
  # bootstrap results
  ## bres <- t(simplify2array(bout))
  ## # bootstrap variance matrix 
  ## # V <- cov(bres)
  ## # bootstrap standard error
  ## bSigma <- apply(bres, 2,
  ##                 function(b) (quantile(b, .75, type=1, na.rm = T) -
  ##                                quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
  ## # critical value for uniform confidence band
  ## bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
  crit.val <- quantile(bout, 1-alp, type=1)

  list(bres=bout, crit.val=crit.val)
}



test.empboot <- function(inf.func, DIDparams) {

  # setup needed variables
  data <- DIDparams$data
  idname <- DIDparams$idname
  clustervars <- DIDparams$clustervars
  biters <- DIDparams$biters
  tname <- DIDparams$tname
  tlist <- unique(data[,tname])[order(unique(data[,tname]))]
  alp <- DIDparams$alp
  
  # just get n obsevations (for clustering below...)
  dta <- data[ data[,tname]==tlist[1], ]
  n <- nrow(dta)
  
  # bootstrap
  bout <- pbapply::pbsapply(1:biters, FUN=function(b) {
    if (length(clustervars) > 0) {
      stop("no clustering allowed")
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
      bdata <- BMisc::blockBootSample(data, idname)
    }
    # multiply weights onto influence function
    Jb <- t(apply(Ub*inf.func, c(2,3), mean))
    CvMb <- n*sum(apply(Jb^2, 2, mean))
    # return bootstrap draw
    CvMb
  })
  # bootstrap results
  ## bres <- t(simplify2array(bout))
  ## # bootstrap variance matrix 
  ## # V <- cov(bres)
  ## # bootstrap standard error
  ## bSigma <- apply(bres, 2,
  ##                 function(b) (quantile(b, .75, type=1, na.rm = T) -
  ##                                quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
  ## # critical value for uniform confidence band
  ## bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
  crit.val <- quantile(bout, 1-alp, type=1)

  list(bres=bout, crit.val=crit.val)
}
