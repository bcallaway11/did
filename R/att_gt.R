#' @title att_gt
#'
#' @description \code{att_gt} computes average treatment effects in DID setups  where there are more
#'  than two periods of data and allowing for treatment to occur at different points in time.
#'  See Callaway ans Sant'Anna (2019) for a detailed description.
#'
#' @param outcome The outcome y (in quotations, always!)
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param idname The individual (cross-sectional unit) id name
#' @param first.treat.name The name of the variable in \code{data} that contains the first
#'  period when a particular observation is treated.  This should be a positive
#'  number for all observations in treated groups.  It should be 0 for observations
#'  in the untreated group.
#' @param nevertreated Boolean for using the group which is never treated in the sample as the comparison unit. Default is TRUE.
#' @param aggte boolean for whether or not to compute aggregate treatment effect parameters, default TRUE
#' @param maxe maximum values of periods ahead to be computed in event study. Only used if aggte = T.
#' @param mine minimum values of periods ahead to be computed in event study. Only used if aggte = T.
#' @param w The name of the column containing the sampling weights. If not set, all observations have same weight.
#' @param alp the significance level, default is 0.05
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier boostrap.  If standard errors are clustered, then one
#'  must set \code{bstrap=TRUE}. Default is \code{TRUE}.
#' @param biters The number of boostrap iterations to use.  The default is 1000,
#'  and this is only applicable if \code{bstrap=TRUE}.
#' @param clustervars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level.
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability \code{1-alp}.  The default is \code{TRUE}
#'  and the resulting standard errors will be pointwise.
#' @param printdetails Boolean for showing detailed results or not
#'
#' @param method The method for estimating the propensity score when covariates
#'  are included (not implemented)
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors#' (not implemented)
#' @param pl Boolean for whether or not to use parallel processing (not implemented)
#' @param cores The number of cores to use for parallel processing (not implemented)
#' @param estMethod the method to compute group-time average treatment effects.  The default is "dr" which uses the doubly robust
#' approach in the \code{DRDID} package.  Other built-in methods
#' include "ipw" for inverse probability weighting and "reg" for
#' first step regression estimators.  The user can also pass their
#' own function for estimating group time average treatment
#' effects.  This should be a function
#' \code{f(Y1,Y0,treat,covariates)} where \code{Y1} is an
#' \code{n} x \code{1} vector of outcomes in the post-treatment
#' outcomes, \code{Y0} is an \code{n} x \code{1} vector of
#' pre-treatment outcomes, \code{treat} is a vector indicating
#' whether or not an individual participates in the treatment,
#' and \code{covariates} is an \code{n} x \code{k} matrix of
#' covariates.  The function should return a list that includes
#' \code{ATT} (an estimated average treatment effect), and
#' \code{inf.func} (an \code{n} x \code{1} influence function).
#' The function can return other things as well, but these are
#' the only two that are required. \code{estMethod} is only used
#' if covariates are included.

#' @references Callaway, Brantly and Sant'Anna, Pedro H. C.. "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#'
#' @return \code{MP} object
#'
#' @export
att_gt <- function(yname, 
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
                   seedvec=NULL, pl=FALSE, cores=2,method="logit",
                   estMethod="dr", panel=TRUE) {


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

  #-----------------------------------------------------------------------------
  # Compute all ATT(g,t)
  #-----------------------------------------------------------------------------
  results <- compute.att_gt(nG=nG,
                            nT=nT,
                            glist=glist,
                            tlist=tlist,
                            data=data,
                            n=n,
                            first.treat.name=first.treat.name,
                            yname=yname,
                            tname=tname,
                            w=w,
                            idname=idname,
                            xformla=xformla,
                            method=method,
                            seedvec=seedvec,
                            pl=pl,
                            cores=cores,
                            printdetails=printdetails,
                            control.group=control.group,
                            estMethod=estMethod,
                            panel=panel)

  # extract ATT(g,t) and influence functions
  attgt.list <- results$attgt.list
  inffunc <- results$inffunc

  attgt.results <- process.attgt(results)
  group <- attgt.results$group
  att <- attgt.results$att
  tt <- attgt.results$tt
  inffunc1 <- attgt.results$inf.func
  

  # estimate variance
  # this is analogous to cluster robust standard errors that
  # are clustered at the unit level
  V <- t(inffunc1)%*%inffunc1/n

  # if clustering along another dimension...we require using the
  # bootstrap (in principle, could come up with analytical standard
  # errors here though)
  if ( (length(clustervars) > 0) & !bstrap) {
    warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  }


  # bootstrap variance matrix
  if (bstrap) {

    bout <- mboot(inffunc1, DIDparams=dp)
    bres <- bout$bres
    V <- bout$V
  }


  
  #-----------------------------------------------------------------------------
  # compute Wald pre-test
  #-----------------------------------------------------------------------------

  # select which periods are pre-treatment
  pre <- which(group > tt)

  # pseudo-atts in pre-treatment periods
  preatt <- as.matrix(att[pre])

  # covariance matrix of pre-treatment atts
  preV <- as.matrix(V[pre,pre])

  # check if there are actually any pre-treatment periods
  if (length(preV) == 0) {
    message("No pre-treatment periods to test")
    W  <- NULL
    Wpval <- NULL
  } else if (det(preV) == 0) {
    # singluar covariance matrix for pre-treatment periods
    warning("Not returning pre-test Wald statistic due to singular covariance matrix")
    W <- NULL
    Wpval <- NULL
  } else {
    # everything is working...
    W <- n*t(preatt)%*%solve(preV)%*%preatt
    q <- length(pre) # number of restrictions
    Wpval <- round(1-pchisq(W,q),5)
  }


  #-----------------------------------------------------------------------------
  # compute confidence intervals / bands
  #-----------------------------------------------------------------------------

  # critical value from N(0,1), for pointwise
  cval <- qnorm(1-alp/2)

  # in order to get uniform confidencs bands
  # HAVE to use the bootstrap
  if (bstrap){
    if (cband) {
      # for uniform confidence band
      # compute new critical value
      # see paper for details
      bSigma <- apply(bres, 2,
                      function(b) (quantile(b, .75, type=1, na.rm = T) -
                                     quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
      bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
      cval <- quantile(bT, 1-alp, type=1, na.rm = T)
      ## this is the appropriate matrix for constructing confidence bands - BUT NOT FOR WALD TEST!!!
      # V <- diag(bSigma^2)
    }
  }

  ## #-----------------------------------------------------------------------------
  ## # Compute all summaries of the ATT(g,t)
  ## #-----------------------------------------------------------------------------
  ## aggeffects <- NULL
  ## if (aggte & (nT > 2)) {
  ##   aggeffects <- compute.aggte(glist,
  ##                               tlist,
  ##                               group,
  ##                               tt,
  ##                               att,
  ##                               first.treat.name,
  ##                               inffunc1,
  ##                               n,
  ##                               clustervars,
  ##                               dta,
  ##                               idname,
  ##                               bstrap,
  ##                               biters,
  ##                               alp,
  ##                               cband,
  ##                               maxe,
  ##                               mine)
  ## }

  # Return this list
  return(MP(group=group, t=tt, att=att, V=V, c=cval, inffunc=inffunc1, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp))

}


process.attgt <- function(attgt.results.list) {
  attgt.list <- attgt.results.list$attgt
  inffunc <- attgt.results.list$inffunc
  nG <- length(unique(unlist(getListElement(attgt.list, "group"))))
  nT <- length(unique(unlist(getListElement(attgt.list, "year"))))
  # pick up number of observations from the influence function
  # this will be flexible for panel or repeated cross sections
  n <- length(inffunc[1,1,]) 
  # create vectors to hold the results
  group <- c()
  att <- c()
  tt <- c()
  i <- 1

  # matrix to hold influence function
  # (note: this is relying on having a balanced panel,
  # which we do currently enforce)
  inffunc1 <- matrix(0, ncol=nG*nT, nrow=n) 

  # populate result vectors and matrices
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      att[i] <- attgt.list[[i]]$att
      inffunc1[,i] <- inffunc[f,s,]
      i <- i+1
    }
  }

  list(group=group, att=att, tt=tt, inf.func=inffunc1)
}
