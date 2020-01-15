#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#' @title att_gt_het
#'
#' @description \code{att_gt_het} computes the difference of ATT across two subpopulations in the case where there are more
#'  than two periods of data and allowing for treatment to occur at different points in time
#'  extending the method of Abadie (2005).  This method relies on once individuals are treated
#'  they remain in the treated state for the duration.
#'
#' @param outcome The outcome y (in quotations, always!)
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param aggte boolean for whether or not to compute aggregate treatment effect parameters, default TRUE
#' @param w The name of the column containing the weights
#' @param idname The individual (cross-sectional unit) id name
#' @param first.treat.name The name of the variable in \code{data} that contains the first
#'  period when a particular observation is treated.  This should be a positive
#'  number for all observations in treated groups.  It should be 0 for observations
#'  in the untreated group.
#' @param alp the significance level, default is 0.05
#' @param method The method for estimating the propensity score when covariates
#'  are included
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier boostrap.  If standard errors are clustered, then one
#'  must set \code{bstrap=TRUE}.
#' @param biters The number of boostrap iterations to use.  The default is 100,
#'  and this is only applicable if \code{bstrap=TRUE}.
#' @param clustervars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level.
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors#' (not implemented)
#' @param pl Boolean for whether or not to use parallel processing
#' @param cores The number of cores to use for parallel processing
#' @param printdetails Boolean for showing detailed results or not
#' @param maxe maximum values of periods ahead to be computed in event study
#' @param mine minimum values of periods ahead to be computed in event study
#' @param nevertreated Boolean for using the group which is never treated in the sample as the comparison unit. Default is TRUE.
#' @param het The name of the column containing the (binary) categories for heterogeneity
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#' @return \code{MP} object
#'
#' @export
att_gt_het <- function(outcome, data, tname,
                       aggte=TRUE, w=NULL,
                       idname=NULL, first.treat.name, alp=0.05,
                       method="logit",
                       bstrap=FALSE, biters=1000, clustervars=NULL,

                       seedvec=NULL, pl=FALSE, cores=2,
                       printdetails=TRUE,
                       maxe = NULL,
                       mine = NULL,
                       nevertreated = T,
                       het) {

  ## make sure that data is a data.frame
  ## this gets around RStudio's default of reading data as tibble
  if (!all( class(data) == "data.frame")) {
    #warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  # weights if null
  if(is.character(w))  w <- data[, as.character(w)]
  if(is.null(w)) {
    w <- as.vector(rep(1, nrow(data)))
  } else if(min(w) < 0) stop("'w' must be non-negative")
  # het if null
  if(is.null(het)) {
    stop("Please specifiy 'het'. If het=NULL, use 'att_gt' instead of 'att_gt_het'.")
  }
  if(is.character(het))  het <- data[, as.character(het)]

  het.dim <- length(unique(het))
  if(het.dim!=2)  {
    stop("'het' must be a binary variable.")
  }

  data$w1 <- w * (het==1)
  data$w0 <- w * (het==0)

  data$y <- data[, as.character(outcome)] ##data[,as.character(formula.tools::lhs(formla))]
  ##figure out the dates and make balanced panel
  tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

  flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
  if ( length(flist[flist==0]) == 0) {
    warning("dataset does not have any observations in the control group.  make sure to set data[,first.treat.name] = 0 for observations in the control group.")
  }
  # First treated groups
  flist <- flist[flist>0]

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
  tlen <- length(tlist)
  # How many treated groups
  flen <- length(flist)

  data <- BMisc::makeBalancedPanel(data, idname, tname)
  #dta is used to get a matrix of size n (like in cross sectional data)
  dta <- data[ data[,tname]==tlist[1], ]  ## use this for the influence function



  #################################################################
  #################################################################

  #################################################################
  #----------------------------------------------------------------------------
  # Results for het==1
  #----------------------------------------------------------------------------
  results_het1 <- compute.att_gt_het(flen, tlen, flist, tlist, data, dta, first.treat.name,
                                     outcome, tname, idname, method, seedvec,
                                     pl, cores, printdetails, nevertreated, het=1)
  fatt_het1 <- results_het1$fatt
  inffunc_het1 <- results_het1$inffunc
  #----------------------------------------------------------------------------
  # Results for het==0
  #----------------------------------------------------------------------------
  results_het0 <- compute.att_gt_het(flen, tlen, flist, tlist, data, dta, first.treat.name,
                                     outcome, tname, idname, method, seedvec,
                                     pl, cores, printdetails, nevertreated, het=0)
  fatt_het0 <- results_het0$fatt
  inffunc_het0 <- results_het0$inffunc
  #----------------------------------------------------------------------------

  ## process the results from computing the spatt
  group <- c()
  t    <- c()
  att_het1 <- c()
  att_het0 <- c()
  i <- 1

  inffunc1_het1 <- matrix(0, ncol=flen*(tlen), nrow=nrow(dta)) ## note, this might not work in unbalanced case
  inffunc1_het0 <- matrix(0, ncol=flen*(tlen), nrow=nrow(dta))


  for (f in 1:length(flist)) {
    for (s in 1:(length(tlist))) {
      group[i] <- fatt_het1[[i]]$group
      t[i] <- fatt_het1[[i]]$year

      att_het1[i] <- fatt_het1[[i]]$att
      att_het0[i] <- fatt_het0[[i]]$att

      inffunc1_het1[,i] <- inffunc_het1[f,s,]
      inffunc1_het0[,i] <- inffunc_het0[f,s,]

      i <- i+1
    }
  }

  # THIS IS ANALOGOUS TO CLUSTER ROBUST STD ERRORS (in our specific setup)
  n <- nrow(dta)
  V <- NULL

  aggeffects <- NULL
  aggeffects_het1 <- NULL
  aggeffects_het0 <- NULL

  if (aggte) {
    aggeffects_het1 <- compute.aggte_het(flist, tlist, group, t, att_het1, first.treat.name, inffunc1_het1,
                                         n, clustervars, dta, idname, bstrap, biters, alp, maxe, mine, het=1)

    aggeffects_het0 <- compute.aggte_het(flist, tlist, group, t, att_het0, first.treat.name, inffunc1_het0,
                                         n, clustervars, dta, idname, bstrap, biters, alp, maxe, mine, het=0)


    aggeffects <- list(simple.att = aggeffects_het1$simple.att - aggeffects_het0$simple.att,
                       simple.att.inf.func = aggeffects_het1$simple.att.inf.func - aggeffects_het0$simple.att.inf.func,

                       dynamic.att = aggeffects_het1$dynamic.att - aggeffects_het0$dynamic.att,
                       dynamic.att.inf.func = aggeffects_het1$dynamic.att.inf.func - aggeffects_het0$dynamic.att.inf.func,

                       dynamic.att.e = aggeffects_het1$dynamic.att.e - aggeffects_het0$dynamic.att.e,
                       dyn.inf.func.e = aggeffects_het1$dyn.inf.func.e - aggeffects_het0$dyn.inf.func.e,
                       e = aggeffects_het1$e
    )


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
        bres <- as.vector(simplify2array(bout))

        bSigma <- (quantile(bres, .75, type=1, na.rm = T) - quantile(bres, .25, type=1, na.rm = T)) /
          (qnorm(.75) - qnorm(.25))
        #bT <- abs(bres/bSigma))
        return(bSigma)
        #return(sqrt( mean( bres^2)) /sqrt(n))
      } else {
        return(sqrt( mean( (thisinffunc)^2 ) /n ))
      }
    }



    aggeffects$simple.se <- getSE_inf(as.matrix(aggeffects$simple.att.inf.func))
    aggeffects$dynamic.se <- getSE_inf(as.matrix(aggeffects$dynamic.att.inf.func))

    aggeffects$dynamic.se.e <- sqrt(colMeans((aggeffects$dyn.inf.func.e)^2)/n)

    aggeffects$c.dynamic <- qnorm(1 - alp/2)

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
        Rb <- (base::colMeans(Ub*(aggeffects$dyn.inf.func.e), na.rm = T))
        Rb
      })
      bres <- t(simplify2array(bout))
      # Non-degenerate dimensions
      ndg.dim <- base::colSums(bres)!=0
      bres <- bres[,ndg.dim]

      #V.dynamic <- cov(bres)
      bSigma <- apply(bres, 2, function(b) (quantile(b, .75, type=1,na.rm = T) - quantile(b, .25, type=1,na.rm = T))/(qnorm(.75) - qnorm(.25)))
      bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
      aggeffects$c.dynamic <- quantile(bT, 1-alp, type=1,na.rm = T)
      #aggeffects$dynamic.se.e <- bSigma
      aggeffects$dynamic.se.e <- rep(0,length(ndg.dim))
      aggeffects$dynamic.se.e[ndg.dim] <- bSigma
    }


  }

  out <- list(group=group,
              t=t,
              att_het1=att_het1, att_het0=att_het0,
              inffunc_het1=inffunc_het1, inffunc_het0=inffunc_het0,
              n=n,
              aggte_het1=aggeffects_het1, aggte_het0=aggeffects_het0,
              aggte = aggeffects,
              alp = alp)

  return(out)
}



#' @title compute.att_gt_het
#'
#' @description \code{compute.att_gt_het} does the main work for computing
#'  mutliperiod group-time average treatment effects with heterogeneity
#'
#' @inheritParams att_gt_het
#'
#' @return a list with length equal to the number of groups times the
#'  number of time periods; each element of the list contains a \code{QTE}
#'  object that contains group-time average treamtent effect as well
#'  as which group it is for and which time period it is for and
#'  the influence function which is used externally to compute
#'  standard errors.
#'
#' @keywords internal
#'
#' @export
compute.att_gt_het <- function(flen, tlen, flist, tlist, data, dta,
                               first.treat.name, outcome, tname, idname,
                               method, seedvec,
                               pl, cores, printdetails, nevertreated, het) {

  yname <- outcome ##as.character(formula.tools::lhs(formla))

  fatt <- list()
  counter <- 1
  inffunc <- array(data=0, dim=c(flen,tlen,nrow(dta)))

  het.val <- het

  for (f in 1:flen) {
    ##satt <- list()
    for (t in 1:(tlen)) {
      #pret <- t
      pret <- utils::tail(which(tlist < flist[f]),1)
      if (flist[f]<=tlist[(t)]) {
        ## set an index for the pretreatment period
        #pret <- utils::tail(which(tlist < flist[f]),1)

        ## print a warning message if there are no pre-treatment
        ##  periods
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", flist[f]))
        }

        ## print the details of which iteration we are on
        if (printdetails) {
          cat(paste("current period:", tlist[t]), "\n")
          cat(paste("current group:", flist[f]), "\n")
          cat(paste("set pretreatment period to be", tlist[pret]), "\n")
        }
      }

      ## --------------------------------------------------------
      post.treat <- 1*(flist[f]<=tlist[(t)])

      if (tlist[pret] == tlist[t]){
        fatt[[counter]] <- list(att=0, group=flist[f], year=tlist[(t)], post=post.treat)
        inffunc[f,t,] <- 0
      } else {
        ## results for the case with panel data
        ## get dataset with current period and pre-treatment period
        disdat <- data[(data[,tname]==tlist[t] | data[,tname]==tlist[pret]),]
        ## transform it into "cross-sectional" data where
        ## one of the columns contains the change in the outcome
        ## over time
        ########################## TOASK
        # SHOULD we keep the data on t=t_pre or t=t_post? right now, it is t=t_pre
        disdat <- BMisc::panel2cs(disdat, yname, idname, tname)

        #THIS IS THE PART WE CAN CHANGE FOR THE NOT YET TREATED!!
        ## set up control group
        if(nevertreated ==T){
          disdat$C <- 1*(disdat[,first.treat.name] == 0)
        }

        if(nevertreated ==F){
          disdat$C <- 1*( (disdat[,first.treat.name] == 0) +
                            (disdat[,first.treat.name] > max(disdat[, tname])) )
        }

        ## set up for particular treated group
        disdat$G <- 1*(disdat[,first.treat.name] == flist[f])

        ## drop missing factors
        disdat <- droplevels(disdat)

        ## give short names for data in this iteration
        G <- disdat$G
        C <- disdat$C
        dy <- disdat$dy * ((-1)^(1+post.treat))
        #n <- nrow(disdat)
        w <- (disdat$w1) * het.val + (disdat$w0) * (1 - het.val)


        ## set up weights
        attw <- w * G/mean(w * G)
        attw2a <- w * C
        attw2 <- attw2a/mean(attw2a)
        att <- mean((attw - attw2)*dy)

        if(is.na(att)) att <- 0

        ## save results for this iteration
        fatt[[counter]] <- list(att=att, group=flist[f], year=tlist[(t)], post=post.treat)

        ## --------------------------------------------
        ## get the influence function

        ## weigts for het==1
        wg <- w * G/mean(w * G)
        wc1 <- w * C
        wc <- wc1 / mean(wc1)

        ## influence function for treated group
        psig <- wg*(dy - mean(wg*dy))
        # influence function for the control group
        psic <- wc*(dy - mean(wc*dy))

        ## save the influnce function as the difference between
        ## the treated and control influence functions;
        ## we save this as a 3-dimensional array
        ## and then process afterwards
        infl.att <- (psig - psic)
        if(base::anyNA(infl.att)) infl.att <- 0

        inffunc[f,t,] <- infl.att
      }

      counter <- counter+1
    }

  }

  list(fatt=fatt, inffunc=inffunc)
}






#' @title compute.aggte_het
#'
#' @description does the heavy lifting on computing aggregated group-time
#'  average treatment effects
#'
#' @inheritParams att_gt_het
#'
#' @return \code{AGGTE_het} object
#'
#' @keywords internal
#'
#' @export
compute.aggte_het <- function(flist, tlist, group, t, att, first.treat.name, inffunc1, n,
                              clustervars, dta, idname, bstrap, biters, alp, maxe, mine, het) {

  if ( (length(clustervars) > 0) & !bstrap) {
    warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  }
  getSE <- function(whichones, weights, wif=NULL) {
    weights <- as.matrix(weights) ## just in case pass vector
    thisinffunc <- inffunc1[,whichones]%*%weights  ##multiplies influence function times weights and sums to get vector of weighted IF (of length n)
    if (!is.null(wif)) {
      thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
    }

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
      bres <- as.vector(simplify2array(bout))

      bSigma <- (quantile(bres, .75, type=1, na.rm = T) - quantile(bres, .25, type=1, na.rm = T)) /
        (qnorm(.75) - qnorm(.25))
      #bT <- abs(bres/bSigma))
      return(bSigma)
      #return(sqrt( mean( bres^2)) /sqrt(n))
    } else {
      return(sqrt( mean( (thisinffunc)^2 ) /n ))
    }
  }



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
      bres <- as.vector(simplify2array(bout))

      bSigma <- (quantile(bres, .75, type=1, na.rm = T) - quantile(bres, .25, type=1, na.rm = T)) /
        (qnorm(.75) - qnorm(.25))
      #bT <- abs(bres/bSigma))
      return(bSigma)
      #return(sqrt( mean( bres^2)) /sqrt(n))
    } else {
      return(sqrt( mean( (thisinffunc)^2 ) /n ))
    }
  }

  ## do some recoding to make sure time periods are 1 unit apart
  ## and then put these back together at the end
  originalt <- t
  originalgroup <- group
  originalflist <- flist
  originaltlist <- tlist
  uniquet <- seq(1,length(unique(t)))
  ## function to switch from "new" t values to
  ##  original t values
  t2orig <- function(t) {
    unique(c(originalt,0))[which(c(uniquet,0)==t)]
  }
  ## function to switch between "original"
  ##  t values and new t values
  orig2t <- function(orig) {
    c(uniquet,0)[which(unique(c(originalt,0))==orig)]
  }
  t <- sapply(originalt, orig2t)
  group <- sapply(originalgroup, orig2t)
  flist <- sapply(originalflist, orig2t)

  weights.agg = (dta$w1) * (het) + (dta$w0) * (1 - het)
  #weights.agg = 1
  ## some variables used throughout
  ever.treated <- 1 * (dta[,first.treat.name]>0)
  mean.w.ever.treated <- mean(weights.agg * ever.treated)

  # Probability of being in group g (among ever-treated only!)
  pg <- sapply(originalflist,
               function(g) mean(weights.agg * ever.treated * (dta[,first.treat.name]==g))/
                 mean(weights.agg * ever.treated))
  pgg <- pg
  pg <- pg[match(group, flist)] ## make it have the same length as att
  attg <- split(att, group)
  tg <- split(t, group)
  keepers <- which(group <= t)
  G <-  unlist(lapply(dta[,first.treat.name], orig2t))

  ## simple att
  simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))
  # Estimation effect coming from P(G=g| Ever treated)
  # Part 1: est effect from P(G=g) treating P(ever treated) as known
  simple.oif1 <- sapply(keepers,
                        function(k) ( (weights.agg *(G==group[k]) - mean(weights.agg * (G==group[k]))) /
                                        mean.w.ever.treated
                        )
  )
  # Part 2: est effect from  P(ever treated) treating P(G=g) as known
  #simple.oif2 <- sapply(keepers,
  #                      function(j) mean(weights.agg * (G==group[j])) *
  #                        apply(sapply(keepers, function(k) (weights.agg*(G==group[k]) - mean(weights.agg*(G==group[k])))),1,sum))
  simple.oif2 <- sapply(keepers,
                        function(j) ((mean(weights.agg * (G==group[j]))/(mean.w.ever.treated^2)) *
                                       (weights.agg * ever.treated - mean.w.ever.treated)
                        )
  )
  # Estimation effect from numerator
  simple.oif <- (simple.oif1 - simple.oif2)/(sum(pg[keepers]))
  #Estimation effect from denominator of the weights (normalization)
  simple.oif3 <- rowSums(simple.oif) %*%  t(matrix(pg[keepers]/sum(pg[keepers])))
  #Estimation effect from estimated weights (in total)
  simple.oif <- simple.oif - simple.oif3

  simple.se <- getSE(keepers, pg[keepers]/sum(pg[keepers]), simple.oif)

  simple.att.inf.func <- inffunc1[,keepers]%*% as.matrix(pg[keepers]/sum(pg[keepers]))
  simple.att.inf.func <- simple.att.inf.func + simple.oif %*% as.matrix(att[keepers])


  ######################
  # IF WE WANT TO PLOT PRE-TREATMENT, JUST NEED TO START eseq FROM NEGATIVE NUMBERS!
  ###################
  ###################
  ## Dynamic Treatment Effects
  if(is.null(maxe)) maxe <- max(t-group)+1
  if (maxe > (max(t-group)+1)) maxe <- max(t-group)+1

  eseq <- seq(1,maxe)
  dynamic.att.e <- sapply(eseq, function(e) {
    whiche <- which(t - group + 1 == e)
    atte <- att[whiche]
    pge <- pg[whiche]/(sum(pg[whiche]))
    sum(atte*pge)
  })


  time.treated <- ((max(originalt)) - dta[,first.treat.name] +1) * (dta[,first.treat.name]>0)

  dynamic.e.inf.f <- sapply(eseq, function(e) {
    whiche <- which(t - group + 1 == e)
    pge <- pg[whiche]/sum(pg[whiche])

    ## some variables used
    atleast.e.treated <- 1 * (time.treated >=e)
    mean.w.atleast.e.treated <- mean(weights.agg * atleast.e.treated)

    # Estimation effect coming from P(G=g| treated for at least e periods)
    # Part 1: est effect from P(G=g) treating P(treated for at least e periods) as known
    dynamic.oif1 <- sapply(whiche,
                           function(k) ( (weights.agg *(G==group[k]) - mean(weights.agg * (G==group[k]))) /
                                           mean.w.atleast.e.treated
                           )
    )

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


  n.pre.treat <- (time.treated - max(t)+1) * (dta[,first.treat.name]>0)
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
      Rb <- (base::colMeans(Ub*(cbind(dynamic.pre.e.inf.f, dynamic.e.inf.f)), na.rm = T))
      Rb
    })
    bres <- t(simplify2array(bout))
    # Non-degenerate dimensions
    ndg.dim <- base::colSums(bres)!=0
    #V.dynamic <- cov(bres)
    bres <- bres[,ndg.dim]


    # Part of the code for simulataneous confidence bands for event studies
    bSigma <- apply(bres, 2,
                    function(b) (quantile(b, .75, type=1,na.rm = T) -
                                   quantile(b, .25, type=1,na.rm = T))/(qnorm(.75) - qnorm(.25)))
    bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
    c.dynamic <- quantile(bT, 1-alp, type=1,na.rm = T)
    dynamic.se.e <- rep(0,length(ndg.dim))
    dynamic.se.e[ndg.dim] <- bSigma
  }



  AGGTE(simple.att=simple.att, simple.se=simple.se,
        dynamic.att=dynamic.att, dynamic.se=dynamic.se,
        dynamic.att.e=c(dynamic.att.pre.e, dynamic.att.e), dynamic.se.e=dynamic.se.e,  c.dynamic=c.dynamic,
        # Influence functions
        dyn.inf.func.e = cbind(dynamic.pre.e.inf.f, dynamic.e.inf.f),
        simple.att.inf.func = simple.att.inf.func,
        dynamic.att.inf.func = dynamic.if,

        group=originalflist, times=originaltlist,
        e = c(eseq.pre, eseq))
}





#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#' @title att_gt_het2
#'
#' @description \code{att_gt_het} computes the difference of ATT across two subpopulations in the case where there are more
#'  than two periods of data and allowing for treatment to occur at different points in time
#'  extending the method of Abadie (2005).  This method relies on once individuals are treated
#'  they remain in the treated state for the duration.
#'
#' @param outcome The outcome y (in quotations, always!)
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param aggte boolean for whether or not to compute aggregate treatment effect parameters, default TRUE
#' @param w The name of the column containing the weights
#' @param idname The individual (cross-sectional unit) id name
#' @param first.treat.name The name of the variable in \code{data} that contains the first
#'  period when a particular observation is treated.  This should be a positive
#'  number for all observations in treated groups.  It should be 0 for observations
#'  in the untreated group.
#' @param alp the significance level, default is 0.05
#' @param method The method for estimating the propensity score when covariates
#'  are included
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier boostrap.  If standard errors are clustered, then one
#'  must set \code{bstrap=TRUE}.
#' @param biters The number of boostrap iterations to use.  The default is 100,
#'  and this is only applicable if \code{bstrap=TRUE}.
#' @param clustervars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level.
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors#' (not implemented)
#' @param pl Boolean for whether or not to use parallel processing
#' @param cores The number of cores to use for parallel processing
#' @param printdetails Boolean for showing detailed results or not
#' @param maxe maximum values of periods ahead to be computed in event study
#' @param mine minimum values of periods before treatment to be computed in event study
#' @param nevertreated Boolean for using the group which is never treated in the sample as the comparison unit. Default is TRUE.
#' @param het The name of the column containing the (binary) categories for heterogeneity
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#' @return \code{MP} object
#'
#' @export
att_gt_het2 <- function(outcome, data, tname,
                        aggte=TRUE, w=NULL,
                        idname=NULL, first.treat.name, alp=0.05,
                        method="logit",
                        bstrap=FALSE, biters=1000, clustervars=NULL,

                        seedvec=NULL, pl=FALSE, cores=2,
                        printdetails=TRUE,
                        maxe = NULL,
                        mine = NULL,
                        nevertreated = T,
                        het) {

  ## make sure that data is a data.frame
  ## this gets around RStudio's default of reading data as tibble
  if (!all( class(data) == "data.frame")) {
    #warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  # weights if null
  if(is.character(w))  w <- data[, as.character(w)]
  if(is.null(w)) {
    w <- as.vector(rep(1, nrow(data)))
  } else if(min(w) < 0) stop("'w' must be non-negative")
  # het if null
  if(is.null(het)) {
    stop("Please specifiy 'het'. If het=NULL, use 'att_gt' instead of 'att_gt_het'.")
  }
  if(is.character(het))  het <- data[, as.character(het)]

  het.dim <- length(unique(het))
  if(het.dim!=2)  {
    stop("'het' must be a binary variable.")
  }
  data$w <- w
  data$w1 <- w * (het==1)
  data$w0 <- w * (het==0)

  data$y <- data[, as.character(outcome)] ##data[,as.character(formula.tools::lhs(formla))]
  ##figure out the dates and make balanced panel
  tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

  flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
  if ( length(flist[flist==0]) == 0) {
    warning("dataset does not have any observations in the control group.  make sure to set data[,first.treat.name] = 0 for observations in the control group.")
  }
  # First treated groups
  flist <- flist[flist>0]

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
  tlen <- length(tlist)
  # How many treated groups
  flen <- length(flist)

  data <- BMisc::makeBalancedPanel(data, idname, tname)
  #dta is used to get a matrix of size n (like in cross sectional data)
  dta <- data[ data[,tname]==tlist[1], ]  ## use this for the influence function



  #################################################################
  #################################################################

  #################################################################
  #----------------------------------------------------------------------------
  # Results for het==1
  #----------------------------------------------------------------------------
  results_het1 <- compute.att_gt_het2(flen, tlen, flist, tlist, data, dta, first.treat.name,
                                      outcome, tname, idname, method, seedvec,
                                      pl, cores, printdetails, nevertreated, het=1)
  fatt_het1 <- results_het1$fatt
  inffunc_het1 <- results_het1$inffunc
  #----------------------------------------------------------------------------
  # Results for het==0
  #----------------------------------------------------------------------------
  results_het0 <- compute.att_gt_het2(flen, tlen, flist, tlist, data, dta, first.treat.name,
                                      outcome, tname, idname, method, seedvec,
                                      pl, cores, printdetails, nevertreated, het=0)
  fatt_het0 <- results_het0$fatt
  inffunc_het0 <- results_het0$inffunc
  #----------------------------------------------------------------------------

  ## process the results from computing the spatt
  group <- c()
  t    <- c()
  att_het1 <- c()
  att_het0 <- c()
  i <- 1

  inffunc1_het1 <- matrix(0, ncol=flen*(tlen), nrow=nrow(dta)) ## note, this might not work in unbalanced case
  inffunc1_het0 <- matrix(0, ncol=flen*(tlen), nrow=nrow(dta))


  for (f in 1:length(flist)) {
    for (s in 1:(length(tlist))) {
      group[i] <- fatt_het1[[i]]$group
      t[i] <- fatt_het1[[i]]$year

      att_het1[i] <- fatt_het1[[i]]$att
      att_het0[i] <- fatt_het0[[i]]$att

      inffunc1_het1[,i] <- inffunc_het1[f,s,]
      inffunc1_het0[,i] <- inffunc_het0[f,s,]

      i <- i+1
    }
  }

  # THIS IS ANALOGOUS TO CLUSTER ROBUST STD ERRORS (in our specific setup)
  n <- nrow(dta)
  V <- NULL

  aggeffects <- NULL
  aggeffects_het1 <- NULL
  aggeffects_het0 <- NULL

  if (aggte) {
    aggeffects_het1 <- compute.aggte_het(flist, tlist, group, t, att_het1, first.treat.name, inffunc1_het1,
                                         n, clustervars, dta, idname, bstrap, biters, alp, maxe, mine, het=1)

    aggeffects_het0 <- compute.aggte_het(flist, tlist, group, t, att_het0, first.treat.name, inffunc1_het0,
                                         n, clustervars, dta, idname, bstrap, biters, alp, maxe, mine,het=0)


    aggeffects <- list(simple.att = aggeffects_het1$simple.att - aggeffects_het0$simple.att,
                       simple.att.inf.func = aggeffects_het1$simple.att.inf.func - aggeffects_het0$simple.att.inf.func,

                       dynamic.att = aggeffects_het1$dynamic.att - aggeffects_het0$dynamic.att,
                       dynamic.att.inf.func = aggeffects_het1$dynamic.att.inf.func - aggeffects_het0$dynamic.att.inf.func,

                       dynamic.att.e = aggeffects_het1$dynamic.att.e - aggeffects_het0$dynamic.att.e,
                       dyn.inf.func.e = aggeffects_het1$dyn.inf.func.e - aggeffects_het0$dyn.inf.func.e,
                       e = aggeffects_het1$e
    )


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
        bres <- as.vector(simplify2array(bout))

        bSigma <- (quantile(bres, .75, type=1, na.rm = T) - quantile(bres, .25, type=1, na.rm = T)) /
          (qnorm(.75) - qnorm(.25))
        #bT <- abs(bres/bSigma))
        return(bSigma)
        #return(sqrt( mean( bres^2)) /sqrt(n))
      } else {
        return(sqrt( mean( (thisinffunc)^2 ) /n ))
      }
    }



    aggeffects$simple.se <- getSE_inf(as.matrix(aggeffects$simple.att.inf.func))
    aggeffects$dynamic.se <- getSE_inf(as.matrix(aggeffects$dynamic.att.inf.func))

    aggeffects$dynamic.se.e <- sqrt(colMeans((aggeffects$dyn.inf.func.e)^2)/n)
    aggeffects$c.dynamic <- qnorm(1 - alp/2)

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
        Rb <- (base::colMeans(Ub*(aggeffects$dyn.inf.func.e), na.rm = T))
        Rb
      })
      bres <- t(simplify2array(bout))
      # Non-degenerate dimensions
      ndg.dim <- base::colSums(bres)!=0
      #V.dynamic <- cov(bres)
      bres <- bres[,ndg.dim]

      #V.dynamic <- cov(bres)
      bSigma <- apply(bres, 2, function(b) (quantile(b, .75, type=1,na.rm = T) - quantile(b, .25, type=1,na.rm = T))/(qnorm(.75) - qnorm(.25)))
      bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
      aggeffects$c.dynamic <- quantile(bT, 1-alp, type=1,na.rm = T)
      aggeffects$dynamic.se.e <- rep(0,length(ndg.dim))
      aggeffects$dynamic.se.e[ndg.dim] <- bSigma

    }


  }


  out <- list(group=group,
              t=t,
              att_het1=att_het1, att_het0=att_het0,
              inffunc_het1=inffunc_het1, inffunc_het0=inffunc_het0,
              n=n,
              aggte_het1=aggeffects_het1, aggte_het0=aggeffects_het0,
              aggte = aggeffects,
              alp = alp)

  return(out)
}




#' @title compute.att_gt_het2
#'
#' @description \code{compute.att_gt_het} does the main work for computing
#'  mutliperiod group-time average treatment effects with heterogeneity
#'
#' @inheritParams att_gt_het2
#'
#' @return a list with length equal to the number of groups times the
#'  number of time periods; each element of the list contains a \code{QTE}
#'  object that contains group-time average treamtent effect as well
#'  as which group it is for and which time period it is for and
#'  the influence function which is used externally to compute
#'  standard errors.
#'
#' @keywords internal
#'
#' @export
compute.att_gt_het2 <- function(flen, tlen, flist, tlist, data, dta,
                                first.treat.name, outcome, tname, idname,
                                method, seedvec,
                                pl, cores, printdetails, nevertreated, het) {

  yname <- outcome ##as.character(formula.tools::lhs(formla))

  fatt <- list()
  counter <- 1
  inffunc <- array(data=0, dim=c(flen,tlen,nrow(dta)))

  het.val <- het

  for (f in 1:flen) {
    ##satt <- list()
    for (t in 1:(tlen)) {
      #pret <- t
      pret <- utils::tail(which(tlist < flist[f]),1)
      if (flist[f]<=tlist[(t)]) {
        ## set an index for the pretreatment period
        #pret <- utils::tail(which(tlist < flist[f]),1)

        ## print a warning message if there are no pre-treatment
        ##  periods
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", flist[f]))
        }

        ## print the details of which iteration we are on
        if (printdetails) {
          cat(paste("current period:", tlist[t]), "\n")
          cat(paste("current group:", flist[f]), "\n")
          cat(paste("set pretreatment period to be", tlist[pret]), "\n")
        }
      }

      ## --------------------------------------------------------
      post.treat <- 1*(flist[f]<=tlist[(t)])
      if (tlist[pret] == tlist[t]){
        fatt[[counter]] <- list(att=0, group=flist[f], year=tlist[(t)], post=post.treat)
        inffunc[f,t,] <- 0
      } else {

      ## results for the case with panel data
      ## get dataset with current period and pre-treatment period
      disdat <- data[(data[,tname]==tlist[t] | data[,tname]==tlist[pret]),]
      ## transform it into "cross-sectional" data where
      ## one of the columns contains the change in the outcome
      ## over time
      ########################## TOASK
      # SHOULD we keep the data on t=t_pre or t=t_post? right now, it is t=t_pre
      disdat <- BMisc::panel2cs(disdat, yname, idname, tname)

      #THIS IS THE PART WE CAN CHANGE FOR THE NOT YET TREATED!!
      ## set up control group
      if(nevertreated ==T){
        disdat$C <- 1*(disdat[,first.treat.name] == 0)
      }

      if(nevertreated ==F){
        disdat$C <- 1*( (disdat[,first.treat.name] == 0) +
                          (disdat[,first.treat.name] > max(disdat[, tname])) )
      }

      ## set up for particular treated group
      disdat$G <- 1*(disdat[,first.treat.name] == flist[f])

      ## drop missing factors
      disdat <- droplevels(disdat)

      ## give short names for data in this iteration
      G <- disdat$G
      C <- disdat$C
      dy <- disdat$dy * ((-1)^(1+post.treat))
      #n <- nrow(disdat)
      whet <- (disdat$w1) * het.val + (disdat$w0) * (1 - het.val)
      w <- disdat$w


      ## set up weights
      # Treated get the weights based on whet
      attw <- whet * G/mean(whet * G)
      # Comparison group gets weights that does not vary with het
      attw2a <- w * C
      attw2 <- attw2a/mean(attw2a)
      att <- mean((attw - attw2)*dy)

      if(is.na(att)) att <- 0

      ## save results for this iteration
      fatt[[counter]] <- list(att=att, group=flist[f], year=tlist[(t)], post=post.treat)

      ## --------------------------------------------
      ## get the influence function

      ## weigts for het==1
      wg <- whet * G/mean(whet * G)
      wc1 <- w * C
      wc <- wc1 / mean(wc1)

      ## influence function for treated group
      psig <- wg*(dy - mean(wg*dy))
      # influence function for the control group
      psic <- wc*(dy - mean(wc*dy))

      ## save the influnce function as the difference between
      ## the treated and control influence functions;
      ## we save this as a 3-dimensional array
      ## and then process afterwards
      infl.att <- (psig - psic)
      if(base::anyNA(infl.att)) infl.att <- 0

      inffunc[f,t,] <- infl.att
      }

      counter <- counter+1
    }

  }

  list(fatt=fatt, inffunc=inffunc)
}
