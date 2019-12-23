#' @title mp.spatt
#'
#' @description \code{mp.spatt} computes the ATT in the case where there are more
#'  than two periods of data and allowing for treatment to occur at different points in time
#'  extending the method of Abadie (2005).  This method relies on once individuals are treated
#'  they remain in the treated state for the duration.
#'
#' @param outcome The outcome y (in quotations, always!)
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param aggte boolean for whether or not to compute aggregate treatment effect parameters, default TRUE
#' @param w A vector of weights for each observation (not implemented)
#' @param idname The individual (cross-sectional unit) id name
#' @param first.treat.name The name of the variable in \code{data} that contains the first
#'  period when a particular observation is treated.  This should be a positive
#'  number for all observations in treated groups.  It should be 0 for observations
#'  in the untreated group.
#' @param alp the significance level, default is 0.05
#' @param method The method for estimating the propensity score when covariates
#'  are included
#' @param se Boolean whether or not to compute standard errors
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier boostrap.  If standard errors are clustered, then one
#'  must set \code{bstrap=TRUE}.
#' @param biters The number of boostrap iterations to use.  The default is 100,
#'  and this is only applicable if \code{bstrap=TRUE}.
#' @param clustervars A vector of variables to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level.
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability \code{1-alp}.  The default is \code{FALSE}
#'  and the resulting standard errors will be pointwise.
#' @param citers Computing uniform confidence bands requires the bootstrap,
#'  if \code{cband = TRUE}, then this is the number of boostrap iterations
#'  to compute the conidence band.  The default is 100.
#' @param seedvec Optional value to set random seed; can possibly be used
#'  in conjunction with bootstrapping standard errors#' (not implemented)
#' @param pl Boolean for whether or not to use parallel processing
#' @param cores The number of cores to use for parallel processing
#' @param printdetails Boolean for showing detailed results or not
#' @param maxe maximum values of periods ahead to be computed in event study
#' @param nevertreated Boolean for using the group which is never treated in the sample as the comparison unit. Default is TRUE.
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#'
#' @return \code{MP} object
#'
#' @export
mp.spatt <- function(outcome, data, tname,
                     aggte=TRUE, w=NULL,
                     idname=NULL, first.treat.name, alp=0.05,
                     method="logit", se=TRUE,
                     bstrap=FALSE, biters=1000, clustervars=NULL,
                     cband=FALSE, citers=1000,
                     seedvec=NULL, pl=FALSE, cores=2,
                     printdetails=TRUE,
                     maxe = NULL,
                     nevertreated = T) {

  ## make sure that data is a data.frame
  ## this gets around RStudio's default of reading data as tibble
  if (!all( class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  # weights if null
  if(is.character(w))  w <- data[, as.character(w)]
  if(is.null(w)) {
    w <- as.vector(rep(1, nrow(data)))
  } else if(min(w) < 0) stop("'w' must be non-negative")
  data$w <- w

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
  ## more error handling after we have balanced the panel

  # Size of each group (divide by lenth(tlist because we are working with long data))
  #gsize2 <- aggregate(data[,first.treat.name], by=list(data[,first.treat.name]),
  #                   function(x) length(x)/length(tlist))
  gsize <- aggregate(data[,"w"], by=list(data[,first.treat.name]),
                     function(x) sum(x)/length(tlist))


  #################################################################

  #################################################################
  results <- compute.mp.spatt(flen, tlen, flist, tlist, data, dta, first.treat.name,
                              outcome, tname, w, idname, method, seedvec, se,
                              pl, cores, printdetails, nevertreated)

  fatt <- results$fatt
  inffunc <- results$inffunc

  ## process the results from computing the spatt
  group <- c()
  t    <- c()
  att <- c()
  i <- 1

  inffunc1 <- matrix(0, ncol=flen*(tlen-1), nrow=nrow(dta)) ## note, this might not work in unbalanced case
  ## for (f in 1:length(fatt)) {
  ##     for (s in 1:(length(fatt[[f]])-1)) {
  ##         group[i] <- fatt[[f]]$group
  ##         t[i] <- fatt[[f]][[s]]$year
  ##         att[i] <- fatt[[f]][[s]]$att
  ##         inffunc1[,i] <- inffunc[f,s,]
  ##         i <- i + 1
  ##     }
  ## }

  for (f in 1:length(flist)) {
    for (s in 1:(length(tlist)-1)) {
      group[i] <- fatt[[i]]$group
      t[i] <- fatt[[i]]$year
      att[i] <- fatt[[i]]$att
      inffunc1[,i] <- inffunc[f,s,]
      i <- i+1
    }
  }

  # THIS IS ANALOGOUS TO CLUSTER ROBUST STD ERRORS (in our specific setup)
  n <- nrow(dta)
  V <- t(inffunc1)%*%inffunc1/n

  if ( (length(clustervars) > 0) & !bstrap) {
    warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  }

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
      Rb <- sqrt(n)*(apply(Ub*(inffunc1), 2, mean))
      Rb
    })
    bres <- t(simplify2array(bout))
    V <- cov(bres)
  }


  ## TODO: handle case with repeated cross sections; this part is conceptually easier because many
  ##  off-diagonal (though not all) will be 0.

  ## get the actual estimates

  ## wald test for pre-treatment periods:
  #THIS NEED TO COME BEFORE THE OTHER CODE OF THE CONFIDENCE BANDS!!
  pre <- which(t < group)
  preatt <- as.matrix(att[pre])
  preV <- V[pre,pre]


  ## new code
  cval <- qnorm(1-alp/2)
  if (cband) {
    bSigma <- apply(bres, 2, function(b) (quantile(b, .75, type=1, na.rm = T) - quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
    bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
    cval <- quantile(bT, 1-alp, type=1, na.rm = T)
    ##bT1 <- apply(bres, 1, function(b) max( abs(b)*diag(V)^(-.5) ))
    ##cval1 <- quantile(bT1, 1-alp, type=1)
    V <- diag(bSigma^2) ## this is the appropriate matrix for
    ## constructing confidence bands - BUT NOT FOR WALD TEST!!!
  }

  aggeffects <- NULL
  if (aggte) {
    aggeffects <- compute.aggte(flist, tlist, group, t, att, first.treat.name, inffunc1,
                                n, clustervars, dta, idname, bstrap, biters, alp, cband, maxe)
  }


  if (length(preV) == 0) {
    message("No pre-treatment periods to test")
    return(MP(group=group, t=t, att=att, V=V, c=cval, inffunc=inffunc1, n=n, aggte=aggeffects, alp = alp))
  }

  if (det(preV) == 0) { ##matrix not invertible
    warning("Covariance matrix is singular. We do not report pre-test Wald statistic for pre-trends.")
    return(MP(group=group, t=t, att=att, V=V, c=cval, inffunc=inffunc1, n=n, aggte=aggeffects, alp = alp))
  } else Vinv <- solve(preV)

  W <- n*t(preatt)%*%Vinv%*%preatt
  q <- length(pre)##sum(1-as.numeric(as.character(results$post))) ## number of restrictions
  Wpval <- round(1-pchisq(W,q),5)


  return(MP(group=group, t=t, att=att, V=V, c=cval, inffunc=inffunc1, n=n, W=W, Wpval=Wpval, aggte=aggeffects,
            alp = alp))
}



#' @title compute.mp.spatt
#'
#' @description \code{compute.mp.spatt} does the main work for computing
#'  mutliperiod group-time average treatment effects
#'
#' @inheritParams mp.spatt
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
compute.mp.spatt <- function(flen, tlen, flist, tlist, data, dta,
                             first.treat.name, outcome, tname, w, idname,
                             method, seedvec, se,
                             pl, cores, printdetails, nevertreated) {

  yname <- outcome ##as.character(formula.tools::lhs(formla))

  fatt <- list()
  counter <- 1
  inffunc <- array(data=0, dim=c(flen,tlen,nrow(dta)))


  for (f in 1:flen) {
    ##satt <- list()
    for (t in 1:(tlen-1)) {
      pret <- t
      if (flist[f]<=tlist[(t+1)]) {
        ## set an index for the pretreatment period
        pret <- utils::tail(which(tlist < flist[f]),1)

        ## print a warning message if there are no pre-treatment
        ##  periods
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", flist[f]))
        }

        ## print the details of which iteration we are on
        if (printdetails) {
          cat(paste("current period:", tlist[t+1]), "\n")
          cat(paste("current group:", flist[f]), "\n")
          cat(paste("set pretreatment period to be", tlist[pret]), "\n")
        }
      }

      ## --------------------------------------------------------
      ## results for the case with panel data
      ## get dataset with current period and pre-treatment period
      disdat <- data[(data[,tname]==tlist[t+1] | data[,tname]==tlist[pret]),]
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
      dy <- disdat$dy
      n <- nrow(disdat)
      w <- disdat$w

      ## set up weights
      attw <- w * G/mean(w * G)
      attw2a <- w * C
      attw2 <- attw2a/mean(attw2a)
      att <- mean((attw - attw2)*dy)
      ## save results for this iteration
      fatt[[counter]] <- list(att=att, group=flist[f], year=tlist[(t+1)], post=1*(flist[f]<=tlist[(t+1)]))

      ## --------------------------------------------
      ## get the influence function

      ## weigts
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
      inffunc[f,t,] <- psig - psic

      counter <- counter+1
    }

  }

  list(fatt=fatt, inffunc=inffunc)
}

#' @title MP
#'
#' @description multi-period object
#'
#' @param group which group (defined by period first treated) an group-time average treatment effect is for
#' @param t which time period a group-time average treatment effect is for
#' @param att the group-average treatment effect for group \code{group} and time period \code{t}
#' @param c critical value if one is obtaining uniform confidence bands
#' @param V the variance matrix for group-time average treatment effects
#' @param inffunc the influence function for estimating group-time average treatment effects
#' @param n the number of observations
#' @param W the Wald statistic for pre-testing the common trends assumption
#' @param Wpval the p-value of the Wald statistic for pre-testing the
#'  common trends assumption
#' @param aggte an aggregate treatment effects object
#' @param alp the significance level, default is 0.05
#'
#' @return MP object
#' @export
MP <- function(group, t, att, V, c, inffunc, n=NULL, W=NULL, Wpval=NULL, aggte=NULL, alp = 0.05) {
  out <- list(group=group, t=t, att=att, V=V, c=c, inffunc=inffunc, n=n, W=W, Wpval=Wpval, aggte=aggte, alp = alp)
  class(out) <- "MP"
  out
}

#' @title summary.MP
#'
#' @description prints a summary of a \code{MP} object
#'
#' @param object an \code{MP} object
#' @param ... extra arguments
#'
#' @export
summary.MP <- function(object, ...) {
  mpobj <- object
  out <- cbind(mpobj$group, mpobj$t, mpobj$att, sqrt(diag(mpobj$V)/mpobj$n))
  citation()
  colnames(out) <- c("group", "time", "att","se")
  cat("\n")
  print(kable(out))
  cat("\n\n")
  cat("P-value for pre-test of DID assumption:  ")
  cat(as.character(mpobj$Wpval))
  cat("\n\n")
}


#' @title gplot
#'
#' @description does the heavy lifting for making a plot of an group-time
#'  average treatment effect
#'
#' @inheritParams ggdid
#'
#' @return a \code{ggplot2} object
#'
#' @keywords internal
#'
#' @export
gplot <- function(ssresults, ylim=NULL, xlab=NULL, ylab=NULL, title="Group", xgap=1) {
  dabreaks <- ssresults$year[seq(1, length(ssresults$year), xgap)]

  c.point <-  stats::qnorm(1 - ssresults$alp/2)

  p <- ggplot(ssresults,
              aes(x=year, y=att, ymin=(att-c*ate.se),
                  ymax=att+c*ate.se, post=post)) +

    geom_point(aes(colour=post), size=1.5) +
    geom_errorbar(aes(colour=post), width=0.1) +
    #geom_ribbon(aes(ymin= (att-c.point*ate.se), ymax=  (att+c.point*ate.se), alpha=0.35))+
    #geom_ribbon(aes(ymin=  (att-c*ate.se), ymax =  (att+c*ate.se), alpha=0.25))+

    scale_y_continuous(limits=ylim) +
    scale_x_discrete(breaks=dabreaks, labels=as.character(dabreaks)) +
    scale_colour_hue(drop=FALSE) +
    ylab("") +
    xlab("") +
    ggtitle(paste(title, unique(ssresults$group))) +
    theme_bw() +
    theme(plot.title = element_text(color="darkgray", face="bold", size=8)) +
    theme(axis.title = element_text(color="darkgray", face="bold", size=8))
  p
}

#' @title ggdid
#'
#' @description Function to plot \code{MP} objects
#'
#' @param mpobj an \code{MP} object
#' @param type the type of plot, should be one of "attgt", "dynamic",
#'  "selective", "calendar", "dynsel".  "attgt" is the default and plots
#'  all group-time average treatment effects separately by group (including
#'  pre-treatment time periods); "dynamic" plots dynamic treatment effects --
#'  these are the same as event studies; "selective" plots average effects
#'  of the treatment separately by group (which allows for selective treatment
#'  timing); "calendar" plots average treatment effects by time period; and
#'  "dynsel" plots dynamic effects allowing for selective treatment timing
#'  (this also requires setting the additional paramater e1)
#' @param ylim optional y limits for the plot; settng here makes the y limits
#'  the same across different plots
#' @param xlab optional x-axis label
#' @param ylab optional y-axis label
#' @param title optional plot title
#' @param xgap optional gap between the labels on the x-axis.  For example,
#'  \code{xgap=3} indicates that the labels should show up for every third
#'  value on the x-axis.  The default is 1.
#' @param ncol The number of columns to include in the resulting plot.  The
#'  default is 1.
#' @param e1 only used when plot type is "dynsel", this specifies the number
#'  of post-treatment periods that need to be available for particular groups
#'  to be included in the resulting plot when there are dynamic treatment
#'  effects and selective treatment timing
#'
#'
#' @export
ggdid <- function(mpobj, type=c("attgt", "dynamic"), ylim=NULL,
                  xlab=NULL, ylab=NULL, title="Group", xgap=1, ncol=1, e1=1) {

  type <- type[1]

  G <- length(unique(mpobj$group))
  Y <- length(unique(mpobj$t))## drop 1 period bc DID
  g <- unique(mpobj$group)[order(unique(mpobj$group))] ## -1 to drop control group
  y <- unique(mpobj$t)

  if (type=="attgt") {
    results <- data.frame(year=rep(y,G))
    results$group <- unlist(lapply(g, function(x) { rep(x, Y) }))##c(rep(2004,G),rep(2006,G),rep(2007,G))
    results$att <- mpobj$att
    n <- mpobj$n
    results$ate.se <- sqrt(diag(mpobj$V)/n)
    results$post <- as.factor(1*(results$year >= results$group))
    results$year <- as.factor(results$year)
    results$c <- mpobj$c
    vcovatt <- mpobj$V/n
    alp <- mpobj$alp

    ##results <- mp2ATT(results, vcovatt)

    mplots <- lapply(g, function(g) {
      thisdta <- subset(results, group==g)
      gplot(thisdta, ylim, xlab, ylab, title, xgap)
    })

    do.call("grid.arrange", c(mplots))
  } else if (type=="dynamic") {
    aggte <- mpobj$aggte
    #if (mpobj$c > 2) warning("uniform bands not implemented yet for this plot...")
    elen <- length(aggte$dynamic.att.e)
    results <- cbind.data.frame(year=as.factor(seq(1:elen)),
                                att=aggte$dynamic.att.e,
                                ate.se=aggte$dynamic.se.e,
                                post=as.factor(1),
                                c=aggte$c.dynamic,
                                alp = mpobj$alp)
    #qnorm(.975))
    p <- gplot(results, ylim, xlab, ylab, title, xgap)
    p
  }
}


#' @title compute.aggte
#'
#' @description does the heavy lifting on computing aggregated group-time
#'  average treatment effects
#'
#' @inheritParams mp.spatt
#'
#' @return \code{AGGTE} object
#'
#' @keywords internal
#'
#' @export
compute.aggte <- function(flist, tlist, group, t, att, first.treat.name, inffunc1, n,
                          clustervars, dta, idname, bstrap, biters, alp, cband, maxe) {

  if ( (length(clustervars) > 0) & !bstrap) {
    warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
  }

  ## internal function for computing standard errors
  ##  this method is used across different types of
  ##  aggregate treatment effect parameters and is just
  ##  based on using the right influence function and weights
  ##  -- these are specific to which aggregate treatment
  ##  effect parameter is being considered.
  ## @param wif is the influence function for the weights
  getSE <- function(whichones, weights, wif=NULL) {
    weights <- as.matrix(weights) ## just in case pass vector
    thisinffunc <- inffunc1[,whichones]%*%weights  ##multiplies influence function times weights and sums to get vector of weighted IF (of length n)
    if (!is.null(wif)) {
      thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
    }

    if (bstrap) {
      bout <- lapply(1:biters, FUN=function(b) {
        sercor <- idname %in% clustervars ## boolean for whether or not to account for serial correlation
        clustervars <- clustervars[-which(clustervars==idname)]
        if (length(clustervars) > 1) {
          stop("can't handle that many cluster variables")
        }
        if (length(clustervars) > 0) {
          n1 <- length(unique(dta[,clustervars]))
          Vb <- matrix(sample(c(-1,1), n1, replace=TRUE),
                       nrow=n1)
          Vb <- cbind.data.frame(unique(dta[,clustervars]), Vb)
          Ub <- data.frame(dta[,clustervars])
          Ub <- Vb[match(Ub[,1], Vb[,1]),]
          Ub <- Ub[,-1]
          Ub <- as.matrix(Ub)
          ## n1 <- length(unique(dta[,clustervars]))
          ## Vb <- matrix(sample(c(-1,1), n1*ncol(inffunc1), replace=T),
          ##              nrow=n1)
          ## Vb <- cbind.data.frame(unique(dta[,clustervars]), Vb)
          ## colnames(Vb)[1] <- "clvar"
          ## Ub <- data.frame(dta[,clustervars])
          ## colnames(Ub)[1] <- "clvar"
          ## Ub <- merge(Ub, Vb, by="clvar")
          ## Ub <- Ub[,-1]
        } else {
          Ub <- matrix(sample(c(-1,1), nrow(thisinffunc), replace=TRUE), ncol=1)

        }
        ## allow for serial correlation
        ##if (sercor) {
        ##    Ub[,-1] <- Ub[,1]
        ##} ## this doesn't matter here because there is only one influence function
        ## drop cluster for serial correlation

        ##ift <- do.call(magic::adiag, psiitout)
        ##ifunc <- rbind(ift, psiiu)

        ##Ub <- sample(c(-1,1), n, replace=T)
        mb <- Ub*(thisinffunc)
        apply(mb,2,sum)/sqrt(nrow(dta))
      })
      bres <- simplify2array(bout)
      return(sqrt( mean( bres^2)) /sqrt(n))
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

  weights.agg = dta$w
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

  dynamic.se.e <- sqrt(base::colMeans(dynamic.e.inf.f^2)/n)

  # Average over all dynamic att(e)
  dynamic.att <- mean(dynamic.att.e)
  # Influence function of the average over all dynamic att(e)
  dynamic.if <- rowMeans(dynamic.e.inf.f)
  # std. error
  dynamic.se <- sqrt((mean(dynamic.if^2)/n))

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
      Rb <- sqrt(n)*(base::colMeans(Ub*(dynamic.e.inf.f), na.rm = T))
      Rb
    })
    bres <- t(simplify2array(bout))
    #V.dynamic <- cov(bres)
  }


  ## new code
  c.dynamic <- qnorm(1 - alp/2)

  if (cband) {
    bSigma <- apply(bres, 2, function(b) (quantile(b, .75, type=1,na.rm = T) - quantile(b, .25, type=1,na.rm = T))/(qnorm(.75) - qnorm(.25)))
    bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
    c.dynamic <- quantile(bT, 1-alp, type=1,na.rm = T)
    ##bT1 <- apply(bres, 1, function(b) max( abs(b)*diag(V)^(-.5) ))
    ##cval1 <- quantile(bT1, 1-alp, type=1)
    dynamic.se.e <- bSigma/sqrt(n)
    ## constructing confidence bands
  }

  AGGTE(simple.att=simple.att, simple.se=simple.se,
        dynamic.att=dynamic.att, dynamic.se=dynamic.se,
        dynamic.att.e=dynamic.att.e, dynamic.se.e=dynamic.se.e,  c.dynamic=c.dynamic,
        groups=originalflist,times=originaltlist)
}


#' @title AGGTE
#'
#' @description \code{AGGTE} class for aggregate treatment effects
#'
#' @param simple.att simple weighted average of group-time average treatment
#'  effects
#' @param simple.se the standard error for \code{simple.att}
#' @param dynamic.att aggregated group-time average treatment effects when
#'  there are dynamic treatment effects
#' @param dynamic.se the standard error for \code{dynamic.att}
#' @param dynamic.att.e aggregated group-time average treatment effects
#'  when there are dynamic treatment effects for each length of exposure
#'  to treatment
#' @param dynamic.se.e the standard error for \code{dynamic.att.e}
#' @param c.dynamic the (simultaneous) critical value \code{dynamic.att.e}
#' @param groups vector of all groups
#' @param times vector of all times
AGGTE <- function(simple.att=NULL, simple.se=NULL,
                  dynamic.att=NULL, dynamic.se=NULL,
                  dynamic.att.e=NULL, dynamic.se.e=NULL, c.dynamic=NULL,
                  groups=NULL, times=NULL) {

  out <- list(simple.att=simple.att, simple.se=simple.se,
              dynamic.att=dynamic.att, dynamic.se=dynamic.se,
              dynamic.att.e=dynamic.att.e, dynamic.se.e=dynamic.se.e, c.dynamic=c.dynamic,
              groups=groups,times=times)
  class(out) <- "AGGTE"
  out
}


#' @title summary.AGGTE
#'
#' @description print a summary of an AGGTE object
#'
#' @param object an AGGTE object
#' @param type which type of summary to print, options are "dynamic", "selective", "calendar", and "dynsel"
#' @param e1 if the type is "dynsel", this is the number of post-treatment periods required in order for a group to be used to construct aggregated parameters with selective treatment timing and dynamic effects; otherwise not used
#' @param ... other variables
#'
#' @export
summary.AGGTE <- function(object, type=c("dynamic"), e1=1, ...) {
  citation()
  sep <- "          "
  cat("Overall Summary Measures", "\n")
  cat("------------------------", "\n")
  cat("Simple ATT    : ", object$simple.att, "\n")
  cat("  SE          : ", object$simple.se, "\n")
  cat("Dynamic ATT   : ", object$dynamic.att, "\n")
  cat("  SE          : ", object$dynamic.se, "\n")
  type <- type[1]
  if (type == "dynamic") {
    cat("Dynamic Treatment Effects", "\n")
    cat("-------------------------")
    elen <- length(object$dynamic.att.e)
    printmat <- cbind(seq(1:elen), object$dynamic.att.e, object$dynamic.se.e)
    colnames(printmat) <- c("e","att","se")
    print(kable(printmat))
  }
}
