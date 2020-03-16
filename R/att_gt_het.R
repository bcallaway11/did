#' @title att_gt_het
#'
#' @description \code{att_gt_het} computes the difference of average treatment effects across two subpopulations
#' in DID setups where there are more than two periods of data and
#' allowing for treatment to occur at different points in time. Here, we allow for different trends
#' between the two supopulations. See Marcus ans Sant'Anna (2020) for a detailed description.
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
#' @param het The name of the column containing the (binary) categories for heterogeneity
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
#' @param pl Boolean for whether or not to use parallel processing (not implemented) is TRUE.
#' @param cores The number of cores to use for parallel processing (not implemented)
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#' @return \code{MP} object
#'
#' @export
att_gt_het <-function(outcome, data,
                      tname, idname=NULL,first.treat.name,
                      nevertreated = T,
                      het,
                      aggte=TRUE,
                      maxe = NULL,
                      mine = NULL,
                      w=NULL,
                      alp=0.05,
                      bstrap=T, biters=1000, clustervars=NULL,
                      cband=T,
                      printdetails=TRUE,
                      seedvec=NULL, pl=FALSE, cores=2,method="logit")
{
  #-------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------
  #                             Data pre-processing and error checking
  #-------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------
  ## make sure that data is a data.frame
  df <- data
  ## make sure that data is a data.frame
  ## this gets around RStudio's default of reading data as tibble
  if (!all( class(df) == "data.frame")) {
    #warning("class of data object was not data.frame; converting...")
    df <- as.data.frame(df)
  }
  # weights if null
  if(is.character(w))  w <- df[, as.character(w)]
  if(is.null(w)) {
    w <- as.vector(rep(1, nrow(df)))
  } else if(min(w) < 0) stop("'w' must be non-negative")
  # het if null
  if(is.null(het)) {
    stop("Please specifiy 'het'. If het=NULL, use 'att_gt' instead of 'att_gt_het'.")
  }
  if(is.character(het))  het <- df[, as.character(het)]

  het.dim <- length(unique(het))
  if(het.dim!=2)  {
    stop("'het' must be a binary variable.")
  }

  df$w1 <- w * (het==1)
  df$w0 <- w * (het==0)

  df$y <- df[, as.character(outcome)] ##data[,as.character(formula.tools::lhs(formla))]
  ##figure out the dates and make balanced panel
  tlist <- unique(df[,tname])[order(unique(df[,tname]))] ## this is going to be from smallest to largest

  flist <- unique(df[,first.treat.name])[order(unique(df[,first.treat.name]))]

  # Check if there is a never treated grup
  if ( length(flist[flist==0]) == 0) {
    if(nevertreated){
      stop("It seems you do not have a never-treated group in the data. If you do have a never-treated group in the data, make sure to set data[,first.treat.name] = 0 for the observation in this group. Otherwise, select nevertreated = F so you can use the not-yet treated units as a comparison group.")
    } else {
      warning("It seems like that there is not a never-treated group in the data. In this case, we cannot identity the ATT(g,t) for the group that is treated las, nor any ATT(g,t) for t higher than or equal to the largest g.\n \nIf you do have a never-treated group in the data, make sure to set data[,first.treat.name] = 0 for the observation in this group.")
      # Drop all time periods with time periods >= latest treated
      df <- base::subset(df,(df[,tname] < max(flist)))
      # Replace last treated time with zero
      lines.gmax = df[,first.treat.name]==max(flist)
      df[lines.gmax,first.treat.name] <- 0

      ##figure out the dates
      tlist <- unique(df[,tname])[order(unique(df[,tname]))] ## this is going to be from smallest to largest
      # Figure out the groups
      flist <- unique(df[,first.treat.name])[order(unique(df[,first.treat.name]))]
    }
  }

  # First treated groups
  flist <- flist[flist>0]

  ##################################
  ## do some error checking
  if (!is.numeric(tlist)) {
    warning("not guaranteed to order time periods correclty if they are not numeric")
  }

  ## check that first.treat doesn't change across periods for particular individuals
  if (!all(sapply( split(df, df[,idname]), function(df) {
    length(unique(df[,first.treat.name]))==1
  }))) {
    stop("Error: the value of first.treat must be the same across all periods for each particular individual.")
  }
  ####################################
  # How many time periods
  tlen <- length(tlist)
  # How many treated groups
  flen <- length(flist)

  df <- BMisc::makeBalancedPanel(df, idname, tname)
  #dta is used to get a matrix of size n (like in cross sectional data)
  dta <- df[ df[,tname]==tlist[1], ]  ## use this for the influence function

  #-------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------
  #                         Compute all ATT(g,t) for each het group
  #-------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------

   # Results for het==1
  #----------------------------------------------------------------------------
  results_het1 <- compute.att_gt_het(flen, tlen, flist, tlist, df, dta, first.treat.name,
                                     outcome, tname, idname, method, seedvec,
                                     pl, cores, printdetails, nevertreated, het=1)
  fatt_het1 <- results_het1$fatt
  inffunc_het1 <- results_het1$inffunc
  #----------------------------------------------------------------------------
  # Results for het==0
  #----------------------------------------------------------------------------
  results_het0 <- compute.att_gt_het(flen, tlen, flist, tlist, df, dta, first.treat.name,
                                     outcome, tname, idname, method, seedvec,
                                     pl, cores, printdetails, nevertreated, het=0)
  fatt_het0 <- results_het0$fatt
  inffunc_het0 <- results_het0$inffunc
  #----------------------------------------------------------------------------

  ## process the results from computing the att
  group <- c()
  tt    <- c()
  att_het1 <- c()
  att_het0 <- c()
  i <- 1

  inffunc1_het1 <- matrix(0, ncol=flen*(tlen), nrow=nrow(dta)) ## note, this might not work in unbalanced case
  inffunc1_het0 <- matrix(0, ncol=flen*(tlen), nrow=nrow(dta))


  for (f in 1:length(flist)) {
    for (s in 1:(length(tlist))) {
      group[i] <- fatt_het1[[i]]$group
      tt[i] <- fatt_het1[[i]]$year

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

  #-------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------
  #                         Compute all summaries of the ATT(g,t)
  #-------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------
  aggeffects <- NULL
  aggeffects_het1 <- NULL
  aggeffects_het0 <- NULL

  if (aggte) {
    aggeffects_het1 <- compute.aggte_het(flist, tlist, group, tt, att_het1, first.treat.name, inffunc1_het1,
                                         n, clustervars, dta, idname, bstrap, biters, alp, maxe, mine, het=1)

    aggeffects_het0 <- compute.aggte_het(flist, tlist, group, tt, att_het0, first.treat.name, inffunc1_het0,
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
        return(as.numeric(bSigma))
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
      bSigma <- apply(bres, 2, function(b) (quantile(b, .75, type=1,na.rm = T) -
                                              quantile(b, .25, type=1,na.rm = T))/(qnorm(.75) - qnorm(.25)))
      bT <- apply(bres, 1, function(b) max( abs(b/bSigma)))
      aggeffects$c.dynamic <- quantile(bT, 1-alp, type=1,na.rm = T)
      aggeffects$dynamic.se.e <- rep(0,length(ndg.dim))
      aggeffects$dynamic.se.e[ndg.dim] <- as.numeric(bSigma)
    }
  }

  # Return this list

  out <- list(group=group,
              t=tt,
              att_het1=att_het1, att_het0=att_het0,
              inffunc_het1=inffunc_het1, inffunc_het0=inffunc_het0,
              n=n,
              aggte_het1=aggeffects_het1, aggte_het0=aggeffects_het0,
              aggte = aggeffects,
              alp = alp)

  return(out)
}
