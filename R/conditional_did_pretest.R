#' @title Pre-Test of Conditional Parallel Trends Assumption
#'
#' @description An integrated moments test for the conditional parallel trends
#'  assumption holding in all pre-treatment time periods for all groups
#'
#' @inheritParams att_gt
#'
#' @examples
#' \dontrun{
#' data(mpdta)
#' pre.test <- conditional_did_pretest(yname="lemp",
#'                                     tname="year",
#'                                     idname="countyreal",
#'                                     gname="first.treat",
#'                                     xformla=~lpop,
#'                                     data=mpdta)
#' summary(pre.test)
#' }
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro H. C.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://arxiv.org/abs/1803.09015v2> (2018).
#'
#' @return an \code{\link{MP.TEST}} object
#' @export
conditional_did_pretest <- function(yname,
                                    tname,
                                    idname=NULL,
                                    gname,
                                    xformla=NULL,
                                    data,
                                    panel=TRUE,
                                    allow_unbalanced_panel=FALSE,
                                    control_group=c("nevertreated","notyettreated"),
                                    weightsname=NULL,
                                    alp=0.05,
                                    bstrap=TRUE,
                                    cband=TRUE,
                                    biters=1000,
                                    clustervars=NULL,
                                    est_method="ipw",
                                    print_details=FALSE,
                                    pl=FALSE,
                                    cores=1) {

  message("We are no longer updating this function.  It should continue to work, but most users find the pre-tests already reported by the `att_gt` function to be sufficient for most empirical applications.")

  # this is a DIDparams object
  dp <- pre_process_did(yname=yname,
                        tname=tname,
                        idname=idname,
                        gname=gname,
                        xformla=xformla,
                        data=data,
                        panel=panel,
                        allow_unbalanced_panel=allow_unbalanced_panel,
                        control_group=control_group,
                        weightsname=weightsname,
                        alp=alp,
                        bstrap=bstrap,
                        cband=cband,
                        biters=biters,
                        clustervars=clustervars,
                        est_method=est_method,
                        print_details=print_details,
                        pl=pl,
                        cores=cores
  )


  data <- dp$data
  tlist <- dp$tlist
  glist <- dp$glist
  n <- dp$n

  # check if possible to do test
  # note: tlist[2] contains the 2nd time period
  # (which is the first period where able to calculate ATT(g,t)'s)
  if ( max(glist) <= tlist[2] ) {
    stop("There are no pre-treatment periods to use to conduct test.")
  }


  # set which weight function to use
  # the only option that will work with current setup
  # is indicator so hard-code it here
  weightfun <- indicator


  if (allow_unbalanced_panel) {
    stop("Conditional pre-test not currently supported for unbalanced panel.")
  }

  # create dataset with n observations;
  # recover covariates from this dataset
  ifelse(panel & (allow_unbalanced_panel==FALSE),
         dta <- data[ data[,tname]==tlist[1], ],
         dta <- data
         )

  # need X twice, once to loop over, and once to calculate weighting function
  X <- model.matrix(xformla, dta)
  X1 <- X

  #for debugging:
  # X <- as.matrix(X[1:100,])

  cat("Step 1 of 2: Computing test statistic....\n")
  out <- pbapply::pblapply(1:nrow(X), function(i) {
    # these are the weights for the conditional moment test
    # indicator weights
    www <- as.numeric(weightfun(X1, X[i,]))
    # for indicator weights, just choose rows where weights = 1
    rightids <- dta[,idname][www==1]

    # create a new dataset and set the outcome to be the outcomes multiplied by the
    # weighting function (***this *trick* is only going to work for indicator
    # weights*** (otherwise you will multiply by www twice))
    thisdata <- data
    thisdata[,yname] <- 0
    thisdata[ thisdata[,idname] %in% rightids,yname] <- data[ thisdata[,idname] %in% rightids,yname]
    thisdata$y <- thisdata[,yname]

    # set new parameters to pass to call to compute.att_gt
    thisdp <- dp
    thisdp$data <- thisdata
    thisdp$print_details <- FALSE

    # compute the test statistic with call to compute.att_gt
    Jres <- compute.att_gt(thisdp)

    # turn in into more usable format
    J.results <- process_attgt(Jres$attgt.list)
    group <- J.results$group
    att <- J.results$att
    tt <- J.results$tt
    inf.func <- Jres$inffunc

    # return the results for this weighting function
    list(J=att, group=group, t=tt, inf.func=inf.func)
  }, cl=cores)

  # drop post-treatment (g,t); ***TODO: this is an obvious place
  # to make the code faster -- instead of dropping these, never
  # compute them***
  out <- lapply(out, function(Js) {
    # which elements of results to keep
    keepers <- Js$group > Js$t
    this.group <- Js$group[keepers]
    this.t <- Js$t[keepers]
    this.J <- as.matrix(Js$J[keepers])
    this.inf.func <- as.matrix(Js$inf.func[,keepers])
    list(J=this.J, group=this.group, t=this.t, inf.func=this.inf.func)
  })

  # grab an array of the influence function for the test statistic
  Jinf.func <- simplify2array(BMisc::getListElement(out, "inf.func"))

  # get the test statistic for all values of g,t,x
  J.inner <- sapply(out, function(o) o$J)

  # handle case with 1 pre-treatment period differently from multiple periods
  ifelse(class(J.inner)=="matrix", J <- t(J.inner), J <- as.matrix(J.inner))

  # compute CvM test statistic by averaging over X, and summing over g and t
  CvM <- n*sum(apply(J^2, 2, mean))


  #-----------------------------------------------------------------------------
  # use the multiplier bootstrap to simulate limiting distribution of
  # test statistic
  cat("Step 2 of 2: Simulating limiting distribution of test statistic....\n")
  boot.res <- test.mboot(Jinf.func, dp, cores=cores)


  #-----------------------------------------------------------------------------
  ## # some debugging code
  ## # keeping in case helpful later on
  ## ddd <- 1
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
  ## #-----------------------------------------------------------------------------

  # bootstrap results
  CvMb <- boot.res$bres
  # bootstrap critical value
  CvM.crit.val <- boot.res$crit.val
  # bootstrap p-value
  CvMpval <- 1-ecdf(CvMb)(CvM)

  #-----------------------------------------------------------------------------
  # KS Test - not reporting any more
  #
  # KSb <- unlist(lapply(bout, function(b1) b1$KSb))
  # KSocval <- quantile(KSb, probs=(1-alp), type=1)
  # KSpval <- 1-ecdf(KSb)(KS)
  #-----------------------------------------------------------------------------

  # return test results
  out <- MP.TEST(CvM=CvM, CvMb=CvMb, CvMcval=CvM.crit.val, CvMpval=CvMpval, clustervars=clustervars, xformla=xformla)
  out
}



## #' @title expf
## #'
## #' @description exponential weighting function
## #'
## #' @param X matrix of X's from the data
## #' @param u a particular value to multiply times the X's
## #'
## #' @return numeric vector
## #' @examples
## #' data(mpdta)
## #' dta <- subset(mpdta, year==2007)
## #' X <- model.matrix(~lpop, data=dta)
## #' X <- expf(X, X[1,])
## #'
## #' @export
## expf <- function(X, u) {
##     exp(X%*%u)
## }

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
}



#' @title Multiplier Bootstrap for Conditional Moment Test
#'
#' @description A slightly modified multiplier bootstrap procedure
#'  for the pre-test of the conditional parallel trends assumption
#'
#' @inheritParams mboot
#' @param cores The number of cores to use to bootstrap the test
#'  statistic in parallel.  Default is \code{cores=1} which
#'  corresponds to not running parallel.
#'
#' @return list
#' \item{bres}{CvM test statistics for each bootstrap iteration}
#' \item{crit.val}{critical value for CvM test statistic}
#'
#' @export
test.mboot <- function(inf.func, DIDparams, cores=1) {

  # setup needed variables
  data <- DIDparams$data
  idname <- DIDparams$idname
  clustervars <- DIDparams$clustervars
  biters <- DIDparams$biters
  tname <- DIDparams$tname
  tlist <- unique(data[,tname])[order(unique(data[,tname]))]
  alp <- DIDparams$alp
  panel <- DIDparams$panel

  # just get n obsevations (for clustering below...)
  ifelse(panel,
         dta <- data[ data[,tname]==tlist[1], ],
         dta <- data)
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

  crit.val <- quantile(bout, 1-alp, type=1)

  list(bres=bout, crit.val=crit.val)
}
