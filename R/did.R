#' @title mp.spatt
#'
#' @description \code{mp.spatt} computes the ATT in the case where there are more
#'  than two periods of data and allowing for treatment to occur at different points in time
#'  extending the method of Abadie (2005).  This method relies on once individuals are treated
#'  they remain in the treated state for the duration.
#'
#' @param formla The formula y ~ d where y is the outcome and d is the
#'  treatment indicator (d should be binary)
#' @param xformla A optional one sided formula for additional covariates that
#'  will be adjusted for.  E.g ~ age + education.  Additional covariates can
#'  also be passed by name using the x paramater.
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param aggte boolean for whether or not to compute aggregate treatment effect parameters, default TRUE
#' @param w A vector of weights for each observation (not implemented)
#' @param panel Boolean indicating whether the data is panel or repeated cross
#'  sections
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
#'
#' @examples
#' data(mpdta)
#'
#' ## with covariates
#' out1 <- mp.spatt(lemp ~ treat, xformla=~lpop, data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year",
#'                 bstrap=FALSE, se=TRUE, cband=FALSE)
#' ## summarize the group-time average treatment effects
#' summary(out1)
#' ## summarize the aggregated treatment effect parameters
#' summary(out1$aggte)
#'
#' ## without any covariates
#' out2 <- mp.spatt(lemp ~ treat, xformla=NULL, data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year",
#'                 bstrap=FALSE, se=TRUE, cband=FALSE)
#' summary(out2)
#'
#' @references Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
#'
#' @return \code{MP} object
#'
#' @export
mp.spatt <- function(formla, xformla=NULL, data, tname,
                     aggte=TRUE, w=NULL, panel=FALSE,
                     idname=NULL, first.treat.name, alp=0.05,
                     method="logit", se=TRUE,
                     bstrap=FALSE, biters=100, clustervars=NULL,
                     cband=FALSE, citers=100,
                     seedvec=NULL, pl=FALSE, cores=2,
                     printdetails=TRUE) {

    ## make sure that data is a data.frame
    ## this gets around RStudio's default of reading data as tibble
    if (!all( class(data) == "data.frame")) {
        warning("class of data object was not data.frame; converting...")
        data <- as.data.frame(data)
    }

    data$y <- data[,BMisc::lhs.vars(formla)] ##data[,as.character(formula.tools::lhs(formla))]
    ##figure out the dates and make balanced panel
    tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

    flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
    if ( length(flist[flist==0]) == 0) {
        warning("dataset does not have any observations in the control group.  make sure to set data[,first.treat.name] = 0 for observations in the control group.")
    }
    flist <- flist[flist>0]

    ##################################
    ## do some error checking
    if (!is.numeric(tlist)) {
        warning("not guaranteed to order time periods correclty if they are not numeric")
    }

    ####################################


    tlen <- length(tlist)
    flen <- length(flist)
    if (panel) {
        data <- makeBalancedPanel(data, idname, tname)
        dta <- data[ data[,tname]==tlist[1], ]  ## use this for the influence function
        
        ## check that first.treat doesn't change across periods for particular individuals
        if (!all(sapply( split(data, data[,idname]), function(df) {
            length(unique(df[,first.treat.name]))==1
        }))) {
            stop("Error: the value of first.treat must be the same across all periods for each particular individual.")
        }
    } else {

        dta <- data ## this is for repeated cross sections case though
        ## i'm not sure it's working correctly overall
    }

    if (is.null(xformla)) {
        xformla <- ~1
    }


    #################################################################
    ## more error handling after we have balanced the panel
    gsize <- aggregate(data[,first.treat.name], by=list(data[,first.treat.name]), function(x) length(x)/length(tlist))
    reqsize <- length(rhs.vars(xformla)) + 5 ## 5 is just to give a buffer, could increase or decrease
    gsize <- subset(gsize, x < reqsize) ## x is name of column from aggregate
    if (nrow(gsize) > 0) {
        gpaste <-  paste(gsize[,1], collapse=",")
        warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ", gpaste, "\n  and consider dropping these..."))
    }

    #################################################################

    results <- compute.mp.spatt(flen, tlen, flist, tlist, data, dta, first.treat.name,
                                formla, xformla, tname, w, panel, idname, method, seedvec, se,
                                pl, cores, printdetails)



    ## if (!panel) { ## if not panel use empirical bootstrap
    ##     fatt <- results$fatt
    ##     warning("only reporting point estimates for data with repeated cross sections")
    ##     return(fatt)
    ## }

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


    ## new code
    cval <- qnorm(1-alp/2)
    if (cband) {
        bSigma <- apply(bres, 2, function(b) (quantile(b, .75, type=1) - quantile(b, .25, type=1))/(qnorm(.75) - qnorm(.25)))
        bT <- apply(bres, 1, function(b) max( abs(b)*bSigma^(-.5) ))
        cval <- quantile(bT, 1-alp, type=1)
        ##bT1 <- apply(bres, 1, function(b) max( abs(b)*diag(V)^(-.5) ))
        ##cval1 <- quantile(bT1, 1-alp, type=1)
        V <- diag(bSigma) ## this is the appropriate matrix for
         ## constructing confidence bands
    }

    aggeffects <- NULL
    if (aggte) {
        aggeffects <- compute.aggte(flist, tlist, group, t, att, first.treat.name, inffunc1, n, clustervars, dta, idname, bstrap, biters)
    }

    ## wald test for pre-treatment periods
    pre <- which(t < group)
    preatt <- as.matrix(att[pre])
    preV <- V[pre,pre]

    if (length(preV) == 0) {
        message("No pre-treatment periods to test")
        return(MP(group=group, t=t, att=att, V=V, c=cval, inffunc=inffunc1, n=n, aggte=aggeffects))
    }

    if (det(preV) == 0) { ##matrix not invertible
        warning("Not returning pre-test Wald statistic due to singular covariance matrix")
        return(MP(group=group, t=t, att=att, V=V, c=cval, inffunc=inffunc1, n=n, aggte=aggeffects))
    }

    W <- n*t(preatt)%*%solve(preV)%*%preatt
    q <- length(pre)##sum(1-as.numeric(as.character(results$post))) ## number of restrictions
    Wpval <- round(1-pchisq(W,q),5)


    return(MP(group=group, t=t, att=att, V=V, c=cval, inffunc=inffunc1, n=n, W=W, Wpval=Wpval, aggte=aggeffects))
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
                             first.treat.name, formla,
                             xformla, tname, w, panel, idname,
                             method, seedvec, se,
                             pl, cores, printdetails) {

    yname <- BMisc::lhs.vars(formla) ##as.character(formula.tools::lhs(formla))

    fatt <- list()
    counter <- 1
    inffunc <- array(data=0, dim=c(flen,tlen,nrow(dta)))
    for (f in 1:flen) {
            ##satt <- list()
        for (t in 1:(tlen-1)) {
            pret <- t
            if (flist[f]<=tlist[(t+1)]) {
                ## set an index for the pretreatment period
                pret <- tail(which(tlist < flist[f]),1)

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
            if (panel) {
                ## get dataset with current period and pre-treatment period
                disdat <- data[(data[,tname]==tlist[t+1] | data[,tname]==tlist[pret]),]
                ## transform it into "cross-sectional" data where
                ## one of the columns contains the change in the outcome
                ## over time
                disdat <- panel2cs(disdat, yname, idname, tname)

                ## set up control group
                disdat$C <- 1*(disdat[,first.treat.name] == 0)

                ## set up for particular treated group
                disdat$G <- 1*(disdat[,first.treat.name] == flist[f])

                ## drop missing factors
                disdat <- droplevels(disdat)

                ## set up xformla in no covariates case
                if (is.null(xformla)) {
                    xformla <- ~1
                }

                ## set up formula for propensity score, estimate it,
                ## get coefficients and get propensity score
                pformla <- xformla

                pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))

                pscore.reg <- glm(pformla, family=binomial(link="logit"),
                                  data=subset(disdat, C+G==1))
                thet <- coef(pscore.reg)

                ## error handling for too many covariates
                if (any(is.na(thet))) {
                    warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ", flist[f], " at time period: ", tlist[t+1]))
                }

                ## estimate propensity score
                pscore <- predict(pscore.reg, newdata=disdat, type="response")


                ## give short names for data in this iteration
                G <- disdat$G
                C <- disdat$C
                dy <- disdat$dy
                x <- model.matrix(xformla, data=disdat)
                n <- nrow(disdat)

                ## set up weights
                attw1 <- G/mean(G)
                attw2a <- pscore*C/(1-pscore)
                attw2 <- attw2a/mean(attw2a)
                att <- mean((attw1 - attw2)*dy)

                ## save results for this iteration
                fatt[[counter]] <- list(att=att, group=flist[f], year=tlist[(t+1)], post=1*(flist[f]<=tlist[(t+1)]))

                ## --------------------------------------------
                ## get the influence function

                ## weigts
                wg <- G/mean(G)
                wc1 <- C*pscore / (1-pscore)
                wc <- wc1 / mean(wc1)

                ## influence function for treated group
                psig <- wg*(dy - mean(wg*dy))

                ## influence function for control group (see paper)
                M <- as.matrix(apply(as.matrix((C/(1-pscore))^2 * g(x,thet) * (dy - mean(wc*dy)) * x), 2, mean) / mean(wc1))
                A1 <- (G + C)*g(x,thet)^2/(pscore*(1-pscore))
                A1 <- (t(A1*x)%*%x/n)
                A2 <- ((G + C)*(G-pscore)*g(x,thet)/(pscore*(1-pscore)))*x
                A <- A2%*%MASS::ginv(A1)
                psic <- wc*(dy - mean(wc*dy)) + A%*%M

                ## save the influnce function as the difference between
                ## the treated and control influence functions;
                ## we save this as a 3-dimensional array
                ## and then process afterwards
                inffunc[f,t,] <- psig - psic
            } else {
                ## --------------------------------------------
                ## this is repeated cross sections case
                ## unlike the panel data case, here we calculate averages
                ## over the entire data set to estimate some things

                ## set up short variables for the control group and
                ## treated groups, and a dummy variable for an
                ## observation in the current period
                data$C <- 1*(data[,first.treat.name] == 0)
                data$G <- 1*(data[,first.treat.name] == flist[f])
                data$T <- 1*(data[,tname]==tlist[t+1])
                data$preT <- 1*(data[,tname]==tlist[pret])

                ## estimate the propensity score
                ## unlike the panel case, here we can use all time periods
                ## though we drop observations that are not in group G
                ## or the control group
                if (is.null(xformla)) {
                    xformla <- ~1
                }
                pformla <- xformla
                pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))
                pscore.reg <- glm(pformla, family=binomial(link="logit"),
                                  data=subset(data, C+G==1))

                thet <- coef(pscore.reg)
                ## error handling for too many covariates
                if (any(is.na(thet))) {
                    warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ", flist[f], " at time period: ", tlist[t+1]))
                }
                pscore <- predict(pscore.reg, newdata=data, type="response")
                data$pscore <- pscore

                ## set up the weights

                ## first set of weights for group G in period T
                wt1 <- data$G * data$T
                nwt1 <- wt1/mean(wt1)

                ## weights for group G in period G-1
                wt2 <- data$G * data$preT
                nwt2 <- wt2/mean(wt2)

                ## weights for group C in period T
                wc1 <- data$T * pscore * data$C / (1-pscore)
                nwc1 <- wc1/mean(wc1)

                ## weights for group C in period G-1
                wc2 <- data$preT * pscore * data$C / (1-pscore)
                nwc2 <- wc2/mean(wc2)

                ## average treatment effect
                att <- mean( ((nwt1 - nwt2) - (nwc1 - nwc2))*data$y)

                ## save results
                fatt[[counter]] <- list(att=att, group=flist[f], year=tlist[(t+1)], post=1*(flist[f]<=tlist[(t+1)]))

                ## ------------------------------------------
                ## get the influence function

                ## for the treated group
                y <- data$y
                psit1 <- nwt1*(y - mean(nwt1*y))
                psit2 <- nwt2*(y - mean(nwt2*y))

                ## for the untreated group
                x <- model.matrix(xformla, data=data)
                M1 <- as.matrix(apply( (data$T * data$C/(1-pscore))^2 * g(x,thet)*(y - mean(nwc1*y))*x, 2, mean) / mean(wc1))
                M2 <- as.matrix(apply( (data$preT * data$C/(1-pscore))^2 * g(x,thet)*(y - mean(nwc2*y))*x, 2, mean) / mean(wc2))
                A1 <- ((data$G + data$C)*g(x,thet)^2 / (pscore*(1-pscore)))*x
                A1 <- t(A1)%*%x
                A2 <- (((data$G + data$C)*(data$G - pscore)*g(x,thet))/(pscore*(1-pscore)))*x
                xi <- t(solve(A1)%*%t(A2))
                psic1 <- nwc1*(y - mean(nwc1*y)) + xi%*%M1
                psic2 <- nwc2*(y - mean(nwc2*y)) + xi%*%M2

                inffunc[f,t,] <- psit1 - psit2 - (psic1 - psic2)
            }

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
#'
#' @return MP object
#' @export
MP <- function(group, t, att, V, c, inffunc, n=NULL, W=NULL, Wpval=NULL, aggte=NULL) {
    out <- list(group=group, t=t, att=att, V=V, c=c, inffunc=inffunc, n=n, W=W, Wpval=Wpval, aggte=aggte)
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
mp.spatt.test <- function(formla, xformlalist=NULL, data, tname,
                          weightfun=NULL, w=NULL, panel=FALSE,
                          idname=NULL, first.treat.name,
                          alp=0.05, method="logit",
                          biters=100, clustervarlist=NULL,
                          pl=FALSE, cores=2) {


    data$y <- data[,BMisc::lhs.vars(formla)] ##data[,as.character(formula.tools::lhs(formla))]
    ##figure out the dates and make balanced panel
    tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

    flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
    flist <- flist[flist>0]


    ##################################
    ## eventually can put this in its own functions as it is duplicate
    ## code for the estimations
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

    tlen <- length(tlist)
    flen <- length(flist)
    if (panel) {
        data <- makeBalancedPanel(data, idname, tname)
        dta <- data[ data[,tname]==tlist[1], ]  ## use this for the influence function
    } else {
        warning("not guaranteed to work correctly for repeated cross sections")
        dta <- data ## this is for repeated cross sections case though
        ## i'm not sure it's working correctly overall
    }

    n <- nrow(dta)

    if (is.null(weightfun)) {
        weightfun <- expf
    }

    xformlalist <- lapply(xformlalist, function(ff) if (is.null(ff)) ~1 else ff)

    thecount <- 1
    innercount <- 1

    outlist <- lapply(xformlalist, function(xformla) {

        ##
        ##X1 <- apply(model.matrix(xformla, data), 2, function(col) {
        ##    pnorm(col)
        ##})## for the entire dataset
        X <- apply(model.matrix(xformla, dta), 2, function(col) {
            ##pnorm(col)
            ecdf(col)(col)
        })

        ## for debugging: X <- as.matrix(X[1:100,])

        X <- unique(X)

        X1 <- apply(model.matrix(xformla, dta), 2, function(col) {
            ecdf(col)(col)
            ##pnorm(col)
        })

        thetlist <- list()
        pscorelist <- list()
        for (f in 1:flen) {

            disdat <- data[(data[,tname]==tlist[1]),]

            disdat$C <- 1*(disdat[,first.treat.name] == 0)

            disdat$G <- 1*(disdat[,first.treat.name] == flist[f])

            disdat <- droplevels(disdat)

            pformla <- xformla
            pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))##formula.tools::lhs(pformla) <- as.name("G")
            pscore.reg <- glm(pformla, family=binomial(link="logit"),
                              data=subset(disdat, C+G==1))
            thetlist[[f]] <- coef(pscore.reg)
            pscorelist[[f]] <- predict(pscore.reg, newdata=disdat, type="response")
        }



        cat("\n Step", thecount, "of", length(xformlalist), ":.....................\n")
        thecount <<- thecount+1
        out <- pbapply::pblapply(1:nrow(X), function(i) {
            www <- as.numeric(weightfun(X1, X[i,]))##exp(X1%*%X[i,])##plogis(X1%*%X[i,]) ##(1*(apply((X1 <= X[i,]), 1, all)))
            yname <- BMisc::lhs.vars(formla) ##as.character(formula.tools::lhs(formla))

            fatt <- list()
            counter <- 1
            inffunc <- array(data=0, dim=c(sum(sapply(flist, function(g) 1*(g>tlist[-1]))), nrow(dta)))
            Jout <- c()
            groupout <- c()
            tout <- c()
            for (f in 1:flen) {
                for (t in 1:(tlen-1)) {
                    pret <- t
                    ## if (flist[f]<=tlist[(t+1)]) {
                    ##     pret <- tail(which(tlist < flist[f]),1) ## remember, this is just an index
                    ## }

                    if (flist[f] <= tlist[t+1]) break

                    disdat <- data[(data[,tname]==tlist[t+1] | data[,tname]==tlist[pret]),]
                    disdat <- panel2cs(disdat, yname, idname, tname)

                    disdat$C <- 1*(disdat[,first.treat.name] == 0)

                    disdat$G <- 1*(disdat[,first.treat.name] == flist[f])

                    disdat <- droplevels(disdat)

                    ## try to do tthis his outside of the loop
                    ## pformla <- xformla
                    ## formula.tools::lhs(pformla) <- as.name("G")
                    ## pscore.reg <- glm(pformla, family=binomial(link="logit"),
                    ##                   data=subset(disdat, C+G==1))
                    ## thet <- coef(pscore.reg)
                    ## pscore <- predict(pscore.reg, newdata=disdat, type="response")
                    thet <-thetlist[[f]]
                    pscore <- pscorelist[[f]]

                    G <- disdat$G
                    C <- disdat$C
                    dy <- disdat$dy
                    x <- model.matrix(xformla, data=disdat)
                    n <- nrow(disdat)

                    Jwg1 <- G/pscore
                    Jwg <- Jwg1/mean(Jwg1)
                    Jwc1 <- pscore*C/(1-pscore)
                    Jwc <- Jwc1/mean(Jwc1)
                    ##Jwc1 <- C/(1-pscore)
                    ##Jwc <- Jwc1/mean(Jwc1)

                    J <- mean( (Jwg - Jwc)*www*dy )

                    fatt[[counter]] <- list(J=J, group=flist[f], year=tlist[(t+1)], post=1*(flist[f]<=tlist[(t+1)]))

                    ## get the influence function

                    psig <- Jwg*(www*dy - mean(Jwg*www*dy))

                    M <- as.matrix(apply(as.matrix((C/(1-pscore))^2 * g(x,thet) * (www*dy - mean(Jwc*www*dy)) * x), 2, mean) / mean(Jwc1))
                    A1 <- (G + C)*g(x,thet)^2/(pscore*(1-pscore))
                    A1 <- (t(A1*x)%*%x/n)
                    A2 <- ((G + C)*(G-pscore)*g(x,thet)/(pscore*(1-pscore)))*x
                    A <- A2%*%MASS::ginv(A1)
                    psic <- Jwc*(www*dy - mean(Jwc*www*dy)) + A%*%M

                    ## Jg <- mean(Jwg*dy*www)
                    ## Jc <- mean(Jwc*dy*www)

                    ## Mg <- as.matrix(apply(as.matrix((G/(pscore))^2 * g(x,thet) * ( (dy*www - Jg) )* x), 2, mean) / mean( G/pscore  ))

                    ## Mc <- -as.matrix(apply(as.matrix((C/(1-pscore))^2 * g(x,thet) * ( (dy*www - Jc)) * x), 2, mean) / mean( C / (1-pscore) ) )

                    ## A1 <- (G + C)*g(x,thet)^2/(pscore*(1-pscore))
                    ## A1 <- (t(A1*x)%*%x/n)
                    ## A2 <- ((G + C)*(G-pscore)*g(x,thet)/(pscore*(1-pscore)))*x
                    ## A <- A2%*%MASS::ginv(A1)

                    ## psig <- Jwg*dy*www - A%*%Mg
                    ## psic <- Jwc*dy*www - A%*%Mc

                    inffunc[counter,] <- psig - psic
                    Jout[counter] <- J
                    groupout[counter] <- flist[f]
                    tout[counter] <- tlist[t+1]

                    counter <- counter+1
                }
            }

            list(J=Jout, group=groupout, t=tout, inffunc=inffunc)

        })


        outinffunc <- lapply(out, function(o) t(o$inffunc))

        J <- t(sapply(out, function(o) o$J))
        KS <- sqrt(n) * sum(apply(J,2,function(j) max(abs(j))))
        CvM <- n*sum(apply(J, 2, function(j) mean( j^2 )))

        if (!pl) {
            cores <- 1
        }

        innercount <<- 1
        lapply(clustervarlist, function(clustervars) {

            cat("\n >>> Inner Step", innercount, "of", length(clustervarlist), ":.....................\n")
            innercount <<- innercount+1
            bout <- pbapply::pblapply(1:biters, cl=cores, FUN=function(b) {
                Jb <- t(sapply(outinffunc, function(inffunc1) {
                    ## new version
                    if (idname %in% clustervars) {
                        clustervars <- clustervars[-which(clustervars==idname)]
                    }
                    ##clustervars <- clustervars[-which(clustervars==idname)]
                    if (length(clustervars) > 1) {
                        stop("can't handle that many cluster variables")
                    }
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
                    ##Rb <- sqrt(n)*(apply(Ub*(inffunc1), 2, mean))
                    ##Rb
                    Jb1 <- apply(Ub*inffunc1, 2, mean)
                    Jb1
                }))

                ##bres <- t(simplify2array(bout))
                ##V <- cov(bres)
                CvMb <- n*sum(apply(Jb, 2, function(j) mean( j^2 )))
                KSb <- sqrt(n) * sum(apply(Jb,2,function(j) max(abs(j))))
                list(CvMb=CvMb, KSb=KSb)
            })

            CvMb <- unlist(lapply(bout, function(b1) b1$CvMb))
            CvMocval <- quantile(CvMb, probs=(1-alp), type=1)
            CvMpval <- 1-ecdf(CvMb)(CvM)

            KSb <- unlist(lapply(bout, function(b1) b1$KSb))
            KSocval <- quantile(KSb, probs=(1-alp), type=1)
            KSpval <- 1-ecdf(KSb)(KS)


            ## bres <- lapply(bout, simplify2array)
            ## CvMb <- sapply(bres, function(b) apply(b, 1, function(bb) mean(bb^2)))
            ## CvMb <- n*apply(CvMb, 2, sum)
            ## CvMocval <- quantile(CvMb, probs=(1-alp), type=1)
            ## CvMpval <- 1-ecdf(CvMb)(CvM)
            ## bres <- sapply(bres, function(b) apply(b, 1, function(bb) max(abs(bb))))
            ## KSb <- sqrt(n)*apply(bres, 2, sum)
            ## KSocval <- quantile(KSb, probs=(1-alp), type=1)
            ## KSpval=1-ecdf(KSb)(KS)
            out <- MP.TEST(CvM=CvM, CvMb=CvMb, CvMcval=CvMocval, CvMpval=CvMpval, KS=KS, KSb=KSb, KScval=KSocval, KSpval=KSpval, clustervars=clustervars, xformla=xformla)
            out
        })
    })

    return(unlist(outlist, recursive=FALSE))
}


#' @title MP.TEST
#'
#' @description MP.TEST objects
#'
#' @param CvM Cramer von Mises test statistic
#' @param CvMb a vector of boostrapped Cramer von Mises test statistics
#' @param CvMcval CvM critical value
#' @param CvMpval p-value for CvM test
#' @param KS Kolmogorov-Smirnov test statistic
#' @param KSb a vector of boostrapped KS test statistics
#' @param KScval KS critical value
#' @param KSpval p-value for KS test
#' @param clustervars vector of which variables were clustered on for the test
#' @param xformla formla for the X variables used in the test
#'
#' @export
MP.TEST <- function(CvM, CvMb, CvMcval, CvMpval, KS, KSb, KScval, KSpval, clustervars, xformla) {
    out <- list(CvM=CvM, CvMb=CvMb, CvMcval=CvMcval, CvMpval=CvMpval, KS=KS, KSb=KSb, KScval=KScval, KSpval=KSpval, clustervars=clustervars, xformla=xformla)
    class(out) <- "MP.TEST"
    out
}

#' @title summary.MP.TEST
#'
#' @description print a summary of test results
#'
#' @param object an MP.TEST object
#' @param ... other variables
#'
#' @export
summary.MP.TEST <- function(object, ... ) {
    CvM <- object$CvM
    CvMcval <- object$CvMcval
    CvMpval <- object$CvMpval
    citation()
    cat("Cramer von Mises: \n")
    cat("  Test Statistic: ", CvM, "\n")
    cat("  Critical Value: ", CvMcval, "\n")
    cat("  P-value       : ", CvMpval, "\n \n")
    cat("Clustering on   : ", paste0(object$clustervars,sep=","), "\n")
    cat("X formula       : ", as.character(object$xformla), "\n")
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
    p <- ggplot(ssresults,
                aes(x=year, y=att, ymin=(att-c*ate.se),
                    ymax=att+c*ate.se, post=post)) +
        geom_point(aes(colour=post), size=1.5) +
        geom_errorbar(aes(colour=post), width=0.1) +
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
#' @examples
#' \dontrun{
#' data(mpdta)
#' out <- mp.spatt(lemp ~ treat, xformla=~lpop, data=mpdta,
#'                 panel=TRUE, first.treat.name="first.treat",
#'                 idname="countyreal", tname="year",
#'                 bstrap=FALSE, se=TRUE, cband=FALSE)
#' ggdid(out)
#' }
#'
#' @export
ggdid <- function(mpobj, type=c("attgt", "dynamic", "selective", "calendar", "dynsel"), ylim=NULL, xlab=NULL, ylab=NULL, title="Group", xgap=1, ncol=1, e1=1) {

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

        ##results <- mp2ATT(results, vcovatt)

        mplots <- lapply(g, function(g) {
            thisdta <- subset(results, group==g)
            gplot(thisdta, ylim, xlab, ylab, title, xgap)
        })

        do.call("grid.arrange", c(mplots))
    } else if (type=="dynamic") {
        aggte <- mpobj$aggte
        if (mpobj$c > 2) warning("uniform bands not implemented yet for this plot...")
        elen <- length(aggte$dynamic.att.e)
        results <- cbind.data.frame(year=as.factor(seq(1:elen)), att=aggte$dynamic.att.e, ate.se=aggte$dynamic.se.e, post=as.factor(1), c=qnorm(.975))
        p <- gplot(results, ylim, xlab, ylab, title, xgap)
        p
    } else if (type=="selective") {
        aggte <- mpobj$aggte
        if (mpobj$c > 2) warning("uniform bands not implemented yet for this plot...")
        results <- cbind.data.frame(year=as.factor(aggte$groups), att=aggte$selective.att.g, ate.se=aggte$selective.se.g, post=as.factor(1), c=qnorm(.975))
        p <- gplot(results, ylim, xlab, ylab, title, xgap)
        p
    } else if (type=="calendar") {
        aggte <- mpobj$aggte
        if (mpobj$c > 2) warning("uniform bands not implemented yet for this plot...")
        results <- cbind.data.frame(year=as.factor(aggte$times[aggte$times >= min(aggte$groups)]), att=aggte$calendar.att.t, ate.se=aggte$calendar.se.t, post=as.factor(1), c=qnorm(.975))
        p <- gplot(results, ylim, xlab, ylab, title, xgap)
        p
    } else if (type=="dynsel") {
        if (is.null(e1)) {
            stop("must provide value of e1 for reporting dynamic effects with selective treatment timing")
        }
        aggte <- mpobj$aggte
        if (mpobj$c > 2) warning("uniform bands not implemented yet for this plot...")
        results <- cbind.data.frame(year=as.factor(1:e1), att=aggte$dynsel.att.ee1[[e1]]$dte[1:e1], ate.se=aggte$dynsel.se.ee1[[e1]]$se[1:e1], post=as.factor(1), c=qnorm(.975))
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
compute.aggte <- function(flist, tlist, group, t, att, first.treat.name, inffunc1, n, clustervars, dta, idname, bstrap, biters) {

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
            return(sqrt( mean( (thisinffunc)^2 ) ) / sqrt(n))
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

    ## some variables used throughout
    pg <- sapply(originalflist, function(g) sum(1*(dta[,first.treat.name]==g & dta[,first.treat.name]>0))/ sum(1*(dta[,first.treat.name]>0)))
    pgg <- pg
    pg <- pg[match(group, flist)] ## make it have the same length as att
    attg <- split(att, group)
    tg <- split(t, group)
    keepers <- which(group <= t)
    G <-  unlist(lapply(dta[,first.treat.name], orig2t))
    p0 <- mean(1*G==0)
    pT <- 1-p0
    T <- 1*(dta[first.treat.name]>0)


    ## simple att
    simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))
    simple.oif1 <- sapply(keepers, function(k) (1*(G==group[k]) - mean(1*(G==group[k]))) / sum(pg[keepers]))
    simple.oif2 <- sapply(keepers, function(j) mean(1*(G==group[j])) * apply(sapply(keepers, function(k) (1*(G==group[k]) - mean(1*(G==group[k])))),1,sum))
    simple.se <- getSE(keepers, pg[keepers]/sum(pg[keepers]), simple.oif1-simple.oif2)

    ## Selective Treatment Timing
    ## Note: for selective.att.g, don't need to adjust standard
    ##  errors for estimating weights because they are known
    selective.att.g <- sapply(flist, function(g) {
        whichg <- which( (group == g) & (g <= t))
        attg <- att[whichg]
        mean(attg)
    })
    selective.se.g <- sapply(flist, function(g) {
        whichg <- which( (group == g) & (g <= t))
        getSE(whichg, pg[whichg]/sum(pg[whichg]))
    })
    ## Note: for estimating selective.att do need to adjust
    ##  standard errors because need to estimate P(G=g)
    selective.att <- sum(selective.att.g * pgg)
    keepers <- which(group <= t)
    selective.weights <- pg[keepers]/(max(t) - group[keepers] + 1)  ## note could just use these directly by mulitiplying att[keepers] and taking sum
    selective.oif1 <- sapply(keepers, function(k) 1*(G==group[k])/p0)##T/p0
    selective.oif2 <- sapply(keepers, function(k) (1-T)*pg[k]/p0^2 )##(1-T)/p0^2##apply(sapply(pgg, function(p) p*(1-T)),1,sum)
    selective.se <- getSE(keepers, selective.weights, selective.oif1 - selective.oif2)


    ## Dynamic Treatment Effects
    eseq <- seq(1,max(t-group)+1)
    dynamic.att.e <- sapply(eseq, function(e) {
        whiche <- which(t - group + 1 == e)
        atte <- att[whiche]
        pge <- pg[whiche]/(sum(pg[whiche]))
        sum(atte*pge)
    })
    dynamic.se.e <- sapply(eseq, function(e) {
        whiche <- which(t - group + 1 == e)
        pge <- pg[whiche]/(sum(pg[whiche]))
        dynamic.oif1 <- sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k]))) / sum(pg[whiche]))
        dynamic.oif2 <- sapply(whiche, function(j) mean(1*(G==group[j])) * apply(sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k])))),1,sum))
        getSE(whiche, pge, dynamic.oif1 - dynamic.oif2)
    })
    dynamic.att <- mean(dynamic.att.e)
    keepers <- which(group <= t)
    dynamic.weights <- lapply(eseq, function(e) {
        whiche <- which(t - group + 1 == e)
        pge <- pg[whiche]/(sum(pg[whiche]))
        dynamic.oif1 <- sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k]))) / sum(pg[whiche]))
        dynamic.oif2 <- sapply(whiche, function(j) mean(1*(G==group[j])) * apply(sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k])))),1,sum))
        list(whiche=whiche, pge=pge, oif=dynamic.oif1-dynamic.oif2)
    })
    which.dynamic.weights <- unlist(lapply(dynamic.weights, function(d) d$whiche))
    dynamic.oif <- do.call(cbind,lapply(dynamic.weights, function(d) d$oif))
    dynamic.weights <- unlist(lapply(dynamic.weights, function(d) d$pge))[order(which.dynamic.weights)] / length(unique(eseq))
    dynamic.se <- getSE(keepers, dynamic.weights, wif=dynamic.oif)


    ## Calendar Time Effects
    tseq <- unique(t)[unique(t) >= min(group)]
    calendar.att.t <- sapply(tseq, function(t1) {
        whicht <- which((group <= t1) & (t==t1))
        attt <- att[whicht]
        pgt <- pg[whicht]/(sum(pg[whicht]))
        sum(attt*pgt)
    })
    calendar.se.t <- sapply(tseq, function(t1) {
        whicht <- which((group <= t1) & (t==t1))
        pgt <- pg[whicht]/(sum(pg[whicht]))
        calendar.oif1 <- sapply(whicht, function(k) (1*(G==group[k]) - mean(1*(G==group[k]))) / sum(pg[whicht]))
        calendar.oif2 <- sapply(whicht, function(j) mean(1*(G==group[j])) * apply(sapply(whicht, function(k) (1*(G==group[k]) - mean(1*(G==group[k])))),1,sum))
        getSE(whicht, pgt, calendar.oif1 - calendar.oif2)
    })
    calendar.att <- mean(calendar.att.t)
    keepers <- which(group <= t)
    calendar.weights <- lapply(tseq, function(t1) {
        whicht <- which((group <= t1) & (t==t1))
        pgt <- pg[whicht]/(sum(pg[whicht]))
        calendar.oif1 <- sapply(whicht, function(k) (1*(G==group[k]) - mean(1*(G==group[k]))) / sum(pg[whicht]))
        calendar.oif2 <- sapply(whicht, function(j) mean(1*(G==group[j])) * apply(sapply(whicht, function(k) (1*(G==group[k]) - mean(1*(G==group[k])))),1,sum))
        list(whicht=whicht, pgt=pgt, oif=(calendar.oif1-calendar.oif2))
    })
    which.calendar.weights <- unlist(lapply(calendar.weights, function(t1) t1$whicht))
    calendar.oif <- do.call(cbind,lapply(calendar.weights, function(t1) t1$oif))
    calendar.weights <- unlist(lapply(calendar.weights, function(t1) t1$pgt))[order(which.calendar.weights)] / length(unique(tseq))
    calendar.se <- getSE(keepers, calendar.weights, calendar.oif)


    ## Selective Treatment Timing and Dynamic Treatment Effects
    eseq <- seq(1,max(t-group)+1)
    e1seq <- eseq
    dynsel.att.ee1 <- lapply(e1seq, function(e1) {
        list(dte=sapply(eseq, function(e) {
            whiche <- which( (t - group + 1 == e) &
                             ( max(t) - group + 1 >= e1) &
                             ( e <= e1 ))
            atte <- att[whiche]
            pge <- pg[whiche]/sum(pg[whiche])
            sum(atte*pge)
        }), e1=e1)
    })

    dynsel.se.ee1 <- lapply(e1seq, function(e1) {
        list(se=sapply(eseq[eseq <= e1], function(e) {
            whiche <- which( (t - group + 1 == e) &
                             ( max(t) - group + 1 >= e1) &
                             ( e <= e1 ))
            pge <- pg[whiche]/sum(pg[whiche])
            dynsel.oif1 <- sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k]))) / sum(pg[whiche]))
            dynsel.oif2 <- sapply(whiche, function(j) mean(1*(G==group[j])) * apply(sapply(whiche, function(k) (1*(G==group[k]) - mean(1*(G==group[k])))),1,sum))
            getSE(whiche, pge, dynsel.oif1-dynsel.oif2)
        }), e1=e1)
    })
    dynsel.att.e1 <- sapply(dynsel.att.ee1, function(d) sum(d$dte)/d$e1)
    keepers <- lapply(e1seq, function(e1) { which( (group <= t) & (max(t) - group + 1 >= e1) ) })
    dynsel.weights <- lapply(e1seq, function(e1) {
        list(whiche=unlist(lapply(eseq, function(e) {
            whiche <- which( (t - group + 1 == e) &
                             ( max(t) - group + 1 >= e1) &
                             ( e <= e1 ))
            whiche
        })), pge=unlist(lapply(eseq, function(e) {
            whiche <- which( (t - group + 1 == e) &
                             ( max(t) - group + 1 >= e1) &
                             ( e <= e1 ))
            pg[whiche]/sum(pg[whiche])
        })), e1=e1)
    })

    ## TODO:  accounting for estimation of weights (these
    ##  don't seem to have any effect
    ## kind of hack, but just post-process the weights
    dynsel.weights <- lapply(dynsel.weights, function(d) {
        d$pge <- d$pge/sum(d$pge)
        d
    })
    dynsel.se.e1 <- sapply(dynsel.weights, function(d) {
        getSE(d$whiche, d$pge)
    })

    AGGTE(simple.att=simple.att, simple.se=simple.se, selective.att=selective.att, selective.se=selective.se, selective.att.g=selective.att.g, selective.se.g=selective.se.g, dynamic.att=dynamic.att, dynamic.se=dynamic.se, dynamic.att.e=dynamic.att.e, dynamic.se.e=dynamic.se.e, calendar.att=calendar.att, calendar.se=calendar.se, calendar.att.t=calendar.att.t, calendar.se.t=calendar.se.t, dynsel.att.e1=dynsel.att.e1, dynsel.se.e1=dynsel.se.e1, dynsel.att.ee1=dynsel.att.ee1, dynsel.se.ee1=dynsel.se.ee1,groups=originalflist,times=originaltlist)
}


#' @title AGGTE
#'
#' @description \code{AGGTE} class for aggregate treatment effects
#'
#' @param simple.att simple weighted average of group-time average treatment
#'  effects
#' @param simple.se the standard error for \code{simple.att}
#' @param selective.att aggregated group-time average treament effects when
#'  there is selective treatment timing
#' @param selective.se the standard error for \code{selective.att}
#' @param selective.att.g aggregated group-time average treatment effects
#'  when there is selective treatment timing for each particular group
#' @param selective.se.g the standard error for \code{selective.att.g}
#' @param dynamic.att aggregated group-time average treatment effects when
#'  there are dynamic treatment effects
#' @param dynamic.se the standard error for \code{dynamic.att}
#' @param dynamic.att.e aggregated group-time average treatment effects
#'  when there are dynamic treatment effects for each length of exposure
#'  to treatment
#' @param dynamic.se.e the standard error for \code{dynamic.att.e}
#' @param calendar.att the aggregated group-time average treatment effects
#'  when there are calendar time effects
#' @param calendar.se the standard error for \code{calendar.att}
#' @param calendar.att.t the aggregated group-time average treatment effects
#'  when there are calendar time effects for each time period
#' @param calendar.se.t the standard error for \code{calendar.att.t}
#' @param dynsel.att.e1 aggregated group-time average treatment effects when
#'  there are dynamic treatment effects and selective treatment timing.
#'  Here, e1 is the number of periods that a group is required to be treated
#'  in order to be included in the results.
#' @param dynsel.se.e1 the standard error for \code{dynsel.att.e1}
#' @param dynsel.att.ee1 aggregated group-time average treatment effects when
#'  there are dynamic treatment effects and selective treatment timing.
#'  Here, e1 is the number of periods that a group is required to be treated
#'  in order to be included in the results and for each length of exposure
#'  to treatment
#' @param dynsel.se.ee1 the standard error for \code{dynsel.att.ee1}
#' @param groups vector of all groups
#' @param times vector of all times
AGGTE <- function(simple.att=NULL, simple.se=NULL, selective.att=NULL, selective.se=NULL, selective.att.g=NULL, selective.se.g=NULL, dynamic.att=NULL, dynamic.se=NULL, dynamic.att.e=NULL, dynamic.se.e=NULL, calendar.att=NULL, calendar.se=NULL, calendar.att.t=NULL, calendar.se.t=NULL, dynsel.att.e1=NULL, dynsel.se.e1=NULL, dynsel.att.ee1=NULL, dynsel.se.ee1=NULL, groups=NULL, times=NULL) {
    out <- list(simple.att=simple.att, simple.se=simple.se, selective.att=selective.att, selective.se=selective.se, selective.att.g=selective.att.g, selective.se.g=selective.se.g, dynamic.att=dynamic.att, dynamic.se=dynamic.se, dynamic.att.e=dynamic.att.e, dynamic.se.e=dynamic.se.e, calendar.att=calendar.att, calendar.se=calendar.se, calendar.att.t=calendar.att.t, calendar.se.t=calendar.se.t, dynsel.att.e1=dynsel.att.e1, dynsel.se.e1=dynsel.se.e1, dynsel.att.ee1=dynsel.att.ee1, dynsel.se.ee1=dynsel.se.ee1,groups=groups,times=times)
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
summary.AGGTE <- function(object, type=c("dynamic","selective","calendar","dynsel"), e1=1, ...) {
    citation()
    sep <- "          "
    cat("Overall Summary Measures", "\n")
    cat("------------------------", "\n")
    cat("Simple ATT    : ", object$simple.att, "\n")
    cat("  SE          : ", object$simple.se, "\n")
    cat("Selective ATT : ", object$selective.att, "\n")
    cat("  SE          : ", object$simple.se, "\n")
    cat("Dynamic ATT   : ", object$dynamic.att, "\n")
    cat("  SE          : ", object$dynamic.se, "\n")
    cat("Calendar ATT  : ", object$calendar.att, "\n")
    cat("  SE          : ", object$calendar.se, "\n\n")
    type <- type[1]
    if (type == "dynamic") {
        cat("Dynamic Treatment Effects", "\n")
        cat("-------------------------")
        elen <- length(object$dynamic.att.e)
        printmat <- cbind(seq(1:elen), object$dynamic.att.e, object$dynamic.se.e)
        colnames(printmat) <- c("e","att","se")
        print(kable(printmat))
    }

    ## issue is that groups and times are a bit off..., they are getting set though
    if (type=="selective") {
        cat("Selective Treatment Timing", "\n")
        cat("--------------------------", "\n")
        printmat <- cbind(object$groups, object$selective.att.g, object$selective.se.g)
        colnames(printmat) <- c("g","att","se")
        print(kable(printmat))
    }

    if (type=="calendar") {
        cat("Calendar Time Effects", "\n")
        cat("---------------------")
        printmat <- cbind(object$times[object$times >= min(object$groups)], object$calendar.att.t, object$calendar.se.t)
        colnames(printmat) <- c("t","att","se")
        print(kable(printmat))
    }

    if (type=="dynsel") {
        if (is.null(e1)) {
            stop("must provide value of e1 for reporting dynamic effects with selective treatment timing")
        }
        cat("Selective Treatment Timing and Dynamic Effects with e=", e1, "\n")
        cat("-------------------------------------------------------")
        printmat <- cbind(seq(1:e1), object$dynsel.att.ee1[[e1]]$dte[1:e1],
                          object$dynsel.se.ee1[[e1]]$se[1:e1])
        colnames(printmat) <- c("e","att","se")
        print(kable(printmat))
    }
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
    apply(X <= u, 1, function(b) 1*all(b))
}

#' @title onefun
#'
#' @description just return the value 1
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
#' X <- onefun(X, X[1,])
#'
#' @export
onefun <- function(X, u) {
    rep(1, nrow(X))
}

#' @title g
#'
#' @description Logit pdf
#'
#' @param x nxk data matrix
#' @param thet kx1 vector of parameters
#'
#' @return nx1 vector
#'
#' @keywords internal
g <- function(x,thet) {
    x <- as.matrix(x)
    thet <- as.matrix(thet)
    gval <- 1/((1+exp(x%*%thet))^2)
    as.numeric(gval)
}

## ## x nxk matrix
## ## thet kx1 vector
## ## return nx1 vector
## #' @title G
## #'
## #' @description Logit cdf
## #'
## #' @param x nxk data matrix
## #' @param thet kx1 vector of parameters
## #'
## #' @return nx1 vector
## #' @keywords internal
G <- function(x,thet) {
    x <- as.matrix(x)
    thet <- as.matrix(thet)
    Gval <- exp(x%*%thet)/(1+exp(x%*%thet))
    as.numeric(Gval)
}
