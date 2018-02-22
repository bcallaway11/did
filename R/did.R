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
#'  period when a particular observation is treated
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
#' @references Callaway and Sant'Anna (2018)
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

    data$y <- data[,as.character(formula.tools::lhs(formla))]
    ##figure out the dates and make balanced panel
    tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

    flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
    flist <- flist[flist>0]

    if (!is.numeric(tlist)) {
        warning("not guaranteed to order time periods correclty if they are not numeric")
    }
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


    results <- compute.mp.spatt(flen, tlen, flist, tlist, data, dta, first.treat.name,
                                formla, xformla, tname, w, panel, idname, method, seedvec, se,
                                pl, cores, printdetails)
    fatt <- results$fatt
    inffunc <- results$inffunc

    ## process the results from computing the spatt
    group <- c()
    t    <- c()
    att <- c()
    i <- 1
    
    inffunc1 <- matrix(0, ncol=flen*(tlen-1), nrow=nrow(dta)) ## note, this might not work in unbalanced case 
    for (f in 1:length(fatt)) {
        for (s in 1:(length(fatt[[f]])-1)) {
            group[i] <- fatt[[f]]$group
            t[i] <- fatt[[f]][[s]]$year
            att[i] <- fatt[[f]][[s]]$ate
            inffunc1[,i] <- inffunc[f,s,]
            i <- i + 1
        }
    }


    n <- nrow(dta)
    V <- t(inffunc1)%*%inffunc1/n

    if ( (length(clustervars) > 0) & !bstrap) {
        warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
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
                Vb <- matrix(sample(c(-1,1), n1*ncol(inffunc1), replace=T),
                             nrow=n1)
                Vb <- cbind.data.frame(unique(dta[,clustervars]), Vb)
                Ub <- data.frame(dta[,clustervars])
                Ub <- Vb[match(Ub[,1], Vb[,1]),]
                Ub <- Ub[,-1]
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
                Ub <- matrix(sample(c(-1,1), n*ncol(inffunc1), replace=T),
                             nrow=nrow(inffunc1))
            }
            ## allow for serial correlation
            if (sercor) {
                Ub[,-1] <- Ub[,1]
            }
            ## drop cluster for serial correlation
            
            ##ift <- do.call(magic::adiag, psiitout)
            ##ifunc <- rbind(ift, psiiu)
            
            ##Ub <- sample(c(-1,1), n, replace=T)
            mb <- Ub*(inffunc1)
            apply(mb,2,sum)/sqrt(n)
        })
        bres <- t(simplify2array(bout))
        V <- cov(bres)
    }


    ## TODO: handle case with repeated cross sections; this part is conceptually easier because many
    ##  off-diagonal (though not all) will be 0.

    ## get the actual estimates


    ## compute critical value for uniform confidence band
    c <- qnorm(1-alp/2)
    if (cband) {
        Sig <- matrix(0, nrow=nrow(V), ncol=ncol(V))
        diag(Sig) <- diag(V)
        Siggy <- expm::sqrtm(solve(Sig))
        ZB <- MASS::mvrnorm(citers, rep(0,nrow(V)), V)
        KSB <- apply(ZB %*% Siggy, 1, max)
        c <- quantile(KSB, 1-alp, type=7)
    }

    
    ## wald test for pre-treatment periods
    pre <- which(t < group)
    preatt <- as.matrix(att[pre])
    preV <- V[pre,pre]

    W <- NULL ##n*t(preatt)%*%solve(preV)%*%preatt

    q <- length(pre)##sum(1-as.numeric(as.character(results$post))) ## number of restrictions
    Wpval <- NULL ##round(1-pchisq(W,q),5)

    aggeffects <- NULL
    if (aggte) {
        aggeffects <- compute.aggte(flist, group, t, att, first.treat.name, inffunc1, n, clustervars, dta, idname, bstrap, biters)
    }

    
    return(MP(group=group, t=t, att=att, V=V, c=c, inffunc=inffunc1, n=n, W=W, Wpval=Wpval, aggte=aggeffects))
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


    fatt <- list()
    inffunc <- array(data=0, dim=c(flen,tlen,nrow(dta)))
    for (f in 1:flen) {
        satt <- list()
        for (t in 1:(tlen-1)) {
            pret <- t
            if (flist[f]<=tlist[(t+1)]) {
                pret <- tail(which(tlist < flist[f]),1) ## remember, this is just an index
                if (printdetails) {
                    cat(paste("current period:", tlist[t+1]), "\n")
                    cat(paste("current group:", flist[f]), "\n")
                    cat(paste("set pretreatment period to be", tlist[pret]), "\n")
                }
            }
            S <- (dta[,first.treat.name]==0 |
                    dta[,first.treat.name]==flist[f])
            r <- mean(S)
            disdat <- data[(data[,tname]==tlist[t+1] | data[,tname]==tlist[pret]) &
                           (data[,first.treat.name]==0 | data[,first.treat.name]==flist[f]),]
            disdat <- droplevels(disdat)
            satt[[t]] <- c(spatt(formla, xformla, t=tlist[t+1], tmin1=tlist[pret],
                      tname=tname, data=disdat, w=w, panel=panel,
                      idname=idname, 
                      iters=NULL, alp=NULL, method=method, plot=NULL, se=se,
                      retEachIter=NULL, seedvec=seedvec, pl=pl, cores=cores),
                      group=flist[f], year=tlist[(t+1)], post=1*(flist[f]<=tlist[(t+1)]))
            inffunc[f,t,S] <- satt[[t]]$inffunc
            ## ifu <- satt[[t]]$inffuncu
            ## t(ifu)%*%as.matrix(ifu)/n
        }
        fatt[[f]] <- c(satt, group=flist[f])
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
#' @param mpobj an \code{MP} object
#'
#' @export
summary.MP <- function(mpobj) {
    out <- cbind(mpobj$group, mpobj$t, mpobj$att, sqrt(diag(mpobj$V)/mpobj$n))
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
#'
#' @return list containing test results
#' @export
mp.spatt.test <- function(formla, xformlalist=NULL, data, tname, w=NULL, panel=FALSE,
                     idname=NULL, first.treat.name,
                     alp=0.05, method="logit", plot=FALSE, se=TRUE,
                     bstrap=FALSE, biters=100, clustervarlist=NULL,
                     cband=FALSE, citers=100,
                     retEachIter=FALSE, seedvec=NULL, pl=FALSE, cores=2) {


    data$y <- data[,as.character(formula.tools::lhs(formla))]
    ##figure out the dates and make balanced panel
    tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

    flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
    flist <- flist[flist>0]

    if (!is.numeric(tlist)) {
        warning("not guaranteed to order time periods correclty if they are not numeric")
    }
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

    thecount <- 1

    outlist <- lapply(xformlalist, function(xformla) {

        ## 
        X1 <- apply(model.matrix(xformla, data), 2, function(col) {
            pnorm(col)
        })## for the entire dataset
        X <- apply(model.matrix(xformla, dta), 2, function(col) {
            pnorm(col)
        })

        cat("\n Step", thecount, "of", length(xformlalist), ":.....................\n")
        thecount <<- thecount+1
        out <- pbapply::pblapply(1:nrow(X), function(i) {
            www <- exp(X1%*%X[i,])##plogis(X1%*%X[i,]) ##(1*(apply((X1 <= X[i,]), 1, all)))
            data$lhs <- data$y * www
            formula.tools::lhs(formla) <- as.name("lhs")
            out1 <- mp.spatt(formla=formla, xformla=xformla, data=data, tname=tname, panel=panel, idname=idname, first.treat.name=first.treat.name, iters=iters, alp=alp, se=TRUE, bstrap=FALSE, printdetails=FALSE, aggte=FALSE, cband=FALSE)
            out1
        }, cl=8)


        keepers <- which(out[[1]]$group > out[[1]]$t)
        out <- lapply(out, function(o) {
            list(group=o$group[keepers],
                 t=o$t[keepers],
                 att=o$att[keepers],
                 V=o$V[keepers,keepers],
                 inffunc=o$inffunc[,keepers])
        })

        outinffunc <- lapply(out, function(o) o$inffunc)

        J <- t(sapply(out, function(o) o$att))
        KS <- sqrt(n) * sum(apply(J,2,function(j) max(abs(j))))

        othercount <<- 1
        lapply(clustervarlist, function(clustervars) {

            cat("\n >>> Inner Step", othercount, "of", length(clustervarlist), ":.....................\n")
            othercount <<- othercount+1
            bout <- pbapply::pblapply(1:biters, cl=8, FUN=function(b) {
                lapply(outinffunc, function(inffunc1) {
                    sercor <- idname %in% clustervars ## boolean for whether or not to account for serial correlation
                    clustervars <- clustervars[-which(clustervars==idname)]
                    if (length(clustervars) > 1) {
                        stop("can't handle that many cluster variables")
                    }
                    if (length(clustervars) > 0) {
                        n1 <- length(unique(dta[,clustervars]))
                        Vb <- matrix(sample(c(-1,1), n1*ncol(inffunc1), replace=T),
                                     nrow=n1)
                        Vb <- cbind.data.frame(unique(dta[,clustervars]), Vb)
                        ##colnames(Vb)[1] <- "clvar"
                        Ub <- data.frame(dta[,clustervars])
                        ##colnames(Ub)[1] <- "clvar"
                        Ub <- Vb[match(Ub[,1], Vb[,1]),]
                        ##Ub <- merge(Ub, Vb, by="clvar")
                        Ub <- Ub[,-1]
                    } else {
                        Ub <- matrix(sample(c(-1,1), n*ncol(inffunc1), replace=T),
                                     nrow=nrow(inffunc1))
                    }
                    ## allow for serial correlation
                    if (sercor) {
                        Ub[,-1] <- Ub[,1]
                    }
                    ## drop cluster for serial correlation
                    
                    ##ift <- do.call(magic::adiag, psiitout)
                    ##ifunc <- rbind(ift, psiiu)
                    
                    ##Ub <- sample(c(-1,1), n, replace=T)
                    mb <- Ub*(inffunc1)
                    apply(mb,2,mean)
                })
            })
            
            bres <- lapply(bout, simplify2array)
            bres <- sapply(bres, function(b) apply(b, 1, function(bb) max(abs(bb))))
            KSb <- sqrt(n)*apply(bres, 2, sum)
            ocval <- quantile(KSb, probs=(1-alp), type=1)
            pval=1-ecdf(KSb)(KS)
            out <- list(KS=KS, KSb=KSb, cval=ocval, pval=pval, clustervar=clustervars, xformla=xformla)
            out
        })
    })
    
    return(unlist(outlist, recursive=FALSE))
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
#'
#' @export
ggdid <- function(mpobj, ylim=NULL, xlab=NULL, ylab=NULL, title="Group", xgap=1, ncol=1) {
    G <- length(unique(mpobj$group))
    Y <- length(unique(mpobj$t))## drop 1 period bc DID
    g <- unique(mpobj$group)[order(unique(mpobj$group))] ## -1 to drop control group
    y <- unique(mpobj$t)
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
compute.aggte <- function(flist, group, t, att, first.treat.name, inffunc1, n, clustervars, dta, idname, bstrap, biters) {

    if ( (length(clustervars) > 0) & !bstrap) {
        warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
    }
    
    getSE <- function(whichones, weights) {
        weights <- as.matrix(weights) ## just in case pass vector
        thisinffunc <- inffunc1[,whichones]%*%weights  ##multiplies influence function times weights and sums to get vector of weighted IF (of length n)
        
        if (bstrap) {
            bout <- lapply(1:biters, FUN=function(b) {
                sercor <- idname %in% clustervars ## boolean for whether or not to account for serial correlation
                clustervars <- clustervars[-which(clustervars==idname)]
                if (length(clustervars) > 1) {
                    stop("can't handle that many cluster variables")
                }
                if (length(clustervars) > 0) {
                    n1 <- length(unique(dta[,clustervars]))
                    Vb <- matrix(sample(c(-1,1), n1, replace=T),
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
                    Ub <- matrix(sample(c(-1,1), nrow(thisinffunc), replace=T), ncol=1)
                                 
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
    uniquet <- seq(1,length(unique(t)))
    t2orig <- function(t) {
        unique(originalt)[which(uniquet==t)]
    }
    orig2t <- function(orig) {
        uniquet[which(unique(originalt)==orig)]
    }
    t <- sapply(originalt, orig2t)
    group <- sapply(originalgroup, orig2t)
    flist <- sapply(originalflist, orig2t)

##### aggregated treatment effect parameters
    pg <- sapply(originalflist, function(g) sum(1*(dta[,first.treat.name]==g & dta[,first.treat.name]>0))/ sum(1*(dta[,first.treat.name]>0)))
    pgg <- pg
    pg <- pg[match(group, flist)] ## make it have the same length as att
    attg <- split(att, group)
    tg <- split(t, group)
    keepers <- which(group <= t)
    simple.att <- sum(att[keepers]*pg[keepers])/(sum(pg[keepers]))
    simple.se <- getSE(keepers, pg[keepers]/sum(pg[keepers]))

    selective.att.g <- sapply(flist, function(g) {
        whichg <- which( (group == g) & (g <= t))
        attg <- att[whichg]
        mean(attg)
    })
    selective.se.g <- sapply(flist, function(g) {
        whichg <- which( (group == g) & (g <= t))
        getSE(whichg, pg[whichg]/sum(pg[whichg]))
    })
    selective.att <- sum(selective.att.g * pgg)
    keepers <- which(group <= t)
    selective.weights <- pg[keepers]/(max(t) - group[keepers] + 1)  ## note could just use these directly by mulitiplying att[keepers] and taking sum
    selective.se <- getSE(keepers, selective.weights)

    eseq <- seq(1,max(group-t)-1)
    ##eseq <- seq(1, length(unique(t[(group-t) <= 0])))
    dynamic.att.e <- sapply(eseq, function(e) {
        whiche <- which(t - group + 1 == e)
        atte <- att[whiche]
        pge <- pg[whiche]/(sum(pg[whiche]))
        sum(atte*pge)
    })
    dynamic.se.e <- sapply(eseq, function(e) {
        whiche <- which(t - group + 1 == e)
        pge <- pg[whiche]/(sum(pg[whiche]))
        getSE(whiche, pge)
    })
    dynamic.att <- mean(dynamic.att.e)
    keepers <- which(group <= t)
    dynamic.weights <- lapply(eseq, function(e) {
        whiche <- which(t - group + 1 == e)
        pge <- pg[whiche]/(sum(pg[whiche]))
        list(whiche=whiche, pge=pge)
    })
    which.dynamic.weights <- unlist(lapply(dynamic.weights, function(d) d$whiche))
    dynamic.weights <- unlist(lapply(dynamic.weights, function(d) d$pge))[order(which.dynamic.weights)] / length(unique(eseq))
    dynamic.se <- getSE(keepers, dynamic.weights)

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
        getSE(whicht, pgt)
    })
    calendar.att <- mean(calendar.att.t)
    keepers <- which(group <= t)
    calendar.weights <- lapply(tseq, function(t1) {
        whicht <- which((group <= t1) & (t==t1))
        pgt <- pg[whicht]/(sum(pg[whicht]))
        list(whicht=whicht, pgt=pgt)
    })
    which.calendar.weights <- unlist(lapply(calendar.weights, function(t1) t1$whicht))
    calendar.weights <- unlist(lapply(calendar.weights, function(t1) t1$pgt))[order(which.calendar.weights)] / length(unique(tseq))
    calendar.se <- getSE(keepers, calendar.weights)

    eseq <- seq(1,max(group-t)-1)
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
        list(se=sapply(eseq, function(e) {
            whiche <- which( (t - group + 1 == e) &
                             ( max(t) - group + 1 >= e1) &
                             ( e <= e1 ))
            pge <- pg[whiche]/sum(pg[whiche])
            getSE(whiche, pge)
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
    ## kind of hack, but just post-process the weights
    dynsel.weights <- lapply(dynsel.weights, function(d) {
        d$pge <- d$pge/sum(d$pge)
        d
    })
    dynsel.se.e1 <- sapply(dynsel.weights, function(d) {
        getSE(d$whiche, d$pge)
    })

    AGGTE(simple.att=simple.att, simple.se=simple.se, selective.att=selective.att, selective.se=selective.se, selective.att.g=selective.att.g, selective.se.g=selective.se.g, dynamic.att=dynamic.att, dynamic.se=dynamic.se, dynamic.att.e=dynamic.att.e, dynamic.se.e=dynamic.se.e, calendar.att=calendar.att, calendar.se=calendar.se, calendar.att.t=calendar.att.t, calendar.se.t=calendar.se.t, dynsel.att.e1=dynsel.att.e1, dynsel.se.e1=dynsel.se.e1, dynsel.att.ee1=dynsel.att.ee1, dynsel.se.ee1=dynsel.se.ee1)
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
#' @param selective.att.g the standard error for \code{selective.att.g}
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
AGGTE <- function(simple.att=NULL, simple.se=NULL, selective.att=NULL, selective.se=NULL, selective.att.g=NULL, selective.se.g=NULL, dynamic.att=NULL, dynamic.se=NULL, dynamic.att.e=NULL, dynamic.se.e=NULL, calendar.att=NULL, calendar.se=NULL, calendar.att.t=NULL, calendar.se.t=NULL, dynsel.att.e1=NULL, dynsel.se.e1=NULL, dynsel.att.ee1=NULL, dynsel.se.ee1=NULL) {
    out <- list(simple.att=simple.att, simple.se=simple.se, selective.att=selective.att, selective.se=selective.se, selective.att.g=selective.att.g, selective.se.g=selective.se.g, dynamic.att=dynamic.att, dynamic.se=dynamic.se, dynamic.att.e=dynamic.att.e, dynamic.se.e=dynamic.se.e, calendar.att=calendar.att, calendar.se=calendar.se, calendar.att.t=calendar.att.t, calendar.se.t=calendar.se.t, dynsel.att.e1=dynsel.att.e1, dynsel.se.e1=dynsel.se.e1, dynsel.att.ee1=dynsel.att.ee1, dynsel.se.ee1=dynsel.se.ee1)
    class(out) <- "AGGTE"
    out
}
