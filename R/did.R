#' @title mp.spatt
#'
#' @description \code{mp.spatt} computes the ATT in the case where there are more
#'  than two periods of data and allowing for treatment to occur at different points in time
#'  extending the method of Abadie (2005).  This method relies on once individuals are treated
#'  they remain in the treated state for the duration.
#'
#' @param first.treat.name Give the column name of the variable that forms groups based on when an observation is first treated
#'
#' @inheritParams spatt
#'
#' @return \code{QTE} object
#' 
#' @export
mp.spatt <- function(formla, xformla=NULL, data, tname, w=NULL, panel=FALSE,
                     idname=NULL, first.treat.name,
                     iters=100, alp=0.05, method="logit", plot=FALSE, se=TRUE,
                     bstrap=FALSE, biters=100, clustervars=NULL,
                     cband=FALSE, citers=100,
                     retEachIter=FALSE, seedvec=NULL, pl=FALSE, cores=2,
                     printdetails=TRUE) {


    ##TODO: make this handle passing in treatment indicators in a more "natural" way
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

    ## get all the results; importantly this now returns the influence functions
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
                      iters=iters, alp=alp, method=method, plot=plot, se=se,
                      retEachIter=retEachIter, seedvec=seedvec, pl=pl, cores=cores),
                      group=flist[f], year=tlist[(t+1)], post=1*(flist[f]<=tlist[(t+1)]))
            inffunc[f,t,S] <- satt[[t]]$inffunc
            ## browser()
            ## ifu <- satt[[t]]$inffuncu
            ## t(ifu)%*%as.matrix(ifu)/n
        }
        fatt[[f]] <- c(satt, group=flist[f])
    }

    
    
    
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
                colnames(Vb)[1] <- "clvar"
                Ub <- data.frame(dta[,clustervars])
                colnames(Ub)[1] <- "clvar"
                Ub <- merge(Ub, Vb, by="clvar")
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

    return(list(group=group, t=t, att=att, V=V, c=c, inffunc=inffunc1))
}


## The idea here is to combine the weighting function with Y and run the previous
## code for computing group-time average treatment effects
mp.spatt.test <- function(formla, xformla=NULL, data, tname, w=NULL, panel=FALSE,
                     idname=NULL, first.treat.name,
                     iters=100, alp=0.05, method="logit", plot=FALSE, se=TRUE,
                     bstrap=FALSE, biters=100, clustervars=NULL,
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

    ## 
    X1 <- apply(model.matrix(xformla, data), 2, function(col) {
        pnorm(col)
    })## for the entire dataset
    X <- apply(model.matrix(xformla, dta), 2, function(col) {
        pnorm(col)
    })

    ## Extra code for debugging:
    ## for (i in 1:nrow(X)) {
    ##     www <- exp(X1%*%X[i,])##plogis(X1%*%X[i,]) ##(1*(apply((X1 <= X[i,]), 1, all))) -- this definitely won't work because of the dummy variables.
    ##     print(max(www))
    ## }

     ##TODO: just get results for just untreated g-t
    cat("\nStep 1 of 2:.....................\n")
    out <- pbapply::pblapply(1:nrow(X), function(i) {
        www <- exp(X1%*%X[i,])##plogis(X1%*%X[i,]) ##(1*(apply((X1 <= X[i,]), 1, all)))
        data$lhs <- data$y * www
        formula.tools::lhs(formla) <- as.name("lhs")
        out1 <- mp.spatt(formla=formla, xformla=xformla, data=data, tname=tname, panel=panel, idname=idname, first.treat.name=first.treat.name, iters=iters, alp=alp, se=TRUE, bstrap=FALSE, printdetails=FALSE)
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

    cat("\nStep 2 of 2:.....................\n")
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
                colnames(Vb)[1] <- "clvar"
                Ub <- data.frame(dta[,clustervars])
                colnames(Ub)[1] <- "clvar"
                Ub <- merge(Ub, Vb, by="clvar")
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

    return(list(KS=KS, KSb=KSb, cval=ocval, pval=1-ecdf(KSb)(KS)))
}
