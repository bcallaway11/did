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
