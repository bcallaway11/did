#' @title Pre-Test of Conditional Parallel Trends Assumption
#'
#' @description An integrated moments test for the conditional parallel trends
#'  assumption holding in all pre-treatment time periods for all groups
#'
#' @inheritParams att_gt
#' @param cores The number of cores to use for parallel processing. This
#'  parallelizes Step 1 (computing the test statistic); Step 2's multiplier
#'  bootstrap is vectorized and runs in a single process.
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
    stop("There are no pre-treatment periods available to conduct the conditional pre-test. Ensure your data has sufficient pre-treatment periods.")
  }


  # set which weight function to use
  # the only option that will work with current setup
  # is indicator so hard-code it here
  weightfun <- indicator


  if (allow_unbalanced_panel) {
    stop("The conditional pre-test is not supported for unbalanced panel data. Set `allow_unbalanced_panel = FALSE` to balance the panel, or use the pre-test reported by att_gt() instead.")
  }

  # create dataset with n observations;
  # recover covariates from this dataset
  if (panel & (allow_unbalanced_panel == FALSE)) {
    dta <- data[ data[,tname]==tlist[1], ]
  } else {
    dta <- data
  }

  # need X twice, once to loop over, and once to calculate weighting function
  X <- model.matrix(xformla, dta)
  X1 <- X

  #for debugging:
  # X <- as.matrix(X[1:100,])

  # Internal compute.att_gt speedups for the per-unit loop below:
  #  (1) pretreatment_cells_only -- post-treatment (g,t) cells are dropped
  #      right after the loop (the keepers filter retains only group > t), so
  #      compute.att_gt never computes them (this was the TODO below).
  #  (2) setup_precomp / y_override -- compute.att_gt's per-call setup (data
  #      copy, design matrix, per-period slicing, alignment check) depends on
  #      the outcome only through the per-period outcome slices, and the loop
  #      below only modifies the outcome column; build the y-invariant bundle
  #      once and pass each unit's modified outcome as a plain vector instead
  #      of copying the whole data.frame n times. Used only when the balanced,
  #      id-aligned panel precompute applies; otherwise the legacy data-copy
  #      path inside the loop is kept.
  # Both are bit-identical to the legacy loop (same kept cells, same inputs to
  # every 2x2 estimator, and compute.att_gt consumes no RNG).
  basedp <- dp
  basedp$print_details <- FALSE
  basedp$pretreatment_cells_only <- TRUE
  use_y_override <- FALSE
  if (dp$panel) {
    setupdp <- basedp
    setupdp$return_setup <- TRUE
    setup_precomp <- compute.att_gt(setupdp)
    if (!is.null(setup_precomp) && isTRUE(setup_precomp$use_precompute_panel)) {
      basedp$setup_precomp <- setup_precomp
      use_y_override <- TRUE
      id_col <- data[, idname]
      y_col <- data[, yname]
    }
  }

  cat("Step 1 of 2: Computing test statistic....\n")
  out <- pbapply::pblapply(1:nrow(X), function(i) {
    # these are the weights for the conditional moment test
    # indicator weights
    www <- as.numeric(weightfun(X1, X[i,]))
    # for indicator weights, just choose rows where weights = 1
    rightids <- dta[,idname][www==1]

    # set new parameters to pass to call to compute.att_gt
    thisdp <- basedp

    # set the outcome to be the outcomes multiplied by the weighting function
    # (***this *trick* is only going to work for indicator weights***
    # (otherwise you will multiply by www twice))
    if (use_y_override) {
      # pass only the modified outcome vector; the y-invariant setup bundle
      # built above is reused by compute.att_gt
      ynew <- numeric(length(y_col))
      sel <- id_col %in% rightids
      ynew[sel] <- y_col[sel]
      thisdp$y_override <- ynew
    } else {
      # legacy path: copy the whole dataset and overwrite the outcome column
      thisdata <- data
      thisdata[,yname] <- 0
      thisdata[ thisdata[,idname] %in% rightids,yname] <- data[ thisdata[,idname] %in% rightids,yname]
      thisdata$y <- thisdata[,yname]
      thisdp$data <- thisdata
    }

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

  # drop post-treatment (g,t). compute.att_gt now skips these cells entirely
  # (pretreatment_cells_only above), so this filter is a no-op retained as a
  # safety net in case a code path reintroduces post-treatment cells.
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

  # Orient J as (nX rows = covariate threshold values) x (n_gt cols = pre-treatment
  # (g,t) cells) so the observed CvM statistic below uses the SAME (1/nX)
  # normalization as the bootstrap null distribution in test.mboot() (which forms
  # Jb = t(apply(Ub*inf.func, c(2,3), mean)) and is therefore (nX x n_gt)).
  #
  # sapply() returns an (n_gt x nX) MATRIX when there are multiple pre-treatment
  # cells, but a length-nX VECTOR when there is only one -- so transpose the matrix
  # case and column-ify the vector case.
  #
  # BUG FIX (gh issue: pre-test spuriously rejects): this was previously
  #   ifelse(class(J.inner) == "matrix", J <- t(J.inner), J <- as.matrix(J.inner))
  # which was correct under R < 4.0 but silently broke in R >= 4.0, where
  # class(<matrix>) became c("matrix","array") (length 2). ifelse() then evaluated
  # BOTH assignment branches and the second (as.matrix, no transpose) always won, so
  # for the multi-cell case J stayed (n_gt x nX). Because sum(apply(J^2, 2, mean)) =
  # (1/nrow(J)) * sum(J^2), that mis-orientation scaled the observed CvM by
  # nX/n_gt = n/n_gt relative to the bootstrap, driving the p-value to ~0 and making
  # the pre-test over-reject whenever there was more than one pre-treatment (g,t)
  # cell. Use is.matrix() to restore the intended orientation.
  if (is.matrix(J.inner)) {
    J <- t(J.inner)
  } else {
    J <- as.matrix(J.inner)
  }

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
  # Row i is 1 iff X[i, j] <= u[j] for EVERY column j. Vectorized: rep(u, each =
  # nrow(X)) lines u up with X's column-major layout, so (X <= .) is the elementwise
  # comparison and a row passes iff none of its entries exceed u (rowSums == ncol).
  # Bit-identical to the prior t(apply(X, 1, `<=`)) + apply(all) for finite X
  # (model.matrix output is finite), but O(n*p) without the per-row closure.
  1 * (rowSums(X <= rep(u, each = nrow(X))) == ncol(X))
}



#' @title Multiplier Bootstrap for Conditional Moment Test
#'
#' @description A slightly modified multiplier bootstrap procedure
#'  for the pre-test of the conditional parallel trends assumption
#'
#' @inheritParams mboot
#' @param cores Unused; retained for backward compatibility. The
#'  multiplier bootstrap is computed with vectorized matrix operations
#'  in a single process, so this argument has no effect here. In
#'  \code{conditional_did_pretest}, \code{cores} parallelizes only
#'  Step 1 (computing the test statistic).
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
  if (panel) {
    dta <- data[ data[,tname]==tlist[1], ]
  } else {
    dta <- data
  }
  n <- nrow(dta)

  # if include id as variable to cluster on
  # drop it as we do this automatically
  if (idname %in% clustervars) {
    clustervars <- clustervars[-which(clustervars==idname)]
  }

  # we can only handle up to 2-way clustering
  # (in principle could do more, but not high priority now)
  if (length(clustervars) > 1) {
    stop("At most one cluster variable (beyond 'idname') is supported. Please reduce to one.")
  }

  # ---------------------------------------------------------------------------
  # Vectorized multiplier bootstrap (replaces a biters-long pbsapply loop that
  # multiplied a full n x k x nX influence array by fresh Rademacher weights and
  # ran apply(., c(2,3), mean) on every draw -- O(n*k*nX) compute AND an O(n*k*nX)
  # transient array allocation per draw, ~1.2 GB x biters at n in the thousands).
  #
  # The per-draw statistic is, writing inf[i,j,x] for the influence array,
  #   Jb[x,j] = (1/n) sum_i U_b[i] inf[i,j,x]
  #   CvMb_b  = n * mean_x( sum_j Jb[x,j]^2 ) = (n/nX) * sum_{x,j} Jb[x,j]^2.
  # Stacking the weights of all draws into V (n x biters) and reshaping inf to
  # M (n x k*nX), Jb for every draw is crossprod(V, M)/n in one BLAS call, and
  # CvMb is (n/nX) * rowSums(.^2). Drawing all biters*n (or biters*n1) Rademacher
  # values in a single sample() call consumes the RNG stream identically to the
  # per-draw loop (each value is one R_unif_index(2)), so results match the old
  # loop up to floating-point summation order (~1e-14).
  # ---------------------------------------------------------------------------
  dms <- dim(inf.func)
  k <- dms[2]; nX <- dms[3]

  # Build the n x biters multiplier matrix; columns are independent draws.
  if (length(clustervars) > 0) {
    # Rademacher weights are constant within cluster; map each unit to its cluster's
    # weight via match() exactly as the original loop did.
    uc <- unique(dta[, clustervars])
    n1 <- length(uc)
    Vc <- matrix(sample(c(-1, 1), as.numeric(n1) * biters, replace = TRUE), nrow = n1)
    V <- Vc[match(dta[, clustervars], uc), , drop = FALSE]
  } else {
    V <- matrix(sample(c(-1, 1), as.numeric(n) * biters, replace = TRUE), nrow = n)
  }

  # Tile over the X dimension so peak memory is n*k*chunk doubles instead of the
  # full n*k*nX reshape. CvMb is a plain sum over X, so chunking is exact (the only
  # difference is summation order). Chunk targets ~5e6 doubles (~40 MB) per slice.
  chunk <- max(1L, as.integer(5e6 %/% (as.numeric(n) * k)))
  bout <- numeric(biters)
  for (start in seq.int(1L, nX, by = chunk)) {
    xs <- start:min(start + chunk - 1L, nX)
    Mc <- matrix(inf.func[, , xs, drop = FALSE], nrow = n)   # n x (k*length(xs))
    Jf <- crossprod(V, Mc) / n                                # biters x (k*length(xs))
    bout <- bout + rowSums(Jf^2)
  }
  bout <- (n / nX) * bout

  crit.val <- quantile(bout, 1-alp, type=1)

  list(bres=bout, crit.val=crit.val)
}
