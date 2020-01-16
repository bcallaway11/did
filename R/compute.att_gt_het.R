#' @title compute.att_gt_het
#'
#' @description \code{compute.att_gt_het} does the main work for computing
#'  mutliperiod group-time average treatment effects with heterogeneity
#'
#' @inheritParams att_gt_het
#'
#' @return a list with length equal to the number of groups times the
#'  number of time periods; each element of the list contains an
#'  object that contains group-time average treamtent effect as well
#'  as which group it is for and which time period it is for. It also exports
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
    for (t in 1:(tlen)) {
      pret <- utils::tail(which(tlist < flist[f]),1)
      if (flist[f]<=tlist[(t)]) {
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
        ## get dataset with current period and pre-treatment period
        disdat <- data[(data[,tname]==tlist[t] | data[,tname]==tlist[pret]),]

        ## set up control group
        if(nevertreated ==T){
          disdat$C <- 1*(disdat[,first.treat.name] == 0)
        }
        if(nevertreated ==F){
          disdat$C <- 1*((disdat[,first.treat.name] == 0) +
                           (disdat[,first.treat.name] > max(disdat[, tname]))) *
                            (disdat[,first.treat.name] != flist[f])
        }

        ## set up dummy for particular treated group
        disdat$G <- 1*(disdat[,first.treat.name] == flist[f])

        ## transform it into "cross-sectional" data where one of the columns contains the change in the outcome
        ## over time
        disdat <- BMisc::panel2cs(disdat, yname, idname, tname)

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
