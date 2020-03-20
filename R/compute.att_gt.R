#' @title compute.att_gt
#'
#' @description \code{compute.att_gt} does the main work for computing
#'  mutliperiod group-time average treatment effects
#'
#' @inheritParams att_gt
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
compute.att_gt <- function(nG,
                           nT,
                           glist,
                           tlist,
                           data,
                           dta,
                           first.treat.name,
                           yname,
                           tname,
                           w,
                           idname,
                           xformla,
                           method,
                           seedvec,
                           pl,
                           cores,
                           printdetails,
                           nevertreated,
                           estMethod,
                           panel) {


  # will populate with all att(g,t)
  attgt.list <- list()

  # place holder in lists
  counter <- 1

  # number of time periods
  tlist.length <- length(tlist)

  # 3-dimensional array which will store influence function
  # across groups and times
  inffunc <- array(data=0, dim=c(nG,nT,nrow(dta)))

  # loop over groups
  for (g in 1:nG) {

    # loop over time periods
    for (t in 1:(tlist.length-1)) {


      # set pre-treatment time period (this is updated later
      # if g <= t (i.e. for "already treated" groups)
      pret <- t
      
      
      # code to update pre-treatment time periods
      if (glist[g]<=tlist[(t+1)]) {

        # set an index for the pretreatment period
        # this recovers the right pre-treatment period for this group
        # it is the most recent pre-treatment period (g-1)
        pret <- utils::tail(which(tlist < glist[g]),1)
        
        # print a warning message if there are no pre-treatment period
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the group first treated at ", glist[g], "\nUnits from this group are dropped"))

          # if there are not pre-treatment periods, code will
          # break, jump out of this loop
          break
        }
        
      }

      # print the details of which iteration we are on
      if (printdetails) {
        cat(paste("current period:", tlist[(t+1)]), "\n")
        cat(paste("current group:", glist[g]), "\n")
        cat(paste("set pretreatment period to be", tlist[pret]), "\n")
      }

      #-----------------------------------------------------------------------------
      # results for the case with panel data
      #-----------------------------------------------------------------------------
      if (panel) {
        
        # post treatment dummy variable
        post.treat <- 1*(glist[g]<=tlist[t+1])

        # get dataset with current period and pre-treatment period
        disdat <- data[(data[,tname]==tlist[t+1] | data[,tname]==tlist[pret]),]

        # sete up control group
        if(nevertreated ==T){
          # use the "never treated" group as the control group
          disdat$C <- 1*(disdat[,first.treat.name] == 0)
        }
        if(nevertreated ==F){
          # use "not yet treated as control"
          # that is, never treated + units that are eventually treated,
          # but not treated by the current period
          disdat$C <- 1*((disdat[,first.treat.name] == 0) |
                           (disdat[,first.treat.name] > tlist[t+1]))
          ## disdat$C <- 1*((disdat[,first.treat.name] == 0) +
          ##                  (disdat[,first.treat.name] > max(disdat[, tname]))) *
          ##   (disdat[,first.treat.name] != glist[g])
        }

        # set up dummy for particular treated group
        disdat$G <- 1*(disdat[,first.treat.name] == glist[g])

        # transform  disdat it into "cross-sectional" data where one of the columns
        # contains the change in the outcome over time.
        # dy is computed as latest year - earliest year. "Y" is outcome
        # in the pre period, "yt1" is outcome in the post period
        disdat <- BMisc::panel2cs(disdat, yname, idname, tname)

        # still total number of units (not just included in G or C)
        n <- nrow(disdat)

        # pick up the indices for units that will be used to compute ATT(g,t)
        disidx <- disdat$G==1 | disdat$C==1

        # pick up the data that will be used to compute ATT(g,t)
        disdat <- disdat[disidx,]

        # drop missing factors
        disdat <- base::droplevels(disdat)

        # give short names for data in this iteration
        G <- disdat$G
        C <- disdat$C
        Ypre <- disdat$Y
        Ypost <- disdat$yt1
        dy <- disdat$dy 
        n1 <- nrow(disdat) # num obs. for computing ATT(g,t)
        w <- disdat$w


        
        # Currently, we still pass to functions if there are no covariates
        # but perhaps
        # if there are no covariates, just manually compute
        # att_gt

        # matrix of covariates
        covariates <- model.matrix(xformla, data=disdat)


        #-----------------------------------------------------------------------------
        # code for actually computing att(g,t)
        #-----------------------------------------------------------------------------

        if (class(estMethod) == "function") {
          # user-specified function
          attgt <- estMethod(Y1=Ypost, Y0=Ypre,
                             treat=G,
                             covariates=covariates)
        } else if (estMethod == "ipw") {
          # inverse-probability weights
          attgt <- DRDID::ipw_did_panel(Ypost, Ypre, G,
                                        covariates=covariates,
                                        boot=FALSE)
        } else if (estMethod == "reg") {
          # regression
          attgt <- DRDID::reg_did_panel(Ypost, Ypre, G,
                                        covariates=covariates,
                                        boot=FALSE)
        } else {
          # doubly robust, this is default
          attgt <- DRDID::drdid_panel(Ypost, Ypre, G,
                                      covariates=covariates,
                                      boot=FALSE)
        }
        
        # save results for this att(g,t)
        attgt.list[[counter]] <- list(att=attgt$ATT, group=glist[g], year=tlist[(t+1)], post=post.treat)

        # recover the influence function
        # start with vector of 0s because influence function
        # for units that are not in G or C will be equal to 0
        inf.func <- rep(0, n)

        # n/n1 adjusts for estimating the
        # att_gt only using observations from groups
        # G and C
        
        inf.func[disidx] <- (n/n1)*attgt$inf.func

        # save it in influence function matrix
        inffunc[g,t,] <- inf.func

        # update counter
        counter <- counter+1
      } # end panel 
    } # end looping over t

  }

  return(list(attgt.list=attgt.list, inffunc=inffunc))
}
