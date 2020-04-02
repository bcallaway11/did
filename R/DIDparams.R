#' @title DIDparams
#'
#' @description object to hold did parameters
#'
#' @inheritParams att_gt
#' @param n The number of observations.  This is equal to the
#'  number of units (which may be different from the number
#'  of rows in a panel dataset).
#' @param nG The number of groups
#' @param nT The number of time periods
#' @param tlist a vector containing each time period
#' @param glist a vector containing each group
#'
#' @export
DIDparams <- function(yname, 
                   tname,
                   idname=NULL,
                   first.treat.name,
                   xformla=NULL,
                   data,
                   control.group,
                   weightsname=NULL,
                   alp=0.05,
                   bstrap=T,
                   biters=1000,
                   clustervars=NULL,
                   cband=T,
                   printdetails=TRUE,
                   pl=FALSE,
                   cores=1,
                   estMethod="dr",
                   panel=TRUE,
                   n=NULL,
                   nG=NULL,
                   nT=NULL,
                   tlist=NULL,
                   glist=NULL) {

  out <- list(yname=yname,
              tname=tname,
              idname=idname,
              first.treat.name=first.treat.name,
              xformla=xformla,
              data=data,
              control.group=control.group,
              weightsname=weightsname,
              alp=alp,
              bstrap=bstrap,
              biters=biters,
              clustervars=clustervars,
              cband=cband,
              printdetails=printdetails,
              pl=pl,
              cores=cores,
              estMethod=estMethod,
              panel=panel,
              n=n,
              nG=nG,
              nT=nT,
              tlist=tlist,
              glist=glist)
  class(out) <- "DIDparams"
  out
}
