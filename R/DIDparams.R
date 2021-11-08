#' @title DIDparams
#'
#' @description Object to hold did parameters that are passed across functions
#'
#' @inheritParams att_gt
#' @inheritParams pre_process_did
#' @param n The number of observations.  This is equal to the
#'  number of units (which may be different from the number
#'  of rows in a panel dataset).
#' @param nG The number of groups
#' @param nT The number of time periods
#' @param tlist a vector containing each time period
#' @param glist a vector containing each group
#' @param true_repeated_cross_sections Whether or not the data really
#'  is repeated cross sections.  (We include this because unbalanced
#'  panel code runs through the repeated cross sections code)
#'
#' @export
DIDparams <- function(yname,
                   tname,
                   idname=NULL,
                   gname,
                   xformla=NULL,
                   data,
                   control_group,
                   anticipation=0,
                   weightsname=NULL,
                   alp=0.05,
                   bstrap=TRUE,
                   biters=1000,
                   clustervars=NULL,
                   cband=TRUE,
                   print_details=TRUE,
                   pl=FALSE,
                   cores=1,
                   est_method="dr",
                   base_period="varying",
                   panel=TRUE,
                   true_repeated_cross_sections,
                   n=NULL,
                   nG=NULL,
                   nT=NULL,
                   tlist=NULL,
                   glist=NULL,
                   call=NULL) {

  out <- list(yname=yname,
              tname=tname,
              idname=idname,
              gname=gname,
              xformla=xformla,
              data=data,
              control_group=control_group,
              anticipation=anticipation,
              weightsname=weightsname,
              alp=alp,
              bstrap=bstrap,
              biters=biters,
              clustervars=clustervars,
              cband=cband,
              print_details=print_details,
              pl=pl,
              cores=cores,
              est_method=est_method,
              base_period=base_period,
              panel=panel,
              true_repeated_cross_sections=true_repeated_cross_sections,
              n=n,
              nG=nG,
              nT=nT,
              tlist=tlist,
              glist=glist,
              call=call)
  class(out) <- "DIDparams"
  out
}
