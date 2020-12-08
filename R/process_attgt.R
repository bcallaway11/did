#' @title Process Results from \code{\link{compute.att_gt}}
#'
#' @param attgt.results.list list of results from \code{\link{compute.att_gt}}
#'
#' @return list with elements:
#' \item{group}{which group a set of results belongs to}
#' \item{tt}{which time period a set of results belongs to}
#' \item{att}{the group time average treatment effect}
#' \item{inf.func}{the influence function for that group time average treatment effect}
#'
#' @export
process_attgt <- function(attgt.results.list) {
  attgt.list <- attgt.results.list$attgt
  inffunc <- attgt.results.list$inffunc
  nG <- length(unique(unlist(BMisc::getListElement(attgt.list, "group"))))
  nT <- length(unique(unlist(BMisc::getListElement(attgt.list, "year"))))
  # pick up number of observations from the influence function
  # this will be flexible for panel or repeated cross sections

  #n <- length(inffunc[1,1,]) 

  n <- length(inffunc[,1])

  # create vectors to hold the results
  group <- c()
  att <- c()
  tt <- c()
  i <- 1

  # matrix to hold influence function
  # (note: this is relying on having a balanced panel,
  # which we do currently enforce)

  #inffunc1 <- matrix(0, ncol=nG*nT, nrow=n) 

  # populate result vectors and matrices
  for (f in 1:nG) {
    for (s in 1:nT) {
      group[i] <- attgt.list[[i]]$group
      tt[i] <- attgt.list[[i]]$year
      att[i] <- attgt.list[[i]]$att
      #inffunc1[,i] <- inffunc[f,s,]
      i <- i+1
    }
  }

  list(group=group, att=att, tt=tt, inf.func=inffunc)
}
