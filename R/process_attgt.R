#' @title Process Results from [compute.att_gt()]
#'
#' @param attgt.list list of results from [compute.att_gt()]
#'
#' @return list with elements:
#' \item{group}{which group a set of results belongs to}
#' \item{tt}{which time period a set of results belongs to}
#' \item{att}{the group time average treatment effect}
#'
#' @export
process_attgt <- function(attgt.list) {
  nG <- length(unique(vapply(attgt.list, `[[`, numeric(1), "group")))
  nT <- length(unique(vapply(attgt.list, `[[`, numeric(1), "year")))

  if (length(attgt.list) != nG * nT) {
    stop("Number of results does not match expected group-time combinations")
  }

  group <- vapply(attgt.list, `[[`, numeric(1), "group")
  att   <- vapply(attgt.list, `[[`, numeric(1), "att")
  tt    <- vapply(attgt.list, `[[`, numeric(1), "year")

  list(group=group, att=att, tt=tt)
}
