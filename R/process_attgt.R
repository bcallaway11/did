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
  group <- vapply(attgt.list, `[[`, numeric(1), "group")
  att   <- vapply(attgt.list, `[[`, numeric(1), "att")
  tt    <- vapply(attgt.list, `[[`, numeric(1), "year")

  list(group=group, att=att, tt=tt)
}
