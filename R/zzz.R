utils::globalVariables(c("year","att","att.se","post","group","x"))

#' @title citation
#'
#' @description print the citation for the relevant paper
#'
#' @keywords internal
citation <- function() {
  cat("Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  \"Difference-in-Differences with Multiple Time Periods.\" Journal of Econometrics, Vol. 225, No. 2, pp. 200-230, 2021. <https://doi.org/10.1016/j.jeconom.2020.12.001>, <https://arxiv.org/abs/1803.09015> \n")
}
