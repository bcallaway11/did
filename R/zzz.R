utils::globalVariables(c("year","att","ate.se","post","group"))

#' @title citation
#'
#' @description print the citation for the relevant paper
#'
#' @keywords internal
citation <- function() {
    cat("\n Callaway and Sant'Anna (2018) Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment <https://ssrn.com/abstract=3148250> \n")
}
