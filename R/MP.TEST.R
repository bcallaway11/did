#' @title MP.TEST
#'
#' @description An object that holds results from computing pre-test of the
#'  conditional parallel trends assumption
#'
#' @param CvM Cramer von Mises test statistic
#' @param CvMb a vector of bootstrapped Cramer von Mises test statistics
#' @param CvMcval CvM critical value
#' @param CvMpval p-value for CvM test
#' @param KS Kolmogorov-Smirnov test statistic
#' @param KSb a vector of bootstrapped KS test statistics
#' @param KScval KS critical value
#' @param KSpval p-value for KS test
#' @param clustervars vector of which variables were clustered on for the test
#' @param xformla formla for the X variables used in the test
#'
#' @export
MP.TEST <- function(CvM=NULL, CvMb=NULL, CvMcval=NULL, CvMpval=NULL, KS=NULL, KSb=NULL, KScval=NULL, KSpval=NULL, clustervars=NULL, xformla=NULL) {
    out <- list(CvM=CvM, CvMb=CvMb, CvMcval=CvMcval, CvMpval=CvMpval, KS=KS, KSb=KSb, KScval=KScval, KSpval=KSpval, clustervars=clustervars, xformla=xformla)
    class(out) <- "MP.TEST"
    out
}

#' @title summary.MP.TEST
#'
#' @description print a summary of test results
#'
#' @param object an MP.TEST object
#' @param ... other variables
#'
#' @export
summary.MP.TEST <- function(object, ... ) {
    CvM <- object$CvM
    CvMcval <- object$CvMcval
    CvMpval <- object$CvMpval
    citation()
    cat("Cramer von Mises: \n")
    cat("  Test Statistic: ", CvM, "\n")
    cat("  Critical Value: ", CvMcval, "\n")
    cat("  P-value       : ", CvMpval, "\n \n")
    cat("Clustering on   : ", paste0(object$clustervars,sep=","), "\n")
    cat("X formula       : ", as.character(object$xformla), "\n")
}

