#' @title summary.AGGTE
#'
#' @description print a summary of an AGGTE object
#'
#' @param object an AGGTE object
#' @param type which type of summary to print, options are "dynamic", "selective", "calendar", and "dynsel"
#' @param e1 if the type is "dynsel", this is the number of post-treatment periods required in order for a group to be used to construct aggregated parameters with selective treatment timing and dynamic effects; otherwise not used
#' @param ... other variables
#'
#' @export
summary.AGGTE <- function(object, type=c("dynamic"), e1=1, ...) {
  citation()
  sep <- "          "
  cat("Overall Summary Measures", "\n")
  cat("------------------------", "\n")
  cat("Simple ATT    : ", object$simple.att, "\n")
  cat("  SE          : ", object$simple.se, "\n")
  cat("Dynamic ATT   : ", object$dynamic.att, "\n")
  cat("  SE          : ", object$dynamic.se, "\n")
  type <- type[1]
  if (type == "dynamic") {
    cat("Dynamic Treatment Effects", "\n")
    cat("-------------------------")
    elen <- length(object$dynamic.att.e)
    printmat <- cbind(object$dynamic.att.e, object$dynamic.att.e, object$dynamic.se.e)
    colnames(printmat) <- c("e","att","se")
    print(kable(printmat))
  }
}

