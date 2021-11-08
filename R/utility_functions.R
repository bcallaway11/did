#' @title trimmer
#'
#' @description A utility function to find observations that appear
#'  to violate support conditions.  This function is not called
#'  anywhere in the code, but it is just useful for debugging some
#'  common issues that users run into.
#'
#' @param g is a particular group (below I pass in 2009)
#' @inheritParams att_gt
#' @param threshold the cutoff for which observations are flagged as
#'  likely violators of the support condition.
#'
#' @return list of ids of observations that likely violate support conditions
#' 
#' @export
trimmer <- function(g, tname, idname, gname, xformla, data, control_group="notyettreated", threshold=.999) {

  time.period <- data[,tname]
  this.data <- data[time.period == (g-1),]
  if (control_group == "notyettreated") {
    # not yet treated
    this.data <- this.data[(this.data[,gname] >= g) |
                           (this.data[,gname] == 0), ]
  } else {
    # never treated
    this.data <- this.data[(this.data[,gname] == g) |
                           (this.data[,gname] == 0), ]
  }
  this.data$D <- 1*this.data[,gname]==g
  this.pscore_reg <- glm(BMisc::toformula("D", BMisc::rhs.vars(xformla)),
                         data=this.data,
                         family=binomial(link=logit))
  this.pscore <- predict(this.pscore_reg, type="response")
  dropper <- (this.pscore >  threshold) & (this.data$D==1)
  if (sum(dropper) > 0) {
    print("hard to match treated observations: ")
    print(this.data[dropper,idname])
    return(this.data[dropper,idname])
  }
}

