#' @title Plot \code{did} objects using \code{ggplot2}
#'
#' @description Function to plot objects from the \code{did} package
#'
#' @param object either a \code{MP} object or \code{AGGTEobj} object
#' @param ... other arguments
#'
#' @export
ggdid <- function(object, ...) {
  UseMethod("ggdid", object)
}


## #' @param type the type of plot, should be one of "attgt", "dynamic",
## #'  "selective", "calendar", "dynsel".  "attgt" is the default and plots
## #'  all group-time average treatment effects separately by group (including
## #'  pre-treatment time periods); "dynamic" plots dynamic treatment effects --
## #'  these are the same as event studies; "selective" plots average effects
## #'  of the treatment separately by group (which allows for selective treatment
## #'  timing); "calendar" plots average treatment effects by time period; and
## #'  "dynsel" plots dynamic effects allowing for selective treatment timing
## #'  (this also requires setting the additional paramater e1)


#' @title Plot \code{MP} objects using \code{ggplot2}
#'
#' @description A function to plot \code{MP} objects
#'
#' @inheritParams ggdid
#' @param ylim optional y limits for the plot; settng here makes the y limits
#'  the same across different plots
#' @param xlab optional x-axis label
#' @param ylab optional y-axis label
#' @param title optional plot title
#' @param xgap optional gap between the labels on the x-axis.  For example,
#'  \code{xgap=3} indicates that the labels should show up for every third
#'  value on the x-axis.  The default is 1.
#' @param ncol The number of columns to include in the resulting plot.  The
#'  default is 1.
#'
#' @export
ggdid.MP <- function(object,
                     ylim=NULL,
                     xlab=NULL,
                     ylab=NULL,
                     title="Group",
                     xgap=1,
                     ncol=1,
                     ...) {

  mpobj <- object

  G <- length(unique(mpobj$group))
  Y <- length(unique(mpobj$t))## drop 1 period bc DID
  g <- unique(mpobj$group)[order(unique(mpobj$group))] ## -1 to drop control group
  y <- unique(mpobj$t)

  results <- data.frame(year=rep(y,G))
  results$group <- unlist(lapply(g, function(x) { rep(x, Y) }))
  results$att <- mpobj$att
  n <- mpobj$n
  results$att.se <- sqrt(diag(mpobj$V)/n)
  results$post <- as.factor(1*(results$year >= results$group))
  results$year <- as.factor(results$year)
  results$c <- mpobj$c
  vcovatt <- mpobj$V/n
  alp <- mpobj$alp

  mplots <- lapply(g, function(g) {
    thisdta <- subset(results, group==g)
    gplot(thisdta, ylim, xlab, ylab, title, xgap)
  })

  do.call("grid.arrange", c(mplots))
}


#' @title Plot \code{AGGTEobj} objects
#'
#' @description A function to plot \code{AGGTEobj} objects
#'
#' @inheritParams ggdid.MP
#'
#' @export
ggdid.AGGTEobj <- function(object,
                           ylim=NULL,
                           xlab=NULL,
                           ylab=NULL,
                           title="",
                           xgap=1,
                           ...) {
  
  post.treat <- 1*(object$egt >= 0)
  results <- cbind.data.frame(year=as.factor(object$egt),
                              att=object$att.egt,
                              att.se=object$se.egt,
                              post=as.factor(post.treat))
  results$c <- object$crit.val.egt
  
  p <- gplot(results, ylim, xlab, ylab, title, xgap)
  p
}
