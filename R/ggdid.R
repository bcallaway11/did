#' @title ggdid
#'
#' @description Function to plot \code{MP} objects
#'
#' @param mpobj an \code{MP} object
#' @param type the type of plot, should be one of "attgt", "dynamic",
#'  "selective", "calendar", "dynsel".  "attgt" is the default and plots
#'  all group-time average treatment effects separately by group (including
#'  pre-treatment time periods); "dynamic" plots dynamic treatment effects --
#'  these are the same as event studies; "selective" plots average effects
#'  of the treatment separately by group (which allows for selective treatment
#'  timing); "calendar" plots average treatment effects by time period; and
#'  "dynsel" plots dynamic effects allowing for selective treatment timing
#'  (this also requires setting the additional paramater e1)
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
#' @param e1 only used when plot type is "dynsel", this specifies the number
#'  of post-treatment periods that need to be available for particular groups
#'  to be included in the resulting plot when there are dynamic treatment
#'  effects and selective treatment timing
#'
#'
#' @export
ggdid <- function(object, ...) {
  UseMethod("ggdid", object)
}
                 

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

  #if (type=="attgt") {
  results <- data.frame(year=rep(y,G))
  results$group <- unlist(lapply(g, function(x) { rep(x, Y) }))##c(rep(2004,G),rep(2006,G),rep(2007,G))
  results$att <- mpobj$att
  n <- mpobj$n
  results$att.se <- sqrt(diag(mpobj$V)/n)
  results$post <- as.factor(1*(results$year >= results$group))
  results$year <- as.factor(results$year)
  results$c <- mpobj$c
  vcovatt <- mpobj$V/n
  alp <- mpobj$alp

  ##results <- mp2ATT(results, vcovatt)

  mplots <- lapply(g, function(g) {
    thisdta <- subset(results, group==g)
    gplot(thisdta, ylim, xlab, ylab, title, xgap)
  })

  do.call("grid.arrange", c(mplots))
  #} else if (type=="dynamic") {
  ##   aggte <- mpobj$aggte
  ##   #if (mpobj$c > 2) warning("uniform bands not implemented yet for this plot...")
  ##   elen <- length(aggte$dynamic.att.e)
  ##   results <- cbind.data.frame(year=as.factor(seq(1:elen)),
  ##                               att=aggte$dynamic.att.e,
  ##                               att.se=aggte$dynamic.se.e,
  ##                               post=as.factor(1),
  ##                               c=aggte$c.dynamic,
  ##                               alp = mpobj$alp)
  ##   #qnorm(.975))
  ##   p <- gplot(results, ylim, xlab, ylab, title, xgap)
  ##   p
  ## }
}

ggdid.AGGTEobj <- function(object,
                           ylim=NULL,
                           xlab=NULL,
                           ylab=NULL,
                           title="Group",
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
