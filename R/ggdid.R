#' @title Plot `did` objects using `ggplot2`
#'
#' @description Function to plot objects from the `did` package
#'
#' @param object either a `MP` object or `AGGTEobj` object. See `help(ggdid.MP)` and `help(ggdid.AGGTEobj)`.
#' @param ... other arguments
#'
#' @export
ggdid <- function(object, ...) {
  UseMethod("ggdid", object)
}


## #' @param type the type of plot, should be one of "attgt", "dynamic",
## #'  "group", "calendar", "dynsel".  "attgt" is the default and plots
## #'  all group-time average treatment effects separately by group (including
## #'  pre-treatment time periods); "dynamic" plots dynamic treatment effects --
## #'  these are the same as event studies; "group" plots average effects
## #'  of the treatment separately by group (which allows for selective treatment
## #'  timing); "calendar" plots average treatment effects by time period; and
## #'  "dynsel" plots dynamic effects allowing for selective treatment timing
## #'  (this also requires setting the additional paramater e1)


#' @title Plot `MP` objects using `ggplot2`
#'
#' @description A function to plot `MP` objects
#'
#' @inheritParams ggdid
#' @param ylim optional y limits for the plot; setting here makes the y limits
#'  the same across different plots
#' @param xlab optional x-axis label
#' @param ylab optional y-axis label
#' @param title optional plot title
#' @param xgap optional gap between the labels on the x-axis.  For example,
#'  `xgap=3` indicates that the labels should show up for every third
#'  value on the x-axis.  The default is 1.
#' @param ncol The number of columns to include in the resulting plot.  The
#'  default is 1.
#' @param legend Whether or not to include a legend (which will indicate color
#'  of pre- and post-treatment estimates).  Default is `TRUE`.
#' @param group Vector for which groups to include in the plots of ATT(g,t).
#'   Default is NULL, and, in this case, plots for all groups will be included (`ggdid.MP` only).
#' @param ref_line A reference line at this value, usually to compare confidence
#'   intervals to 0. Set to NULL to omit.
#' @param theming Set to FALSE to skip all theming so you can do it yourself.
#' @param grtitle Title to append before each group name (`ggdid.MP` only).
#'
#' @export
ggdid.MP <- function(object,
                     ylim=NULL,
                     xlab=NULL,
                     ylab=NULL,
                     title="Group",
                     xgap=1,
                     ncol=1,
                     legend=TRUE,
                     group=NULL,
                     ref_line = 0,
                     theming = TRUE,
                     grtitle = "Group",
                     ...) {

  mpobj <- object

  G <- length(unique(mpobj$group))
  Y <- length(unique(mpobj$t))## drop 1 period bc DID
  g <- unique(mpobj$group)[order(unique(mpobj$group))] ## -1 to drop control group
  y <- unique(mpobj$t)

  results <- data.frame(year=rep(y,G))
  results$group <- unlist(lapply(g, function(x) { rep(x, Y) }))
  results$grtitle <- paste(grtitle,results$group)
  results$att <- mpobj$att
  n <- mpobj$n
  results$att.se <- mpobj$se #sqrt(diag(mpobj$V)/n)
  results$post <- as.factor(1*(results$year >= results$group))
  results$year <- as.factor(results$year)
  results$c <- mpobj$c
  #vcovatt <- mpobj$V/n
  alp <- mpobj$alp


  if(is.null(group)) group <- g
  if(any(!(group %in% g))){
    warning('Some of the specified groups do not exist in the data. Reporting all available groups.')
    group <- g
  }

  mplots <- gplot(results[results$group %in% group,],
                  ylim, xlab, ylab, title, xgap,
                  legend, ref_line, theming) +
    facet_wrap(~grtitle, ncol = ncol, scales = 'free')

  return(mplots)
}


#' @title Plot `AGGTEobj` objects
#'
#' @description A function to plot `AGGTEobj` objects
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
                           legend=TRUE,
                           ref_line = 0,
                           theming = TRUE,
                           ...) {

  if ( !(object$type %in% c("dynamic","group","calendar")) ) {
    stop(paste0("Plot method not available for this type of aggregation"))
  }

  post.treat <- 1*(object$egt >= 0)
  results <- cbind.data.frame(year=object$egt,
                              att=object$att.egt,
                              att.se=object$se.egt,
                              post=as.factor(post.treat))
  results$c <- ifelse(is.null(object$crit.val.egt), abs(qnorm(.025)), object$crit.val.egt)

  if (title == "") {
    # get title right depending on which aggregation
    title <- ifelse(object$type=="group", "Average Effect by Group", ifelse(object$type=="dynamic", "Average Effect by Length of Exposure", "Average Effect by Time Period"))
  }

  if (object$type == "group") {
    # alternative plot if selective/group treatment timing plot
    p <- splot(results, ylim, xlab, ylab, title, legend, ref_line, theming)
  } else {
    p <- gplot(results, ylim, xlab, ylab, title, xgap, legend, ref_line, theming)
  }

  p
}
