
#' @title gplot
#'
#' @description does the heavy lifting for making a plot of an group-time
#'  average treatment effect
#'
#' @inheritParams ggdid
#'
#' @return a \code{ggplot2} object
#'
#' @keywords internal
#'
#' @export
gplot <- function(ssresults, ylim=NULL, xlab=NULL, ylab=NULL, title="Group", xgap=1,
                  legend=TRUE) {
  dabreaks <- ssresults$year[seq(1, length(ssresults$year), xgap)]

  c.point <-  qnorm(1 - ssresults$alp/2)

  p <- ggplot(ssresults,
              aes(x=as.numeric(year), y=att, ymin=(att-c*att.se),
                  ymax=(att+c*att.se))) +
    
    geom_point(aes(colour=post), size=1.5) +
    #geom_ribbon(aes(x=as.numeric(year)), alpha=0.2) + 
    geom_errorbar(aes(colour=post), width=0.1) +   
    scale_y_continuous(limits=ylim) +
    #scale_x_discrete(breaks=dabreaks, labels=as.character(dabreaks)) +
    scale_x_continuous(breaks=as.numeric(dabreaks), labels=as.character(dabreaks)) + 
    scale_colour_hue(drop=FALSE) +
    ylab("") +
    xlab("") +
    ggtitle(paste(title, unique(ssresults$group))) +
    theme_bw() +
    theme(plot.title = element_text(color="darkgray", face="bold", size=12)) +
    theme(axis.title = element_text(color="darkgray", face="bold", size=12))

  if (!legend) {
    p  <- p + ggpubr::rremove("legend")
  }
  
  p
}


#' @title splot
#'
#' @description alternative plotting function for the case with selective treatment
#'  timing
#'
#' @inheritParams ggdid
#' @inheritParams gplot
#'
#' @return a \code{ggplot2} object
#'
#' @keywords internal
#'
#' @export
splot <- function(ssresults, ylim=NULL, xlab=NULL, ylab=NULL, title="Group",
                  legend=TRUE) {

  # names of variables are "weird" for this function because this code builds
  # on the same infrastructure as for plotting group-time average treatment
  # effects and aggregations using event time or calendar time

  # the group variable is saved in year
  
  p <- ggplot(ssresults,
              aes(y=as.factor(year), x=att, xmin=(att-c*att.se),
                  xmax=(att+c*att.se))) +    
    geom_point(aes(colour=post), size=1.5) +
    #geom_ribbon(aes(x=as.numeric(year)), alpha=0.2) + 
    geom_errorbarh(aes(colour=post), height=0.1)  +   
    scale_y_discrete(breaks=as.factor(ssresults$year)) +
    #scale_x_discrete(breaks=dabreaks, labels=as.character(dabreaks)) +
    scale_x_continuous(limits=ylim) + 
    scale_colour_hue(drop=FALSE) +
    ylab("group") +
    xlab("att") +
    ggtitle(paste(title, unique(ssresults$group))) +
    theme_bw() +
    theme(plot.title = element_text(color="darkgray", face="bold", size=8)) +
    theme(axis.title = element_text(color="darkgray", face="bold", size=8)) +
    theme(legend.position="none")

  
  if (!legend) {
    p  <- p + ggpubr::rremove("legend")
  }
  
  p
}
