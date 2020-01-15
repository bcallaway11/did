
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
gplot <- function(ssresults, ylim=NULL, xlab=NULL, ylab=NULL, title="Group", xgap=1) {
  dabreaks <- ssresults$year[seq(1, length(ssresults$year), xgap)]

  c.point <-  stats::qnorm(1 - ssresults$alp/2)

  p <- ggplot(ssresults,
              aes(x=year, y=att, ymin=(att-c*ssresults$att.se),
                  ymax=(att+c*ssresults$att.se), post=post)) +

    geom_point(aes(colour=post), size=1.5) +
    geom_errorbar(aes(colour=post), width=0.1) +
    #geom_ribbon(aes(ymin= (att-c.point*att.se), ymax=  (att+c.point*att.se), alpha=0.35))+
    #geom_ribbon(aes(ymin=  (att-c*att.se), ymax =  (att+c*att.se), alpha=0.25))+

    scale_y_continuous(limits=ylim) +
    scale_x_discrete(breaks=dabreaks, labels=as.character(dabreaks)) +
    scale_colour_hue(drop=FALSE) +
    ylab("") +
    xlab("") +
    ggtitle(paste(title, unique(ssresults$group))) +
    theme_bw() +
    theme(plot.title = element_text(color="darkgray", face="bold", size=8)) +
    theme(axis.title = element_text(color="darkgray", face="bold", size=8))
  p
}
