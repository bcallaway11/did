
#' @title gplot
#'
#' @description does the heavy lifting for making a plot of an group-time
#'  average treatment effect
#'
#' @inheritParams ggdid
#'
#' @return a `ggplot2` object
#'
#' @keywords internal
#'
#' @export
gplot <- function(ssresults, ylim=NULL, xlab=NULL, ylab=NULL, title="Group", xgap=1,
                  legend=TRUE, ref_line = 0, theming = TRUE) {
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
    scale_color_manual(drop=FALSE, values=c("#e87d72","#56bcc2"), breaks = c(0, 1), labels = c('Pre','Post')) +
    labs(x = xlab, y = ylab, title = title, color = NULL)

  if (!is.null(ref_line)) {
    p <- p + geom_hline(aes(yintercept = ref_line), linetype = 'dashed')
  }
  if (theming) {
    p <- p + ggpubr::theme_pubr() +
      theme(plot.title = element_text(color="darkgray", face="bold", size=12),
            axis.title = element_text(color="darkgray", face="bold", size=12),
            strip.background = element_rect(fill = 'white', color = 'white'),
            strip.text = element_text(color = 'darkgray', face = 'bold', size = 12, hjust = 0),
            legend.position = 'bottom')
  }
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
#' @return a `ggplot2` object
#'
#' @keywords internal
#'
#' @export
splot <- function(ssresults, ylim=NULL, xlab=NULL, ylab=NULL, title="Group",
                  legend=TRUE, ref_line = 0, theming = TRUE) {

  # names of variables are "weird" for this function because this code builds
  # on the same infrastructure as for plotting group-time average treatment
  # effects and aggregations using event time or calendar time

  # the group variable is saved in year

  if (is.null(xlab)) {
    xlab <- 'ATT'
  }
  if (is.null(ylab)) {
    ylab <- 'Group'
  }

  p <- ggplot(ssresults,
              aes(y=as.factor(year), x=att, xmin=(att-c*att.se),
                  xmax=(att+c*att.se))) +
    geom_point(aes(colour=post), size=1.5) +
    #geom_ribbon(aes(x=as.numeric(year)), alpha=0.2) +
    geom_errorbarh(aes(colour=post), height=0.1)  +
    scale_y_discrete(breaks=as.factor(ssresults$year)) +
    #scale_x_discrete(breaks=dabreaks, labels=as.character(dabreaks)) +
    scale_x_continuous(limits=ylim) +
    scale_color_manual(drop=FALSE, values=c("#e87d72","#56bcc2"), breaks = c(0, 1), labels = c('Pre','Post')) +
    labs(x = xlab, y = ylab, title = title)

  if (!is.null(ref_line)) {
    p <- p + geom_vline(aes(xintercept = ref_line), linetype = 'dashed')
  }
  if (theming) {
    p <- p + ggpubr::theme_pubr() +
      theme(plot.title = element_text(color="darkgray", face="bold", size=12),
            axis.title = element_text(color="darkgray", face="bold", size=12),
            legend.position = 'none')
  }

  if (!legend) {
    p  <- p + ggpubr::rremove("legend")
  }

  p
}
