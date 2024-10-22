#' Difference in Differences
#'
#' `did` implements Difference in Differences with multiple periods and variation in
#' treatment timing.
#'
#' @keywords internal
"_PACKAGE"

#' @import stats
#' @import utils
#' @import ggplot2
#' @import ggpubr
#' @import BMisc
#' @import data.table
#' @importFrom tidyr gather
#' @importFrom methods is
NULL
utils::globalVariables(c('.','.G','.y', 'asif_never_treated', 'treated_first_period', 'count', 'constant', '.rowid',
                         'V1', 'control_group', 'cohort', 'cohort_size', 'period', 'period_size', 'y1', 'y0',
                         'i.weights', 'y', 'cluster', 'id', '..cols_to_keep', '..g'))
