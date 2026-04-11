# edid-imports.R
# @importFrom declarations for symbols used in edid-*.R files that are not
# already declared in imports.R.
#
# Note: stats::pnorm, stats::qnorm, stats::quantile, stats::sd, stats::setNames
# are already declared in imports.R. We add only symbols NOT already present.
# data.table is @imported (not just @importFrom) in imports.R, so dcast is
# available without an additional entry, but we add explicit entries for clarity.

#' @importFrom stats sd quantile as.formula
NULL
