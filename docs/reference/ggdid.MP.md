# Plot `MP` objects using `ggplot2`

A function to plot `MP` objects

## Usage

``` r
# S3 method for class 'MP'
ggdid(
  object,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  title = "Group",
  xgap = 1,
  ncol = 1,
  legend = TRUE,
  group = NULL,
  ref_line = 0,
  theming = TRUE,
  grtitle = "Group",
  ...
)
```

## Arguments

- object:

  either a `MP` object or `AGGTEobj` object. See `help(ggdid.MP)` and
  [`help(ggdid.AGGTEobj)`](https://bcallaway11.github.io/did/reference/ggdid.AGGTEobj.md).

- ylim:

  optional y limits for the plot; setting here makes the y limits the
  same across different plots

- xlab:

  optional x-axis label

- ylab:

  optional y-axis label

- title:

  optional plot title

- xgap:

  optional gap between the labels on the x-axis. For example, `xgap=3`
  indicates that the labels should show up for every third value on the
  x-axis. The default is 1.

- ncol:

  The number of columns to include in the resulting plot. The default is
  1.

- legend:

  Whether or not to include a legend (which will indicate color of pre-
  and post-treatment estimates). Default is `TRUE`.

- group:

  Vector for which groups to include in the plots of ATT(g,t). Default
  is NULL, and, in this case, plots for all groups will be included
  (`ggdid.MP` only).

- ref_line:

  A reference line at this value, usually to compare confidence
  intervals to 0. Set to NULL to omit.

- theming:

  Set to FALSE to skip all theming so you can do it yourself.

- grtitle:

  Title to append before each group name (`ggdid.MP` only).

- ...:

  other arguments
