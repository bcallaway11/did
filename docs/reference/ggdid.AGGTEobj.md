# Plot `AGGTEobj` objects

A function to plot `AGGTEobj` objects

## Usage

``` r
# S3 method for class 'AGGTEobj'
ggdid(
  object,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  title = "",
  xgap = 1,
  legend = TRUE,
  ref_line = 0,
  theming = TRUE,
  ...
)
```

## Arguments

- object:

  either a `MP` object or `AGGTEobj` object. See
  [`help(ggdid.MP)`](https://bcallaway11.github.io/did/reference/ggdid.MP.md)
  and `help(ggdid.AGGTEobj)`.

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

- legend:

  Whether or not to include a legend (which will indicate color of pre-
  and post-treatment estimates). Default is `TRUE`.

- ref_line:

  A reference line at this value, usually to compare confidence
  intervals to 0. Set to NULL to omit.

- theming:

  Set to FALSE to skip all theming so you can do it yourself.

- ...:

  other arguments
