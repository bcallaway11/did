# gplot

does the heavy lifting for making a plot of an group-time average
treatment effect

## Usage

``` r
gplot(
  ssresults,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  title = "Group",
  xgap = 1,
  legend = TRUE,
  ref_line = 0,
  theming = TRUE
)
```

## Value

a `ggplot2` object
