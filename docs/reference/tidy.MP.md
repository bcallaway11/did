# Tidy an MP object into a data frame

Returns a tidy data frame of group-time average treatment effect
estimates from an
[`att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.md)
result.

## Usage

``` r
# S3 method for class 'MP'
tidy(x, ...)
```

## Arguments

- x:

  a model of class MP produced by the
  [`att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.md)
  function

- ...:

  Additional arguments to tidying method.

## Value

A data frame with one row per ATT(g,t) estimate and columns:

- term:

  ATT(g,t) label

- group:

  the treatment cohort g

- time:

  the time period t

- estimate:

  the ATT(g,t) point estimate

- std.error:

  standard error

- statistic:

  t-statistic (`estimate / std.error`)

- p.value:

  two-sided pointwise p-value (`2 * (1 - pnorm(abs(statistic)))`).
  Marginal per-estimate; does **not** account for multiple testing
  across ATT(g,t) cells.

- conf.low, conf.high:

  simultaneous confidence band limits, using the bootstrap uniform
  critical value when `bstrap=TRUE` and `cband=TRUE`, otherwise
  pointwise

- point.conf.low, point.conf.high:

  pointwise confidence interval limits using `qnorm(1 - alp/2)`, always
