# Tidy an AGGTEobj into a data frame

Returns a tidy data frame of aggregated treatment effect estimates from
an [`aggte()`](https://bcallaway11.github.io/did/reference/aggte.md)
result.

## Usage

``` r
# S3 method for class 'AGGTEobj'
tidy(x, ...)
```

## Arguments

- x:

  a model of class AGGTEobj produced by the
  [`aggte()`](https://bcallaway11.github.io/did/reference/aggte.md)
  function

- ...:

  Additional arguments to tidying method.

## Value

A data frame whose columns depend on `type`:

- type:

  the aggregation type: `"simple"`, `"dynamic"`, `"group"`, or
  `"calendar"`

- term:

  label for each estimate

- estimate:

  point estimate

- std.error:

  standard error

- statistic:

  t-statistic (`estimate / std.error`)

- p.value:

  two-sided pointwise p-value (`2 * (1 - pnorm(abs(statistic)))`).
  Marginal per-estimate; does **not** account for multiple testing
  across event times or groups.

- conf.low, conf.high:

  simultaneous confidence band limits. When `bstrap=TRUE` and
  `cband=TRUE` these use the bootstrap uniform critical value
  (`crit.val.egt`); otherwise they equal the pointwise intervals. For
  `type="simple"` and the overall average row of `type="group"`, a
  single scalar is returned so simultaneous and pointwise coincide.

- point.conf.low, point.conf.high:

  pointwise confidence interval limits always using `qnorm(1 - alp/2)`.

## Details

The key distinction between `conf.low`/`conf.high` and
`point.conf.low`/`point.conf.high` is that the former accounts for
multiple testing across all estimates (simultaneous coverage), while the
latter provides marginal (per-estimate) coverage only. Use the
simultaneous bands when you want to make joint inferences across all
event times or groups.
