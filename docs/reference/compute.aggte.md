# Compute Aggregated Treatment Effect Parameters

Does the heavy lifting on computing aggregated group-time average
treatment effects

## Usage

``` r
compute.aggte(
  MP,
  type = "group",
  balance_e = NULL,
  min_e = -Inf,
  max_e = Inf,
  na.rm = FALSE,
  bstrap = NULL,
  biters = NULL,
  cband = NULL,
  alp = NULL,
  clustervars = NULL,
  call = NULL
)
```

## Arguments

- MP:

  an MP object (i.e., the results of the
  [`att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.md)
  method)

- type:

  Which type of aggregated treatment effect parameter to compute. One
  option is "simple" (this just computes a weighted average of all
  group-time average treatment effects with weights proportional to
  group size). Other options are "dynamic" (this computes average
  effects across different lengths of exposure to the treatment and is
  similar to an "event study"; here the overall effect averages the
  effect of the treatment across all positive lengths of exposure);
  "group" (this is the default option and computes average treatment
  effects across different groups; here the overall effect averages the
  effect across different groups); and "calendar" (this computes average
  treatment effects across different time periods; here the overall
  effect averages the effect across each time period).

- balance_e:

  If set (and if one computes dynamic effects), it balances the sample
  with respect to event time. For example, if `balance.e=2`, `aggte`
  will drop groups that are not exposed to treatment for at least three
  periods. (the initial period when `e=0` as well as the next two
  periods when `e=1` and the `e=2`). This ensures that the composition
  of groups does not change when event time changes.

- min_e:

  For event studies, this is the smallest event time to compute dynamic
  effects for. By default, `min_e = -Inf` so that effects at all lengths
  of exposure are computed.

- max_e:

  For event studies, this is the largest event time to compute dynamic
  effects for. By default, `max_e = Inf` so that effects at all lengths
  of exposure are computed.

- na.rm:

  Logical value if we are to remove missing Values from analyses.
  Defaults is FALSE.

- bstrap:

  Boolean for whether or not to compute standard errors using the
  multiplier bootstrap. Default is `TRUE` (in addition, cband is also by
  default `TRUE` indicating that uniform confidence bands will be
  returned). If `bstrap=FALSE`, analytical standard errors are reported;
  these are cluster-robust when `clustervars` is supplied.

- biters:

  The number of bootstrap iterations to use. The default is 1000, and
  this is only applicable if `bstrap=TRUE`.

- cband:

  Boolean for whether or not to compute a uniform confidence band that
  covers all of the group-time average treatment effects with fixed
  probability `1-alp`. In order to compute uniform confidence bands,
  `bstrap` must also be set to `TRUE`. The default is `TRUE`.

- alp:

  the significance level, default is 0.05

- clustervars:

  A vector of variables names to cluster on. At most, there can be two
  variables (otherwise will throw an error) and one of these must be the
  same as idname which allows for clustering at the individual level.
  Clustered standard errors are available with the multiplier bootstrap
  (`bstrap=TRUE`) or analytically (`bstrap=FALSE`).

- call:

  The function call to aggte

## Value

[`AGGTEobj`](https://bcallaway11.github.io/did/reference/AGGTEobj.md)
object
