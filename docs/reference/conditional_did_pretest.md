# Pre-Test of Conditional Parallel Trends Assumption

An integrated moments test for the conditional parallel trends
assumption holding in all pre-treatment time periods for all groups

## Usage

``` r
conditional_did_pretest(
  yname,
  tname,
  idname = NULL,
  gname,
  xformla = NULL,
  data,
  panel = TRUE,
  allow_unbalanced_panel = FALSE,
  control_group = c("nevertreated", "notyettreated"),
  weightsname = NULL,
  alp = 0.05,
  bstrap = TRUE,
  cband = TRUE,
  biters = 1000,
  clustervars = NULL,
  est_method = "ipw",
  print_details = FALSE,
  pl = FALSE,
  cores = 1
)
```

## Arguments

- yname:

  The name of the outcome variable

- tname:

  The name of the column containing the time periods

- idname:

  The individual (cross-sectional unit) id name

- gname:

  The name of the variable in `data` that contains the first period when
  a particular observation is treated. This should be a positive number
  for all observations in treated groups. It defines which "group" a
  unit belongs to. It should be 0 for units in the untreated group.

- xformla:

  A formula for the covariates to include in the model. It should be of
  the form `~ X1 + X2`. Default is NULL which is equivalent to
  `xformla=~1`. This is used to create a matrix of covariates which is
  then passed to the 2x2 DID estimator chosen in `est_method`.

  For time-varying covariates: (1) With balanced panel data, in each 2x2
  comparison, the covariates are taken to be the value of the covariates
  in the earlier time period, and all of the underlying computations
  involve changes in Y as a function of those values of covariates. (2)
  With repeated cross sections data and unbalanced panel data, the
  covariates are taken from each time period and computations involve
  Y_post conditional on X_post minus Y_pre conditional on X_pre. A
  byproduct of this is that, with balanced panel data and in the
  presence of time-varying covariates, it is possible to get different
  numerical results according to whether or not
  `allow_unbalanced_panel=TRUE` or `FALSE`.

- data:

  The name of the data.frame that contains the data

- panel:

  Whether or not the data is a panel dataset. The panel dataset should
  be provided in long format – that is, where each row corresponds to a
  unit observed at a particular point in time. The default is TRUE. When
  using a panel dataset, the variable `idname` must be set. When
  `panel=FALSE`, the data is treated as repeated cross sections.

- allow_unbalanced_panel:

  Whether or not function should "balance" the panel with respect to
  time and id. The default value is `FALSE` which means that
  [`att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.md)
  will drop all units where data is not observed in all periods. The
  advantage of this is that the computations are faster (sometimes
  substantially).

- control_group:

  Which units to use as the control group. The default is "nevertreated"
  which sets the control group to be the group of units that never
  participate in the treatment. This group does not change across groups
  or time periods. The other option is to set `group="notyettreated"`.
  In this case, the control group is set to the group of units that have
  not yet participated in the treatment in that time period. This
  includes all never treated units, but it includes additional units
  that eventually participate in the treatment, but have not
  participated yet.

- weightsname:

  The name of the column containing the sampling weights. If not set,
  all observations have same weight.

- alp:

  the significance level, default is 0.05

- bstrap:

  Not used by the pre-test. Critical values and p-values for the Cramer
  von Mises and Kolmogorov-Smirnov test statistics are always computed
  with the multiplier bootstrap (using `biters` iterations), regardless
  of this argument, and no standard errors are reported.

- cband:

  Boolean for whether or not to compute a uniform confidence band that
  covers all of the group-time average treatment effects with fixed
  probability `1-alp`. In order to compute uniform confidence bands,
  `bstrap` must also be set to `TRUE`. The default is `TRUE`.

- biters:

  The number of multiplier bootstrap iterations used to simulate the
  limiting distribution of the test statistics. The default is 1000.

- clustervars:

  A vector of variables names to cluster on (the multiplier bootstrap
  then draws cluster-level multipliers). At most, there can be two
  variables (otherwise will throw an error) and one of these must be the
  same as idname which allows for clustering at the individual level.

- est_method:

  the method to compute group-time average treatment effects. The
  default for `conditional_did_pretest` is "ipw" for inverse probability
  weighting. Other built-in methods include "dr" which uses the doubly
  robust approach in the `DRDID` package and "reg" for first step
  regression estimators.

- print_details:

  Whether or not to show details/progress of computations. Default is
  `FALSE`.

- pl:

  Whether or not to use parallel processing

- cores:

  The number of cores to use for parallel processing. This parallelizes
  Step 1 (computing the test statistic); Step 2's multiplier bootstrap
  is vectorized and runs in a single process.

## Value

an [`MP.TEST`](https://bcallaway11.github.io/did/reference/MP.TEST.md)
object

## References

Callaway, Brantly and Sant'Anna, Pedro H. C. "Difference-in-Differences
with Multiple Time Periods and an Application on the Minimum Wage and
Employment." Working Paper <https://arxiv.org/abs/1803.09015v2> (2018).

## Examples

``` r
if (FALSE) { # \dontrun{
data(mpdta)
pre.test <- conditional_did_pretest(yname="lemp",
                                    tname="year",
                                    idname="countyreal",
                                    gname="first.treat",
                                    xformla=~lpop,
                                    data=mpdta)
summary(pre.test)
} # }
```
