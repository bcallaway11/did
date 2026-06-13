# MP

Multi-period objects that hold results for group-time average treatment
effects

## Usage

``` r
MP(
  group,
  t,
  att,
  V_analytical,
  se,
  c,
  inffunc,
  n = NULL,
  W = NULL,
  Wpval = NULL,
  aggte = NULL,
  alp = 0.05,
  DIDparams = NULL
)
```

## Arguments

- group:

  which group (defined by period first treated) an group-time average
  treatment effect is for

- t:

  which time period a group-time average treatment effect is for

- att:

  the group-average treatment effect for group `group` and time period
  `t`

- V_analytical:

  Analytical estimator for the asymptotic variance-covariance matrix for
  group-time average treatment effects

- se:

  standard errors for group-time average treatment effects. If bootstrap
  is set to TRUE, this provides bootstrap-based se.

- c:

  simultaneous critical value if one is obtaining simultaneous
  confidence bands. Otherwise it reports the critical value based on
  pointwise normal approximation.

- inffunc:

  the influence function for estimating group-time average treatment
  effects: one column per ATT(g,t) and one row per cross-sectional unit
  (one row per observation with repeated cross sections). The rownames
  hold the unit ids and are the authoritative link between rows and
  units – the row order differs across `faster_mode = TRUE` (internal
  (period, cohort, id) ordering) and `faster_mode = FALSE` (id-sorted),
  so align rows by rowname, never by position.

- n:

  the number of unique cross-sectional units (unique values of idname)

- W:

  the Wald statistic for pre-testing the common trends assumption

- Wpval:

  the p-value of the Wald statistic for pre-testing the common trends
  assumption

- aggte:

  an aggregate treatment effects object

- alp:

  the significance level, default is 0.05

- DIDparams:

  a
  [`DIDparams`](https://bcallaway11.github.io/did/reference/DIDparams.md)
  object. A way to optionally return the parameters of the call to
  [`att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.md) or
  [`conditional_did_pretest()`](https://bcallaway11.github.io/did/reference/conditional_did_pretest.md).

## Value

MP object
