# AGGTEobj

Objects of this class hold results on aggregated group-time average
treatment effects

An object for holding aggregated treatment effect parameters.

## Usage

``` r
AGGTEobj(
  overall.att = NULL,
  overall.se = NULL,
  type = "simple",
  egt = NULL,
  att.egt = NULL,
  se.egt = NULL,
  crit.val.egt = NULL,
  inf.function = NULL,
  min_e = NULL,
  max_e = NULL,
  balance_e = NULL,
  call = NULL,
  DIDparams = NULL
)
```

## Arguments

- overall.att:

  The estimated overall ATT

- overall.se:

  Standard error for overall ATT

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

- egt:

  Holds the length of exposure (for dynamic effects), the group (for
  selective treatment timing), or the time period (for calendar time
  effects)

- att.egt:

  The ATT specific to egt

- se.egt:

  The standard error specific to egt

- crit.val.egt:

  A critical value for computing uniform confidence bands for dynamic
  effects, selective treatment timing, or time period effects.

- inf.function:

  The influence function of the chosen aggregated parameters

- min_e:

  For event studies, this is the smallest event time to compute dynamic
  effects for. By default, `min_e = -Inf` so that effects at all lengths
  of exposure are computed.

- max_e:

  For event studies, this is the largest event time to compute dynamic
  effects for. By default, `max_e = Inf` so that effects at all lengths
  of exposure are computed.

- balance_e:

  If set (and if one computes dynamic effects), it balances the sample
  with respect to event time. For example, if `balance.e=2`, `aggte`
  will drop groups that are not exposed to treatment for at least three
  periods. (the initial period when `e=0` as well as the next two
  periods when `e=1` and the `e=2`). This ensures that the composition
  of groups does not change when event time changes.

- call:

  The function call to aggte

- DIDparams:

  A DIDparams object

## Value

an AGGTEobj
