# Compute Group-Time Average Treatment Effects

`compute.att_gt2` does the (g,t) cell computation and sends it to
estimation, then does all the post-processing after estimation

## Usage

``` r
compute.att_gt2(dp2)
```

## Arguments

- dp2:

  A DIDparams object v2.0

## Value

a list with length equal to the number of groups times the number of
time periods; each element of the list contains an object that contains
group-time average treatment effect as well as which group it is for and
which time period it is for. It also exports the influence function
which is used externally to compute standard errors.
