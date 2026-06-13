# Compute Group-Time Average Treatment Effects

`compute.att_gt` does the main work for computing multiperiod group-time
average treatment effects

## Usage

``` r
compute.att_gt(dp)
```

## Arguments

- dp:

  A DIDparams object

## Value

a list with length equal to the number of groups times the number of
time periods; each element of the list contains an object that contains
group-time average treatment effect as well as which group it is for and
which time period it is for. It also exports the influence function
which is used externally to compute standard errors.
