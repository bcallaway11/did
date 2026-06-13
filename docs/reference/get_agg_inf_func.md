# Get an influence function for particular aggregate parameters

Get an influence function for particular aggregate parameters

This is a generic internal function for combining influence functions
across ATT(g,t)'s to return an influence function for various aggregated
treatment effect parameters.

## Usage

``` r
get_agg_inf_func(att, inffunc1, whichones, weights.agg, wif = NULL)
```

## Arguments

- att:

  vector of group-time average treatment effects

- inffunc1:

  influence function for all group-time average treatment effects
  (matrix)

- whichones:

  which elements of att will be used to compute the aggregated treatment
  effect parameter

- weights.agg:

  the weights to apply to each element of att(whichones); should have
  the same dimension as att(whichones)

- wif:

  extra influence function term coming from estimating the weights;
  should be n x k matrix where k is dimension of whichones

## Value

nx1 influence function
