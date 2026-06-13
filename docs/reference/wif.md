# Compute extra term in influence function due to estimating weights

A function to compute the extra term that shows up in the influence
function for aggregated treatment effect parameters due to estimating
the weights

## Usage

``` r
wif(keepers, pg, weights.ind, G, group)
```

## Arguments

- keepers:

  a vector of indices for which group-time average treatment effects are
  used to compute a particular aggregated parameter

- pg:

  a vector with same length as total number of group-time average
  treatment effects that contains the probability of being in particular
  group

- weights.ind:

  additional sampling weights (nx1)

- G:

  vector containing which group a unit belongs to (nx1)

- group:

  vector of groups

## Value

nxk influence function matrix
