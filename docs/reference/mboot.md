# Multiplier Bootstrap

A function to take an influence function and use the multiplier
bootstrap to compute standard errors and critical values for uniform
confidence bands.

## Usage

``` r
mboot(inf.func, DIDparams, pl = FALSE, cores = 1, return_V = TRUE)
```

## Arguments

- inf.func:

  an influence function

- DIDparams:

  DIDparams object

- pl:

  whether or not to use parallel processing in the multiplier bootstrap,
  default=FALSE

- cores:

  the number of cores to use with parallel processing, default=1

- return_V:

  whether to compute and return the bootstrap variance matrix `V`.
  Default is `TRUE`. Internal callers that only consume `bres`, `se`, or
  `crit.val` set this to `FALSE` to skip the computation (it is the only
  O(biters x k^2) step in the function).

## Value

list with elements

- bres:

  results from each bootstrap iteration

- V:

  variance matrix (`NULL` when `return_V = FALSE`)

- se:

  standard errors

- crit.val:

  a critical value for computing uniform confidence bands
