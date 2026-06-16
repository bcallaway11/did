# Multiplier Bootstrap for Conditional Moment Test

A slightly modified multiplier bootstrap procedure for the pre-test of
the conditional parallel trends assumption

## Usage

``` r
test.mboot(inf.func, DIDparams, cores = 1)
```

## Arguments

- inf.func:

  an influence function

- DIDparams:

  DIDparams object

- cores:

  Unused; retained for backward compatibility. The multiplier bootstrap
  is computed with vectorized matrix operations in a single process, so
  this argument has no effect here. In `conditional_did_pretest`,
  `cores` parallelizes only Step 1 (computing the test statistic).

## Value

list

- bres:

  CvM test statistics for each bootstrap iteration

- crit.val:

  critical value for CvM test statistic
