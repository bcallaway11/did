# reset.sim

a function to create a "reasonable" set of parameters to create
simulated panel data that obeys a parallel trends assumption. In
particular, it provides parameters where the the effect of participating
in the treatment is equal to one in all post-treatment time periods.

After calling this function, the user can change particular values of
the parameters in order to generate dynamics, heterogeneous effects
across groups, etc.

## Usage

``` r
reset.sim(time.periods = 4, n = 5000, ipw = TRUE, reg = TRUE)
```

## Arguments

- time.periods:

  The number of time periods to include

- n:

  The total number of observations

- ipw:

  If TRUE, sets parameters so that DGP is compatible with recovering
  ATT(g,t)'s using IPW (i.e., where logit that just includes a linear
  term in X works). If FALSE, sets parameters that will be incompatible
  with IPW. Either way, these parameters can be specified by the user if
  so desired.

- reg:

  If TRUE, sets parameters so that DGP is compatible with recovering
  ATT(g,t)'s using regressions on untreated untreated potential
  outcomes. If FALSE, sets parameters that will be incompatible with
  using regressions (i.e., regressions that include only linear term in
  X). Either way, these parameters can be specified by the user if so
  desired.

## Value

list of simulation parameters
