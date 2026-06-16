# sim

An internal function that builds simulated data, computes ATT(g,t)'s and
some aggregations. It is useful for testing the inference procedures in
the `did` function.

## Usage

``` r
sim(
  sp_list,
  ret = NULL,
  bstrap = TRUE,
  cband = TRUE,
  control_group = "nevertreated",
  xformla = ~X,
  est_method = "dr",
  clustervars = NULL,
  panel = TRUE
)
```

## Arguments

- sp_list:

  A list of simulation parameters. See `reset.sim` to generate some
  default values for parameters

- ret:

  which type of results to return. The options are `Wpval` (returns 1 if
  the p-value from a Wald test that all pre-treatment ATT(g,t)'s are
  equal is less than .05), `cband` (returns 1 if a uniform confidence
  band covers 0 for groups and times), `simple` (returns 1 if, using the
  simple treatment effect aggregation results in rejecting that this
  aggregated treatment effect parameter is equal to 0), `dynamic`
  (returns 1 if the uniform confidence band from the dynamic treatment
  effect aggregation covers 0 in all pre- and post-treatment periods).
  The default value is NULL, and in this case the function will just
  return the results from the call to `att_gt`.

- bstrap:

  whether or not to use the bootstrap to conduct inference (default is
  TRUE)

- cband:

  whether or not to compute uniform confidence bands in the call to
  `att_gt` (the default is TRUE)

- control_group:

  Whether to use the "nevertreated" comparison group (the default) or
  the "notyettreated" as the comparison group

- xformla:

  Formula for covariates in `att_gt` (default is `~X`)

- est_method:

  Which estimation method to use in `att_gt` (default is "dr")

- clustervars:

  Any additional variables which should be clustered on

- panel:

  whether to simulate panel data (the default) or otherwise repeated
  cross sections data

## Value

When `ret=NULL`, returns the results of the call to `att_gt`, otherwise
returns 1 if the specified test rejects or 0 if not.
