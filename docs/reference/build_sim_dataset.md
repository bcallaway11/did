# build_sim_dataset

A function for building simulated data

## Usage

``` r
build_sim_dataset(sp_list, panel = TRUE)
```

## Arguments

- sp_list:

  A list of simulation parameters. See `reset.sim` to generate some
  default values for parameters

- panel:

  whether to construct panel data (the default) or repeated cross
  sections data

## Value

a data.frame with the following columns

- G observations group

- X value of covariate

- id observation's id

- cluster observation's cluster (by construction there is no
  within-cluster correlation)

- period time period for current observation

- Y outcome

- treat whether or not this unit is ever treated
