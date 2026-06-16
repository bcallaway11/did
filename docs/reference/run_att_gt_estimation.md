# Run ATT estimation for a given group-time pair

`run_att_gt_estimation` does the main work for computing multiperiod
group-time average treatment effects

## Usage

``` r
run_att_gt_estimation(g, t, dp2)
```

## Arguments

- g:

  group of interest (treated group at time t)

- t:

  time period

- dp2:

  A DIDparams object v2.0

## Value

a list with the gt cell and the results after performing estimation
