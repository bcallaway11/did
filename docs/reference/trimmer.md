# trimmer

A utility function to find observations that appear to violate support
conditions. This function is not called anywhere in the code, but it is
just useful for debugging some common issues that users run into.

## Usage

``` r
trimmer(
  g,
  tname,
  idname,
  gname,
  xformla,
  data,
  control_group = "notyettreated",
  threshold = 0.999
)
```

## Arguments

- g:

  is a particular group (below I pass in 2009)

- tname:

  The name of the column containing the time periods

- idname:

  The individual (cross-sectional unit) id name

- gname:

  The name of the variable in `data` that contains the first period when
  a particular observation is treated. This should be a positive number
  for all observations in treated groups. It defines which "group" a
  unit belongs to. It should be 0 for units in the untreated group.

- xformla:

  A formula for the covariates to include in the model. It should be of
  the form `~ X1 + X2`. Default is NULL which is equivalent to
  `xformla=~1`. This is used to create a matrix of covariates which is
  then passed to the 2x2 DID estimator chosen in `est_method`.

  For time-varying covariates: (1) With balanced panel data, in each 2x2
  comparison, the covariates are taken to be the value of the covariates
  in the earlier time period, and all of the underlying computations
  involve changes in Y as a function of those values of covariates. (2)
  With repeated cross sections data and unbalanced panel data, the
  covariates are taken from each time period and computations involve
  Y_post conditional on X_post minus Y_pre conditional on X_pre. A
  byproduct of this is that, with balanced panel data and in the
  presence of time-varying covariates, it is possible to get different
  numerical results according to whether or not
  `allow_unbalanced_panel=TRUE` or `FALSE`.

- data:

  The name of the data.frame that contains the data

- control_group:

  Which units to use as the control group. The default is "nevertreated"
  which sets the control group to be the group of units that never
  participate in the treatment. This group does not change across groups
  or time periods. The other option is to set `group="notyettreated"`.
  In this case, the control group is set to the group of units that have
  not yet participated in the treatment in that time period. This
  includes all never treated units, but it includes additional units
  that eventually participate in the treatment, but have not
  participated yet.

- threshold:

  the cutoff for which observations are flagged as likely violators of
  the support condition.

## Value

list of ids of observations that likely violate support conditions
