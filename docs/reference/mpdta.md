# County Teen Employment Dataset

A dataset containing (the log of) teen employment in 500 counties in the
U.S. from 2003 to 2007. This is a subset of the dataset used in Callaway
and Sant'Anna (2021). See that paper for additional descriptions.

## Usage

``` r
mpdta
```

## Format

A data frame with 2500 rows and 6 variables:

- year:

  the year of the observation

- countyreal:

  a unique identifier for a particular county

- lpop:

  the log of 1000s of population for the county

- lemp:

  the log of teen employment in the county

- first.treat:

  the year that the state where the county is located raised its minimum
  wage, it is set equal to 0 for counties that have minimum wages equal
  to the federal minimum wage over the entire period.

- treat:

  whether or not a particular county is ever treated during the sample
  period. It is time-invariant within each county and equal to 1 exactly
  when first.treat is positive; it is not a year-specific treatment
  indicator.

## Source

Callaway and Sant'Anna (2021)
