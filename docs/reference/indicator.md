# indicator

indicator weighting function

## Usage

``` r
indicator(X, u)
```

## Arguments

- X:

  matrix of X's from the data

- u:

  a particular value to compare X's to

## Value

numeric vector

## Examples

``` r
data(mpdta)
dta <- subset(mpdta, year==2007)
X <- model.matrix(~lpop, data=dta)
X <- indicator(X, X[1,])
```
