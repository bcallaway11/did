# MP.TEST

An object that holds results from computing pre-test of the conditional
parallel trends assumption

## Usage

``` r
MP.TEST(
  CvM = NULL,
  CvMb = NULL,
  CvMcval = NULL,
  CvMpval = NULL,
  KS = NULL,
  KSb = NULL,
  KScval = NULL,
  KSpval = NULL,
  clustervars = NULL,
  xformla = NULL
)
```

## Arguments

- CvM:

  Cramer von Mises test statistic

- CvMb:

  a vector of bootstrapped Cramer von Mises test statistics

- CvMcval:

  CvM critical value

- CvMpval:

  p-value for CvM test

- KS:

  Kolmogorov-Smirnov test statistic

- KSb:

  a vector of bootstrapped KS test statistics

- KScval:

  KS critical value

- KSpval:

  p-value for KS test

- clustervars:

  vector of which variables were clustered on for the test

- xformla:

  formla for the X variables used in the test
