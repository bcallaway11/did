
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Difference-in-Differences <img src="man/figures/logo.png" align="right" alt="" width="155" />

<!-- badges: start -->

[![](http://cranlogs.r-pkg.org/badges/grand-total/did?color=orange)](https://cran.r-project.org/package=did)

[![](http://cranlogs.r-pkg.org/badges/last-month/did?color=green)](https://cran.r-project.org/package=did)

[![](https://www.r-pkg.org/badges/version/did?color=blue)](https://cran.r-project.org/package=did)

[![](https://img.shields.io/badge/devel%20version-2.0.0-blue.svg)](https://github.com/bcallaway11/did)

[![](https://img.shields.io/github/languages/code-size/bcallaway11/did.svg)](https://github.com/bcallaway11/did)

[![CRAN
checks](https://cranchecks.info/badges/summary/did)](https://cran.r-project.org/web/checks/check_results_did.html)

<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

The **did** package contains tools for computing average treatment
effect parameters in a Difference-in-Differences setup allowing for

  - More than two time periods

  - Variation in treatment timing (i.e., units can become treated at
    different points in time)

  - Treatment effect heterogeneity (i.e, the effect of participating in
    the treatment can vary across units and exhibit potentially complex
    dynamics, selection into treatment, or time effects)

  - The parallel trends assumption holds only after conditioning on
    covariates

The main parameters are **group-time average treatment effects**. These
are the average treatment effect for a particular group (group is
defined by treatment timing) in a particular time period. These
parameters are a natural generalization of the average treatment effect
on the treated (ATT) which is identified in the textbook case with two
periods and two groups to the case with multiple periods.

Group-time average treatment effects are also natural building blocks
for more aggregated treatment effect parameters such as overall
treatment effects or event-study-type estimands.

## Getting Started

There has been some recent work on DiD with multiple time periods. The
**did** package implements the framework put forward in

  - [Callaway, Brantly, and Pedro H.C. Sant'Anna.
    Difference-in-differences with multiple time periods. Forthcoming at
    the Journal of Econometrics
    (2020).](https://arxiv.org/abs/1803.09015)

**Higher level discussions of issues are available in**

  - [Our approach to DiD with multiple time
    periods](https://bcallaway11.github.io/did/articles/multi-period-did.html)

## Installation

You can install **did** from CRAN with:

``` r
install.packages("did")
```

or get the latest version from github with:

``` r
# install.packages("devtools")
devtools::install_github("bcallaway11/did")
```

## A short example

The following is a simplified example of the effect of states increasing
their minimum wages on county-level teen employment rates which comes
from [Callaway and Sant’Anna (2020)](https://arxiv.org/abs/1803.09015).

  - [More detailed examples are also
    available](https://bcallaway11.github.io/did/articles/did-basics.html)

A subset of the data is available in the package and can be loaded by

``` r
  library(did)
  data(mpdta)
```

The dataset contains 500 observations of county-level teen employment
rates from 2003-2007. Some states are first treated in 2004, some in
2006, and some in 2007 (see the paper for more details). The important
variables in the dataset are

  - **lemp** This is the log of county-level teen employment. It is the
    outcome variable

  - **first.treat** This is the period when a state first increases its
    minimum wage. It can be 2004, 2006, or 2007. It is the variable that
    defines *group* in this application

  - **year** This is the year and is the *time* variable

  - **countyreal** This is an id number for each county and provides the
    individual identifier in this panel data context

To estimate group-time average treatment effects, use the **att\_gt**
function

``` r
out <- att_gt(yname = "lemp",
              gname = "first.treat",
              idname = "countyreal",
              tname = "year",
              xformla = ~1,
              data = mpdta,
              est_method = "reg"
              )
```

**att\_gt** returns a class **MP** object. This has a lot of
information, but most importantly is has estimates of the group-time
average treatment effects and their standard errors. To see these, we
can call the **summary** function

``` r
summary(out)
#> 
#> Call:
#> att_gt(yname = "lemp", tname = "year", idname = "countyreal", 
#>     gname = "first.treat", xformla = ~1, data = mpdta, est_method = "reg")
#> 
#> Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  "Difference-in-Differences with Multiple Time Periods." Forthcoming at the Journal of Econometrics <https://arxiv.org/abs/1803.09015>, 2020. 
#> 
#> Group-Time Average Treatment Effects:
#>  Group Time ATT(g,t) Std. Error [95% Simult.  Conf. Band]  
#>   2004 2004  -0.0105     0.0244       -0.0742      0.0532  
#>   2004 2005  -0.0704     0.0338       -0.1586      0.0177  
#>   2004 2006  -0.1373     0.0384       -0.2376     -0.0370 *
#>   2004 2007  -0.1008     0.0355       -0.1934     -0.0083 *
#>   2006 2004   0.0065     0.0252       -0.0593      0.0724  
#>   2006 2005  -0.0028     0.0193       -0.0531      0.0476  
#>   2006 2006  -0.0046     0.0179       -0.0512      0.0420  
#>   2006 2007  -0.0412     0.0211       -0.0964      0.0139  
#>   2007 2004   0.0305     0.0154       -0.0097      0.0707  
#>   2007 2005  -0.0027     0.0163       -0.0453      0.0399  
#>   2007 2006  -0.0311     0.0200       -0.0833      0.0211  
#>   2007 2007  -0.0261     0.0174       -0.0715      0.0194  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 
#> P-value for pre-test of parallel trends assumption:  0.16812
#> Control Group:  Never Treated,  Anticipation Periods:  0
#> Estimation Method:  Outcome Regression
```

This provides estimates of group-time average treatment effects for all
groups in all time periods. Group-time average treatment effects are
identified when `t >= g` (these are post-treatment time periods for each
group), but **summary** reports them even in periods when `t < g` –
these can be used to pre-test for the parallel trends assumption. The
`P-value for pre-test of parallel trends assumption` is for a Wald
pre-test of the parallel trends assumption. Here the parallel trends
assumption would not be rejected at conventional significance levels.

It is often also convenient to plot the group-time average treatment
effects. This can be done using the **ggdid** command:

``` r
ggdid(out, ylim = c(-.25,.1))
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="90%" style="display: block; margin: auto;" />

The red dots in the plot are pre-treatment group-time average treatment
effects . Here they are provided with 95% simultaneous confidence
intervals. These are the estimates that can be interpreted as a pre-test
of the parallel trends assumption. The blue dots are post-treatment
group-time average treatment effects. Under the parallel trends
assumption, these can be interpreted as policy effects – here the effect
of the minimum wage on county-level teen employment due to increasing
the minimum wage.

**Event Studies**

Although in the current example it is pretty easy to directly interpret
the group-time average treatment effects, there are many cases where it
is convenient to aggregate the group-time average treatment effects into
a small number of parameters. A main type of aggregation is into an
*event study* plot.

To make an event study plot in the **did** package, one can use the
**aggte** function

``` r
es <- aggte(out, type = "dynamic")
```

Just like for group-time average treatment effects, these can be
summarized and plotted. First, the summary

``` r
summary(es)
#> 
#> Call:
#> aggte(MP = out, type = "dynamic")
#> 
#> Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  "Difference-in-Differences with Multiple Time Periods." Forthcoming at the Journal of Econometrics <https://arxiv.org/abs/1803.09015>, 2020. 
#> 
#> 
#> Overall ATT:  
#>      ATT Std. Error     [95%  Conf. Int.]  
#>  -0.0772     0.0205   -0.1174      -0.037 *
#> 
#> 
#> Dynamic Effects:
#>  event time     ATT Std. Error [95% Simult.  Conf. Band]  
#>          -3  0.0305     0.0157       -0.0079      0.0689  
#>          -2 -0.0006     0.0138       -0.0345      0.0334  
#>          -1 -0.0245     0.0137       -0.0580      0.0091  
#>           0 -0.0199     0.0130       -0.0517      0.0118  
#>           1 -0.0510     0.0177       -0.0942     -0.0077 *
#>           2 -0.1373     0.0385       -0.2317     -0.0428 *
#>           3 -0.1008     0.0376       -0.1930     -0.0086 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 
#> Control Group:  Never Treated,  Anticipation Periods:  0
#> Estimation Method:  Outcome Regression
```

The column `event time` is for each group relative to when they first
participate in the treatment. To give some examples, `event time=0`
corresponds to the *on impact* effect, and `event time=-1` is the
*effect* in the period before a unit becomes treated (checking that this
is equal to 0 is potentially useful as a pre-test).

To plot the event study, use **ggdid**

``` r
ggdid(es)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="90%" style="display: block; margin: auto;" />

The figure here is very similar to the group-time average treatment
effects. Red dots are pre-treatment periods, blue dots are
post-treatment periods. The difference is that the x-axis is in event
time.

**Overall Effect of Participating in the Treatment**

The event study above reported an overall effect of participating in the
treatment. This was computed by averaging the average effects computed
at each length of exposure.

In many cases, a more general purpose overall treatment effect parameter
is give by computing the average treatment effect for each group, and
then averaging across groups. This sort of procedure provides an average
treatment effect parameter with a very similar interpretation to the
Average Treatment Effect on the Treated (ATT) in the two period and two
group case.

To compute this overall average treatment effect parameter, use

``` r
group_effects <- aggte(out, type = "group")
summary(group_effects)
#> 
#> Call:
#> aggte(MP = out, type = "group")
#> 
#> Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  "Difference-in-Differences with Multiple Time Periods." Forthcoming at the Journal of Econometrics <https://arxiv.org/abs/1803.09015>, 2020. 
#> 
#> 
#> Overall ATT:  
#>     ATT Std. Error     [95%  Conf. Int.]  
#>  -0.031     0.0119   -0.0543     -0.0078 *
#> 
#> 
#> Group Effects:
#>  group     ATT Std. Error [95% Simult.  Conf. Band]  
#>   2004 -0.0797     0.0280       -0.1442     -0.0153 *
#>   2006 -0.0229     0.0161       -0.0599      0.0141  
#>   2007 -0.0261     0.0180       -0.0676      0.0154  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 
#> Control Group:  Never Treated,  Anticipation Periods:  0
#> Estimation Method:  Outcome Regression
```

Of particular interest is the `Overall ATT` in the results. Here, we
estimate that increasing the minimum wage decreased teen employment by
3.1% and the effect is marginally statistically significant.

## Additional Resources

We have provided several more vignettes that may be helpful for using
the **did** package

  - [Getting Started with the did
    Package](https://bcallaway11.github.io/did/articles/did-basics.html)

  - [Introduction to DiD with Multiple Time
    Periods](https://bcallaway11.github.io/did/articles/multi-period-did.html)

  - [Pre-Testing in a DiD Setup using the did
    Package](https://bcallaway11.github.io/did/articles/pre-testing.html)

  - [Writing Extensions to the did
    Package](https://bcallaway11.github.io/did/articles/extensions.html)
