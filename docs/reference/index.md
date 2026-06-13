# Package index

## Difference-in-Differences Methods

Tools to compute average treatment effects, particularly in the case
with multiple periods and variation in treatment timing

- [`att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.md) :
  Group-Time Average Treatment Effects
- [`aggte()`](https://bcallaway11.github.io/did/reference/aggte.md) :
  Aggregate Group-Time Average Treatment Effects

## Pre-Testing

Tools to compute “pre-test” the DiD assumption in the case where
multiple periods are available

- [`conditional_did_pretest()`](https://bcallaway11.github.io/did/reference/conditional_did_pretest.md)
  : Pre-Test of Conditional Parallel Trends Assumption

## Plotting and Summarizing

Tools for plotting and summarizing the results of the **att_gt** method
and **aggte** method

- [`ggdid()`](https://bcallaway11.github.io/did/reference/ggdid.md) :

  Plot `did` objects using `ggplot2`

- [`ggdid(`*`<MP>`*`)`](https://bcallaway11.github.io/did/reference/ggdid.MP.md)
  :

  Plot `MP` objects using `ggplot2`

- [`ggdid(`*`<AGGTEobj>`*`)`](https://bcallaway11.github.io/did/reference/ggdid.AGGTEobj.md)
  :

  Plot `AGGTEobj` objects

- [`summary(`*`<MP>`*`)`](https://bcallaway11.github.io/did/reference/summary.MP.md)
  : summary.MP

- [`summary(`*`<AGGTEobj>`*`)`](https://bcallaway11.github.io/did/reference/summary.AGGTEobj.md)
  : Summary Aggregate Treatment Effect Parameter Objects

- [`summary(`*`<MP.TEST>`*`)`](https://bcallaway11.github.io/did/reference/summary.MP.TEST.md)
  : summary.MP.TEST

- [`print(`*`<MP>`*`)`](https://bcallaway11.github.io/did/reference/print.MP.md)
  : print.MP

- [`print(`*`<AGGTEobj>`*`)`](https://bcallaway11.github.io/did/reference/print.AGGTEobj.md)
  : print.AGGTEobj

## Classes

Classes available in the **did** package

- [`MP()`](https://bcallaway11.github.io/did/reference/MP.md) : MP
- [`AGGTEobj()`](https://bcallaway11.github.io/did/reference/AGGTEobj.md)
  : AGGTEobj
- [`MP.TEST()`](https://bcallaway11.github.io/did/reference/MP.TEST.md)
  : MP.TEST
- [`DIDparams()`](https://bcallaway11.github.io/did/reference/DIDparams.md)
  : DIDparams

## Helper Functions

A number of other miscellaneous helper functions

- [`mboot()`](https://bcallaway11.github.io/did/reference/mboot.md) :
  Multiplier Bootstrap

- [`test.mboot()`](https://bcallaway11.github.io/did/reference/test.mboot.md)
  : Multiplier Bootstrap for Conditional Moment Test

- [`pre_process_did()`](https://bcallaway11.github.io/did/reference/pre_process_did.md)
  :

  Process `did` Function Arguments

- [`pre_process_did2()`](https://bcallaway11.github.io/did/reference/pre_process_did2.md)
  :

  Process `did` Function Arguments

- [`process_attgt()`](https://bcallaway11.github.io/did/reference/process_attgt.md)
  :

  Process Results from
  [`compute.att_gt()`](https://bcallaway11.github.io/did/reference/compute.att_gt.md)

- [`indicator()`](https://bcallaway11.github.io/did/reference/indicator.md)
  : indicator

- [`build_sim_dataset()`](https://bcallaway11.github.io/did/reference/build_sim_dataset.md)
  : build_sim_dataset

- [`reset.sim()`](https://bcallaway11.github.io/did/reference/reset.sim.md)
  : reset.sim

- [`sim()`](https://bcallaway11.github.io/did/reference/sim.md) : sim

- [`trimmer()`](https://bcallaway11.github.io/did/reference/trimmer.md)
  : trimmer

- [`glance(`*`<AGGTEobj>`*`)`](https://bcallaway11.github.io/did/reference/glance.AGGTEobj.md)
  : glance model characteristics from AGGTEobj objects

- [`glance(`*`<MP>`*`)`](https://bcallaway11.github.io/did/reference/glance.MP.md)
  : glance model characteristics from MP objects

- [`nobs(`*`<AGGTEobj>`*`)`](https://bcallaway11.github.io/did/reference/nobs.AGGTEobj.md)
  : Number of unique cross-sectional units in an AGGTEobj object

- [`nobs(`*`<MP>`*`)`](https://bcallaway11.github.io/did/reference/nobs.MP.md)
  : Number of unique cross-sectional units in an MP object

- [`tidy(`*`<AGGTEobj>`*`)`](https://bcallaway11.github.io/did/reference/tidy.AGGTEobj.md)
  : Tidy an AGGTEobj into a data frame

- [`tidy(`*`<MP>`*`)`](https://bcallaway11.github.io/did/reference/tidy.MP.md)
  : Tidy an MP object into a data frame

## Data

Available datasets in the package

- [`mpdta`](https://bcallaway11.github.io/did/reference/mpdta.md) :
  County Teen Employment Dataset

## Package Documentation

- [`did`](https://bcallaway11.github.io/did/reference/did-package.md)
  [`did-package`](https://bcallaway11.github.io/did/reference/did-package.md)
  : Difference in Differences
