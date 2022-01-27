* There is one error at https://cloud.r-project.org/web/checks/check_results_did.html.  It's origin is that my package is downstream dependent of `nloptr`.  That package appears to have been updated yesterday, and it appears that is issue is now fixed.

## Test environments

* local Ubuntu 20.04, R 4.0.3
* win-builder (devel and release)
* R-hub (Windows Server, Ubuntu Linux, Fedora Linux)

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

## Downstream dependencies

We checked 1 reverse dependencies (0 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
