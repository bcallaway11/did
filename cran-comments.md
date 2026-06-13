## Test environments

* local Ubuntu 24.04.4 LTS (x86_64), R 4.6.0
* local macOS 26.5 (aarch64-apple-darwin), R 4.6.0
* GitHub Actions (macOS-latest): R release
* GitHub Actions (windows-latest): R release
* GitHub Actions (ubuntu-latest): R devel, R release, R oldrel-1

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

This release requires 'DRDID' (>= 1.3.0), which is available on CRAN
(published 2026-06-10).

## Notable changes

did 2.5.0 is a large release that consolidates all development since 2.3.0
(see NEWS.md) and bumps the required version of 'DRDID' to (>= 1.3.0).

## Package size

The installed size is about 5.3 Mb, of which 4.4 Mb is the `doc` directory
holding the package's five pre-built vignettes (long-form methodological
documentation with many plots). The source tarball is about 4.5 Mb.

## Downstream dependencies

There are 12 reverse dependencies on CRAN: cdid, did2s, etwfe, fastdid, fect,
fetwfe, fixes, modelsummary, NonlinearDiD, optic, parameters, ptetools.

All 12 were checked locally with `revdeplite`. 10 pass with 0 ERRORs, 0
WARNINGs, and 0 NOTEs. The two failures and one warning are pre-existing issues
in those packages unrelated to `did`:

* **fetwfe**: ERROR in examples — `library(bacondecomp)` fails because
  `bacondecomp` is not installed; this is a missing dependency in `fetwfe`'s
  own examples and vignette, not caused by `did`.
* **parameters**: ERROR in tests — requires `clubSandwich` which is not
  installed; this is a missing optional dependency in `parameters` itself.
* **ptetools**: WARNING — calls the deprecated `BMisc::rhs.vars()` in its own
  source code; no connection to changes in `did`.
