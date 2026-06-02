# did 2.3.1.907

  * `edid()` gains an opt-in `estimation_effect` argument (default `FALSE`). When `TRUE`, the influence function is augmented with the first-step nuisance-estimation correction of Ackerberg, Chen & Hahn (2012) for the sieve nuisances (conditional means and propensity ratios) entering the doubly-robust moment. The correction is asymptotically negligible under correct specification (the EIF moments are Neyman orthogonal) and provides finite-sample robustness when a first-step nuisance is misspecified; it propagates automatically to the event-study and overall aggregations. Default behavior is unchanged.

# did 2.3.1.906

  * The cluster-robust multiplier bootstrap (`mboot`) now follows Callaway & Sant'Anna (2021), Remark 10: it draws one multiplier per cluster and aggregates the influence function to cluster sums. For equal-sized clusters the standard errors are identical to before; for unbalanced clusters and repeated cross sections it uses the appropriate cluster-sum aggregation. All clustering input checks (numeric `clustervars`, at most one cluster dimension beyond `idname`, and cluster-variable time-invariance) are preserved

  * Clustered standard errors are now available **without** the bootstrap. With `clustervars` set and `bstrap = FALSE`, `att_gt()` and `aggte()` report cluster-robust standard errors computed analytically from the cluster sums of the influence function (the same cluster-sum aggregation as the bootstrap), at every aggregation level (group-time, simple, group, dynamic, calendar). The pre-test Wald statistic is also reported under clustering in this case

  * Clustered inference (bootstrap and analytical) is supported for panel data, unbalanced panels, and repeated cross sections, whether or not `idname` is supplied, and gives identical standard errors under `faster_mode = TRUE` and `faster_mode = FALSE`. For repeated cross sections without an `idname`, the internal observation id is used to align the cluster identifiers with the influence function

  * Added `edid()`: efficient DiD estimator (Chen, Sant'Anna & Xie 2025) supporting PT-All and PT-Post parallel trends assumptions, analytical EIF-based standard errors (iid and cluster-robust), multiplier bootstrap inference (Rademacher, Mammen, Webb), and overall/event-study/group aggregation with WIF correction.

  * Fixed bug where `faster_mode = TRUE` and `faster_mode = FALSE` produced different ATT estimates when sampling weights (`weightsname`) vary across time. The fast path was always using first-period weights; it now correctly uses the same period's weights as the slow path

  * New `fix_weights` argument in `att_gt()` gives users explicit control over how time-varying sampling weights are resolved in each 2x2 DiD comparison. Options: `NULL` (default, preserves existing behavior), `"varying"` (per-observation weights using RC estimators), `"base_period"` (fix at g-1 for all cells), `"first_period"` (fix at first period). See `?att_gt` for details

  * Runtime message when time-varying weights are detected in balanced panel data, directing users to the `fix_weights` argument

  * Reduced namespace pollution: replaced blanket `import(stats)`, `import(utils)`, and `import(BMisc)` with selective `importFrom()` calls. The `did` package no longer re-exports `stats::filter` or `stats::lag`, which previously masked `dplyr::filter` and `dplyr::lag` when both packages were loaded

  * Fixed `aggte()` crash (`"Error in get(gname): invalid first argument"`) when the user's group column is literally named `gname` and `dreamerr` >= 1.5.0 is installed. The issue was `dreamerr` intercepting `data.table`'s `get()` inside `[.data.table`; replaced with `set()` which is immune to this

  * Expanded `weightsname` documentation explaining how time-varying weights are handled differently for balanced panels vs. repeated cross sections and unbalanced panels

  * Added `nobs()` S3 methods for `MP` and `AGGTEobj` objects, returning the number of unique cross-sectional units as an integer

  * Added `statistic` (t-statistic) and `p.value` (pointwise, two-sided) columns to `tidy()` output for both `MP` and `AGGTEobj` objects, following `broom` conventions

  * Fixed `glance.MP()` returning `NULL` for `ngroup` and `ntime` when `faster_mode = TRUE` (DIDparams2 uses different field names than DIDparams)

  * Fixed influence function aggregation bug in `fix_weights = "varying"` on balanced panel: the slow path now correctly uses `rowsum()` by unit ID instead of assuming stacked ordering

  * Fixed `fix_weights = "base_period"` / `"first_period"` for unbalanced panels: corrected influence function length mismatch after weight-based unit dropping, and fast path now properly excludes units missing from the target period

  * Added validation: `fix_weights = "base_period"` and `"first_period"` are blocked for repeated cross sections (`panel = FALSE`) with a clear error message

  * Added validation: `fix_weights = "varying"` is blocked when using a custom `est_method` function, since the varying path uses RC estimators internally with a different function signature

  * Completed namespace cleanup: added missing `@importFrom` for `stats::ecdf`, `stats::glm`, `stats::predict`; registered `.w`, `D`, `N`, `post`, `weights` as `globalVariables` for data.table column references. `R CMD check` now passes with 0 code-related NOTEs

  * Replaced fragile `exists("use_rc_for_weights")` with direct `dp2$fix_weights` check in `compute.att_gt2()`

  * Fixed `data.table` `.checkTypos` crash in `get_wide_data()` when user column names match local variable names (e.g., column named `tname`)

  * Replaced `get()` with `c()` in time-varying weight detection grouping to avoid potential `dreamerr` interception

  * Inference tests (`test-inference.R`): switched to HTTPS mirror, added `requireNamespace` verification, wrapped install in `skip_on_cran` guard, added proper temp directory cleanup

  * Added GitHub Action to auto-bump dev version in DESCRIPTION on PR merge

  * Substantially expanded test suite covering `glance()`, `ggdid()`, error handling, edge cases, all aggregation types, and systematic `faster_mode` consistency across 36 parameter combinations. Test suite now runs with 0 warnings (previously 66+)

# did 2.3.1.903

  * Added `nobs()` S3 methods and `statistic`/`p.value` columns to `tidy()` output (superseded by 2.3.1.904 entry above)

# did 2.3.1.902

  * Bug fixes, diagnostic improvements, and JEL replication tests

# did 2.3.1.901

  * `att_gt()` now accepts `...` (dots) for passing additional arguments to custom `est_method` functions

  * Fixed bug where groups treated after the last observed period but within the anticipation window were incorrectly coerced to never-treated, contaminating the control group with anticipation effects

  * Fixed data filter inconsistency for always-treated units when `anticipation > 0`

  * Added informative message when `anticipation > 0` clarifying that never-treated units are not affected by anticipation

  * Performance optimizations for `faster_mode = TRUE`: cumulative cohort sizes, cached pre-treatment periods, sparse matrix triplet construction

  * Robustness improvements: replaced `get()` and `:=` with `set()` inside data.table loops to prevent non-standard evaluation scoping bugs

  * Grammar and typo fixes across documentation, vignettes, and error messages

  * Fixed typo in `pre-testing` vignette (`idnam` to `idname`)

  * Fixed `mpdta` data documentation to correctly report 6 variables

  * Removed development path reference from `README.Rmd`

  * Improved error messages throughout the package for better user experience

  * Column name validation: misspelled `yname`, `idname`, `tname`, `gname`, `weightsname`, or `clustervars` now produces a clear message listing the missing columns (fixes #203)

  * `est_method` validation: passing an invalid string or unquoted variable now produces a clear error instead of silently defaulting to `"dr"` (fixes #194)

  * Windows parallel processing: `pl = TRUE` on Windows now warns and falls back to sequential processing instead of crashing (fixes #176)

  * `aggte()` dimension errors: improved guards against matrix subsetting failures when all `att_gt()` estimates are NA (fixes #185, #190)

  * Graceful error handling: when internal 2x2 DiD estimation fails for a specific (g, t) cell (e.g., singular matrix), `att_gt()` now warns and sets that cell's ATT to NA instead of crashing. This applies to both `faster_mode = TRUE` and `faster_mode = FALSE`

# did 2.3.0

  * Code improvements that make the package faster and more memory efficient

  * Improved automated testing and regression testing

  * Check if data is balanced if `panel = TRUE` and `allow_unbalanced_panel = TRUE`. If it is, disable `allow_unbalanced_panel` and proceed with panel data setup. This is different from the previous behavior, which would always proceed as if `panel = FALSE`.

  * Significantly reduced the number of recursive package dependencies, enabling faster installation times and a smaller build footprint.

# did 2.2.0

  * Skipped, as the number of changes was large enough to warrant a minor version bump

# did 2.1.2

  * Added wrapper function for HonestDiD package

  * Fix bug for setups where `gname` is not contained in `tname` (but is in the `tname` range)

  * Fix bug for including too many groups with universal base period in pre-treatment periods

  * Bug fix for anticipation using `notyettreated` comparison group

# did 2.1.1

  * Bug fixes related to unbalanced panel and clustered standard errors

  * Bug fixes for conditional_did_pretest

  * Even faster bootstrap code (thanks Kyle Butts)

  * Updated version requirement for `BMisc` package

  * Bug fix for unbalanced panel and repeated cross sections in pre-treatment periods using universal base period

# did 2.1.0

  * Code is substantially faster/more memory efficient

  * Support for *universal* base period

  * Major improvements to unit testing

  * Completely removed `mp.spatt` and `mp.spatt.test` functions (which were the original names for `att_gt`)

  * Simulation/testing code now exported

  * Removed some slow running checks

  * Multiplier bootstrap code is now written in C++

  * Improvements to error handling, added some additional warning messages, removed some unnecessary warning messages

  * Bug fixes for NA standard errors that occur with very small groups


# did 2.0.1

  * Improved plots

  * Maximum event time for event studies

  * Compute critical value for simultaneous confidence bands even when some standard error is zero (set these to NA)

  * Improved codes for unbalanced panel data: faster and more memory efficient

  * Correct estimates of P(G=g|Eventually treated) with unbalanced panel data. This affects **aggte** objects with unbalanced panel data

  * Bug fixes for summary **aggte** objects

  * Allow clustering for unbalanced panel data

  * Fixed error in calendar-type aggregation within **aggte** function (point estimates were not being weighted by group-size; now they are).

  * Additional error handling


# did 2.0.0
  * Big improvement on code base / functionality / testing

  * Deprecated **mp.spatt** function and replaced it with **att_gt** function

  * Calling **att_gt** is similar to calling **mp.spatt**; instead of formula for outcome of the form `y~treat`, now just pass the name of the outcome variable

  * Deprecated **mp.spatt** function and replaced it with **conditional_did_pretest** function

  * New **est_method** parameter.  Can call any function for 2x2 DID in the **DRDID** package (default is now doubly robust estimation, but inverse probability weights and regression estimators are also supported) as well as provide custom 2x2 DID estimators

  * Bug fixes for including groups that are *already treated* in the first period

  * Allow for user to select control group -- either *never treated* or *not yet treated*

  * Add functionality for uniform confidence bands for all aggregated treatment effect parameters

  * Introduced dynamic effects in pre-treatment periods.  These allow for users to report event study plots that are common that include pre-treatment periods and are common in applied work.  The event study plots in the **did** package are robust to selective treatment timing (unlike standard regression event study plots)

  * Support for using repeated cross sections data instead of panel data is much improved

  * Support for using sampling weights is much improved

  * Big improvement to website, vignettes, and code documentation

  * Code for dealing with unbalanced panels

  * Allow for event studies to be computed over subsets of event times

  * Allow for treatment anticipation via *anticipation* argument

# did 1.2.3
  * Corrected check problems

# did 1.2.2

  * Improved ways to summarize aggregated treatment effect parameters

  * Fixed bug related to needing new version of BMisc

  * Fixed bug related to plotting with no pre-treatment periods

  * Improved ways to easily plot aggregated treatment effect parameters

# did 1.2.1

  * Added some error handling for some cases with small group sizes, and fixed some cryptic error messages

  * Fixes handling for data being in format besides *data.frame* (e.g. *tibble*)

  * Add warnings about small group sizes which are a very common source of trouble

# did 1.2.0

  * Updates for handling repeated cross sections data, both estimation and inference

# did 1.1.2

  * bug fixes for testing without covariates, allowed to pass NULL in addition to ~1

# did 1.1.1

  * fixed issues between BMisc and formula.tools

  * added point estimates for repeated cross sections data

# did 1.1.0

  * bug fixes for the case without any covariates

# did 1.0.0

  * first version of package, functions for computing group-time average treatment effects, combining them into a smaller number of parameters, and pre-testing the common trends assumption



