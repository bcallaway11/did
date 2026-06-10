# did 2.5.0

This is a large release that consolidates all development since 2.3.0. Headline
changes: a substantially faster engine (group-time ATTs are roughly 2.5-3x faster in
common settings, and the conditional pre-test is several times faster and far lighter
on memory), first-class clustered and unbalanced-panel inference, support for
transformation and factor covariates, a point-estimates-only mode, and a long list of
correctness fixes. Numerical results are unchanged up to floating-point precision
except where a bug fix is explicitly noted.

## Dependencies

  * Requires `DRDID (>= 1.3.0)`, which provides faster and more robust 2x2 DiD estimators (used internally for every group-time ATT) and guards against silently-incorrect standard errors on ill-conditioned (near-singular) designs.

## New features

  * New `compute_inffunc` argument in `att_gt()` (default `TRUE`). Set `compute_inffunc = FALSE` for a **point-estimates-only** run: it returns the group-time ATT point estimates (identical to a full run) without influence functions, standard errors, uniform bands, or the pre-test. This is faster and uses much less memory (no \eqn{n \times k} influence-function matrix is ever formed or bootstrapped), which helps for quick exploration or very large datasets. The result cannot be passed to `aggte()` (which now errors with a clear message); `bstrap` and `cband` are set to `FALSE` automatically.

  * Covariate formulas (`xformla`) may now use transformations and other model-matrix terms, e.g. `~ I(X^2)`, `~ log(X)`, `~ poly(X, 2)`, `~ X1 * X2`. These previously errored because pre-processing stored the evaluated model frame (losing the raw variable and creating matrix-valued columns) instead of the raw covariates. Pre-processing now keeps the raw covariates so the design matrix can be rebuilt, and drops rows whose evaluated design is non-finite (e.g. `log()` of a non-positive value).

  * Factor covariates now produce exactly the same estimates, standard errors, and warning messages as adding their dummy columns by hand. Previously the `faster_mode = FALSE` path applied `droplevels()` within each 2x2 comparison, so a factor level absent from a comparison changed the design (and could error with "contrasts can be applied only to factors with 2 or more levels"). The design matrix is now built once over the full sample, with global factor levels, and row-subset per cell.

  * New `fix_weights` argument in `att_gt()` for explicit control over how time-varying sampling weights are resolved in each 2x2 comparison: `NULL` (default, prior behavior), `"varying"` (per-observation weights via the RC estimators), `"base_period"` (fix at g-1), or `"first_period"`. See `?att_gt`. A runtime message points users to it when time-varying weights are detected in balanced panel data.

  * `att_gt()` accepts `...` to forward additional arguments to a custom `est_method` function.

  * Added `nobs()` S3 methods for `MP` and `AGGTEobj` objects (number of unique cross-sectional units), and `statistic` (t-statistic) and `p.value` (pointwise, two-sided) columns to `tidy()` output for both classes, following `broom` conventions.

## Performance and memory

  * `att_gt(faster_mode = TRUE)` (the default) is about 2.5-3x faster on common problems, with identical results. Each (g,t) cell previously built and extracted a `data.table` even though the cohort vectors are already materialized in the pre-computed tensors; the per-cell cohort is now a plain list of vectors passed straight to the DRDID estimators (panel, repeated cross-section, and unbalanced-panel paths), and the unbalanced-panel influence-function aggregation uses `rowsum()` in place of a per-cell `data.table` group-by. Earlier (g,t) work also feeds in: cumulative cohort sizes, cached pre-treatment periods, and sparse-triplet influence-matrix construction.

  * The per-(g,t) overlap-check propensity logit (fit for every `dr`/`ipw` cell to detect propensity-score overlap violations) now uses `fastglm`'s low-level entry point (`fastglmPure`) instead of the `fastglm()` wrapper, skipping its per-call input coercion and family/deviance bookkeeping. The fitted values -- and therefore the overlap decision -- are bit-identical.

  * `att_gt(faster_mode = FALSE)` builds the covariate design matrix once and assembles each 2x2 cell directly from precomputed per-period blocks (outcomes, weights, design) indexed by position, instead of rebuilding `model.matrix()` and reshaping (`get_wide_data()`) the long data for every cell. Bit-identical, with about half the transient allocation; `options(did.disable_precompute = TRUE)` restores the original path. The repeated-cross-section / unbalanced influence-function aggregation now uses `rowsum()` instead of `stats::aggregate()` (about 40x faster on that step). `faster_mode = TRUE` and `faster_mode = FALSE` remain identical to numerical precision for every supported option.

  * The conditional pre-test (`conditional_did_pretest()`) is several times faster end-to-end and far lighter on memory. Its multiplier bootstrap (`test.mboot()`) previously looped over `biters` draws, each multiplying the full `n x k x nX` influence array by fresh weights -- `O(n^2 k)` work **and** an `O(n^2 k)` transient allocation per draw (over 1 GB per draw at a few thousand units). It is now a single tiled matrix contraction (100x+ faster on that step, with the per-draw gigabyte allocations eliminated), numerically identical to the old loop up to floating-point summation order; the `indicator()` weighting function is also vectorized.

  * Internal speedups with identical results: vectorized the multiplier-bootstrap post-processing (`mboot()`), removed a duplicated `n x k` matrix construction in the aggregation estimated-weight influence term (`wif()`), preallocated the sparse influence-function assembly, and removed redundant work in pre-processing and simulation.

## Clustered and unbalanced-panel inference

  * Clustered standard errors are now available **without** the bootstrap. With `clustervars` set and `bstrap = FALSE`, `att_gt()` and `aggte()` report cluster-robust standard errors computed analytically from the cluster sums of the influence function, at every aggregation level (group-time, simple, group, dynamic, calendar), and the pre-test Wald statistic is reported under clustering.

  * The cluster-robust multiplier bootstrap (`mboot`) now follows Callaway & Sant'Anna (2021), Remark 10: it draws one multiplier per cluster and aggregates the influence function to cluster sums. Identical to before for equal-sized clusters; correct cluster-sum aggregation for unbalanced clusters and repeated cross sections.

  * Clustered inference (bootstrap and analytical) is supported for panel data, unbalanced panels, and repeated cross sections, whether or not `idname` is supplied, and is identical under `faster_mode = TRUE` and `faster_mode = FALSE`.

  * `aggte()` no longer silently ignores a `clustervars` request it cannot honor (the aggregation can only use the cluster information `att_gt()` retained). It now warns -- including when an override names a different variable than `att_gt()` clustered on -- and falls back to non-clustered standard errors, instead of silently returning the i.i.d. error or crashing in `mboot()`.

  * Fixed two `faster_mode = TRUE` clustered-standard-error bugs on **unbalanced panels** so the fast path reproduces the slow path: (1) the analytical clustered SE silently fell back to the i.i.d. error because the stored per-unit cluster vector was observation-length and no longer aligned with the influence function; and (2) in `aggte()`, the estimated-weight influence term was added in id-sorted order while the influence function is in first-appearance order, misattributing it and giving a wrong aggregated SE (point estimates were unaffected). Balanced panels and repeated cross sections were unaffected.

## Bug fixes

  * Fixed the conditional parallel-trends pre-test (`conditional_did_pretest()`), which had silently broken under R >= 4.0 and spuriously rejected almost always whenever there was more than one pre-treatment ATT(g,t) cell. The observed Cramér-von Mises statistic was left in `(n_gt x nX)` orientation while its bootstrap null distribution is `(nX x n_gt)`, scaling the observed statistic by `n / n_gt` and driving the p-value to ~0. The root cause was `ifelse(class(J) == "matrix", ...)`: `class()` of a matrix became length-2 (`c("matrix","array")`) in R 4.0, so `ifelse()` evaluated both branches and the no-transpose branch always won. The orientation is now selected with `is.matrix()`.

  * `aggte(type = "group", na.rm = TRUE)` with a finite `max_e` no longer errors ("No valid att_gt() estimates found ...") when a group's only non-missing ATT(g,t) lies past `max_e`; the group filter now applies the same `max_e` window. The default `max_e = Inf` is unchanged.

  * Duplicated `(idname, tname)` rows (the same unit observed more than once in a period, a common long-format data-prep mistake) are now rejected with a clear error in **both** code paths. Previously only `faster_mode = TRUE` caught this; `faster_mode = FALSE` silently produced incorrect estimates.

  * A `weightsname` column with negative values or a non-positive mean is now rejected with a clear, identical error in both code paths, instead of silently producing `NA`/`NaN` estimates.

  * Fixed `faster_mode = TRUE` vs `FALSE` ATT disagreement when sampling weights (`weightsname`) vary across time: the fast path was always using first-period weights and now uses the same period's weights as the slow path.

  * Fixed influence-function aggregation for `fix_weights = "varying"` on balanced panels (now aggregates by unit id with `rowsum()` rather than assuming stacked order), and a length mismatch for `fix_weights = "base_period"`/`"first_period"` on unbalanced panels after weight-based unit dropping.

  * Fixed `glance.MP()` returning `NULL` for `ngroup`/`ntime` under `faster_mode = TRUE`.

  * Fixed an `aggte()` crash ("Error in get(gname): invalid first argument") when the group column is literally named `gname` and `dreamerr >= 1.5.0` is installed (data.table's `get()` was intercepted; replaced with `set()`).

  * Fixed groups treated after the last observed period but within the anticipation window being coerced to never-treated (contaminating the control group with anticipation effects), and a data-filter inconsistency for always-treated units when `anticipation > 0`. Added an informative message clarifying that never-treated units are unaffected by `anticipation`.

  * When internal 2x2 estimation fails for a specific (g,t) cell (e.g. a singular design), `att_gt()` now warns and sets that cell's ATT to `NA` instead of crashing, in both `faster_mode = TRUE` and `FALSE` (#185, #190).

  * `pl = TRUE` on Windows now warns and falls back to sequential processing instead of crashing (#176).

## Validation and clearer errors

  * Misspelled `yname`, `idname`, `tname`, `gname`, `weightsname`, or `clustervars` now produce a clear message listing the missing columns (#203).

  * An invalid `est_method` (an unrecognized string or an unquoted variable) now errors clearly instead of silently defaulting to `"dr"` (#194).

  * `fix_weights = "base_period"`/`"first_period"` are blocked for repeated cross sections (`panel = FALSE`); `fix_weights = "varying"` is blocked with a custom `est_method` function (whose signature differs from the internal RC path). Both with clear messages.

  * `panel = TRUE` (the default) without `idname` now errors with "Must provide idname when panel = TRUE. Set panel = FALSE for repeated cross sections." Previously this failed with a cryptic internal `data.table` error (`faster_mode = TRUE`) or a misleading "All observations dropped while converting data to balanced panel" message (`faster_mode = FALSE`).

  * A non-numeric outcome variable (`yname`) is now rejected up front with a clear message in both code paths (logical 0/1 outcomes remain allowed). Previously a character or factor outcome "ran" to completion with all-`NA` ATTs and misleading per-cell warnings, and a list-column outcome failed with a cryptic `complete.cases()` error.

  * `alp` must now be a single number strictly between 0 and 1 (e.g. `alp = 1.5` previously inverted the confidence bands silently or errored deep inside `quantile()`), and `biters` must be a single positive whole number when `bstrap = TRUE` (a negative value previously crashed inside the bootstrap's linear-algebra code with no hint about the cause).

## Documentation, namespace, and internals

  * Reduced namespace pollution: replaced blanket `import(stats)`, `import(utils)`, and `import(BMisc)` with selective `importFrom()` calls. `did` no longer re-exports `stats::filter`/`stats::lag` (which previously masked `dplyr::filter`/`dplyr::lag`). `R CMD check` passes with 0 code-related NOTEs.

  * Replaced fragile `ifelse(cond, x <- a, x <- b)` side-effect idioms (which relied on R's branch-evaluation order) with `if/else`, and `get()`/`:=` with `set()` inside `data.table` loops, throughout; behavior is unchanged.

  * Substantially expanded the test suite: `glance()`, `ggdid()`, error handling, edge cases, all aggregation types, systematic `faster_mode` consistency across dozens of parameter combinations, and JEL replication tests. The suite now runs with 0 warnings (previously 66+). Added a GitHub Action to auto-bump the dev version in `DESCRIPTION` on PR merge.

  * Expanded `weightsname` documentation (how time-varying weights are handled for balanced panels vs. repeated cross sections / unbalanced panels); grammar and typo fixes across docs, vignettes, and error messages; corrected the `mpdta` data documentation.

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



