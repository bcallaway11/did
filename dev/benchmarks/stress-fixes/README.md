# Performance checks for the stress-test fixes

The main benchmark compared clean `origin/master` commit `74b88eb07f1faa644df1271055a0f77848f6550c` with the combined fixes at `1fbfd29`, using the same installed DRDID 1.4.0. It detected no material end-to-end slowdown in the tested workloads. Three end-to-end medians increased by 0.4–0.9% (0.5–1 ms), with overlapping interquartile ranges. Across the 28 end-to-end cases, the geometric mean final/baseline time ratio was 0.950.

The isolated weighted rank guard costs more: 100 calls took 4 ms on master and 9 ms with the correction (about 50 microseconds extra per call). The corrected guard evaluates the sampling-weighted design; master checked the wrong unweighted matrix. The corresponding weighted panel workloads showed no measured end-to-end regression. This does not mean every internal check is faster.

## Reproduce

The examples below run from the repository root, with R and all package dependencies installed. The script may also be invoked by absolute path from another working directory: it finds the repository and worker from its own location. Copied source trees stay outside the tracked benchmark directory. Both workers assert DRDID 1.4.0.

```sh
Rscript dev/benchmarks/stress-fixes/performance.R 7 2 /tmp/did-main-performance main
Rscript dev/benchmarks/stress-fixes/performance.R 7 2 /tmp/did-supplemental-performance supplemental
```

The driver archives the fixed baseline and snapshots the current working source. Its metadata records the source commit, dirty-file summary, and R-source hashes. Re-running after further code changes measures those changes. Output includes copied source trees, complete worker sessions, raw timings, medians, interquartile ranges, and final/baseline ratios. Do not commit the output source trees.

## Measurement design

- R 4.6.0 on aarch64 macOS Tahoe 26.5.1, Accelerate BLAS; DRDID 1.4.0. Full versions are in `session-info.txt`.
- Data seed: 20261028. T = 6 with n = 1,000 and 5,000 panel units; repeated cross sections have one independent observation ID per row. The supplement adds n = 3,000, T = 12.
- Separate sequential R processes in baseline/final/final/baseline order. Each path is warmed twice in each process, then measured seven times: 14 measurements per source/case.
- Cases are randomly ordered within each timed round. Loading, data generation, and explicit garbage collection occur outside timings. Automatic allocation/collection costs during the function call remain included.
- Other agents and the parent process stopped heavy R work during timing.
- Bootstrap cases use 499 iterations. The small conditional pre-test uses n = 60 and 99 iterations.
- Millisecond elapsed-time resolution and local machine noise limit interpretation of small differences. IQR overlap is descriptive and is not a formal equivalence test. These workloads do not establish a universal performance bound.

## Main measurements

All values below are elapsed milliseconds. Guard microbenchmarks batch 100 calls, or 1,000 calls for intercept-only checks.

| Case | Master median [Q25, Q75] | Fixed median [Q25, Q75] | Fixed / master |
|---|---:|---:|---:|
| `aggte_1k_dynamic` | 4.00 [4.00, 5.00] | 4.00 [4.00, 4.75] | 1.000 |
| `aggte_1k_dynamic_boot` | 18.00 [17.00, 18.75] | 17.00 [17.00, 18.00] | 0.944 |
| `aggte_5k_dynamic` | 6.00 [6.00, 6.00] | 6.00 [6.00, 6.00] | 1.000 |
| `aggte_5k_simple` | 2.00 [1.25, 2.00] | 1.50 [1.00, 2.00] | 0.750 |
| `aggte_unbalanced_dynamic` | 30.50 [30.00, 32.75] | 29.50 [29.00, 30.00] | 0.967 |
| `conditional_pretest_60` | 88.00 [82.25, 98.50] | 81.50 [81.00, 82.00] | 0.926 |
| `guard_overlap_5k` | 75.00 [72.50, 77.00] | 73.00 [72.00, 73.00] | 0.973 |
| `guard_overlap_intercept_5k` | 12.50 [12.00, 13.00] | 12.00 [12.00, 13.00] | 0.960 |
| `guard_rank_5k` | 4.00 [4.00, 4.00] | 4.00 [3.00, 4.00] | 1.000 |
| `guard_rank_intercept_5k` | 7.00 [7.00, 7.00] | 7.00 [7.00, 7.00] | 1.000 |
| `guard_rank_weighted_5k` | 4.00 [4.00, 4.00] | 9.00 [9.00, 9.75] | 2.250 |
| `panel_1k_fast_dr_boot` | 30.00 [29.00, 30.75] | 26.00 [26.00, 26.75] | 0.867 |
| `panel_1k_fast_dr_nox` | 19.50 [19.00, 20.00] | 17.50 [17.00, 18.00] | 0.897 |
| `panel_1k_fast_dr_x` | 25.00 [24.00, 26.00] | 22.00 [22.00, 22.00] | 0.880 |
| `panel_1k_fast_reg_x` | 18.00 [17.00, 18.00] | 16.00 [16.00, 16.00] | 0.889 |
| `panel_1k_slow_dr_x` | 16.00 [16.00, 17.00] | 16.00 [15.25, 16.00] | 1.000 |
| `panel_5k_fast_dr_boot` | 81.50 [81.00, 84.50] | 78.00 [77.00, 78.00] | 0.957 |
| `panel_5k_fast_dr_nox` | 40.00 [39.25, 42.50] | 37.00 [37.00, 38.00] | 0.925 |
| `panel_5k_fast_dr_weighted` | 71.50 [71.00, 74.50] | 67.50 [67.00, 68.00] | 0.944 |
| `panel_5k_fast_dr_x` | 65.00 [64.00, 68.25] | 63.00 [62.00, 64.00] | 0.969 |
| `panel_5k_fast_reg_x` | 44.50 [43.00, 45.75] | 41.00 [41.00, 41.75] | 0.921 |
| `panel_5k_slow_dr_nox` | 29.50 [29.00, 30.75] | 29.00 [28.25, 29.00] | 0.983 |
| `panel_5k_slow_dr_weighted` | 50.00 [49.25, 52.75] | 50.00 [49.00, 50.75] | 1.000 |
| `panel_5k_slow_dr_x` | 48.00 [46.00, 49.00] | 46.50 [45.25, 47.00] | 0.969 |
| `panel_5k_slow_reg_x` | 24.00 [23.00, 25.75] | 24.00 [23.00, 24.00] | 1.000 |
| `rc_1k_fast_dr_boot` | 69.00 [68.00, 69.75] | 66.50 [64.25, 68.00] | 0.964 |
| `rc_1k_fast_dr_x` | 49.00 [49.00, 49.00] | 48.00 [47.25, 49.75] | 0.980 |
| `rc_5k_fast_dr_x` | 148.00 [147.00, 149.00] | 149.00 [148.00, 153.00] | 1.007 |
| `rc_5k_slow_dr_x` | 135.50 [134.00, 139.00] | 136.00 [136.00, 137.75] | 1.004 |
| `unbalanced_1k_fast_dr_boot` | 54.50 [53.25, 56.00] | 52.00 [50.25, 52.75] | 0.954 |
| `unbalanced_1k_fast_dr_x` | 48.50 [47.00, 49.75] | 46.00 [46.00, 48.00] | 0.948 |
| `unbalanced_5k_fast_dr_x` | 137.50 [135.25, 144.00] | 137.00 [136.00, 138.00] | 0.996 |
| `unbalanced_5k_slow_dr_x` | 113.50 [112.00, 118.00] | 114.50 [113.25, 118.25] | 1.009 |

Raw observations are in `main-raw-timings.csv`; the complete summary, including CPU medians, is in `main-comparison.csv`.

## Supplemental measurements

The first supplemental run at `1a16c46` detected two regressions: n = 3,000, T = 12 dynamic aggregation with `balance_e = 2` increased from 8 to 9 ms; n = 5,000 slow varying-weight DR increased from 169 to 175.5 ms. Both directions repeated across the two process rounds with non-overlapping IQRs. Those observations remain in `supplemental-before-optimization-{comparison,raw-timings}.csv`.

The code was then changed to replace repeated cohort/event scans with vectorized support counting and to avoid allocating normalized RC weights/support subsets when all weights remain positive after normalization. Independent differential checks preserved selector and support decisions. The same seeded supplemental design was measured again in a new quiet ABBA run.

In the final supplemental run, no material end-to-end slowdown was detected. Slow varying-weight DR measured 169 versus 170 ms (0.6%, overlapping IQRs); the T = 12 balance_e = 2 case measured 9 versus 9 ms. The baseline balance measurement itself moved from 8 to 9 ms between runs, so the full 1 ms change should not be attributed solely to optimization. The geometric mean final/baseline time ratio across the 12 supplemental cases was 0.968.

The frozen optimized source was a working tree above `1a16c46`; exact R-file hashes and the dirty-file summary are in `optimized-source-manifest.txt`. The main 33-case measurements preceded these last two targeted optimizations; the final run remeasured all 12 supplemental cases.

| Case | Master median [Q25, Q75] | Fixed median [Q25, Q75] | Fixed / master |
|---|---:|---:|---:|
| `aggte_1k_dynamic_balance_e2` | 3.00 [3.00, 3.00] | 3.00 [3.00, 3.00] | 1.000 |
| `aggte_3k_T12_dynamic_balance_e2` | 9.00 [8.25, 9.00] | 9.00 [8.00, 9.00] | 1.000 |
| `aggte_3k_T12_dynamic_balance_e5` | 9.00 [9.00, 9.00] | 9.00 [9.00, 9.00] | 1.000 |
| `aggte_5k_dynamic_balance_e2` | 4.00 [4.00, 4.00] | 4.00 [4.00, 4.00] | 1.000 |
| `panel_1k_fast_dr_varying` | 49.50 [49.00, 50.75] | 47.00 [47.00, 48.00] | 0.949 |
| `panel_1k_fast_reg_weighted` | 19.00 [18.00, 19.00] | 17.00 [17.00, 17.00] | 0.895 |
| `panel_5k_fast_dr_varying` | 142.00 [140.50, 142.75] | 140.00 [138.50, 141.00] | 0.986 |
| `panel_5k_fast_reg_weighted` | 48.00 [48.00, 48.00] | 45.00 [45.00, 45.00] | 0.937 |
| `panel_5k_slow_dr_varying` | 169.00 [168.00, 170.75] | 170.00 [169.25, 171.00] | 1.006 |
| `panel_5k_slow_reg_weighted` | 26.50 [26.00, 27.00] | 26.00 [26.00, 26.00] | 0.981 |
| `rc_1k_fast_reg_weighted` | 29.00 [28.25, 30.00] | 27.00 [27.00, 27.75] | 0.931 |
| `unbalanced_1k_fast_reg_weighted` | 30.00 [30.00, 30.75] | 28.00 [28.00, 28.00] | 0.933 |

Final raw observations are in `supplemental-raw-timings.csv`; the complete final summary is in `supplemental-comparison.csv`.

### RC preflight microbenchmark

`support-performance.R` compares the original correct RC support helper from `1a16c46` with the optimized helper: 500 calls per timed batch, 20 batches per source/case, alternating before/after order. Positive support measured 80 → 46 microseconds per call. The zero-weight and normalization-underflow fallbacks measured 78 → 80 and 70 → 74 microseconds, respectively; the additional minimum scan has a small cost when the fast path is unavailable. All support semantics remain unchanged. These are internal-check measurements, separate from end-to-end runtime.

| Support | Before median [Q25, Q75], μs/call | After median [Q25, Q75], μs/call |
|---|---:|---:|
| positive | 80.0 [80.0, 80.0] | 46.0 [44.0, 46.0] |
| zero | 78.0 [76.0, 78.0] | 80.0 [79.5, 82.0] |
| underflow | 70.0 [68.0, 70.0] | 74.0 [74.0, 74.0] |

To reproduce the helper checks and microbenchmark:

```sh
Rscript dev/benchmarks/stress-fixes/guard-support-differential.R
Rscript dev/benchmarks/stress-fixes/support-performance.R /tmp/did-rc-support-performance.csv
```

## Guard optimization check

A fixed-seed 2,000-case comparison against the first audit implementation preserved all supplied-weight overlap/rank guard decisions, including zero weights, normalization underflow, and designs close to the rank cutoff. The RC support optimization additionally passed 4,000 decisions against an independent four-subgroup summation reference; its script is included. A proposed symmetric-Gram shortcut was rejected because its floating-point accumulation changed a near-cutoff rank decision. The final weighted path retains the original two-matrix weighted Gram calculation; the no-user-weights path retains the released unweighted guard.
