# Guard for the parallel cell-loop lever, the `cores` argument (formerly only options(edid_mc_cores)). The
# (g,t) cells are independent, so the forked path (parallel::mclapply, edid-fit.R) MUST be numerically
# identical to the serial path -- a future cache/RNG regression inside a forked worker would otherwise pass
# CI undetected. Fork-based, so this is a no-op on Windows (cores is ignored there) and is skipped.
#
# Worker error propagation: when a forked worker fails, parallel::mclapply returns a "try-error" object and
# the reducer at edid-fit.R re-raises it (`if (inherits(.r, "try-error")) stop(...)`), rather than silently
# dropping the cell. That reducer is verified by inspection; it is not exercised here because forcing a fork
# crash deterministically (the try-error branch is unreachable on the serial lapply path) would make the
# test brittle without adding real coverage of the production arithmetic.

test_that("cores > 1 is bit-identical to the serial path (att / se / EIF)", {
  skip_on_cran()
  skip_on_os("windows")          # mclapply forking is unavailable on Windows; cores is ignored there
  data(mpdta, package = "did")
  run <- function(k) {
    set.seed(7)
    edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
         gname = "first.treat", xformla = ~ lpop, aggregate = "none", cores = k)
  }
  f1 <- run(1L)
  f2 <- run(2L)
  expect_equal(f2$att_gt$att, f1$att_gt$att, tolerance = 1e-12)
  expect_equal(f2$att_gt$se,  f1$att_gt$se,  tolerance = 1e-12)
  expect_equal(f2$eif,        f1$eif,        tolerance = 1e-12)
})

test_that("the edid_mc_cores option still works as a session-wide default for cores", {
  skip_on_cran()
  skip_on_os("windows")
  data(mpdta, package = "did")
  old <- options(edid_mc_cores = 2L)
  on.exit(options(old))
  set.seed(7)
  fo <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
             gname = "first.treat", xformla = ~ lpop, aggregate = "none")   # cores defaults to the option
  options(old)
  set.seed(7)
  f1 <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
             gname = "first.treat", xformla = ~ lpop, aggregate = "none", cores = 1L)
  expect_equal(fo$att_gt$att, f1$att_gt$att, tolerance = 1e-12)
  expect_equal(fo$att_gt$se,  f1$att_gt$se,  tolerance = 1e-12)
})
