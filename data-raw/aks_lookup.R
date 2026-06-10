# data-raw/aks_lookup.R
# ===========================================================================
# Build inst/extdata/aks_lookup/aks_lookup.rds from the vendored MissAdapt
# lookup tables (the .mat files shipped alongside it for provenance).
#
# This is the canonical conversion script for the data consumed by
# edid_adaptive() / .edid_aks_core() (R/edid-adaptive.R). Re-running it must
# be a byte-stable no-op once the .rds is current (verified at the end).
#
# PROVENANCE ----------------------------------------------------------------
# Source repository : https://github.com/lsun20/MissAdapt
# Source commit     : 98d823a0818eebbec37ce7d1acf9ca0b78aee46b
#                     (obtained via `git -C <clone> rev-parse HEAD`,
#                      2024-10-21; the vendored .mat files in
#                      inst/extdata/aks_lookup/ are byte-identical, md5:
#                      policy.mat     ef296a9b370e7a7d1b5be3e3e4304b82,
#                      thresholds.mat 27807a651e1a4f0114585cf35fe46f8f,
#                      emse_corr.mat  e7e1f51ec5685edc750dbf1c0b364584)
# Also archived as  : Zenodo record 16890198 (replication package of
#                     Armstrong, Kline & Sun, "Adapting to Misspecification",
#                     Econometrica 93(6), 2025, 1981-2005)
# License           : MIT, Copyright (c) 2023 Sophie Sun (see the LICENSE
#                     file of the MissAdapt repository, reproduced in this
#                     package's inst/COPYRIGHTS). The MIT notice is also
#                     embedded as an attribute of the .rds.
#
# GRID CONVENTIONS ----------------------------------------------------------
# corr_grid : abs(tanh(seq(-3, -0.05, 0.05))), length 60, DECREASING from
#             0.99505 to 0.04996. This is the |correlation| grid indexing the
#             COLUMNS of psi_mat and the entries of st / ht / mse_lambda.
#             It is exactly the `Sigma_UO_grid` hard-coded in the authors'
#             R/calculate_adaptive_estimates.R; the SIGNED grid
#             tanh(seq(-3, -0.05, 0.05)) is stored by the authors themselves
#             as `Sigma_UO_grid` inside emse_corr.mat (checked below).
# y_grid    : the t_O (over-identification statistic) grid, 481 points on
#             [-12, 12] in steps of 0.05, stored as `y_grid` in policy.mat.
# psi_mat   : 481 x 60; psi_mat[i, j] = delta*(y_grid[i]; corr_grid[j]^2),
#             the minimax shrinkage function of AKS Theorem 1(ii).
#             ORIENTATION: ROWS = y-grid points, COLUMNS = corr-grid points.
#             Verified against the authors' own usage in
#             R/calculate_adaptive_estimates.R:
#                 psi.function <- splinefun(Sigma_UO_grid, policy$psi.mat[i,])
#             i.e. row i (a 60-vector across the corr grid) is splined for
#             each of the Ky = length(policy$y.grid) = 481 y-grid points;
#             and dimensionally (481 x 60, matching length(y_grid) = 481 and
#             length(Sigma_UO_grid) = 60, so the transpose cannot be splined
#             this way).
# st, ht    : length-60 soft- / hard-threshold lookups (thresholds.mat),
#             indexed by corr_grid.
# mse_lambda: length-60 ERM lambda lookup (`MSE_lambda_mat` in
#             emse_corr.mat), indexed by corr_grid.
#
# USAGE ---------------------------------------------------------------------
#   Rscript data-raw/aks_lookup.R        # from the package root
# Requires the R.matlab package (in Suggests). Stops loudly if any sanity
# check fails; writes inst/extdata/aks_lookup/aks_lookup.rds otherwise.
# ===========================================================================

if (!requireNamespace("R.matlab", quietly = TRUE)) {
  stop("data-raw/aks_lookup.R requires the R.matlab package (install.packages('R.matlab')).")
}

src_dir <- file.path("inst", "extdata", "aks_lookup")
out_rds <- file.path(src_dir, "aks_lookup.rds")
for (f in c("policy.mat", "thresholds.mat", "emse_corr.mat")) {
  if (!file.exists(file.path(src_dir, f))) {
    stop("Vendored MissAdapt file not found: ", file.path(src_dir, f),
         " -- run this script from the package root.")
  }
}

SOURCE_REPO   <- "https://github.com/lsun20/MissAdapt"
SOURCE_COMMIT <- "98d823a0818eebbec37ce7d1acf9ca0b78aee46b"
MIT_NOTICE <- paste(
  "MIT License",
  "",
  "Copyright (c) 2023 Sophie Sun",
  "",
  "Permission is hereby granted, free of charge, to any person obtaining a copy",
  "of this software and associated documentation files (the \"Software\"), to deal",
  "in the Software without restriction, including without limitation the rights",
  "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell",
  "copies of the Software, and to permit persons to whom the Software is",
  "furnished to do so, subject to the following conditions:",
  "",
  "The above copyright notice and this permission notice shall be included in all",
  "copies or substantial portions of the Software.",
  "",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR",
  "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,",
  "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE",
  "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER",
  "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,",
  "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE",
  "SOFTWARE.",
  sep = "\n")

# ---- (a) read the vendored .mat files -------------------------------------
policy     <- R.matlab::readMat(file.path(src_dir, "policy.mat"))
thresholds <- R.matlab::readMat(file.path(src_dir, "thresholds.mat"))
mse        <- R.matlab::readMat(file.path(src_dir, "emse_corr.mat"))

# ---- (b) convert to the structure .edid_aks_core() consumes ----------------
corr_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
tab <- list(
  y_grid     = as.numeric(policy$y.grid),          # 481 t_O grid points
  psi_mat    = unname(policy$psi.mat),             # 481 x 60, rows = y, cols = corr
  st         = as.numeric(thresholds$st.mat),      # 60 soft thresholds
  ht         = as.numeric(thresholds$ht.mat),      # 60 hard thresholds
  mse_lambda = as.numeric(mse$MSE.lambda.mat),     # 60 ERM lambdas
  corr_grid  = corr_grid                           # 60 |corr| grid (decreasing)
)

attr(tab, "provenance") <- paste0(
  "Converted from policy.mat / thresholds.mat / emse_corr.mat of the MissAdapt ",
  "replication package of Armstrong, Kline & Sun, 'Adapting to Misspecification', ",
  "Econometrica 93(6), 2025, 1981-2005 (Zenodo record 16890198). ",
  "The source .mat files are vendored byte-identically in inst/extdata/aks_lookup/. ",
  "See data-raw/aks_lookup.R for the conversion and its sanity checks.")
attr(tab, "source_repo")   <- SOURCE_REPO
attr(tab, "source_commit") <- SOURCE_COMMIT
attr(tab, "license")       <- MIT_NOTICE
attr(tab, "grid_conventions") <- c(
  corr_grid  = "abs(tanh(seq(-3, -0.05, 0.05))); length 60, decreasing 0.99505 -> 0.04996; indexes psi_mat columns and st/ht/mse_lambda entries",
  y_grid     = "t_O grid; 481 points, -12 to 12 in steps of 0.05",
  psi_mat    = "481 x 60; rows = y_grid points, columns = corr_grid points; psi_mat[i, j] = delta*(y_grid[i]; corr_grid[j]^2)")

# ---- (c) sanity checks ------------------------------------------------------
fail <- function(...) stop("aks_lookup sanity check failed: ", ..., call. = FALSE)

## dimensions and orientation (rows = y, cols = corr)
if (length(tab$y_grid) != 481L)             fail("y_grid length ", length(tab$y_grid), " != 481")
if (!identical(dim(tab$psi_mat), c(481L, 60L)))
  fail("psi_mat is ", paste(dim(tab$psi_mat), collapse = "x"),
       ", expected 481 x 60 (rows = y_grid, cols = corr_grid)")
if (length(tab$corr_grid) != 60L)           fail("corr_grid length != 60")
if (length(tab$st) != 60L || length(tab$ht) != 60L || length(tab$mse_lambda) != 60L)
  fail("st/ht/mse_lambda are not all length 60")

## the y grid is what the docs say: [-12, 12] in steps of 0.05, symmetric
if (!isTRUE(all.equal(tab$y_grid, seq(-12, 12, by = 0.05), tolerance = 1e-12)))
  fail("y_grid is not seq(-12, 12, 0.05)")
if (max(abs(tab$y_grid + rev(tab$y_grid))) != 0) fail("y_grid is not exactly symmetric")

## the corr grid convention, cross-checked against the authors' OWN stored
## grid: emse_corr.mat carries the signed grid tanh(seq(-3, -0.05, 0.05)) as
## `Sigma_UO_grid`; our corr_grid must be its absolute value.
if (!isTRUE(all.equal(abs(as.numeric(mse$Sigma.UO.grid)), tab$corr_grid, tolerance = 1e-12)))
  fail("corr_grid does not match abs(Sigma_UO_grid) stored in emse_corr.mat")
if (any(diff(tab$corr_grid) >= 0)) fail("corr_grid is not strictly decreasing")

## monotonicity: delta*(y; rho^2) is nondecreasing in y for every fixed corr
## column (it is a monotone shrinkage rule). Tolerance only for FP roundoff.
mono <- vapply(seq_len(ncol(tab$psi_mat)),
               function(j) all(diff(tab$psi_mat[, j]) >= -1e-12), logical(1L))
if (!all(mono)) fail(sum(!mono), " psi_mat columns are not monotone nondecreasing in y")

## sign convention: sign(psi) = sign(y) away from y = 0, |psi| <= |y|
## (shrinkage toward 0), and psi(0; .) = 0 up to the tabulation error.
nz <- abs(tab$y_grid) > 0.025
if (!all(sign(tab$psi_mat[nz, ]) == sign(tab$y_grid[nz])))
  fail("sign(psi_mat) != sign(y_grid) somewhere away from y = 0")
if (max(abs(tab$psi_mat) - abs(tab$y_grid)) > 1e-4)
  fail("|psi| > |y| somewhere: not a shrinkage rule?")
if (max(abs(tab$psi_mat[which.min(abs(tab$y_grid)), ])) > 1e-4)
  fail("psi(0; .) is not ~0")

## odd symmetry: delta* is odd in y *in theory*; the tabulated solver output
## satisfies it only approximately (max residual ~3.63e-4 at this commit), so
## this is a tolerance check, NOT exact equality -- do not tighten it.
odd_resid <- max(abs(tab$psi_mat + tab$psi_mat[rev(seq_along(tab$y_grid)), ]))
if (odd_resid > 1e-3) fail("odd-symmetry residual ", format(odd_resid), " > 1e-3")
message(sprintf("odd-symmetry residual max|psi(y)+psi(-y)| = %.3e (approximate, as expected)", odd_resid))

## thresholds / lambdas are positive and in their documented ranges
if (any(tab$st <= 0) || any(tab$ht <= 0) || any(tab$mse_lambda <= 0))
  fail("st/ht/mse_lambda not all positive")
if (any(tab$ht < tab$st)) fail("hard threshold below soft threshold somewhere")

## end-to-end regression against the authors' published vignette example
## (MissAdapt README; de Chaisemartin & D'Haultfoeuille 2020, Table 3 inputs):
## reproduce the interpolation inline (same spline calls as the authors'
## calculate_adaptive_estimates.R and as .edid_aks_core) and require the
## adaptive estimate the authors report (0.36 per 100, i.e. 0.0036).
YR <- 0.0026; VR <- 0.0009^2; YU <- 0.0043; VU <- 0.0014^2
VUR <- 0.7236 * sqrt(VR * VU)
YO <- YR - YU; VO <- VR - 2 * VUR + VU; VUO <- VUR - VU
tO <- YO / sqrt(VO); corr <- VUO / sqrt(VO) / sqrt(VU)
GMM <- YU - VUO / VO * YO
psi_grid <- vapply(seq_along(tab$y_grid), function(i) {
  stats::splinefun(tab$corr_grid, tab$psi_mat[i, ], method = "fmm", ties = mean)(abs(corr))
}, numeric(1L))
adaptive <- VUO / sqrt(VO) *
  stats::splinefun(tab$y_grid, psi_grid, method = "natural")(tO) + GMM
if (abs(tO - (-1.747359190350)) > 1e-9)   fail("vignette t_O mismatch: ", format(tO, digits = 12))
if (abs(corr - (-0.769619216098)) > 1e-9) fail("vignette corr mismatch: ", format(corr, digits = 12))
if (abs(adaptive - 0.003565247561) > 1e-9)
  fail("vignette adaptive estimate mismatch: ", format(adaptive, digits = 12))
if (round(100 * adaptive, 2) != 0.36) fail("vignette adaptive != 0.36 per 100")
message(sprintf("vignette check: t_O = %.4f, corr = %.4f, adaptive = %.6f (README: -1.75, -0.77, 0.0036)",
                tO, corr, adaptive))

## exact equality with the currently shipped .rds (data components; the
## attributes may legitimately differ across script revisions). Skipped on a
## first build where no .rds exists yet.
if (file.exists(out_rds)) {
  old <- readRDS(out_rds)
  for (nm in names(tab)) {
    if (!identical(unname(tab[[nm]]), unname(old[[nm]])))
      fail("component `", nm, "` differs from the currently shipped aks_lookup.rds")
  }
  message("all data components identical to the currently shipped aks_lookup.rds")
} else {
  message("no shipped aks_lookup.rds found; writing a fresh one")
}

# ---- (d) write the .rds (deterministic settings => byte-stable reruns) -----
saveRDS(tab, out_rds, version = 3L, compress = "gzip")
message("wrote ", out_rds, " (", file.size(out_rds), " bytes, md5 ",
        tools::md5sum(out_rds), ")")
