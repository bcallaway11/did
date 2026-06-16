#!/usr/bin/env Rscript
# Byte-identity self-check: compare pre/post fingerprints cell-by-cell.
# Unweighted configs must be byte-identical (max|delta| == 0, NOT just ~tol).
# Weighted ridge/LW configs are expected to MOVE (reported, not asserted == 0).
pre  <- readRDS("/tmp/gate_runs/ridgefix/baseline.rds")$results
post <- readRDS("/tmp/gate_runs/ridgefix/postfix.rds")$results

is_weighted <- function(k) startsWith(k, "wnocov_") || startsWith(k, "wcov_")

cell_delta <- function(a, b) {
  if (is.null(a$cells) || is.null(b$cells)) {
    # fall back to overall att/se
    return(max(abs(c(a$att - b$att, a$se - b$se)), na.rm = TRUE))
  }
  ca <- a$cells; cb <- b$cells
  if (!identical(dim(ca), dim(cb))) return(Inf)
  cols <- intersect(c("att", "se", "lambda", "cond"), names(ca))
  d <- 0
  for (cc in cols) {
    va <- ca[[cc]]; vb <- cb[[cc]]
    dd <- abs(va - vb)
    dd[is.na(va) & is.na(vb)] <- 0          # NA==NA -> 0
    d <- max(d, suppressWarnings(max(dd, na.rm = TRUE)))
  }
  d
}

cat("=== Byte-identity self-check (cell-level max|delta| over att/se/lambda/cond) ===\n")
unw_max <- 0
for (k in names(pre)) {
  if (!is.null(post[[k]])) {
    d <- cell_delta(pre[[k]], post[[k]])
    tag <- if (is_weighted(k)) "WEIGHTED (expected to move)" else "unweighted (must be 0)"
    if (!is_weighted(k)) unw_max <- max(unw_max, d)
    cat(sprintf("  %-26s max|delta|=%.3e   [%s]\n", k, d, tag))
  }
}
cat(sprintf("\n>>> UNWEIGHTED max|delta| across all cells = %.3e  (PASS iff == 0)\n", unw_max))
cat(sprintf(">>> BYTE-IDENTICAL UNWEIGHTED: %s\n", if (unw_max == 0) "YES" else "NO -- STOP"))
