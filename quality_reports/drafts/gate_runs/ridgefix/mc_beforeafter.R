#!/usr/bin/env Rscript
# Before/after MC: same 400 weighted-ridge reps run with the n_eff fix ON (after) and
# with n_eff forced back to the RAW count (before), by overriding n_eff_edid in the
# package namespace. Confirms the fix improves weighted calibration (raw-n under-
# regularizes -> noisier weights -> the plug-in SE under-covers more).
suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE)))

R    <- as.integer(Sys.getenv("MCREPS", "400"))
disp <- 1.0
n_g3 <- 40L; n_g5 <- 40L; n_never <- 60L; np <- 6L
n    <- n_g3 + n_g5 + n_never

gen <- function(seed) {
  set.seed(seed)
  uid <- rep(seq_len(n), each = np); tid <- rep(seq_len(np), times = n)
  ft  <- c(rep(3L, n_g3 * np), rep(5L, n_g5 * np), rep(Inf, n_never * np))
  ufe <- rep(rnorm(n, 0, 1), each = np); tfe <- rep(seq(0, 1, length.out = np), times = n)
  e   <- rnorm(n * np, 0, 0.7)
  t3  <- (uid <= n_g3) & (tid >= 3L); t5 <- (uid > n_g3 & uid <= n_g3 + n_g5) & (tid >= 5L)
  data.frame(unit = uid, time = tid, outcome = ufe + tfe + e + 1.5 * t3 + 2.5 * t5,
             first_treat = ft, w = rep(exp(rnorm(n, 0, disp)), each = np))
}
fit_overall <- function(df) {
  fit <- suppressWarnings(suppressMessages(
    edid(yname = "outcome", tname = "time", idname = "unit", gname = "first_treat",
         data = df, weightsname = "w", omega_cov_shrink = "ridge")))
  c(att = unname(fit$overall$overall.att), se = unname(fit$overall$overall.se))
}
run_mc <- function() {
  res <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("att", "se")))
  for (r in seq_len(R)) res[r, ] <- tryCatch(fit_overall(gen(1000L + r)), error = function(e) c(NA, NA))
  ok <- is.finite(res[,"att"]) & is.finite(res[,"se"]) & res[,"se"] > 0
  res[ok, , drop = FALSE]
}

# anchor estimand (n_eff fix doesn't move the point estimand's plim materially; use the after-fix plim)
set.seed(1)
anchor <- mean(replicate(40, fit_overall(gen(sample.int(1e6,1)))["att"]), na.rm = TRUE)

# AFTER (n_eff fix as implemented)
after <- run_mc()

# BEFORE: override n_eff_edid to return the raw count (length(active_mask) or n_full),
# i.e. ignore weights -> reproduces the pre-fix raw-n intensity.
ns <- asNamespace("did")
orig <- get("n_eff_edid", ns)
raw_neff <- function(w, active_mask, n_full = length(active_mask)) n_full
assignInNamespace("n_eff_edid", raw_neff, ns)
before <- run_mc()
assignInNamespace("n_eff_edid", orig, ns)   # restore

summ <- function(M, lab) {
  att <- M[,"att"]; se <- M[,"se"]; mcsd <- sd(att)
  cov <- mean(abs(att - anchor) <= qnorm(0.975) * se)
  cat(sprintf("  %-8s reps=%d  MCsd=%.4f  meanSE=%.4f  meanSE/MCsd=%.3f  cover95=%.3f\n",
              lab, nrow(M), mcsd, mean(se), mean(se)/mcsd, cov))
  c(mcsd = mcsd, meanSE = mean(se), ratio = mean(se)/mcsd, cover = cov)
}
cat(sprintf("\n=== WEIGHTED RIDGE before/after (disp=%.1f, n=%d, anchor=%.4f) ===\n", disp, n, anchor))
b <- summ(before, "BEFORE"); a <- summ(after, "AFTER")
cat(sprintf("\n  delta meanSE = %+.4f   delta cover = %+.3f   (fix should RAISE both toward calibration)\n",
            a["meanSE"] - b["meanSE"], a["cover"] - b["cover"]))
saveRDS(list(anchor = anchor, before = b, after = a), "/tmp/gate_runs/ridgefix/mc_beforeafter.rds")
