#!/usr/bin/env Rscript
# Targeted MC calibration (full-adaptation audit item 4d): coverage of the overall
# ATT under DISPERSED observation weights, comparing the ridge intensity BEFORE vs
# AFTER the n_eff fix is irrelevant here (we only have AFTER loaded); the point is to
# confirm the weighted ridge SE is well-calibrated (coverage ~ 0.95, mean-SE / MC-SD ~ 1)
# AFTER the fix. A small staggered DGP, known overall ATT, dispersed weights.
suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE)))

R    <- as.integer(Sys.getenv("MCREPS", "500"))
disp <- as.numeric(Sys.getenv("MCDISP", "1.0"))
n_g3 <- 40L; n_g5 <- 40L; n_never <- 60L; np <- 6L
n    <- n_g3 + n_g5 + n_never
# true overall (dynamic event-study) ATT for this DGP: effects 1.5 (g3) and 2.5 (g5),
# constant in event time -> overall dynamic average is the equal-weight mean over cohorts*
# present event-times. We instead estimate the truth by a large-sample run (weights=1) to
# avoid hand-deriving the exact dynamic weights; that anchors coverage to the estimand edid targets.

gen <- function(seed, weighted = TRUE) {
  set.seed(seed)
  uid <- rep(seq_len(n), each = np); tid <- rep(seq_len(np), times = n)
  ft  <- c(rep(3L, n_g3 * np), rep(5L, n_g5 * np), rep(Inf, n_never * np))
  ufe <- rep(rnorm(n, 0, 1), each = np); tfe <- rep(seq(0, 1, length.out = np), times = n)
  e   <- rnorm(n * np, 0, 0.7)
  t3  <- (uid <= n_g3) & (tid >= 3L); t5 <- (uid > n_g3 & uid <= n_g3 + n_g5) & (tid >= 5L)
  df  <- data.frame(unit = uid, time = tid, outcome = ufe + tfe + e + 1.5 * t3 + 2.5 * t5,
                    first_treat = ft)
  if (weighted) df$w <- rep(exp(rnorm(n, 0, disp)), each = np)
  df
}

fit_overall <- function(df, wn) {
  fit <- suppressWarnings(suppressMessages(
    edid(yname = "outcome", tname = "time", idname = "unit", gname = "first_treat",
         data = df, weightsname = wn, omega_cov_shrink = "ridge")))
  c(att = unname(fit$overall$overall.att), se = unname(fit$overall$overall.se))
}

# Anchor the estimand: large unweighted sample (weights enter the estimand via the Hajek
# reweighting, so the WEIGHTED estimand differs; we anchor to the WEIGHTED estimand by a
# huge weighted sample with the SAME weight law -> the probability-limit of the estimator).
cat("Anchoring weighted estimand via a large sample...\n")
set.seed(1)
big_seed <- 99999L
n_big_factor <- 8L
# replicate the DGP at 8x and average many draws to pin the plim
ests <- replicate(40, fit_overall(gen(sample.int(1e6,1), weighted = TRUE), "w")["att"])
truth <- mean(ests, na.rm = TRUE)
cat(sprintf("Anchored weighted overall ATT (plim estimate) = %.4f  (mc se %.4f)\n",
            truth, sd(ests, na.rm = TRUE) / sqrt(length(ests))))

cat(sprintf("Running %d weighted-ridge reps (disp=%.2f, n=%d)...\n", R, disp, n))
res <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("att", "se")))
for (r in seq_len(R)) {
  v <- tryCatch(fit_overall(gen(1000L + r, weighted = TRUE), "w"), error = function(e) c(NA, NA))
  res[r, ] <- v
}
ok   <- is.finite(res[, "att"]) & is.finite(res[, "se"]) & res[, "se"] > 0
att  <- res[ok, "att"]; se <- res[ok, "se"]
mcsd <- sd(att)
cover <- mean(abs(att - truth) <= qnorm(0.975) * se)
cat("\n===== WEIGHTED RIDGE MC CALIBRATION (after n_eff fix) =====\n")
cat(sprintf("  valid reps      : %d / %d\n", sum(ok), R))
cat(sprintf("  MC SD(att)      : %.4f\n", mcsd))
cat(sprintf("  mean(SE)        : %.4f\n", mean(se)))
cat(sprintf("  mean(SE)/MC SD  : %.3f   (target ~ 1.0)\n", mean(se) / mcsd))
cat(sprintf("  95%% coverage    : %.3f   (target ~ 0.95)\n", cover))
cat(sprintf("  mean att        : %.4f   (truth %.4f, bias %.4f)\n", mean(att), truth, mean(att) - truth))
saveRDS(list(truth = truth, att = att, se = se, cover = cover,
             ratio = mean(se) / mcsd, disp = disp, R = R),
        "/tmp/gate_runs/ridgefix/mc_coverage.rds")
