#!/usr/bin/env Rscript
# ============================================================================
# Targeted MC calibration (full-adaptation audit item 4d) for the effective-n
# ridge fix -- EVENT-STUDY estimands.
#
# Coverage of the 95% CI for:
#   (i)  ES_avg = the dynamic event-study average ATT  (fit$overall$overall.*)
#   (ii) ES(0)  = the on-impact effect, event time 0    (fit$event_study att.egt/se.egt at egt==0)
# under DISPERSED time-invariant (per-unit) weights, in two regimes:
#   - Schmitt-like : log-normal weights s=0.80 -> per-cell Kish ESS ~ 54%
#   - Bailey-like  : log-normal weights s=0.90 -> per-cell Kish ESS ~ 48%
# comparing the ridge/LW intensity denominator:
#   - BEFORE (raw n) : n_eff_edid overridden to return the raw active count (pre-fix behavior)
#   - AFTER  (n_eff) : the implemented per-cell Kish ESS fix
#
# PASS criteria (FAIL loudly on regression):
#   * AFTER coverage >= BEFORE coverage - tol_mc (no coverage regression), for BOTH
#     estimands in BOTH regimes;
#   * AFTER mean(SE) >= BEFORE mean(SE) - tol_mc (ridge over-shrinks weights LESS under
#     raw-n -> the fix should not shrink the reported SE), i.e. the fix moves SE in the
#     calibration-improving direction or holds it;
#   * the gain is retained: AFTER mean(SE)/MC-SD is no worse than BEFORE by more than tol_mc.
# tol_mc accounts for Monte-Carlo noise at this rep count (2 * binomial se on coverage).
# ============================================================================
suppressWarnings(suppressMessages(pkgload::load_all(
  "/Users/pcostag/Documents/GitHub/did", quiet = TRUE)))

R   <- as.integer(Sys.getenv("MCREPS", "500"))
n_g3 <- 40L; n_g5 <- 40L; n_never <- 60L; np <- 6L
n    <- n_g3 + n_g5 + n_never
kish <- function(w) (sum(w)^2) / (length(w) * sum(w^2))

# Two dispersion regimes calibrated to Kish targets (see calib_kish.R):
regimes <- list(
  Schmitt = list(s = 0.82, target_kish = 0.54),
  Bailey  = list(s = 0.93, target_kish = 0.48)
)

gen <- function(seed, s) {
  set.seed(seed)
  uid <- rep(seq_len(n), each = np); tid <- rep(seq_len(np), times = n)
  ft  <- c(rep(3L, n_g3 * np), rep(5L, n_g5 * np), rep(Inf, n_never * np))
  ufe <- rep(rnorm(n, 0, 1), each = np); tfe <- rep(seq(0, 1, length.out = np), times = n)
  e   <- rnorm(n * np, 0, 0.7)
  t3  <- (uid <= n_g3) & (tid >= 3L)
  t5  <- (uid > n_g3 & uid <= n_g3 + n_g5) & (tid >= 5L)
  wunit <- exp(rnorm(n, 0, s))
  list(df = data.frame(unit = uid, time = tid,
                       outcome = ufe + tfe + e + 1.5 * t3 + 2.5 * t5,
                       first_treat = ft, w = rep(wunit, each = np)),
       kish = kish(wunit))
}

# Returns c(es_avg_att, es_avg_se, es0_att, es0_se) for one fit.
fit_es <- function(df) {
  fit <- suppressWarnings(suppressMessages(
    edid(yname = "outcome", tname = "time", idname = "unit", gname = "first_treat",
         data = df, weightsname = "w", omega_cov_shrink = "ridge", aggregate = "all")))
  es     <- fit$event_study
  i0     <- which(es$egt == 0)
  c(es_avg_att = unname(fit$overall$overall.att),
    es_avg_se  = unname(fit$overall$overall.se),
    es0_att    = unname(es$att.egt[i0]),
    es0_se     = unname(es$se.egt[i0]))
}

run_mc <- function(s, seeds) {
  M <- matrix(NA_real_, length(seeds), 4,
              dimnames = list(NULL, c("es_avg_att","es_avg_se","es0_att","es0_se")))
  kvec <- numeric(length(seeds))
  for (i in seq_along(seeds)) {
    g <- gen(seeds[i], s)
    kvec[i] <- g$kish
    M[i, ] <- tryCatch(fit_es(g$df), error = function(e) rep(NA_real_, 4))
  }
  attr(M, "kish") <- mean(kvec, na.rm = TRUE)
  M
}

# ---- Anchor the weighted estimands per regime (large-sample plim, after-fix) -------------
anchor_regime <- function(s, nanchor = 60L) {
  set.seed(1)
  A <- vapply(seq_len(nanchor),
              function(j) fit_es(gen(sample.int(1e6, 1), s)$df)[c("es_avg_att","es0_att")],
              numeric(2))
  c(es_avg = mean(A["es_avg_att", ], na.rm = TRUE),
    es0    = mean(A["es0_att", ], na.rm = TRUE),
    es_avg_mcse = sd(A["es_avg_att", ], na.rm = TRUE) / sqrt(ncol(A)),
    es0_mcse    = sd(A["es0_att", ], na.rm = TRUE)   / sqrt(ncol(A)))
}

# ---- n_eff override machinery (BEFORE = raw count) ---------------------------------------
ns   <- asNamespace("did")
orig <- get("n_eff_edid", ns)
raw_neff <- function(w, active_mask, n_full = length(active_mask)) n_full

summarize <- function(M, anchor_avg, anchor_0) {
  ok_a <- is.finite(M[,"es_avg_att"]) & is.finite(M[,"es_avg_se"]) & M[,"es_avg_se"] > 0
  ok_0 <- is.finite(M[,"es0_att"])    & is.finite(M[,"es0_se"])    & M[,"es0_se"]    > 0
  f <- function(att, se, anch) {
    mcsd <- sd(att)
    cov  <- mean(abs(att - anch) <= qnorm(0.975) * se)
    c(reps = length(att), mcsd = mcsd, meanSE = mean(se),
      ratio = mean(se) / mcsd, cover = cov, bias = unname(mean(att) - anch))
  }
  list(es_avg = f(M[ok_a,"es_avg_att"], M[ok_a,"es_avg_se"], anchor_avg),
       es0    = f(M[ok_0,"es0_att"],    M[ok_0,"es0_se"],    anchor_0))
}

OUT <- list()
for (rn in names(regimes)) {
  s <- regimes[[rn]]$s
  seeds <- 20000L + seq_len(R)
  cat(sprintf("\n########## REGIME %s (s=%.2f, Kish target ~%.0f%%) ##########\n",
              rn, s, 100 * regimes[[rn]]$target_kish))
  anch <- anchor_regime(s)
  cat(sprintf("  anchored weighted estimands: ES_avg=%.4f (mcse %.4f)  ES(0)=%.4f (mcse %.4f)\n",
              anch["es_avg"], anch["es_avg_mcse"], anch["es0"], anch["es0_mcse"]))

  # AFTER (fix as implemented)
  assignInNamespace("n_eff_edid", orig, ns)
  M_after  <- run_mc(s, seeds)
  # BEFORE (raw-n override)
  assignInNamespace("n_eff_edid", raw_neff, ns)
  M_before <- run_mc(s, seeds)
  assignInNamespace("n_eff_edid", orig, ns)   # restore

  realized_kish <- attr(M_after, "kish")
  s_before <- summarize(M_before, anch["es_avg"], anch["es0"])
  s_after  <- summarize(M_after,  anch["es_avg"], anch["es0"])

  pr <- function(tag, x) cat(sprintf(
    "    %-7s %-7s reps=%d MCsd=%.4f meanSE=%.4f meanSE/MCsd=%.3f cover95=%.3f bias=%+.4f\n",
    tag[1], tag[2], x["reps"], x["mcsd"], x["meanSE"], x["ratio"], x["cover"], x["bias"]))
  cat(sprintf("  realized mean per-rep panel Kish = %.3f\n", realized_kish))
  cat("  --- ES_avg ---\n")
  pr(c("BEFORE","ES_avg"), s_before$es_avg); pr(c("AFTER","ES_avg"), s_after$es_avg)
  cat("  --- ES(0) ---\n")
  pr(c("BEFORE","ES(0)"),  s_before$es0);    pr(c("AFTER","ES(0)"),  s_after$es0)

  OUT[[rn]] <- list(s = s, kish = realized_kish, anchor = anch,
                    before = s_before, after = s_after)
}

saveRDS(OUT, "/tmp/gate_runs/ridgefix/mc_es_coverage.rds")
cat("\nSaved /tmp/gate_runs/ridgefix/mc_es_coverage.rds\n")
