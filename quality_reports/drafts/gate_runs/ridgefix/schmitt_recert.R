## ============================================================================
## SCHMITT RECOVERY -- re-certify weighted (w_base, cluster h) with the FIXED
## (n_eff Kish-ESS) ridge, on the LIVE patched did working tree.
##
## Task structure:
##  (1) PLUG-IN (none) over-id Hausman (U=PT-Post, R=PT-All), gatekeeper that is
##      NEVER ridged -> exact certified window from e=0. Apply effective-n
##      over-id-size caveat (over-id over-rejects at small n_eff).
##  (2) Recovered weighted gain (ARE) on the post-fix ridge default, on the
##      certified window.
##  (3) Schmitt n_eff (Kish ESS) -- the over-id-size context.
##
## Cluster = hospital h (paper level; id == h). Window grown from e=0.
## ============================================================================
source("/tmp/gate_runs/helpers.R")
options(edid_mc_cores = 1L)
SEED <- 20260612
suppressMessages(library(data.table))

p <- as.data.table(readRDS("/tmp/gate_runs/rangel-60/panel_bal.rds"))
p[, w_base := dis_tot[year == min(year)][1], by = id]   # time-invariant 1996 baseline volume
dat <- as.data.frame(p)
CL <- "h"; MINE <- 0L; MAXE <- 5L; WN <- "w_base"
N_units <- length(unique(dat[["id"]]))

cat(sprintf("panel: %d rows, %d units, %d clusters (h); weight=%s\n",
            nrow(dat), N_units, length(unique(dat[["h"]])), WN))

## ---- Schmitt n_eff (global Kish ESS of the mean-1 weights) ----
wv <- p[, .(w = w_base[1]), by = id]$w
wbar <- wv / mean(wv)                       # mean-1 normalization (as edid does)
kish_global <- sum(wbar)^2 / sum(wbar^2)
cat(sprintf("\n[Schmitt n_eff] global Kish ESS = %.1f of %d units (n/n_eff = %.3f); CV(w)=%.3f\n",
            kish_global, N_units, N_units / kish_global, sd(wbar)))

## ============================================================================
## (1) PLUG-IN (none) OVER-ID GATEKEEPER -- window grow from e=0 (NEVER ridge)
## ============================================================================
cat("\n===== (1) PLUG-IN (none) over-id Hausman, WEIGHTED, window grow from e=0 =====\n")
fA <- fit_edid("A_none", data = dat, yname = "y", idname = "id", tname = "year",
               gname = "gvar", pt_assumption = "all", weight_scheme = "efficient",
               weightsname = WN, clustervars = CL, omega_cov_shrink = "none",
               estimation_effect = FALSE, cband = FALSE, seed = SEED)
fp_none <- fit_edid("post_none", data = dat, yname = "y", idname = "id", tname = "year",
                    gname = "gvar", pt_assumption = "post", weightsname = WN,
                    clustervars = CL, omega_cov_shrink = "none",
                    estimation_effect = FALSE, cband = FALSE, seed = SEED)

res0 <- data.table(emax = integer(), df = integer(), chi2 = numeric(),
                   p = numeric(), unstable = character())
for (emax in 0:14) {
  hz <- tryCatch(edid_hausman(fp_none, fA, parameter = "event_study", e_set = 0:emax),
                 error = function(e) NULL)
  if (!is.null(hz)) {
    res0 <- rbind(res0, data.table(emax = emax, df = hz$df, chi2 = hz$statistic,
      p = hz$p_value,
      unstable = ifelse(is.null(hz$leg_unstable), "NA", as.character(hz$leg_unstable))))
    cat(sprintf("  e[0,%2d] df=%2d chi2=%7.2f p=%.4f %s unstable=%s\n",
                emax, hz$df, hz$statistic, hz$p_value,
                ifelse(hz$p_value > 0.05, "PASS", "fail"),
                ifelse(is.null(hz$leg_unstable), "NA", hz$leg_unstable)))
  }
}
cert0 <- -1L
ok0 <- res0[order(emax)]
for (i in seq_len(nrow(ok0))) { if (ok0$p[i] > 0.05) cert0 <- ok0$emax[i] else break }
cat(sprintf("  => PLUG-IN (none) CERTIFIED WINDOW (contiguous from e=0) = e[0,%d]\n", cert0))

## ============================================================================
## (2) Recovered weighted gain -- POST-FIX ridge default + the ladder
## ============================================================================
cat("\n===== (2) Recovered weighted ladder (POST-FIX n_eff ridge), e[0,5], cluster=h =====\n")
ef <- function(label, ocs, ee) fit_edid(label, data = dat, yname = "y", idname = "id",
  tname = "year", gname = "gvar", pt_assumption = "all", weight_scheme = "efficient",
  weightsname = WN, clustervars = CL, omega_cov_shrink = ocs, estimation_effect = ee,
  cband = FALSE, seed = SEED)
edid_win <- function(fit) {
  ag <- aggte_edid(fit, type = "dynamic", min_e = MINE, max_e = MAXE, na.rm = TRUE)
  i0 <- which(abs(ag$egt - 0) < 1e-9)
  list(esavg_att = ag$overall.att, esavg_se = ag$overall.se,
       es0_att = ag$att.egt[i0], es0_se = ag$se.egt[i0],
       egt = ag$egt, att = ag$att.egt, se = ag$se.egt)
}
fAp <- ef("Ap_none_EE", "none", TRUE)
fC  <- ef("C_ridge_EE", "ridge", TRUE)
fB  <- ef("B_lw_EE", "ledoit_wolf", TRUE)
A <- edid_win(fA); Ap <- edid_win(fAp); B <- edid_win(fB); C <- edid_win(fC)

## CS anchors (weighted)
cs_nev <- did::att_gt(yname = "y", tname = "year", idname = "id", gname = "gvar0",
  data = dat, control_group = "nevertreated", weightsname = WN, clustervars = CL,
  bstrap = FALSE, cband = FALSE)
cs_nyt <- did::att_gt(yname = "y", tname = "year", idname = "id", gname = "gvar0",
  data = dat, control_group = "notyettreated", weightsname = WN, clustervars = CL,
  bstrap = FALSE, cband = FALSE)
cw <- function(cs) { ag <- did::aggte(cs, type = "dynamic", na.rm = TRUE, bstrap = FALSE,
  cband = FALSE, min_e = MINE, max_e = MAXE)
  i0 <- which(abs(ag$egt - 0) < 1e-9)
  list(esavg_att = ag$overall.att, esavg_se = ag$overall.se,
       es0_att = ag$att.egt[i0], es0_se = ag$se.egt[i0]) }
cn <- cw(cs_nev); cy <- cw(cs_nyt)

are <- function(anc_se, eff_se) (anc_se / eff_se)^2
lab <- c("A none", "A+ none+EE", "C ridge+EE (default)", "B lw+EE")
V <- list(A, Ap, C, B); names(V) <- lab

cat("\n--- certified-window e[0,5] ladder (EE-corrected SEs) ---\n")
for (k in lab) cat(sprintf("%-22s ES_avg=%8.4f (SE %7.4f) | ES(0)=%8.4f (SE %7.4f)\n",
  k, V[[k]]$esavg_att, V[[k]]$esavg_se, V[[k]]$es0_att, V[[k]]$es0_se))
cat(sprintf("%-22s ES_avg=%8.4f (SE %7.4f) | ES(0)=%8.4f (SE %7.4f)\n",
  "CS-never (anchor)", cn$esavg_att, cn$esavg_se, cn$es0_att, cn$es0_se))
cat(sprintf("%-22s ES_avg=%8.4f (SE %7.4f) | ES(0)=%8.4f (SE %7.4f)\n",
  "CS-notyet", cy$esavg_att, cy$esavg_se, cy$es0_att, cy$es0_se))

cat("\n--- ARE vs CS-never (variance ratio CS/efficient) ---\n")
for (k in lab) cat(sprintf("%-22s ARE_avg=%.2f  ARE_0=%.2f  (vs notyet: avg %.2f e0 %.2f)\n",
  k, are(cn$esavg_se, V[[k]]$esavg_se), are(cn$es0_se, V[[k]]$es0_se),
     are(cy$esavg_se, V[[k]]$esavg_se), are(cy$es0_se, V[[k]]$es0_se)))

## per-cell ridge intensity diagnostics (vanishing-lambda confirm under n_eff)
cat("\n--- ridge lambda / cond# (post-fix) ---\n")
diag <- tryCatch({
  qc <- fC$qc_per_cell
  if (!is.null(qc)) qc else NULL
}, error = function(e) NULL)
condmax <- tryCatch(max(fC$cell_cond, na.rm = TRUE), error = function(e) NA)
cat(sprintf("ridge cell cond# max (if available)= %s\n", as.character(condmax)))

saveRDS(list(kish_global = kish_global, N_units = N_units, n_over_neff = N_units/kish_global,
             res0 = res0, cert0 = cert0,
             A = A, Ap = Ap, B = B, C = C, cn = cn, cy = cy,
             are = list(
               A  = c(avg = are(cn$esavg_se, A$esavg_se),  e0 = are(cn$es0_se, A$es0_se)),
               Ap = c(avg = are(cn$esavg_se, Ap$esavg_se), e0 = are(cn$es0_se, Ap$es0_se)),
               C  = c(avg = are(cn$esavg_se, C$esavg_se),  e0 = are(cn$es0_se, C$es0_se)),
               B  = c(avg = are(cn$esavg_se, B$esavg_se),  e0 = are(cn$es0_se, B$es0_se)))),
        "/tmp/gate_runs/ridgefix/schmitt_recert.rds")
cat("\nDONE schmitt_recert\n")
