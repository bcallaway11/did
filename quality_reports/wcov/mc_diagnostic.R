## =====================================================================
## quality_reports/wcov/mc_diagnostic.R
## Diagnostic: is the weighted-cov calibration consistent with the UNWEIGHTED
## one? Compares coverage / (mean SE / MC SD) / bias across the 2x2 of
## {unweighted, weighted} x {misspec_robust = FALSE, TRUE} on ONE DGP
## (true ATT(g,t) = 1). If weighted ~ unweighted within each misspec level,
## the weighting is calibration-neutral (the FALSE under-coverage is the known
## plug-in-weights property, which misspec_robust=TRUE is meant to repair).
## Free-running: Rscript quality_reports/wcov/mc_diagnostic.R
## =====================================================================
suppressPackageStartupMessages({ library(pkgload); library(parallel) })
pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE)
options(edid_mc_cores = 1L)

NREPS <- 800L
NCORE <- max(1L, min(6L, parallel::detectCores() - 2L))

sim_one <- function(rep_seed) {
  set.seed(rep_seed)
  n <- 320L; Tn <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P   <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  wu   <- exp(0.8 * rnorm(n)); wu <- wu / mean(wu)
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  df <- do.call(rbind, rows)
  fit <- function(wt, ms) tryCatch(suppressWarnings(edid(df, "y","id","t","g", xformla = ~ x1,
            weightsname = wt, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE,
            seed = 1L, misspec_robust = ms, estimation_effect = TRUE)), error = function(e) NULL)
  pick <- function(f, gg, tt) { if (is.null(f)) return(c(att=NA,se=NA)); ag<-f$att_gt
    i <- which(ag$group==gg & ag$time==tt); if (length(i)) c(att=ag$att[i[1]], se=ag$se[i[1]]) else c(att=NA,se=NA) }
  cfgs <- list(uF = fit(NULL,FALSE), uT = fit(NULL,TRUE), wF = fit("w",FALSE), wT = fit("w",TRUE))
  lapply(cfgs, function(f) list(c24 = pick(f,2,4), c34 = pick(f,3,4)))
}

cat(sprintf("MC diagnostic (unw vs wt) x (misspec F/T): %d reps on %d cores...\n", NREPS, NCORE))
res <- mclapply(seq_len(NREPS) + 5000L, sim_one, mc.cores = NCORE)
res <- Filter(function(x) !is.null(x), res)

summ <- function(cfg, cell) {
  M <- do.call(rbind, lapply(res, function(r) r[[cfg]][[cell]]))
  M <- M[is.finite(M[,"att"]) & is.finite(M[,"se"]) & M[,"se"]>0, , drop=FALSE]
  if (!nrow(M)) return(sprintf("%s/%s: no finite reps", cfg, cell))
  att<-M[,"att"]; se<-M[,"se"]
  sprintf("%s/%s n=%4d bias=%+.3f MC_SD=%.3f mean_SE=%.3f ratio=%.2f cover=%.3f",
          cfg, cell, nrow(M), mean(att)-1, sd(att), mean(se), mean(se)/sd(att),
          mean(abs(att-1)/se < qnorm(0.975)))
}
out <- character(0)
for (cfg in c("uF","wF","uT","wT")) for (cell in c("c24","c34")) out <- c(out, summ(cfg, cell))
saveRDS(res, "/Users/pcostag/Documents/GitHub/did/quality_reports/wcov/mc_diagnostic.rds")
cat(paste(out, collapse = "\n"), "\n")
cat("[MC-DIAG] done. Read: uF/wF should match (plug-in-weights under-coverage); uT/wT should both improve.\n")
