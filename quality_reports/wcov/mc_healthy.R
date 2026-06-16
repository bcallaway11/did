## =====================================================================
## quality_reports/wcov/mc_healthy.R
## Calibration of the weighted-cov misspec_robust SE on a HEALTHY-overlap DGP
## (gentle propensity -> little/no trimming), to separate the adversarial-DGP
## conservativeness from any residual weighting miscalibration, AND to validate
## the sieve weight-channel fix. true ATT(g,t) = 1.
## Configs: {kernel, sieve} x {unweighted, weighted}, misspec_robust = TRUE.
## =====================================================================
suppressPackageStartupMessages({ library(pkgload); library(parallel) })
pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE)
options(edid_mc_cores = 1L)

NREPS <- 800L
NCORE <- max(1L, min(6L, parallel::detectCores() - 2L))

sim_one <- function(rep_seed) {
  set.seed(rep_seed)
  n <- 500L; Tn <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 0.45 * x1u                                        # GENTLE propensity -> healthy overlap, ~no trimming
  P   <- exp(cbind(0, eta, 0.55 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.4 * x1u, 1)
  wu   <- exp(0.7 * rnorm(n)); wu <- wu / mean(wu)         # moderate dispersion
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.4 * x1u)                          # gentle linear covariate trend
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  df <- do.call(rbind, rows)
  fit <- function(sm, wt) {
    op <- options(edid_omega_method = sm); on.exit(options(op), add = TRUE)
    tryCatch(suppressWarnings(edid(df, "y","id","t","g", xformla = ~ x1, weightsname = wt,
       weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
       misspec_robust = TRUE, estimation_effect = TRUE)), error = function(e) NULL)
  }
  pick <- function(f, gg, tt) { if (is.null(f)) return(c(att=NA,se=NA)); ag<-f$att_gt
    i <- which(ag$group==gg & ag$time==tt); if (length(i)) c(att=ag$att[i[1]], se=ag$se[i[1]]) else c(att=NA,se=NA) }
  list(ku = pick(fit("kernel", NULL), 2, 4), kw = pick(fit("kernel", "w"), 2, 4),
       su = pick(fit("sieve",  NULL), 2, 4), sw = pick(fit("sieve",  "w"), 2, 4))
}

cat(sprintf("MC healthy-overlap (kernel/sieve x unw/wt), misspec=TRUE: %d reps on %d cores...\n", NREPS, NCORE))
res <- mclapply(seq_len(NREPS) + 7000L, sim_one, mc.cores = NCORE)
res <- Filter(function(x) !is.null(x), res)
summ <- function(cfg) {
  M <- do.call(rbind, lapply(res, function(r) r[[cfg]]))
  M <- M[is.finite(M[,"att"]) & is.finite(M[,"se"]) & M[,"se"]>0, , drop=FALSE]
  if (!nrow(M)) return(sprintf("%s: no finite reps", cfg))
  att<-M[,"att"]; se<-M[,"se"]
  sprintf("%-3s n=%4d bias=%+.3f MC_SD=%.3f mean_SE=%.3f ratio=%.2f cover=%.3f",
          cfg, nrow(M), mean(att)-1, sd(att), mean(se), mean(se)/sd(att),
          mean(abs(att-1)/se < qnorm(0.975)))
}
cat(paste(vapply(c("ku","kw","su","sw"), summ, character(1)), collapse = "\n"), "\n")
cat("[MC-HEALTHY] done. ku/kw = kernel unw/wt; su/sw = sieve unw/wt. ratio~1 & cover~0.95 => well-calibrated.\n")
