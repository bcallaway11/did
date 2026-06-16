## =====================================================================
## quality_reports/wcov/mc_calibration.R
## Phase 5d MC calibration of the WEIGHTED-COVARIATE SE.
## DGP: constant treatment effect tau = 1 => true group-time ATT(g,t) = 1
## for treated post-periods, regardless of the (mean-1) obs weights. Checks
##   mean(analytic SE) ~= MC SD(att_hat)   (calibration ratio ~ 1)
##   coverage of the 95% CI ~= 0.95
## for the FD-validated weighted-cov config (efficient, misspec_robust=FALSE,
## estimation_effect=TRUE), on cells (g=2,t=4) and (g=3,t=4).
## Free-running: Rscript quality_reports/wcov/mc_calibration.R  (writes mc_calibration.rds + log)
## =====================================================================
suppressPackageStartupMessages({ library(pkgload); library(parallel) })
pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE)
options(edid_mc_cores = 1L)

NREPS <- 2000L
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
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)            # constant effect 1 => true ATT(g,t)=1 post
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  df <- do.call(rbind, rows)
  fit <- tryCatch(suppressWarnings(edid(df, "y","id","t","g", xformla = ~ x1, weightsname = "w",
            weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, seed = 1L,
            misspec_robust = FALSE, estimation_effect = TRUE)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ag <- fit$att_gt
  pick <- function(gg, tt) { i <- which(ag$group == gg & ag$time == tt); if (length(i)) c(att = ag$att[i[1]], se = ag$se[i[1]]) else c(att = NA, se = NA) }
  list(c24 = pick(2, 4), c34 = pick(3, 4))
}

cat(sprintf("MC weighted-cov calibration: %d reps on %d cores...\n", NREPS, NCORE))
res <- mclapply(seq_len(NREPS) + 1000L, sim_one, mc.cores = NCORE)
res <- Filter(function(x) !is.null(x), res)

summarize <- function(key) {
  M <- do.call(rbind, lapply(res, function(r) r[[key]]))
  M <- M[is.finite(M[, "att"]) & is.finite(M[, "se"]) & M[, "se"] > 0, , drop = FALSE]
  if (!nrow(M)) return(NULL)
  att <- M[, "att"]; se <- M[, "se"]; truth <- 1
  cov <- mean(abs(att - truth) / se < qnorm(0.975))
  list(key = key, nrep = nrow(M), mean_att = mean(att), bias = mean(att) - truth,
       mc_sd = sd(att), mean_se = mean(se), ratio = mean(se) / sd(att), coverage = cov)
}
out <- Filter(function(x) !is.null(x), lapply(c("c24","c34"), summarize))
saveRDS(out, "/Users/pcostag/Documents/GitHub/did/quality_reports/wcov/mc_calibration.rds")
for (s in out) cat(sprintf("[MC %s] nrep=%d mean_att=%.4f bias=%+.4f MC_SD=%.4f mean_SE=%.4f ratio=%.3f coverage=%.3f\n",
                           s$key, s$nrep, s$mean_att, s$bias, s$mc_sd, s$mean_se, s$ratio, s$coverage))
cat("[MC] done.\n")
