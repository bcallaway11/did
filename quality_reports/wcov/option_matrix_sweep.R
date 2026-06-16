## =====================================================================
## quality_reports/wcov/option_matrix_sweep.R
## Phase 5 option-matrix smoke sweep for the WEIGHTED-COVARIATE path.
## Runs edid(weightsname=, xformla=~x1) across representative combinations
## of the estimation surface and checks each runs without error and returns
## finite, non-garbage output (att finite, se finite & > 0 where applicable).
## Usage: Rscript quality_reports/wcov/option_matrix_sweep.R
## =====================================================================
suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet = TRUE)
options(edid_mc_cores = 1L)

make_panel <- function(n = 320, seed = 21) {
  set.seed(seed)
  Tn  <- 4L
  x1u <- runif(n, -2, 2)
  eta <- 1.1 * x1u + 0.7 * x1u^2 - 0.5
  P   <- exp(cbind(0, eta, 0.6 * eta)); P <- P / rowSums(P)
  gc  <- apply(P, 1L, function(p) sample(c(Inf, 2, 3), 1L, prob = p))
  alph <- rnorm(n, 0.5 * x1u, 1)
  wu   <- exp(0.8 * rnorm(n)); wu <- wu / mean(wu)
  cl   <- sample(1:20, n, replace = TRUE)              # 20 clusters
  rows <- lapply(1:Tn, function(tt) {
    ht  <- (tt - 1) * (0.5 * x1u + 0.45 * x1u^2)
    tau <- ifelse(is.finite(gc) & tt >= gc, 1, 0)
    data.frame(id = 1:n, t = tt, g = ifelse(is.finite(gc), gc, 0),
               x1 = x1u, w = wu, cl = cl, y = alph + 0.3 * tt + ht + tau + rnorm(n))
  })
  do.call(rbind, rows)
}
df <- make_panel()

ok_fit <- function(fit) {
  if (inherits(fit, "error")) return(list(ok = FALSE, msg = conditionMessage(fit)))
  ag <- fit$att_gt
  if (is.null(ag) || !nrow(ag)) return(list(ok = FALSE, msg = "no att_gt"))
  fin_att <- any(is.finite(ag$att)); fin_se <- any(is.finite(ag$se) & ag$se > 0)
  if (!fin_att) return(list(ok = FALSE, msg = "no finite att"))
  if (!fin_se)  return(list(ok = FALSE, msg = "no finite/positive se"))
  list(ok = TRUE, msg = sprintf("att[1]=%.4f se[1]=%.4f", ag$att[1], ag$se[1]))
}

run1 <- function(label, opts = list(), ...) {
  old <- options(); on.exit(options(old), add = TRUE)
  if (length(opts)) do.call(options, opts)
  args <- list(data = df, yname = "y", idname = "id", tname = "t", gname = "g",
               xformla = ~ x1, weightsname = "w", seed = 1L, ...)
  fit <- tryCatch(suppressWarnings(do.call(edid, args)), error = function(e) e)
  r <- ok_fit(fit)
  cat(sprintf("[%s] %-46s %s\n", if (r$ok) "OK " else "ERR", label, r$msg))
  list(label = label, ok = r$ok, fit = if (r$ok) fit else NULL)
}

cat("=== weighted-covariate option-matrix smoke sweep ===\n")
results <- list()
smoothers <- c(kernel = "kernel", kernel_orig = "kernel_orig", sieve = "sieve")
for (sm in names(smoothers)) for (rm in c("direct", "exp")) {
  op <- list(edid_omega_method = smoothers[[sm]])
  results[[length(results)+1]] <- run1(sprintf("eff|%s|%s|plugin",  sm, rm), op, ratio_method = rm, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, misspec_robust = FALSE, estimation_effect = FALSE)
  results[[length(results)+1]] <- run1(sprintf("eff|%s|%s|ee",      sm, rm), op, ratio_method = rm, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, misspec_robust = FALSE, estimation_effect = TRUE)
  results[[length(results)+1]] <- run1(sprintf("eff|%s|%s|misspec", sm, rm), op, ratio_method = rm, weight_scheme = "efficient", aggregate = "none", bstrap = FALSE, misspec_robust = TRUE)
}
# weight_scheme = averaged + gmm
results[[length(results)+1]] <- run1("averaged|kernel|exp",  list(edid_omega_method="kernel"), ratio_method="exp", weight_scheme="averaged",  aggregate="none", bstrap=FALSE, misspec_robust=FALSE)
results[[length(results)+1]] <- run1("gmm|kernel|exp",       list(edid_omega_method="kernel"), ratio_method="exp", weight_scheme="gmm",       aggregate="none", bstrap=FALSE, misspec_robust=FALSE)
# bs_df = "ic"
results[[length(results)+1]] <- run1("eff|sieve|exp|bsdf-ic", list(edid_omega_method="sieve"), ratio_method="exp", weight_scheme="efficient", aggregate="none", bstrap=FALSE, misspec_robust=FALSE, bs_df="ic")
# higher_order
results[[length(results)+1]] <- run1("eff|kernel|exp|higher", list(edid_omega_method="kernel"), ratio_method="exp", weight_scheme="efficient", aggregate="none", bstrap=FALSE, higher_order=TRUE)
# bootstrap (multiplier)
results[[length(results)+1]] <- run1("eff|kernel|exp|bstrap", list(edid_omega_method="kernel"), ratio_method="exp", weight_scheme="efficient", aggregate="none", bstrap=TRUE, biters=199L, misspec_robust=FALSE)
# clustering
results[[length(results)+1]] <- run1("eff|kernel|exp|cluster", list(edid_omega_method="kernel"), ratio_method="exp", weight_scheme="efficient", aggregate="none", bstrap=FALSE, clustervars="cl", misspec_robust=FALSE)
# aggregation schemes (valid: all, overall, event_study, group, calendar, none)
for (agg in c("group","event_study","calendar","overall")) {
  results[[length(results)+1]] <- run1(sprintf("eff|kernel|exp|agg-%s", agg), list(edid_omega_method="kernel"), ratio_method="exp", weight_scheme="efficient", aggregate=agg, bstrap=FALSE, misspec_robust=FALSE)
}

## ---- toolkit on weighted-covariate fits (unrestricted = PT-All, restricted = PT-Post) ----
cat("\n=== toolkit on weighted-covariate fits ===\n")
.mk <- function(...) tryCatch(suppressWarnings(edid(df, "y","id","t","g", xformla=~x1, weightsname="w",
                       weight_scheme="efficient", aggregate="none", bstrap=FALSE, seed=1L, ...)),
                       error=function(e) e)
fit_unr <- .mk(misspec_robust=FALSE)                          # PT-All (over-identified) weighted-cov
fit_res <- .mk(misspec_robust=FALSE, pt_assumption="post")   # PT-Post (just-identified) weighted-cov
tk <- function(label, expr) {
  r <- tryCatch({ force(expr); list(ok=TRUE, msg="ran") }, error=function(e) list(ok=FALSE, msg=conditionMessage(e)))
  cat(sprintf("[%s] %-24s %s\n", if (r$ok) "OK " else "ERR", label, substr(r$msg,1,70)))
  r$ok
}
tkres <- c(
  tk("edid_weights",  edid_weights(fit_unr)),
  tk("edid_sargan",   edid_sargan(fit_unr)),
  tk("edid_hausman",  edid_hausman(fit_unr, fit_res)),
  tk("edid_frontier", edid_frontier(fit_unr, fit_res)),
  tk("edid_adaptive", edid_adaptive(fit_unr, fit_res)),
  tk("summary",       summary(fit_unr)),
  tk("print",         capture.output(print(fit_unr)))
)

n_fit <- length(results); n_ok <- sum(vapply(results, function(r) isTRUE(r$ok), TRUE))
n_tk  <- length(tkres);    n_tkok <- sum(tkres)
cat(sprintf("\n[SWEEP] fits %d/%d OK | toolkit %d/%d OK | %s\n",
            n_ok, n_fit, n_tkok, n_tk,
            if (n_ok == n_fit && n_tkok == n_tk) "ALL CLEAN (PASS)" else "SOME FAILED"))
quit(status = if (n_ok == n_fit && n_tkok == n_tk) 0L else 1L)
