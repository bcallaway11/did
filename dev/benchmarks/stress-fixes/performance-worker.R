# Invoked by performance.R. Keep package loading/data creation outside timings.
args <- commandArgs(trailingOnly = TRUE)
package_path <- normalizePath(args[[1L]])
label <- args[[2L]]
output_path <- args[[3L]]
round <- as.integer(args[[4L]])
repetitions <- as.integer(args[[5L]])
suite <- if (length(args) > 5L) args[[6L]] else "main"
options(warn = 1)
pkgload::load_all(package_path, quiet = TRUE, export_all = TRUE)
stopifnot(utils::packageVersion("DRDID") == package_version("1.4.0"))
set.seed(20261028)
make_panel <- function(n, periods = 6L) {
  X <- rnorm(n); X2 <- rnorm(n); a <- rnorm(n)
  groups <- if (periods <= 6L) c(0, 3, 4, 5, 6) else c(0, seq.int(3L, periods))
  probabilities <- if (periods <= 6L) c(.35, .16, .16, .16, .17) else
    c(.35, rep(.65 / (length(groups) - 1L), length(groups) - 1L))
  G <- sample(groups, n, TRUE, probabilities)
  d <- expand.grid(id = seq_len(n), t = seq_len(periods))
  d$G <- G[d$id]; d$X <- X[d$id]; d$X2 <- X2[d$id]
  d$w <- .5 + (d$id %% 7) / 3
  d$Y <- a[d$id] + .2 * d$t + (.3 + .15 * d$t) * d$X +
    .2 * d$X2 + (d$G > 0 & d$t >= d$G) + rnorm(nrow(d))
  d
}
d1 <- make_panel(1000L); d5 <- make_panel(5000L)
u1 <- d1[-sample(nrow(d1), 800L), ]
u5 <- d5[-sample(nrow(d5), 4000L), ]
r1 <- d1; r1$id <- seq_len(nrow(r1))
r5 <- d5; r5$id <- seq_len(nrow(r5))
pd <- make_panel(60L, 4L)
fit <- function(d, faster = TRUE, covariates = TRUE, method = "dr",
                structure = "panel", bootstrap = FALSE, weighted = FALSE,
                fix_weights = NULL) {
  suppressMessages(suppressWarnings(att_gt(
    "Y", "t", "id", "G", data = d,
    xformla = if (covariates) ~X + X2 else NULL,
    weightsname = if (weighted) "w" else NULL, fix_weights = fix_weights,
    est_method = method, faster_mode = faster,
    panel = structure != "rc", allow_unbalanced_panel = structure == "unbalanced",
    bstrap = bootstrap, cband = bootstrap, biters = 499L,
    anticipation = 0
  )))
}
# Each function exercises the same public entry point on both sources.
cases <- list(
  panel_1k_fast_dr_x = function() fit(d1),
  panel_5k_fast_dr_x = function() fit(d5),
  panel_1k_slow_dr_x = function() fit(d1, faster = FALSE),
  panel_5k_slow_dr_x = function() fit(d5, faster = FALSE),
  panel_1k_fast_dr_nox = function() fit(d1, covariates = FALSE),
  panel_5k_fast_dr_nox = function() fit(d5, covariates = FALSE),
  panel_5k_slow_dr_nox = function() fit(d5, faster = FALSE, covariates = FALSE),
  panel_1k_fast_reg_x = function() fit(d1, method = "reg"),
  panel_5k_fast_reg_x = function() fit(d5, method = "reg"),
  panel_5k_slow_reg_x = function() fit(d5, method = "reg", faster = FALSE),
  panel_5k_fast_dr_weighted = function() fit(d5, weighted = TRUE),
  panel_5k_slow_dr_weighted = function() fit(d5, weighted = TRUE, faster = FALSE),
  unbalanced_1k_fast_dr_x = function() fit(u1, structure = "unbalanced"),
  unbalanced_5k_fast_dr_x = function() fit(u5, structure = "unbalanced"),
  unbalanced_5k_slow_dr_x = function() fit(u5, structure = "unbalanced", faster = FALSE),
  rc_1k_fast_dr_x = function() fit(r1, structure = "rc"),
  rc_5k_fast_dr_x = function() fit(r5, structure = "rc"),
  rc_5k_slow_dr_x = function() fit(r5, structure = "rc", faster = FALSE),
  panel_1k_fast_dr_boot = function() fit(d1, bootstrap = TRUE),
  panel_5k_fast_dr_boot = function() fit(d5, bootstrap = TRUE),
  unbalanced_1k_fast_dr_boot = function() fit(u1, structure = "unbalanced", bootstrap = TRUE),
  rc_1k_fast_dr_boot = function() fit(r1, structure = "rc", bootstrap = TRUE)
)
# Reuse pre-fit objects to time aggregation independently of att_gt().
mp1 <- fit(d1); mp5 <- fit(d5); mpu <- fit(u1, structure = "unbalanced")
cases$aggte_1k_dynamic <- function() suppressWarnings(aggte(mp1, type = "dynamic"))
cases$aggte_5k_dynamic <- function() suppressWarnings(aggte(mp5, type = "dynamic"))
cases$aggte_5k_simple <- function() suppressWarnings(aggte(mp5, type = "simple"))
cases$aggte_1k_dynamic_boot <- function() suppressWarnings(aggte(mp1, type = "dynamic", bstrap = TRUE, cband = TRUE, biters = 499L))
cases$aggte_unbalanced_dynamic <- function() suppressWarnings(aggte(mpu, type = "dynamic"))
cases$conditional_pretest_60 <- function() {
  invisible(capture.output(suppressMessages(suppressWarnings(
    conditional_did_pretest("Y", "t", "id", "G", xformla = ~X,
                           data = pd, biters = 99L)
  ))))
}
# Low-level guard probes isolate the added support/weighted-Gram work.
x_guard <- cbind(1, rnorm(5000L), rnorm(5000L))
y_guard <- rep(c(0, 1), length.out = 5000L)
w_guard <- rep(1, 5000L); w_vary <- rep(c(.5, 1, 2), length.out = 5000L)
overlap_fn <- getFromNamespace("overlap_check_fail", "did")
rank_fn <- getFromNamespace("rcond_check_fail", "did")
supports_weights <- "i.weights" %in% names(formals(overlap_fn))
call_guard <- function(fun, x, y, intercept = FALSE, w = NULL) {
  if (supports_weights) fun(x, y, intercept_only = intercept, i.weights = w)
  else fun(x, y, intercept_only = intercept)
}
cases$guard_overlap_5k <- function() for (i in seq_len(100L)) call_guard(overlap_fn, x_guard, y_guard)
cases$guard_rank_5k <- function() for (i in seq_len(100L)) call_guard(rank_fn, x_guard, y_guard)
cases$guard_overlap_intercept_5k <- function() for (i in seq_len(1000L)) call_guard(overlap_fn, x_guard[, 1L, drop = FALSE], y_guard, TRUE)
cases$guard_rank_intercept_5k <- function() for (i in seq_len(1000L)) call_guard(rank_fn, x_guard[, 1L, drop = FALSE], y_guard, TRUE)
cases$guard_rank_weighted_5k <- function() for (i in seq_len(100L)) call_guard(rank_fn, x_guard, y_guard, w = w_vary)
# Supplemental cases cover the changed support-balancing and forced-RC branches.
if (suite == "supplemental") {
  d12 <- make_panel(3000L, 12L)
  mp12 <- fit(d12)
  vary <- d1; vary$w <- vary$w * exp(.1 * vary$t * sin(vary$id))
  vary5 <- d5; vary5$w <- vary5$w * exp(.1 * vary5$t * sin(vary5$id))
  cases <- list(
    aggte_1k_dynamic_balance_e2 = function() suppressWarnings(aggte(mp1, type = "dynamic", balance_e = 2L)),
    aggte_5k_dynamic_balance_e2 = function() suppressWarnings(aggte(mp5, type = "dynamic", balance_e = 2L)),
    aggte_3k_T12_dynamic_balance_e2 = function() suppressWarnings(aggte(mp12, type = "dynamic", balance_e = 2L)),
    aggte_3k_T12_dynamic_balance_e5 = function() suppressWarnings(aggte(mp12, type = "dynamic", balance_e = 5L)),
    panel_1k_fast_dr_varying = function() fit(vary, weighted = TRUE, fix_weights = "varying"),
    panel_5k_fast_dr_varying = function() fit(vary5, weighted = TRUE, fix_weights = "varying"),
    panel_5k_slow_dr_varying = function() fit(vary5, weighted = TRUE, fix_weights = "varying", faster = FALSE),
    panel_1k_fast_reg_weighted = function() fit(d1, weighted = TRUE, method = "reg"),
    panel_5k_fast_reg_weighted = function() fit(d5, weighted = TRUE, method = "reg"),
    panel_5k_slow_reg_weighted = function() fit(d5, weighted = TRUE, method = "reg", faster = FALSE),
    unbalanced_1k_fast_reg_weighted = function() fit(u1, weighted = TRUE, method = "reg", structure = "unbalanced"),
    rc_1k_fast_reg_weighted = function() fit(r1, weighted = TRUE, method = "reg", structure = "rc")
  )
}
# Warm each path twice. Rotate case order across processes to balance drift.
for (name in names(cases)) for (warm in 1:2) {
  set.seed(20261028 + warm)
  invisible(cases[[name]]())
}
results <- list(); k <- 0L
for (iteration in seq_len(repetitions)) {
  set.seed(20261028 + round * 100L + iteration)
  order <- sample(names(cases))
  for (name in order) {
    gc(verbose = FALSE)
    set.seed(20261028 + iteration)
    timing <- system.time(invisible(cases[[name]]()))
    k <- k + 1L
    results[[k]] <- data.frame(source = label, round = round,
      iteration = iteration, case = name, elapsed = unname(timing[["elapsed"]]),
      user = unname(timing[["user.self"]]), system = unname(timing[["sys.self"]]))
  }
}
out <- list(timings = do.call(rbind, results), session = sessionInfo(),
            source_path = package_path, process_id = Sys.getpid(),
            cases = names(cases), data_seed = 20261028L)
saveRDS(out, output_path)
cat(label, "round", round, "completed", nrow(out$timings), "timings\n")
