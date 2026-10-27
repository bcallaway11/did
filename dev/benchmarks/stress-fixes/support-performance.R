# Compare the optimized RC preflight with its original correct implementation.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else file.path(tempdir(), "did-rc-support-performance.csv")
script <- grep("^--file=", commandArgs(), value = TRUE)[[1L]]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script)), "..", "..", ".."))
setwd(root)
pkgload::load_all(root, quiet = TRUE)
old_source <- tempfile(fileext = ".R")
status <- system2("git", c("show", "1a16c46:R/utility_functions.R"), stdout = old_source)
stopifnot(status == 0L)
old <- new.env(parent = asNamespace("did")); sys.source(old_source, old)
functions <- list(before = old$rc_outcome_has_too_few_obs,
                  after = did:::rc_outcome_has_too_few_obs)
set.seed(20261028)
n <- 5000L
D <- rep(c(0, 0, 1, 1), length.out = n)
post <- rep(c(0, 1, 0, 1), length.out = n)
weights <- list(positive = runif(n, .5, 2),
                zero = sample(c(0, 1), n, TRUE),
                underflow = rep(c(1e300, 1e-300, 1, 1), length.out = n))
for (w in weights) for (f in functions) for (i in seq_len(100L)) invisible(f(D, post, w, 3L))
r <- list(); j <- 0L
for (iteration in seq_len(20L)) for (case in names(weights)) {
  labels <- if (iteration %% 2L) names(functions) else rev(names(functions))
  for (label in labels) {
    gc(FALSE)
    f <- functions[[label]]; w <- weights[[case]]
    tm <- system.time(for (i in seq_len(500L)) invisible(f(D, post, w, 3L)))
    j <- j + 1L
    r[[j]] <- data.frame(iteration, case, source = label,
                         elapsed = unname(tm[["elapsed"]]))
  }
}
r <- do.call(rbind, r)
summary <- do.call(rbind, lapply(split(r, list(r$case, r$source)), function(z) {
  data.frame(case = z$case[[1L]], source = z$source[[1L]],
    median_us_per_call = median(z$elapsed) * 1e6 / 500,
    q25_us_per_call = unname(quantile(z$elapsed, .25)) * 1e6 / 500,
    q75_us_per_call = unname(quantile(z$elapsed, .75)) * 1e6 / 500)
}))
write.csv(summary, output, row.names = FALSE)
write.csv(r, sub("\\.csv$", "-raw.csv", output), row.names = FALSE)
print(summary, row.names = FALSE)
