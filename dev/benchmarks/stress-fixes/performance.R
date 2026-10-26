# Usage after source freeze:
# Rscript dev/benchmarks/stress-fixes/performance.R [repetitions=7] [rounds=2] [output_directory] [suite=main]
# Creates an exact clean-master baseline and a frozen snapshot of the working
# source, then interleaves independent warmed R processes (ABBA across rounds).
# Package loading, input construction, and explicit GC run outside each timing.
args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args)) as.integer(args[[1L]]) else 7L
rounds <- if (length(args) > 1L) as.integer(args[[2L]]) else 2L
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[[1L]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
root <- normalizePath(file.path(script_dir, "..", "..", ".."))
setwd(root)
suite <- if (length(args) > 3L) args[[4L]] else "main"
stopifnot(suite %in% c("main", "supplemental"))
artifact_dir <- if (length(args) > 2L) args[[3L]] else
  file.path(tempdir(), paste0("did-stress-performance-", suite))
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
baseline <- file.path(artifact_dir, "baseline-source")
final <- file.path(artifact_dir, "final-source")
dir.create(baseline, showWarnings = FALSE); dir.create(final, showWarnings = FALSE)
archive <- file.path(artifact_dir, "baseline.tar")
status <- system2("git", c("archive", "--format=tar", "--output", shQuote(archive),
                          "74b88eb", "DESCRIPTION", "NAMESPACE", "R", "data"))
stopifnot(status == 0L)
untar(archive, exdir = baseline)
for (path in c("DESCRIPTION", "NAMESPACE", "R", "data")) {
  destination <- file.path(final, basename(path))
  if (dir.exists(destination)) unlink(destination, recursive = TRUE)
  stopifnot(file.copy(file.path(root, path), final, recursive = TRUE))
}
metadata <- list(baseline = "74b88eb07f1faa644df1271055a0f77848f6550c",
  source_head = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  source_diff = system2("git", c("diff", "--stat"), stdout = TRUE),
  source_hashes = tools::md5sum(list.files(file.path(final, "R"), pattern = "\\.R$",
    full.names = TRUE, recursive = TRUE)),
  suite = suite,
  timestamp = Sys.time(), session = sessionInfo(),
  design = "sequential warmed processes; baseline/final then final/baseline; randomized case order; GC excluded")
saveRDS(metadata, file.path(artifact_dir, "metadata.rds"))
worker <- file.path(script_dir, "performance-worker.R")
files <- character()
for (round in seq_len(rounds)) {
  labels <- if (round %% 2L) c("baseline", "final") else c("final", "baseline")
  for (label in labels) {
    output <- file.path(artifact_dir, paste0(label, "-round-", round, ".rds"))
    log <- file.path(artifact_dir, paste0(label, "-round-", round, ".log"))
    path <- if (label == "baseline") baseline else final
    status <- system2(file.path(R.home("bin"), "Rscript"),
      c(shQuote(worker), shQuote(path), label, shQuote(output), round, repetitions, suite),
      stdout = log, stderr = log)
    if (status != 0L) stop("Benchmark worker failed: ", log)
    files <- c(files, output)
    cat(label, "round", round, "finished\n"); flush.console()
  }
}
raw <- do.call(rbind, lapply(files, function(path) readRDS(path)$timings))
write.csv(raw, file.path(artifact_dir, "raw-timings.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(raw, list(raw$case, raw$source), drop = TRUE), function(z) {
  data.frame(case = z$case[[1L]], source = z$source[[1L]], n = nrow(z),
    median = median(z$elapsed), q25 = unname(quantile(z$elapsed, .25)),
    q75 = unname(quantile(z$elapsed, .75)),
    median_cpu = median(z$user + z$system))
}))
comparison <- merge(summary[summary$source == "baseline", ],
                    summary[summary$source == "final", ], by = "case",
                    suffixes = c("_baseline", "_final"))
comparison$ratio <- comparison$median_final / comparison$median_baseline
comparison$delta_ms <- 1000 * (comparison$median_final - comparison$median_baseline)
comparison$iqr_overlap <- comparison$q25_final <= comparison$q75_baseline &
  comparison$q25_baseline <= comparison$q75_final
write.csv(comparison, file.path(artifact_dir, "comparison.csv"), row.names = FALSE)
saveRDS(list(metadata = metadata, raw = raw, summary = summary, comparison = comparison),
        file.path(artifact_dir, "results.rds"))
cat("Artifacts:", normalizePath(artifact_dir), "\n")
print(comparison[c("case", "median_baseline", "q25_baseline", "q75_baseline",
                   "median_final", "q25_final", "q75_final", "ratio", "delta_ms", "iqr_overlap")], row.names = FALSE)
