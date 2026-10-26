# Independent four-cell support-count reference, including local normalization.
script <- grep("^--file=", commandArgs(), value = TRUE)[[1L]]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script)), "..", "..", ".."))
pkgload::load_all(root, quiet = TRUE)
reference <- function(D, post, w, k, include_treated) {
  mean_w <- mean(w)
  if (!is.finite(mean_w) || mean_w <= 0) return(TRUE)
  normalized <- w / mean_w
  supported <- is.finite(normalized) & normalized > 0
  groups <- if (include_treated) 0:1 else 0L
  any(vapply(groups, function(g) any(vapply(0:1, function(t) {
    sum(supported & D == g & post == t) < k
  }, logical(1))), logical(1)))
}
set.seed(20261028)
checks <- logical()
for (i in seq_len(2000L)) {
  n <- sample(0:400, 1L)
  D <- sample(0:1, n, TRUE); post <- sample(0:1, n, TRUE)
  w <- switch(1L + i %% 8L,
    rep(1, n), runif(n), exp(rnorm(n, 0, 150)),
    sample(c(0, 1), n, TRUE), rep(2^-1074, n),
    sample(c(2^1000, 2^-76, 2^-75, 0), n, TRUE),
    rep(Inf, n), rep(NA_real_, n))
  k <- sample(1:6, 1L)
  for (include_treated in c(TRUE, FALSE)) {
    actual <- did:::rc_outcome_has_too_few_obs(D, post, w, k, include_treated)
    expected <- reference(D, post, w, k, include_treated)
    checks <- c(checks, identical(actual, expected))
  }
}
stopifnot(all(checks))
cat(length(checks), "support decisions match the independent reference\n")
