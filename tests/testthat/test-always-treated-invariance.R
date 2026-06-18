# Regression tests for the latest-control-cohort deletion bug.
#
# Bug (did <= 2.5.0): with control_group = "notyettreated" and NO never-treated
# group, the latest cohort is correctly removed from glist so it can serve as a
# not-yet-treated control. But when some unit is treated in the first period
# (accounting for anticipation), the "drop first-period-treated units" step
# filtered the data by cohort membership in glist (plus the never-treated
# sentinel -- 0 in the slow path, Inf in the fast path) -- and because the latest
# cohort is no longer in glist, that filter ALSO deleted the entire latest
# control cohort, silently corrupting ATT(g,t) for every other group.
#
# Invariant A (MUST): att_gt output must be INVARIANT to the presence of
# effectively-always-treated units (treated <= first period, accounting for
# anticipation). Three trigger pathways are exercised:
#   P1  a literal always-treated cohort (g <= first period), anticipation = 0
#   P2  anticipation >= 1 promotes the EARLIEST cohort to first-period-treated
#       (fires with NO nominal always-treated unit)
#   P3  anticipation un-coerces a late (g = T + 1) cohort from never-treated,
#       removing the sole never-treated group

# --- helpers ----------------------------------------------------------------

# build a balanced panel from a vector of cohorts; `het` controls level spread
.mk_design <- function(cohorts, n_periods, n = 25, seed = 1, het = 1.2) {
  set.seed(seed)
  rows <- data.table::rbindlist(lapply(cohorts, function(g)
    data.table::data.table(g = g,
                           b = exp(stats::rnorm(n, log(1e5), het)),
                           uid = NA_integer_)))
  rows[, uid := .I]
  d <- merge(data.table::CJ(uid = rows$uid, t = 1:n_periods),
             rows[, .(uid, g, b)], by = "uid")
  d[, y := b + 50 * (t - 1) + 200 * (t >= g & g > 0) + stats::rnorm(.N, 0, 5)]
  as.data.frame(d[])
}

.att_keyed <- function(data, drop_uid = integer(0), ...) {
  if (length(drop_uid)) data <- data[!(data$uid %in% drop_uid), ]
  cs <- suppressWarnings(suppressMessages(
    att_gt(yname = "y", tname = "t", idname = "uid", gname = "g",
           data = data, control_group = "notyettreated",
           bstrap = FALSE, ...)))
  stats::setNames(cs$att, paste(cs$group, cs$t, sep = "_"))
}

# assert: ATT(g,t) on cells common to full vs (full minus EAT units) is identical,
# including the NA pattern
.expect_invariant <- function(data, eat_uid, ..., tol = 1e-8) {
  full <- .att_keyed(data, ...)
  drop <- .att_keyed(data, drop_uid = eat_uid, ...)
  common <- intersect(names(full), names(drop))
  expect_true(length(common) > 0)
  expect_equal(is.na(full[common]), is.na(drop[common]))   # NA pattern matches
  ok <- !is.na(full[common])
  expect_equal(unname(full[common][ok]), unname(drop[common][ok]), tolerance = tol)
}

# --- P1: literal always-treated cohort, anticipation = 0 --------------------

test_that("Invariant A holds under P1 (always-treated cohort, no never-treated)", {
  d <- .mk_design(c(1, 2, 3, 4), n_periods = 4, seed = 2024)
  eat <- unique(d$uid[d$g == 1])
  for (em in c("dr", "reg", "ipw")) {
    .expect_invariant(d, eat, est_method = em, faster_mode = TRUE)
    .expect_invariant(d, eat, est_method = em, faster_mode = FALSE)
  }
})

# --- P2: anticipation promotes earliest cohort (no always-treated) ----------

test_that("Invariant A holds under P2 (anticipation>=1 promotes earliest cohort)", {
  d <- .mk_design(c(2, 3, 5), n_periods = 5, seed = 11)
  eat <- unique(d$uid[d$g == 2])          # earliest cohort, promoted at anticipation = 1
  .expect_invariant(d, eat, anticipation = 1, faster_mode = TRUE)
  .expect_invariant(d, eat, anticipation = 1, faster_mode = FALSE)
})

# --- P3: anticipation removes the sole (coerced) never-treated group ---------

test_that("Invariant A holds under P3 (anticipation un-coerces a late cohort)", {
  d <- .mk_design(c(2, 3, 6), n_periods = 5, seed = 7)   # cohort 6 = T+1, coerced-to-never at anti=0
  eat <- unique(d$uid[d$g == 2])
  .expect_invariant(d, eat, anticipation = 1, faster_mode = TRUE)
  .expect_invariant(d, eat, anticipation = 1, faster_mode = FALSE)
})

# --- scaling oracle: effectively-always-treated outcomes are never read ------

test_that("always-treated outcomes are never read (scaling oracle)", {
  d <- .mk_design(c(1, 2, 3, 4), n_periods = 4, seed = 2024)
  d_scaled <- d
  d_scaled$y[d_scaled$g == 1] <- d_scaled$y[d_scaled$g == 1] * 1e6
  a1 <- .att_keyed(d)
  a2 <- .att_keyed(d_scaled)
  expect_equal(is.na(a1), is.na(a2))
  ok <- !is.na(a1)
  expect_equal(unname(a1[ok]), unname(a2[ok]), tolerance = 1e-8)
})

# --- fast/slow parity in the vulnerable regime ------------------------------

test_that("fast and slow paths agree in the no-never + always-treated regime", {
  d <- .mk_design(c(1, 2, 3, 4, 5), n_periods = 5, seed = 99)
  fast <- .att_keyed(d, faster_mode = TRUE)
  slow <- .att_keyed(d, faster_mode = FALSE)
  expect_equal(is.na(fast), is.na(slow))
  ok <- !is.na(fast)
  expect_equal(unname(fast[ok]), unname(slow[ok]), tolerance = 1e-10)
})

# --- structural: latest cohort kept as control, not an estimated group ------

test_that("latest cohort is retained as a control, not deleted, when always-treated present", {
  d <- .mk_design(c(1, 2, 3, 4, 5), n_periods = 5, seed = 5)
  dp <- suppressWarnings(suppressMessages(
    pre_process_did(yname = "y", tname = "t", idname = "uid", gname = "g",
                    data = d, panel = TRUE, allow_unbalanced_panel = FALSE,
                    control_group = "notyettreated", print_details = FALSE)))
  kept <- sort(unique(dp$data$g))
  cs <- suppressWarnings(suppressMessages(
    att_gt(yname = "y", tname = "t", idname = "uid", gname = "g",
           data = d, control_group = "notyettreated", bstrap = FALSE)))
  expect_true(5 %in% kept)                 # latest cohort kept in the data (as control)
  expect_false(1 %in% kept)                # always-treated cohort dropped
  expect_false(5 %in% unique(cs$group))    # latest cohort gets no ATT of its own
})

# --- nevertreated branch was always immune: confirm it stays invariant -------

test_that("control_group='nevertreated' is invariant to always-treated presence", {
  d <- .mk_design(c(1, 2, 3, 4), n_periods = 4, seed = 2024)
  eat <- unique(d$uid[d$g == 1])
  full <- suppressWarnings(suppressMessages(att_gt("y", "t", "uid", "g", data = d,
            control_group = "nevertreated", bstrap = FALSE)))
  drop <- suppressWarnings(suppressMessages(att_gt("y", "t", "uid", "g",
            data = d[!(d$uid %in% eat), ], control_group = "nevertreated", bstrap = FALSE)))
  fk <- stats::setNames(full$att, paste(full$group, full$t, sep = "_"))
  dk <- stats::setNames(drop$att, paste(drop$group, drop$t, sep = "_"))
  common <- intersect(names(fk), names(dk))
  expect_true(length(common) > 0)
  expect_equal(is.na(fk[common]), is.na(dk[common]))
  ok <- !is.na(fk[common])
  expect_equal(unname(fk[common][ok]), unname(dk[common][ok]), tolerance = 1e-8)
})
