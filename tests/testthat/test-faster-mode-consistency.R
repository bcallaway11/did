# =============================================================================
# Systematic faster_mode=TRUE vs FALSE consistency tests
# =============================================================================

# Shared setup
set.seed(20260401)
sp <- did::reset.sim()
data_fm <- did::build_sim_dataset(sp)

# Unbalanced version
data_ub <- data_fm[-c(1, 5, 10), ]

# Helper to compare two att_gt results.
# NOTE: the fast and slow paths order units differently (the exported influence
# function is a ROW PERMUTATION between modes), so a naive expect_equal of the two
# inffunc matrices would falsely fail. We therefore compare permutation-INVARIANT
# summaries -- column sums, the Gram matrix t(IF) %*% IF (which is exactly what the
# analytic variance V = crossprod(IF)/n uses), and the resulting standard errors --
# which DO have to match and which catch any row-permutation / indexing bug.
compare_modes <- function(res_slow, res_fast, label) {
  expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10, label = paste(label, "ATT"))
  expect_equal(res_slow$group, res_fast$group, label = paste(label, "group"))
  expect_equal(res_slow$t, res_fast$t, label = paste(label, "t"))
  if (!is.null(res_slow$inffunc) && !is.null(res_fast$inffunc)) {
    ifs <- as.matrix(res_slow$inffunc); iff <- as.matrix(res_fast$inffunc)
    expect_equal(dim(ifs), dim(iff), label = paste(label, "IF dim"))
    expect_equal(colSums(ifs), colSums(iff), tolerance = 1e-8, label = paste(label, "IF colSums"))
    expect_equal(crossprod(ifs), crossprod(iff), tolerance = 1e-8, label = paste(label, "IF Gram"))
  }
  expect_equal(res_slow$se, res_fast$se, tolerance = 1e-9, label = paste(label, "SE"))
}

# =============================================================================
# Core grid: est_method x panel_type x control_group x base_period
# =============================================================================

# --- Balanced panel ---

for (em in c("dr", "ipw", "reg")) {
  for (cg in c("nevertreated", "notyettreated")) {
    for (bp in c("varying", "universal")) {
      label <- paste(em, "balanced", cg, bp)
      test_that(paste("consistency:", label), {
        res_slow <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 faster_mode = FALSE, bstrap = FALSE)
        ))
        res_fast <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 faster_mode = TRUE, bstrap = FALSE)
        ))
        compare_modes(res_slow, res_fast, label)
      })
    }
  }
}

# --- Unbalanced panel ---

for (em in c("dr", "ipw", "reg")) {
  for (cg in c("nevertreated", "notyettreated")) {
    for (bp in c("varying", "universal")) {
      label <- paste(em, "unbalanced", cg, bp)
      test_that(paste("consistency:", label), {
        res_slow <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_ub, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 allow_unbalanced_panel = TRUE,
                 faster_mode = FALSE, bstrap = FALSE)
        ))
        res_fast <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_ub, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 allow_unbalanced_panel = TRUE,
                 faster_mode = TRUE, bstrap = FALSE)
        ))
        compare_modes(res_slow, res_fast, label)
      })
    }
  }
}

# --- Repeated cross sections ---

for (em in c("dr", "ipw", "reg")) {
  for (cg in c("nevertreated", "notyettreated")) {
    for (bp in c("varying", "universal")) {
      label <- paste(em, "RC", cg, bp)
      test_that(paste("consistency:", label), {
        res_slow <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 panel = FALSE,
                 faster_mode = FALSE, bstrap = FALSE)
        ))
        res_fast <- suppressWarnings(suppressMessages(
          att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
                 idname = "id", gname = "G", est_method = em,
                 control_group = cg, base_period = bp,
                 panel = FALSE,
                 faster_mode = TRUE, bstrap = FALSE)
        ))
        compare_modes(res_slow, res_fast, label)
      })
    }
  }
}

# =============================================================================
# Additional consistency tests beyond the core grid
# =============================================================================

test_that("consistency with anticipation > 0", {
  res_slow <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           anticipation = 1, faster_mode = FALSE, bstrap = FALSE)
  ))
  res_fast <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           anticipation = 1, faster_mode = TRUE, bstrap = FALSE)
  ))
  compare_modes(res_slow, res_fast, "anticipation=1")
})

test_that("consistency without covariates", {
  res_slow <- att_gt(yname = "Y", data = data_fm, tname = "period",
                     idname = "id", gname = "G", est_method = "reg",
                     faster_mode = FALSE, bstrap = FALSE)
  res_fast <- att_gt(yname = "Y", data = data_fm, tname = "period",
                     idname = "id", gname = "G", est_method = "reg",
                     faster_mode = TRUE, bstrap = FALSE)
  compare_modes(res_slow, res_fast, "no covariates")
})

# =============================================================================
# aggte consistency across modes
# =============================================================================

test_that("aggte consistency across faster_mode for all types", {
  mp_slow <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           faster_mode = FALSE, bstrap = FALSE)
  ))
  mp_fast <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_fm, tname = "period",
           idname = "id", gname = "G", est_method = "dr",
           faster_mode = TRUE, bstrap = FALSE)
  ))

  for (tp in c("simple", "dynamic", "group", "calendar")) {
    agg_slow <- suppressWarnings(aggte(mp_slow, type = tp))
    agg_fast <- suppressWarnings(aggte(mp_fast, type = tp))
    expect_equal(agg_slow$overall.att, agg_fast$overall.att,
                 tolerance = 1e-8, label = paste("aggte", tp))
  }
})

# =============================================================================
# NA-propagation semantics in the unbalanced-panel influence aggregation
# =============================================================================

# With allow_unbalanced_panel = TRUE, NA/NaN entries in a 2x2 influence function
# (only reachable via a custom est_method returning a finite ATT alongside a
# partially-NA influence function) must propagate to the unit's aggregated
# influence value -- so that cell's SE is NA -- in BOTH modes. The fast path
# used to zero them silently (sum(na.rm = TRUE)), understating the SE and
# diverging from the slow path; this pins the propagating semantics.
test_that("partial-NA influence functions propagate identically across modes (unbalanced panel)", {
  set.seed(20260609)
  sp_na <- did::reset.sim(time.periods = 3, n = 400)
  data_na <- did::build_sim_dataset(sp_na, panel = TRUE)
  # Drop period 3 for 10 units so the panel is genuinely unbalanced (otherwise
  # pre-processing downgrades allow_unbalanced_panel to a balanced panel).
  drop_ids <- head(unique(data_na$id), 10)
  data_na <- data_na[!(data_na$id %in% drop_ids & data_na$period == 3), ]

  # Finite ATT, one NA injected into the influence function of every 2x2 cell
  na_inj <- function(y, post, D, covariates, i.weights, inffunc, ...) {
    res <- DRDID::reg_did_rc(y = y, post = post, D = D, covariates = covariates,
                             i.weights = i.weights, boot = FALSE, inffunc = TRUE)
    inf <- res$att.inf.func
    inf[1] <- NA_real_
    list(ATT = res$ATT, att.inf.func = inf)
  }

  run_na <- function(fm) suppressWarnings(suppressMessages(
    att_gt(yname = "Y", tname = "period", idname = "id", gname = "G",
           xformla = ~X, data = data_na, panel = TRUE,
           allow_unbalanced_panel = TRUE, est_method = na_inj,
           faster_mode = fm, bstrap = FALSE, cband = FALSE)
  ))

  res_fast <- run_na(TRUE)
  res_slow <- run_na(FALSE)

  # ATT point estimates stay finite and identical across modes
  expect_true(all(is.finite(res_fast$att)))
  expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10)
  # The injected NA propagates: every (g,t) standard error is NA in both modes
  expect_true(all(is.na(res_fast$se)))
  expect_true(all(is.na(res_slow$se)))
  # Exactly one aggregated influence entry per (g,t) cell is NA, in both modes
  na_fast <- sum(is.na(as.matrix(res_fast$inffunc)))
  na_slow <- sum(is.na(as.matrix(res_slow$inffunc)))
  expect_identical(na_fast, na_slow)
  expect_identical(na_fast, length(res_fast$att))
})
