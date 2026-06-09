# =============================================================================
# Comprehensive tests for aggte() aggregation types
# =============================================================================

# Shared setup: known DGP with treatment effect = 1
set.seed(20260401)
sp <- did::reset.sim()
sp$te <- 1  # constant treatment effect
data_agg <- did::build_sim_dataset(sp)

mp_agg <- suppressWarnings(suppressMessages(
  att_gt(yname = "Y", xformla = ~X, data = data_agg, tname = "period",
         idname = "id", gname = "G", est_method = "dr",
         bstrap = FALSE)
))

# =============================================================================
# type = "simple"
# =============================================================================

test_that("aggte simple returns valid overall ATT", {
  agg <- suppressWarnings(aggte(mp_agg, type = "simple"))
  expect_s3_class(agg, "AGGTEobj")
  expect_false(is.na(agg$overall.att))
  expect_equal(agg$overall.att, 1, tolerance = 0.5)
})

test_that("aggte simple returns valid SE", {
  agg <- suppressWarnings(aggte(mp_agg, type = "simple"))
  expect_true(agg$overall.se > 0)
  expect_false(is.na(agg$overall.se))
})

test_that("aggte simple has no egt component", {
  agg <- suppressWarnings(aggte(mp_agg, type = "simple"))
  expect_null(agg$egt)
})

test_that("aggte simple influence function has correct length", {
  agg <- suppressWarnings(aggte(mp_agg, type = "simple"))
  expect_equal(length(agg$inf.function$simple.att), nobs(mp_agg))
})

# =============================================================================
# type = "dynamic"
# =============================================================================

test_that("aggte dynamic returns event-time specific ATTs", {
  agg <- suppressWarnings(aggte(mp_agg, type = "dynamic"))
  expect_false(is.null(agg$egt))
  expect_true(length(agg$egt) > 0)
  # Should have both negative and non-negative event times
  expect_true(any(agg$egt < 0))
  expect_true(any(agg$egt >= 0))
})

test_that("aggte dynamic event times are sorted", {
  agg <- suppressWarnings(aggte(mp_agg, type = "dynamic"))
  expect_equal(agg$egt, sort(agg$egt))
})

test_that("aggte dynamic overall.att averages post-treatment event times", {
  agg <- suppressWarnings(aggte(mp_agg, type = "dynamic"))
  expect_false(is.na(agg$overall.att))
  expect_equal(agg$overall.att, 1, tolerance = 0.5)
})

test_that("aggte dynamic min_e filters event times", {
  agg_full <- suppressWarnings(aggte(mp_agg, type = "dynamic"))
  agg_filt <- suppressWarnings(aggte(mp_agg, type = "dynamic", min_e = -1))
  expect_true(min(agg_filt$egt) >= -1)
  expect_true(length(agg_filt$egt) <= length(agg_full$egt))
})

test_that("aggte dynamic max_e filters event times", {
  agg_full <- suppressWarnings(aggte(mp_agg, type = "dynamic"))
  agg_filt <- suppressWarnings(aggte(mp_agg, type = "dynamic", max_e = 1))
  expect_true(max(agg_filt$egt) <= 1)
  expect_true(length(agg_filt$egt) <= length(agg_full$egt))
})

test_that("aggte dynamic min_e and max_e together", {
  agg <- suppressWarnings(aggte(mp_agg, type = "dynamic", min_e = -1, max_e = 1))
  expect_true(min(agg$egt) >= -1)
  expect_true(max(agg$egt) <= 1)
})

test_that("aggte dynamic balance_e filters groups", {
  agg_unbal <- suppressWarnings(aggte(mp_agg, type = "dynamic"))
  agg_bal <- suppressWarnings(aggte(mp_agg, type = "dynamic", balance_e = 1))
  # Balanced version may have fewer event times or same
  expect_true(length(agg_bal$egt) <= length(agg_unbal$egt))
})

test_that("aggte dynamic SEs are positive where ATT is not NA", {
  agg <- suppressWarnings(aggte(mp_agg, type = "dynamic"))
  non_na <- !is.na(agg$att.egt)
  expect_true(all(agg$se.egt[non_na] > 0))
})

# =============================================================================
# type = "group"
# =============================================================================

test_that("aggte group returns per-group ATTs", {
  agg <- suppressWarnings(aggte(mp_agg, type = "group"))
  expect_false(is.null(agg$egt))
  # egt should contain the group values
  expect_true(all(agg$egt %in% unique(mp_agg$group)))
})

test_that("aggte group overall.att is reasonable", {
  agg <- suppressWarnings(aggte(mp_agg, type = "group"))
  expect_false(is.na(agg$overall.att))
  expect_equal(agg$overall.att, 1, tolerance = 0.5)
})

test_that("aggte group SEs are positive for each group", {
  agg <- suppressWarnings(aggte(mp_agg, type = "group"))
  non_na <- !is.na(agg$att.egt)
  expect_true(all(agg$se.egt[non_na] > 0))
})

# =============================================================================
# type = "calendar"
# =============================================================================

test_that("aggte calendar returns per-period ATTs", {
  agg <- suppressWarnings(aggte(mp_agg, type = "calendar"))
  expect_false(is.null(agg$egt))
  expect_true(length(agg$egt) > 0)
})

test_that("aggte calendar overall.att is reasonable", {
  agg <- suppressWarnings(aggte(mp_agg, type = "calendar"))
  expect_false(is.na(agg$overall.att))
  expect_equal(agg$overall.att, 1, tolerance = 0.5)
})

test_that("aggte calendar SEs are positive for each period", {
  agg <- suppressWarnings(aggte(mp_agg, type = "calendar"))
  non_na <- !is.na(agg$att.egt)
  expect_true(all(agg$se.egt[non_na] > 0))
})

test_that("aggte calendar only includes post-treatment periods", {
  agg <- suppressWarnings(aggte(mp_agg, type = "calendar"))
  min_group <- min(mp_agg$group)
  # All calendar periods should be at or after the earliest treatment
  expect_true(all(agg$egt >= min_group))
})

# =============================================================================
# na.rm behavior
# =============================================================================

test_that("aggte with na.rm=TRUE drops NA ATTs and proceeds", {
  # Create MP with some NA ATTs manually
  mp_tmp <- mp_agg
  mp_tmp$att[1] <- NA
  # Should work with na.rm=TRUE
  agg <- suppressWarnings(aggte(mp_tmp, type = "dynamic", na.rm = TRUE))
  expect_s3_class(agg, "AGGTEobj")
  expect_false(is.na(agg$overall.att))
})

# =============================================================================
# Cross-type consistency
# =============================================================================

test_that("all aggte types return AGGTEobj class", {
  for (tp in c("simple", "dynamic", "group", "calendar")) {
    agg <- suppressWarnings(aggte(mp_agg, type = tp))
    expect_s3_class(agg, "AGGTEobj")
  }
})

test_that("all aggte types have non-NULL DIDparams", {
  for (tp in c("simple", "dynamic", "group", "calendar")) {
    agg <- suppressWarnings(aggte(mp_agg, type = tp))
    expect_false(is.null(agg$DIDparams), label = tp)
  }
})

test_that("aggte preserves overridden bstrap/alp settings", {
  agg <- suppressWarnings(aggte(mp_agg, type = "dynamic", bstrap = FALSE, alp = 0.01))
  expect_equal(agg$DIDparams$bstrap, FALSE)
  expect_equal(agg$DIDparams$alp, 0.01)
})

# =============================================================================
# cband with bootstrap
# =============================================================================

test_that("aggte with cband=TRUE and bstrap=TRUE returns simultaneous critical value", {
  mp_boot <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~X, data = data_agg, tname = "period",
           idname = "id", gname = "G", est_method = "reg",
           bstrap = TRUE, biters = 100, cband = TRUE)
  ))
  agg <- aggte(mp_boot, type = "dynamic", bstrap = TRUE, biters = 100, cband = TRUE)
  # Simultaneous critical value should be at least as large as pointwise z
  expect_true(agg$crit.val.egt >= qnorm(0.975))
})

# =============================================================================
# Regression: type = "group" with na.rm = TRUE and a finite max_e
# =============================================================================

test_that("aggte group with na.rm and finite max_e excludes (not errors on) groups whose only non-NA cell is past max_e", {
  # A group whose only non-NA ATT(g,t) lies PAST max_e used to pass the group filter
  # (which ignored max_e) but then have an all-NA, na.rm-emptied selection in the
  # estimate -> stop("No valid att_gt() estimates ..."). The group should instead be
  # dropped from the group aggregation.
  mp <- mp_agg
  g_first <- sort(unique(mp$group))[1]
  # NA-out that group's earliest post-treatment cells, leaving only a late one (e large)
  post_g <- which(mp$group == g_first & mp$t >= g_first)
  if (length(post_g) >= 2) {
    mp$att[post_g[-length(post_g)]] <- NA   # keep only the last (largest event-time) cell
    last_e <- (mp$t[post_g[length(post_g)]] - g_first)
    me <- max(0, last_e - 1)                # max_e below the only surviving cell's event-time
    expect_error(
      suppressWarnings(suppressMessages(aggte(mp, type = "group", na.rm = TRUE, max_e = me))),
      NA)                                   # must NOT error
    agg <- suppressWarnings(suppressMessages(aggte(mp, type = "group", na.rm = TRUE, max_e = me)))
    expect_s3_class(agg, "AGGTEobj")
    expect_false(g_first %in% agg$egt)      # the offending group is excluded
  }
})

test_that("aggte group default (max_e = Inf) is unchanged by the gnotna max_e filter", {
  a1 <- suppressWarnings(suppressMessages(aggte(mp_agg, type = "group")))
  expect_s3_class(a1, "AGGTEobj")
  expect_false(is.na(a1$overall.att))
})
