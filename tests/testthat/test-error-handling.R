# =============================================================================
# Tests for stop(), warning(), and message() conditions
# =============================================================================

# Shared setup
set.seed(20260401)
sp <- reset.sim()
data_eh <- build_sim_dataset(sp)

# =============================================================================
# att_gt() validation errors
# =============================================================================

test_that("att_gt errors on invalid est_method string", {
  expect_error(
    att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
           gname = "G", est_method = "bad", bstrap = FALSE),
    "must be one of"
  )
})

test_that("att_gt errors on non-character non-function est_method", {
  expect_error(
    att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
           gname = "G", est_method = 42, bstrap = FALSE),
    "must be a character string"
  )
})

test_that("att_gt errors on invalid fix_weights value", {
  expect_error(
    att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
           gname = "G", fix_weights = "bad", bstrap = FALSE),
    "must be NULL or one of"
  )
})

test_that("att_gt rejects non-exact control_group and base_period values in both modes", {
  for (fm in c(FALSE, TRUE)) {
    expect_error(
      att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
             gname = "G", control_group = "NotYetTreated", faster_mode = fm,
             bstrap = FALSE),
      "control_group must be either"
    )
    expect_error(
      att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
             gname = "G", base_period = "Universal", faster_mode = fm,
             bstrap = FALSE),
      "base_period must be either"
    )
  }
})

test_that("att_gt rejects negative or non-numeric anticipation in both modes", {
  for (fm in c(FALSE, TRUE)) {
    expect_error(
      att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
             gname = "G", anticipation = -1, faster_mode = fm, bstrap = FALSE),
      "anticipation must be non-negative"
    )
    expect_error(
      att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
             gname = "G", anticipation = "1", faster_mode = fm, bstrap = FALSE),
      "anticipation must be numeric"
    )
  }
})

test_that("att_gt rejects argument-referenced internal variable names in both modes", {
  reserved_names <- c(".w", ".rowid", ".G", ".C", "post",
                      "asif_never_treated", "treated_first_period")
  reserved_data <- data_eh
  for (nm in reserved_names) {
    reserved_data[[nm]] <- seq_len(nrow(reserved_data))
  }
  msg <- "reserved for internal use"

  for (fm in c(FALSE, TRUE)) {
    for (nm in reserved_names) {
      expect_error(
        att_gt(yname = nm, data = reserved_data, tname = "period",
               idname = "id", gname = "G", faster_mode = fm,
               bstrap = FALSE),
        msg
      )
      expect_error(
        att_gt(yname = "Y", xformla = as.formula(paste("~", nm)),
               data = reserved_data, tname = "period", idname = "id",
               gname = "G", faster_mode = fm, bstrap = FALSE),
        msg
      )
      expect_error(
        att_gt(yname = "Y", data = reserved_data, tname = "period",
               idname = "id", gname = "G", weightsname = nm,
               faster_mode = fm, bstrap = FALSE),
        msg
      )
    }
  }
})

test_that("att_gt errors on panel=TRUE without idname in both modes", {
  for (fm in c(FALSE, TRUE)) {
    expect_error(
      att_gt(yname = "Y", data = data_eh, tname = "period",
             gname = "G", bstrap = FALSE, faster_mode = fm),
      "Must provide idname when panel = TRUE"
    )
  }
})

test_that("att_gt still runs with panel=FALSE and no idname in both modes", {
  res <- list()
  for (fm in c(FALSE, TRUE)) {
    res[[as.character(fm)]] <- suppressWarnings(
      att_gt(yname = "Y", data = data_eh, tname = "period", gname = "G",
             panel = FALSE, bstrap = FALSE, faster_mode = fm)
    )
    expect_true(all(is.finite(res[[as.character(fm)]]$att)))
  }
  expect_equal(res[["FALSE"]]$att, res[["TRUE"]]$att, tolerance = 1e-10)
})

test_that("att_gt errors on invalid alp", {
  for (bad_alp in list(1.5, 0, 1, -0.05, c(0.05, 0.1), "0.05", NA_real_)) {
    expect_error(
      att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
             gname = "G", bstrap = FALSE, alp = bad_alp),
      "alp must be a single number strictly between 0 and 1"
    )
  }
})

test_that("att_gt errors on invalid biters when bootstrapping", {
  for (bad_biters in list(-5, 0, 2.5, c(100, 200), "100", NA_real_)) {
    expect_error(
      att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
             gname = "G", bstrap = TRUE, biters = bad_biters),
      "biters must be a single positive whole number"
    )
  }
  # biters is not used (and so not validated) when bstrap = FALSE
  res <- suppressWarnings(
    att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE, biters = -5)
  )
  expect_s3_class(res, "MP")
})

test_that("att_gt errors on fix_weights with panel=FALSE", {
  expect_error(
    att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
           gname = "G", fix_weights = "base_period", panel = FALSE,
           bstrap = FALSE),
    "not supported for repeated cross sections"
  )
})

test_that("att_gt warns on extra args with built-in est_method", {
  expect_warning(
    att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
           gname = "G", est_method = "reg", bstrap = FALSE, extra_arg = 1),
    "Extra arguments"
  )
})

test_that("att_gt messages about anticipation", {
  expect_message(
    suppressWarnings(
      att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
             gname = "G", anticipation = 1, bstrap = FALSE)
    ),
    "anticipation ="
  )
})

test_that("att_gt computes clustered SEs without the bootstrap (no 'requires bootstrap' warning)", {
  # Clustering at the individual level (idname) without the bootstrap returns the unit-level analytical
  # standard errors -- no "requires bootstrap" warning. A coarser cluster variable yields the analytical
  # cluster-robust SE (covered in test-cluster-analytic.R).
  w <- capture_warnings(
    res <- att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
                  gname = "G", clustervars = "id", bstrap = FALSE)
  )
  expect_false(any(grepl("Clustered standard errors require", w)))
  expect_false(any(grepl("could not be computed analytically", w)))
  expect_true(any(is.finite(res$se)))
})

# =============================================================================
# pre_process_did / pre_process_did2 validation errors
# =============================================================================

test_that("att_gt errors on missing column name (slower mode)", {
  expect_error(
    att_gt(yname = "nonexistent", data = data_eh, tname = "period",
           idname = "id", gname = "G", bstrap = FALSE, faster_mode = FALSE),
    "not found"
  )
})

test_that("att_gt errors on missing column name (faster_mode)", {
  expect_error(
    att_gt(yname = "nonexistent", data = data_eh, tname = "period",
           idname = "id", gname = "G", bstrap = FALSE, faster_mode = TRUE),
    "nonexistent"
  )
})

test_that("att_gt errors on non-numeric tname", {
  bad_data <- data_eh
  bad_data$period <- as.character(bad_data$period)
  expect_error(
    att_gt(yname = "Y", data = bad_data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE),
    "must be numeric"
  )
})

test_that("att_gt errors on non-numeric gname", {
  bad_data <- data_eh
  bad_data$G <- as.character(bad_data$G)
  expect_error(
    att_gt(yname = "Y", data = bad_data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE),
    "must be numeric"
  )
})

test_that("att_gt errors on non-numeric outcome variable in both modes", {
  bad_chr <- data_eh
  bad_chr$Y <- as.character(bad_chr$Y)
  bad_fac <- data_eh
  bad_fac$Y <- factor(ifelse(bad_fac$Y > 0, "hi", "lo"))
  for (fm in c(FALSE, TRUE)) {
    expect_error(
      att_gt(yname = "Y", data = bad_chr, tname = "period", idname = "id",
             gname = "G", bstrap = FALSE, faster_mode = fm),
      "The outcome variable 'Y' must be numeric"
    )
    expect_error(
      att_gt(yname = "Y", data = bad_fac, tname = "period", idname = "id",
             gname = "G", bstrap = FALSE, faster_mode = fm),
      "The outcome variable 'Y' must be numeric"
    )
  }
})

test_that("att_gt still accepts logical and integer outcomes in both modes", {
  log_data <- data_eh
  log_data$Y <- log_data$Y > median(log_data$Y)
  int_data <- data_eh
  int_data$Y <- as.integer(round(int_data$Y))
  for (fm in c(FALSE, TRUE)) {
    res_log <- suppressWarnings(
      att_gt(yname = "Y", data = log_data, tname = "period", idname = "id",
             gname = "G", bstrap = FALSE, faster_mode = fm)
    )
    expect_true(all(is.finite(res_log$att)))
    res_int <- suppressWarnings(
      att_gt(yname = "Y", data = int_data, tname = "period", idname = "id",
             gname = "G", bstrap = FALSE, faster_mode = fm)
    )
    expect_true(all(is.finite(res_int$att)))
  }
})

test_that("att_gt errors on treatment reversals (faster_mode)", {
  bad_data <- data_eh
  # Make one unit switch groups across time
  target_id <- bad_data$id[1]
  periods <- sort(unique(bad_data$period))
  bad_data$G[bad_data$id == target_id & bad_data$period == periods[1]] <- 3
  bad_data$G[bad_data$id == target_id & bad_data$period == periods[2]] <- 4
  expect_error(
    att_gt(yname = "Y", data = bad_data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE, faster_mode = TRUE),
    "same across all periods|time-invariant|must be the same"
  )
})

test_that("att_gt warns on missing data dropped", {
  bad_data <- data_eh
  bad_data$Y[1:5] <- NA
  expect_warning(
    att_gt(yname = "Y", data = bad_data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE),
    "dropped|missing"
  )
})

test_that("att_gt warns on small groups", {
  # Create data with a very small treated group
  small_data <- data_eh
  treated_ids <- unique(small_data$id[small_data$G > 0])
  # Keep only 2 treated units and some controls
  keep_ids <- c(treated_ids[1:2], unique(small_data$id[small_data$G == 0])[1:50])
  small_data <- small_data[small_data$id %in% keep_ids, ]
  expect_warning(
    att_gt(yname = "Y", data = small_data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE),
    "very few observations"
  )
})

test_that("att_gt handles data with .w column (slower mode)", {
  # The .w column gets dropped during column selection in pre_process_did
  # before the collision check, so this should run without error
  bad_data <- data.frame(data_eh)
  bad_data$.w <- 1
  result <- att_gt(yname = "Y", data = bad_data, tname = "period", idname = "id",
                   gname = "G", bstrap = FALSE, faster_mode = FALSE)
  expect_s3_class(result, "MP")
})

test_that("att_gt messages on time-varying weights (panel)", {
  tv_data <- data_eh
  tv_data$tw <- tv_data$period + runif(nrow(tv_data))
  expect_message(
    att_gt(yname = "Y", data = tv_data, tname = "period", idname = "id",
           gname = "G", weightsname = "tw", bstrap = FALSE),
    "Time-varying weights"
  )
})

test_that("clustervars contract enforced in all faster_mode x bstrap combinations", {
  # ?att_gt promises an error for more than two cluster variables (one of which
  # must be idname). Both pre-processors must reject the input with identical
  # wording (also matching mboot()) regardless of bstrap -- the analytical
  # (bstrap = FALSE) slow path used to silently cluster on the first extra
  # variable only.
  d <- data_eh
  d$cl1 <- d$id %% 10  # time-invariant within unit
  d$cl2 <- d$id %% 7
  msg <- "At most one cluster variable \\(beyond 'idname'\\) is supported"
  for (fm in c(FALSE, TRUE)) {
    for (bs in c(FALSE, TRUE)) {
      # three cluster variables (idname + two extras)
      expect_error(
        att_gt(yname = "Y", data = d, tname = "period", idname = "id",
               gname = "G", clustervars = c("id", "cl1", "cl2"),
               bstrap = bs, faster_mode = fm),
        msg
      )
      # two cluster variables, neither of which is idname
      expect_error(
        att_gt(yname = "Y", data = d, tname = "period", idname = "id",
               gname = "G", clustervars = c("cl1", "cl2"),
               bstrap = bs, faster_mode = fm),
        msg
      )
    }
  }
})

test_that("one extra cluster variable still works in all faster_mode x bstrap combinations", {
  # clustered DGP: units within a cluster share a common per-period shock, so
  # the clustered SE genuinely differs from the unclustered one
  set.seed(20260609)
  n_cl <- 40L
  sz <- rep(c(2L, 4L), length.out = n_cl)
  N <- sum(sz)
  cl <- rep(seq_len(n_cl), times = sz)
  gC <- sample(c(2L, 3L, 0L), n_cl, replace = TRUE)[cl]
  eta <- matrix(rnorm(n_cl * 4L, 0, 1.5), n_cl, 4L)
  d <- data.frame(id = rep(seq_len(N), each = 4L), period = rep(1:4, N),
                  cl = rep(cl, each = 4L), G = rep(gC, each = 4L))
  d$Y <- eta[cbind(d$cl, d$period)] + 0.3 * d$period +
    (d$G != 0 & d$period >= d$G) + rnorm(nrow(d))
  for (fm in c(FALSE, TRUE)) {
    se_cl <- NULL
    for (bs in c(FALSE, TRUE)) {
      res <- suppressWarnings(
        att_gt(yname = "Y", data = d, tname = "period", idname = "id",
               gname = "G", clustervars = c("id", "cl"), bstrap = bs,
               biters = 100, cband = FALSE, faster_mode = fm)
      )
      expect_s3_class(res, "MP")
      expect_true(any(is.finite(res$se)))
      if (!bs) se_cl <- res$se
    }
    # the analytical (bstrap = FALSE) path actually clusters: the SE differs
    # from the unclustered analytical SE
    res_iid <- att_gt(yname = "Y", data = d, tname = "period", idname = "id",
                      gname = "G", bstrap = FALSE, faster_mode = fm)
    ok <- is.finite(se_cl) & is.finite(res_iid$se)
    expect_true(any(ok))
    expect_gt(max(abs(se_cl[ok] - res_iid$se[ok])), 1e-8)
  }
})

# =============================================================================
# aggte() validation errors
# =============================================================================

test_that("aggte errors on invalid type", {
  mp_tmp <- att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
                   gname = "G", bstrap = FALSE)
  expect_error(
    aggte(mp_tmp, type = "invalid"),
    "must be one of"
  )
})

test_that("aggte errors when ATTs contain NA and na.rm=FALSE", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  mp_na <- suppressWarnings(suppressMessages(
    att_gt(yname = "Y", data = data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE)
  ))
  # Inject NA directly — no need to rely on estimation failure
  mp_na$att[1] <- NA
  expect_error(
    aggte(mp_na, type = "dynamic", na.rm = FALSE),
    "Missing values"
  )
})

# =============================================================================
# Wald pre-test warnings
# =============================================================================

test_that("att_gt handles singular covariance for Wald test gracefully", {
  # When covariance matrix is singular, W and Wpval should be NULL
  # and the result should still be a valid MP object
  small_sp <- reset.sim()
  small_sp$n <- 50
  small_data <- build_sim_dataset(small_sp)
  result <- suppressWarnings(
    att_gt(yname = "Y", data = small_data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE)
  )
  expect_s3_class(result, "MP")
  # W should be either a valid number or NULL (when singular)
  expect_true(is.null(result$W) || is.numeric(result$W))
})

# =============================================================================
# compute.att_gt warnings
# =============================================================================

test_that("att_gt warns on overlap violations", {
  # Create data where propensity score is near 1 (near-perfect separation)
  sep_data <- data_eh
  # Make covariate perfectly predict treatment for some cells
  sep_data$X_sep <- ifelse(sep_data$G > 0, 100, -100)
  expect_warning(
    att_gt(yname = "Y", xformla = ~X_sep, data = sep_data, tname = "period",
           idname = "id", gname = "G", est_method = "dr", bstrap = FALSE),
    "overlap condition"
  )
})

test_that("slow path warns once per failed cell, matching fast-path wording, with accurate Wald diagnosis", {
  # Each failed (g,t) cell must produce exactly ONE warning with identical text in
  # both modes: the tryCatch wrapper must not re-warn the internal overlap/singular
  # sentinels after the diagnostic warning already fired. And when pre-treatment
  # cells existed but were all dropped for NA/zero variance, the Wald warning must
  # say so instead of claiming no pre-treatment cells exist.
  set.seed(20260610)
  d <- expand.grid(id = 1:200, period = 1:3)
  gvals <- c(0, 2, 3)
  d$G <- gvals[(d$id - 1) %% 3 + 1]
  d$Xsep <- 1 * (d$G > 0)  # separates treated from never-treated: overlap fails in all 4 cells
  d$Y <- 0.1 * d$period + (d$G > 0 & d$period >= d$G) + rnorm(nrow(d), 0, 0.5)

  w_slow <- capture_warnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~Xsep, data = d, tname = "period", idname = "id",
           gname = "G", est_method = "dr", faster_mode = FALSE, bstrap = FALSE)))
  w_fast <- capture_warnings(suppressMessages(
    att_gt(yname = "Y", xformla = ~Xsep, data = d, tname = "period", idname = "id",
           gname = "G", est_method = "dr", faster_mode = TRUE, bstrap = FALSE)))

  # identical warning sets across modes (one per failed cell + the Wald warning)
  expect_identical(w_slow, w_fast)
  # exactly one warning per failed cell, with the fast-path wording
  expect_identical(sum(grepl("overlap condition violated for group", w_slow)), 4L)
  # no leaked internal sentinel via the tryCatch wrapper
  expect_false(any(grepl("Error computing internal 2x2 DiD", w_slow)))
  # accurate Wald diagnosis: pre-treatment cells existed but had NA/zero variance
  expect_true(any(grepl("missing or zero variance", w_slow)))
  expect_false(any(grepl("treated early in the panel", w_slow)))
})

test_that("att_gt warns when no pre-treatment periods for Wald test", {
  # Create data where a group is first treated in period 2 (only 1 pre-treatment period)
  # With base_period="varying", there may be no pre-treatment ATTs for the Wald test
  early_data <- data.frame(
    id = rep(1:100, each = 3),
    period = rep(1:3, 100),
    G = rep(c(rep(2, 50), rep(0, 50)), each = 3),
    Y = rnorm(300)
  )
  expect_warning(
    att_gt(yname = "Y", data = early_data, tname = "period", idname = "id",
           gname = "G", bstrap = FALSE),
    "pre-treatment|Wald"
  )
})

test_that("empty-cell warning names the actually-empty base period under base_period='universal' in both modes", {
  # Group 3's universal base period is 2; group 3 has observations in every
  # OTHER period. The slow path used to name the cell's current period
  # (1/3/4) -- periods where the group does have units -- instead of the empty
  # base period (2), which the fast path reported correctly.
  set.seed(20260610)
  d <- expand.grid(rep = 1:30, period = 1:4, G = c(0, 3))
  d$Y <- rnorm(nrow(d)) + 0.1 * d$period + (d$G > 0 & d$period >= d$G)
  d <- d[!(d$G == 3 & d$period == 2), ]  # group 3 empty only in its base period
  for (fm in c(FALSE, TRUE)) {
    w <- capture_warnings(suppressMessages(
      att_gt(yname = "Y", data = d, tname = "period", gname = "G",
             panel = FALSE, base_period = "universal", bstrap = FALSE,
             faster_mode = fm)
    ))
    # one warning per non-base cell (t = 1, 3, 4), each naming period 2
    expect_identical(
      sum(grepl("No units in group 3 in time period 2", w, fixed = TRUE)), 3L,
      info = paste("faster_mode =", fm)
    )
    expect_false(any(grepl("No units in group 3 in time period [134]", w)),
                 info = paste("faster_mode =", fm))
  }
})

# =============================================================================
# pre-processing warning parity across modes
# =============================================================================

test_that("'no never-treated group' warning is identical across modes and discloses the period filtering", {
  # Both modes drop all periods >= (latest cohort - anticipation); the fast-mode
  # warning used to omit that a chunk of the data was discarded.
  set.seed(20260609)
  d <- data.frame(
    id = rep(1:60, each = 4),
    period = rep(1:4, 60),
    G = rep(c(rep(2, 30), rep(4, 30)), each = 4)
  )
  d$Y <- rnorm(nrow(d))
  msgs <- list()
  for (fm in c(FALSE, TRUE)) {
    w <- capture_warnings(suppressMessages(
      att_gt(yname = "Y", data = d, tname = "period", idname = "id",
             gname = "G", control_group = "nevertreated", bstrap = FALSE,
             faster_mode = fm)
    ))
    m <- grep("No never-treated group", w, value = TRUE)
    expect_length(m, 1L)
    expect_match(m, "filtered out", fixed = TRUE)
    msgs[[as.character(fm)]] <- m
  }
  expect_identical(msgs[["FALSE"]], msgs[["TRUE"]])
})

test_that("balanced-panel coercion warning counts dropped units identically in both modes", {
  # 5 units each miss one period's row: 5 UNITS (20 rows) are dropped while
  # coercing to a balanced panel. The slow path used to report the unit count
  # as 'observations' (a 4x undercount of the rows actually removed).
  drop_ids <- unique(data_eh$id)[1:5]
  d <- data_eh[!(data_eh$id %in% drop_ids & data_eh$period == max(data_eh$period)), ]
  expected <- "5 units are missing in some periods. Converting to balanced panel by dropping them."
  for (fm in c(FALSE, TRUE)) {
    w <- capture_warnings(suppressMessages(
      att_gt(yname = "Y", data = d, tname = "period", idname = "id",
             gname = "G", panel = TRUE, allow_unbalanced_panel = FALSE,
             bstrap = FALSE, faster_mode = fm)
    ))
    expect_true(any(w == expected), info = paste("faster_mode =", fm))
    expect_false(any(grepl("observations while converting", w)),
                 info = paste("faster_mode =", fm))
  }
})
