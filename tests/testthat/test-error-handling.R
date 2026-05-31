# =============================================================================
# Tests for stop(), warning(), and message() conditions
# =============================================================================

# Shared setup
set.seed(20260401)
sp <- did::reset.sim()
data_eh <- did::build_sim_dataset(sp)

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

test_that("att_gt errors on more than 1 extra cluster variable (faster_mode)", {
  expect_error(
    att_gt(yname = "Y", data = data_eh, tname = "period", idname = "id",
           gname = "G", clustervars = c("id", "G", "period"),
           bstrap = FALSE, faster_mode = TRUE),
    "cluster|length 1"
  )
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
  sp <- did::reset.sim()
  data <- did::build_sim_dataset(sp)
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
  small_sp <- did::reset.sim()
  small_sp$n <- 50
  small_data <- did::build_sim_dataset(small_sp)
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
