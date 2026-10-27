# Used data columns may share names with local preprocessing variables.
.name_collision_fit <- function(d, names, faster = TRUE, ...) {
  suppressWarnings(suppressMessages(att_gt(
    yname = names[["lemp"]], tname = names[["year"]],
    idname = names[["countyreal"]], gname = names[["first.treat"]],
    xformla = reformulate(names[["lpop"]]), data = d,
    bstrap = FALSE, cband = FALSE, faster_mode = faster, ...
  )))
}

test_that("fast preprocessing respects columns named args and tlist", {
  data(mpdta, package = "did")
  column_names <- setNames(names(mpdta), names(mpdta))
  ref <- .name_collision_fit(mpdta, column_names)
  for (local_name in c("args", "tlist")) {
    for (column in c("lemp", "year", "countyreal", "first.treat", "lpop")) {
      d <- mpdta
      names(d)[names(d) == column] <- local_name
      renamed <- column_names
      renamed[[column]] <- local_name
      actual <- .name_collision_fit(d, renamed)
      expect_equal(actual$att, ref$att, tolerance = 1e-10)
      expect_equal(actual$se, ref$se, tolerance = 1e-10)
    }
  }
})

test_that("argument name collisions work across preprocessing branches", {
  data(mpdta, package = "did")
  column_names <- setNames(names(mpdta), names(mpdta))
  for (structure in c("unbalanced", "rc", "no_never", "balance_drop", "missing")) {
    d <- mpdta
    extra <- list()
    if (structure == "unbalanced") {
      d <- d[-1L, ]
      extra$allow_unbalanced_panel <- TRUE
    } else if (structure == "rc") {
      d$countyreal <- seq_len(nrow(d))
      extra$panel <- FALSE
    } else if (structure == "no_never") {
      d <- d[d$first.treat > 0, ]
    } else if (structure == "balance_drop") {
      d <- d[!(d$first.treat == 0 & d$year == min(d$year)), ]
    } else {
      d$lemp[1L] <- Inf
    }
    ref <- do.call(.name_collision_fit, c(list(d, column_names), extra))
    for (column in c("year", "countyreal", "first.treat", "lpop")) {
      renamed_data <- d
      names(renamed_data)[names(renamed_data) == column] <- "args"
      renamed <- column_names
      renamed[[column]] <- "args"
      actual <- do.call(.name_collision_fit, c(list(renamed_data, renamed), extra))
      expect_equal(actual$att, ref$att, tolerance = 1e-10)
      expect_equal(actual$se, ref$se, tolerance = 1e-10)
    }
  }
})

test_that("group size checks use counts even if the cohort is named V1", {
  data(mpdta, package = "did")
  ids <- unique(mpdta$countyreal[mpdta$first.treat == 0])
  d <- mpdta[mpdta$first.treat > 0 | mpdta$countyreal %in% ids[1:3], ]
  names(d)[names(d) == "first.treat"] <- "V1"
  expect_error(suppressWarnings(att_gt(
    "lemp", "year", "countyreal", "V1", xformla = ~lpop,
    data = d, bstrap = FALSE, cband = FALSE
  )), "never-treated group is too small")
})


test_that("args can name sampling weights or a cluster column", {
  data(mpdta, package = "did")
  column_names <- setNames(names(mpdta), names(mpdta))
  for (role in c("weightsname", "clustervars")) {
    for (structure in c("balanced", "unbalanced", "rc")) {
      d <- mpdta
      d$sample_var <- if (role == "weightsname")
        1 + d$countyreal %% 5 else d$countyreal %% 31
      extra <- setNames(list("sample_var"), role)
      if (structure == "unbalanced") {
        d <- d[-1L, ]
        extra$allow_unbalanced_panel <- TRUE
      } else if (structure == "rc") {
        d$countyreal <- seq_len(nrow(d))
        extra$panel <- FALSE
      }
      ref <- do.call(.name_collision_fit, c(list(d, column_names), extra))
      names(d)[names(d) == "sample_var"] <- "args"
      extra[[role]] <- "args"
      actual <- do.call(.name_collision_fit, c(list(d, column_names), extra))
      expect_equal(actual$att, ref$att, tolerance = 1e-10)
      expect_equal(actual$se, ref$se, tolerance = 1e-10)
    }
  }
})


test_that("last-cohort filtering uses the cutoff argument rather than a covariate", {
  data(mpdta, package = "did")
  d <- mpdta[mpdta$first.treat > 0, ]
  column_names <- setNames(names(d), names(d))
  ref <- .name_collision_fit(d, column_names)
  names(d)[names(d) == "lpop"] <- "cutoff_t"
  column_names[["lpop"]] <- "cutoff_t"
  actual <- .name_collision_fit(d, column_names)
  expect_equal(actual$att, ref$att, tolerance = 1e-10)
  expect_equal(actual$se, ref$se, tolerance = 1e-10)
})
