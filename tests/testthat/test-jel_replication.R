# Replication tests for the Callaway & Sant'Anna JEL article
# Data source: https://github.com/pedrohcgs/JEL-DiD
# These tests verify that `did` package results match the published article.
# Tests are skipped if the JEL data file is not available.

library(testthat)
library(did)

# Helper: resolve JEL data path, downloading from GitHub if needed (only called after skip_on_cran)
get_jel_data_path <- function() {
  path <- Sys.getenv("JEL_DID_DATA_PATH",
                     unset = file.path(path.expand("~"), "JEL-DiD", "data", "county_mortality_data.csv"))
  if (file.exists(path)) return(path)

  path <- file.path(tempdir(), "county_mortality_data.csv")
  if (file.exists(path)) return(path)

  tryCatch({
    download.file(
      "https://raw.githubusercontent.com/pedrohcgs/JEL-DiD/main/data/county_mortality_data.csv",
      destfile = path,
      quiet = TRUE
    )
  }, error = function(e) {
    # download failed; return path anyway, skip_if_not will handle it
  })
  path
}

# Helper to load and clean JEL data
load_jel_data <- function(path, filter_2xt = TRUE) {
  mydata <- read.csv(path, stringsAsFactors = FALSE)
  mydata$state <- substr(mydata$county, nchar(mydata$county) - 1, nchar(mydata$county))

  # Drop DC and pre-2014 adoption states
  mydata <- mydata[!(mydata$state %in% c("DC", "DE", "MA", "NY", "VT")), ]

  if (filter_2xt) {
    # For 2x2 and 2xT: keep only 2014-adopters and never-adopters
    mydata <- mydata[mydata$yaca == 2014 | is.na(mydata$yaca) | mydata$yaca > 2019, ]
  }

  # Compute covariates
  mydata$perc_white <- mydata$population_20_64_white / mydata$population_20_64 * 100
  mydata$perc_hispanic <- mydata$population_20_64_hispanic / mydata$population_20_64 * 100
  mydata$perc_female <- mydata$population_20_64_female / mydata$population_20_64 * 100
  mydata$unemp_rate <- mydata$unemp_rate * 100
  mydata$median_income <- mydata$median_income / 1000

  # Keep only needed columns
  keep_cols <- c("county_code", "year", "population_20_64", "yaca", "crude_rate_20_64",
                 "perc_female", "perc_white", "perc_hispanic", "unemp_rate", "poverty_rate", "median_income")
  mydata <- mydata[, keep_cols]

  # Drop rows with NA in non-yaca columns
  non_yaca <- setdiff(keep_cols, "yaca")
  mydata <- mydata[complete.cases(mydata[, non_yaca]), ]

  # Keep counties with data in both 2013 and 2014
  county_year_counts <- tapply(mydata$year %in% c(2013, 2014), mydata$county_code, sum)
  valid_counties <- names(county_year_counts[county_year_counts == 2])
  mydata <- mydata[mydata$county_code %in% valid_counties, ]

  # Keep counties with all 11 years of data
  county_counts <- table(mydata$county_code)
  valid_counties <- names(county_counts[county_counts == 11])
  mydata <- mydata[mydata$county_code %in% valid_counties, ]

  mydata
}


# ===========================================================================
# Table 7: 2x2 CS-DiD with Covariates (Sant'Anna-Zhao)
# ===========================================================================
test_that("JEL Table 7: 2x2 CS-DiD point estimates match", {
  skip_on_cran()
  jel_data_path <- get_jel_data_path()
  skip_if_not(file.exists(jel_data_path), "JEL-DiD data not available")

  mydata <- load_jel_data(jel_data_path, filter_2xt = TRUE)

  # Make 2x2 dataset (2013-2014 only)
  short_data <- mydata[mydata$year %in% c(2013, 2014), ]
  short_data$treat_year <- ifelse(!is.na(short_data$yaca) & short_data$yaca == 2014, 2014, 0)
  short_data$county_code <- as.numeric(short_data$county_code)

  # Population weights from 2013
  wt_2013 <- short_data[short_data$year == 2013, c("county_code", "population_20_64")]
  names(wt_2013)[2] <- "set_wt"
  short_data <- merge(short_data, wt_2013, by = "county_code")

  covs_formula <- ~perc_female + perc_white + perc_hispanic + unemp_rate + poverty_rate + median_income

  # Expected point estimates from JEL article (Table 7)
  expected <- list(
    reg_unweighted = -1.6154372119,
    ipw_unweighted = -0.8585625501,
    dr_unweighted  = -1.2256473242,
    reg_weighted   = -3.4592200594,
    ipw_weighted   = -3.8416966846,
    dr_weighted    = -3.7561045985
  )

  for (method in c("reg", "ipw", "dr")) {
    for (wt_info in list(list(name = NULL, label = "unweighted"),
                         list(name = "set_wt", label = "weighted"))) {
      res <- att_gt(
        yname = "crude_rate_20_64", tname = "year", idname = "county_code",
        gname = "treat_year", xformla = covs_formula,
        data = short_data, panel = TRUE, control_group = "nevertreated",
        bstrap = FALSE, est_method = method, weightsname = wt_info$name,
        base_period = "universal"
      )
      agg <- aggte(res, na.rm = TRUE, bstrap = FALSE)

      key <- paste0(method, "_", wt_info$label)
      expect_equal(agg$overall.att, expected[[key]], tolerance = 1e-6,
                   label = paste("Table 7:", key))
    }
  }
})


# ===========================================================================
# 2xT Event Study: No covariates, weighted, reg, nevertreated
# ===========================================================================
test_that("JEL 2xT: event study ATT(g,t) point estimates match", {
  skip_on_cran()
  jel_data_path <- get_jel_data_path()
  skip_if_not(file.exists(jel_data_path), "JEL-DiD data not available")

  mydata <- load_jel_data(jel_data_path, filter_2xt = TRUE)
  mydata$treat_year <- ifelse(!is.na(mydata$yaca) & mydata$yaca == 2014, 2014, 0)
  mydata$county_code <- as.numeric(mydata$county_code)

  wt_2013 <- mydata[mydata$year == 2013, c("county_code", "population_20_64")]
  names(wt_2013)[2] <- "set_wt"
  mydata <- merge(mydata, wt_2013, by = "county_code")

  mod <- att_gt(
    yname = "crude_rate_20_64", tname = "year", idname = "county_code",
    gname = "treat_year", xformla = NULL, data = mydata, panel = TRUE,
    control_group = "nevertreated", bstrap = FALSE,
    est_method = "reg", weightsname = "set_wt", base_period = "universal"
  )

  # Expected ATT(g,t) for g=2014
  expected_att <- c(
    4.1292044043,   # t=2009
    -0.5016807242,  # t=2010
    2.7531791360,   # t=2011
    2.7804626426,   # t=2012
    0.0,            # t=2013 (base period)
    -2.5628745138,  # t=2014
    -1.6973291127,  # t=2015
    0.2189009815,   # t=2016
    -0.8133358354,  # t=2017
    -1.1532954495,  # t=2018
    1.7866564429    # t=2019
  )

  expect_equal(mod$att, expected_att, tolerance = 1e-6)

  # Dynamic aggregation
  es <- aggte(mod, type = "dynamic", bstrap = FALSE)
  expect_equal(es$att.egt, expected_att, tolerance = 1e-6)

  # Overall ATT for e in {0,5}
  agg <- aggte(mod, type = "dynamic", min_e = 0, max_e = 5, bstrap = FALSE)
  expect_equal(agg$overall.att, -0.7035462478, tolerance = 1e-6)
})


# ===========================================================================
# 2xT with covariates: 3 estimation methods
# ===========================================================================
test_that("JEL 2xT: event study with covariates matches across methods", {
  skip_on_cran()
  jel_data_path <- get_jel_data_path()
  skip_if_not(file.exists(jel_data_path), "JEL-DiD data not available")

  mydata <- load_jel_data(jel_data_path, filter_2xt = TRUE)
  mydata$treat_year <- ifelse(!is.na(mydata$yaca) & mydata$yaca == 2014, 2014, 0)
  mydata$county_code <- as.numeric(mydata$county_code)

  wt_2013 <- mydata[mydata$year == 2013, c("county_code", "population_20_64")]
  names(wt_2013)[2] <- "set_wt"
  mydata <- merge(mydata, wt_2013, by = "county_code")

  covs_formula <- ~perc_female + perc_white + perc_hispanic + unemp_rate + poverty_rate + median_income

  for (method in c("reg", "ipw", "dr")) {
    res <- att_gt(
      yname = "crude_rate_20_64", tname = "year", idname = "county_code",
      gname = "treat_year", xformla = covs_formula,
      data = mydata, panel = TRUE, control_group = "nevertreated",
      bstrap = FALSE, est_method = method, weightsname = "set_wt",
      base_period = "universal"
    )

    es <- aggte(res, type = "dynamic", na.rm = TRUE, bstrap = FALSE)

    # Base period (e = -1) should be exactly 0
    base_idx <- which(es$egt == -1)
    expect_equal(es$att.egt[base_idx], 0, label = paste("2xT covs", method, "e=-1"))

    # All estimates should be finite
    expect_true(all(is.finite(es$att.egt[!is.na(es$att.egt)])),
                label = paste("2xT covs", method, "finite"))
  }
})


# ===========================================================================
# GxT: No covariates, weighted, notyettreated
# ===========================================================================
test_that("JEL GxT: staggered event study without covariates matches", {
  skip_on_cran()
  jel_data_path <- get_jel_data_path()
  skip_if_not(file.exists(jel_data_path), "JEL-DiD data not available")

  mydata <- load_jel_data(jel_data_path, filter_2xt = FALSE)  # GxT keeps all groups
  mydata$treat_year <- ifelse(!is.na(mydata$yaca) & mydata$yaca <= 2019, mydata$yaca, 0)
  mydata$county_code <- as.numeric(mydata$county_code)

  wt_2013 <- mydata[mydata$year == 2013, c("county_code", "population_20_64")]
  names(wt_2013)[2] <- "set_wt"
  mydata <- merge(mydata, wt_2013, by = "county_code")

  mod <- att_gt(
    yname = "crude_rate_20_64", tname = "year", idname = "county_code",
    gname = "treat_year", xformla = NULL, data = mydata, panel = TRUE,
    control_group = "notyettreated", bstrap = FALSE,
    weightsname = "set_wt", base_period = "universal"
  )

  # Dynamic aggregation
  es <- aggte(mod, type = "dynamic", bstrap = FALSE)

  # Expected dynamic ATT at key event times
  expected_dynamic <- c(
    `e=-5`  =  2.2186783578,
    `e=-4`  =  0.8579225109,
    `e=-3`  =  1.9161499399,
    `e=-2`  =  2.5644742860,
    `e=-1`  =  0.0,
    `e=0`   = -1.6545648988,
    `e=1`   = -0.2616435324,
    `e=2`   =  1.7055625922,
    `e=3`   = -0.5405232028,
    `e=4`   = -0.5148819184,
    `e=5`   =  1.7866564429
  )

  # Match e = -5 to e = 5
  for (e_val in -5:5) {
    idx <- which(es$egt == e_val)
    key <- paste0("e=", e_val)
    expect_equal(es$att.egt[idx], expected_dynamic[[key]], tolerance = 1e-6,
                 label = paste("GxT no covs", key))
  }

  # Overall ATT for e in {0,5}
  agg <- aggte(mod, type = "dynamic", min_e = 0, max_e = 5, bstrap = FALSE)
  expect_equal(agg$overall.att, 0.0867675805, tolerance = 1e-6,
               label = "GxT no covs overall ATT (e=0:5)")
})


# ===========================================================================
# GxT: With covariates, DR, weighted, notyettreated
# ===========================================================================
test_that("JEL GxT: staggered event study with DR covariates matches", {
  skip_on_cran()
  jel_data_path <- get_jel_data_path()
  skip_if_not(file.exists(jel_data_path), "JEL-DiD data not available")

  mydata <- load_jel_data(jel_data_path, filter_2xt = FALSE)
  mydata$treat_year <- ifelse(!is.na(mydata$yaca) & mydata$yaca <= 2019, mydata$yaca, 0)
  mydata$county_code <- as.numeric(mydata$county_code)

  wt_2013 <- mydata[mydata$year == 2013, c("county_code", "population_20_64")]
  names(wt_2013)[2] <- "set_wt"
  mydata <- merge(mydata, wt_2013, by = "county_code")

  covs_formula <- ~perc_female + perc_white + perc_hispanic + unemp_rate + poverty_rate + median_income

  mod <- att_gt(
    yname = "crude_rate_20_64", tname = "year", idname = "county_code",
    gname = "treat_year", xformla = covs_formula,
    data = mydata, panel = TRUE, control_group = "notyettreated",
    bstrap = FALSE, est_method = "dr", weightsname = "set_wt",
    base_period = "universal"
  )

  # Overall ATT for e in {0,5}
  agg <- aggte(mod, type = "dynamic", min_e = 0, max_e = 5, bstrap = FALSE)
  expect_equal(agg$overall.att, -2.2469982988, tolerance = 1e-6,
               label = "GxT DR covs overall ATT (e=0:5)")

  # Dynamic aggregation at key event times
  es <- aggte(mod, type = "dynamic", bstrap = FALSE)

  expected_dynamic <- c(
    `e=-5`  =  2.6684811691,
    `e=-4`  =  2.1333836537,
    `e=-3`  =  2.8574389179,
    `e=-2`  =  2.9129673584,
    `e=-1`  =  0.0,
    `e=0`   = -1.4929553032,
    `e=1`   = -1.9693922262,
    `e=2`   = -2.7250659699,
    `e=3`   = -5.0556625460,
    `e=4`   = -4.7161630370,
    `e=5`   =  2.4772492894
  )

  for (e_val in -5:5) {
    idx <- which(es$egt == e_val)
    key <- paste0("e=", e_val)
    expect_equal(es$att.egt[idx], expected_dynamic[[key]], tolerance = 1e-6,
                 label = paste("GxT DR covs", key))
  }
})


# ===========================================================================
# Consistency: faster_mode should match regular mode for all JEL analyses
# ===========================================================================
test_that("JEL: faster_mode matches regular mode", {
  skip_on_cran()
  jel_data_path <- get_jel_data_path()
  skip_if_not(file.exists(jel_data_path), "JEL-DiD data not available")

  mydata <- load_jel_data(jel_data_path, filter_2xt = TRUE)

  # 2x2 data
  short_data <- mydata[mydata$year %in% c(2013, 2014), ]
  short_data$treat_year <- ifelse(!is.na(short_data$yaca) & short_data$yaca == 2014, 2014, 0)
  short_data$county_code <- as.numeric(short_data$county_code)

  wt_2013 <- short_data[short_data$year == 2013, c("county_code", "population_20_64")]
  names(wt_2013)[2] <- "set_wt"
  short_data <- merge(short_data, wt_2013, by = "county_code")

  covs_formula <- ~perc_female + perc_white + perc_hispanic + unemp_rate + poverty_rate + median_income

  # Test DR weighted
  res_slow <- att_gt(
    yname = "crude_rate_20_64", tname = "year", idname = "county_code",
    gname = "treat_year", xformla = covs_formula,
    data = short_data, panel = TRUE, control_group = "nevertreated",
    bstrap = FALSE, est_method = "dr", weightsname = "set_wt",
    base_period = "universal", faster_mode = FALSE
  )
  res_fast <- att_gt(
    yname = "crude_rate_20_64", tname = "year", idname = "county_code",
    gname = "treat_year", xformla = covs_formula,
    data = short_data, panel = TRUE, control_group = "nevertreated",
    bstrap = FALSE, est_method = "dr", weightsname = "set_wt",
    base_period = "universal", faster_mode = TRUE
  )

  expect_equal(res_slow$att, res_fast$att, tolerance = 1e-10,
               label = "2x2 DR weighted: ATTs match")

  agg_slow <- aggte(res_slow, na.rm = TRUE, bstrap = FALSE)
  agg_fast <- aggte(res_fast, na.rm = TRUE, bstrap = FALSE)
  expect_equal(agg_slow$overall.att, agg_fast$overall.att, tolerance = 1e-10,
               label = "2x2 DR weighted: aggregate ATT matches")
})
