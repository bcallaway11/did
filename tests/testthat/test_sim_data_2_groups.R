library(DRDID)
library(BMisc)
library(tidyr)
## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test it works without covariates and only two gropups
#-----------------------------------------------------------------------------
test_that("att_gt works with 2 groups", {
  set.seed(20241017)
  n = 5000
  # Generate data with 4 periods
  p3_true = exp(0.5)
  denominator = 1 + p3_true

  #probability of each group: 2 groups in total (G=4)
  p3_true = p3_true / denominator
  pinf_true = 1 / denominator

  # Generate treatment status
  g <- as.numeric(runif(n) <= p3_true)
  g <- ifelse(g==1,3,0)


  #Potential outcomes: Yt_g:=Y_t(g)
  index_unobs_het = g
  nu <- stats::rnorm(n, mean = index_unobs_het, sd = 1)
  index_trend <- 1
  index_att <- 3
  att_rand <- stats::rnorm(n, mean = index_att, sd = 1)
  # Generate untreated potential outcomes
  Yt1_ginf = index_trend + nu + rnorm(n)
  Yt2_ginf = 2*index_trend + nu + rnorm(n)
  Yt3_ginf = 3*index_trend + nu + rnorm(n)
  Yt4_ginf = 4*index_trend + nu + rnorm(n)

  # Generate treated potential outcomes
  Yt1_g3 = index_trend + nu + rnorm(n)
  Yt2_g3 = 2*index_trend + nu + rnorm(n)
  Yt3_g3 = att_rand + 3*index_trend + nu + rnorm(n)
  Yt4_g3 = 1.5*att_rand + 4*index_trend + nu + rnorm(n)

  #Observed data
  y1 <- (g==3) * Yt1_g3 + (g==0) * Yt1_ginf
  y2 <- (g==3) * Yt2_g3 + (g==0) * Yt2_ginf
  y3 <- (g==3) * Yt3_g3 + (g==0) * Yt3_ginf
  y4 <- (g==3) * Yt4_g3 + (g==0) * Yt4_ginf

  data <- data.frame(y1,y2,y3,y4,g)
  data$id <- 1:dim(data)[1]
  # Have data in long format
  data <- data |>
    tidyr::pivot_longer(
      cols = starts_with("y"),
      names_to = "t",
      names_prefix = "y",
      values_to = "y",
    )
  data$t <- as.numeric(data$t)
  data$g <- as.numeric(data$g)

  # Run regression with never-treated and varying base period
  csdid_nt_varying <- did::att_gt(yname = "y",
                          idname = "id",
                          gname = "g",
                          tname = "t",
                          data = data,
                          control_group = "nevertreated",
                          panel = TRUE,
                          xformla = ~1,
                          bstrap = FALSE,
                          cband = FALSE,
                          est_method = "reg",
                          base_period = "varying"
  )

  # Run regression with never-treated and universal base period
  csdid_nt_universal <- did::att_gt(yname = "y",
                                  idname = "id",
                                  gname = "g",
                                  tname = "t",
                                  data = data,
                                  control_group = "nevertreated",
                                  panel = TRUE,
                                  xformla = ~1,
                                  bstrap = FALSE,
                                  cband = FALSE,
                                  est_method = "reg",
                                  base_period = "varying"
  )
  # Run regression with not-yet-treated and varying base period
  csdid_nyt_varying <- did::att_gt(yname = "y",
                                  idname = "id",
                                  gname = "g",
                                  tname = "t",
                                  data = data,
                                  control_group = "notyettreated",
                                  panel = TRUE,
                                  xformla = ~1,
                                  bstrap = FALSE,
                                  cband = FALSE,
                                  est_method = "reg",
                                  base_period = "varying"
  )

  # Run regression with not-yet-treated and universal base period
  csdid_nyt_universal <- did::att_gt(yname = "y",
                                    idname = "id",
                                    gname = "g",
                                    tname = "t",
                                    data = data,
                                    control_group = "notyettreated",
                                    panel = TRUE,
                                    xformla = ~1,
                                    bstrap = FALSE,
                                    cband = FALSE,
                                    est_method = "reg",
                                    base_period = "varying"
  )

  expect_equal(csdid_nyt_universal$att[2], 3, tol=.1)
  expect_equal(csdid_nt_universal$att[2], 3, tol=.1)
  expect_equal(csdid_nyt_varying$att[2], 3, tol=.1)
  expect_equal(csdid_nt_varying$att[2], 3, tol=.1)

  expect_equal(csdid_nyt_universal$att[2], csdid_nt_universal$att[2], tol=.0001)
  expect_equal(csdid_nyt_universal$att[2], csdid_nyt_varying$att[2], tol=.0001)
  expect_equal(csdid_nyt_universal$att[2], csdid_nt_varying$att[2], tol=.0001)


})

