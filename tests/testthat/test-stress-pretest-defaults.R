test_that("conditional pre-test resolves omitted covariates and observation ids", {
  set.seed(401)
  d <- expand.grid(id = seq_len(45), period = seq_len(4))
  d$G <- rep(c(0, 3, 4), each = 15)[d$id]
  d$X <- rnorm(45)[d$id]
  d$Y <- d$period + d$X + rnorm(nrow(d)) +
    (d$G > 0 & d$period >= d$G)

  run <- function(panel, idname, xformla) {
    set.seed(402)
    suppressMessages(capture.output(
      ans <- conditional_did_pretest(
        yname = "Y", tname = "period", gname = "G", data = d,
        panel = panel, idname = idname, xformla = xformla, biters = 31
      )
    ))
    ans
  }

  explicit <- run(TRUE, "id", ~1)
  default <- run(TRUE, "id", NULL)
  expect_equal(default$CvM, explicit$CvM)
  expect_equal(default$CvMb, explicit$CvMb)
  expect_equal(default$CvMpval, explicit$CvMpval)

  # With an intercept-only weighting function every indicator is one, so the
  # conditional statistic reduces to n times the squared pre-treatment ATTs.
  fit <- att_gt(yname = "Y", tname = "period", idname = "id", gname = "G",
                data = d, est_method = "ipw", bstrap = FALSE, cband = FALSE)
  expect_equal(default$CvM,
               fit$n * sum(fit$att[fit$t < fit$group]^2))

  # Treat these observations as independent cross sections; a per-row id and
  # the internally generated id must give the same statistic and bootstrap.
  d$row_id <- seq_len(nrow(d))
  for (formula in list(~1, ~X)) {
    explicit <- run(FALSE, "row_id", formula)
    default <- run(FALSE, NULL, formula)
    expect_gt(explicit$CvM, 0)
    expect_equal(default$CvM, explicit$CvM)
    expect_equal(default$CvMb, explicit$CvMb)
    expect_equal(default$CvMpval, explicit$CvMpval)
  }
})
