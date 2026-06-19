test_that("att_gt does not mutate caller data in either implementation", {
  set.seed(20260619)
  sp <- did::reset.sim(n = 500)
  panel_data <- did::build_sim_dataset(sp)
  rc_data <- did::build_sim_dataset(sp, panel = FALSE)

  cases <- list(
    panel_df = list(data = as.data.frame(panel_data), panel = TRUE),
    panel_dt = list(data = data.table::as.data.table(panel_data), panel = TRUE),
    rc_df = list(data = as.data.frame(rc_data), panel = FALSE),
    rc_dt = list(data = data.table::as.data.table(rc_data), panel = FALSE)
  )

  for (case_name in names(cases)) {
    for (fm in c(FALSE, TRUE)) {
      d <- cases[[case_name]]$data
      before_names <- names(d)
      before_data <- as.data.frame(d)

      suppressWarnings(suppressMessages(
        att_gt(yname = "Y", tname = "period", idname = "id", gname = "G",
               xformla = ~X, data = d, panel = cases[[case_name]]$panel,
               bstrap = FALSE, faster_mode = fm)
      ))

      expect_identical(names(d), before_names,
                       info = paste(case_name, "faster_mode", fm))
      expect_equal(as.data.frame(d), before_data,
                   ignore_attr = TRUE,
                   info = paste(case_name, "faster_mode", fm))
    }
  }
})
