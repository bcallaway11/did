# =============================================================================
# allow_unbalanced_panel = TRUE is silently reset to FALSE when the panel turns
# out to be balanced. The fast path has always announced which branch it took;
# the slow path must emit the identical messages.
# =============================================================================

test_that("both modes announce that a balanced panel resets allow_unbalanced_panel", {
  set.seed(20260814)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)

  for (fm in c(TRUE, FALSE)) {
    expect_message(
      suppressWarnings(att_gt(yname = "Y", xformla = ~X, data = d, tname = "period",
        idname = "id", gname = "G", allow_unbalanced_panel = TRUE,
        faster_mode = fm, bstrap = FALSE)),
      "You have a balanced panel"
    )
  }
})

test_that("both modes announce that an unbalanced panel is used as such", {
  set.seed(20260814)
  sp <- did::reset.sim(time.periods = 4)
  d <- did::build_sim_dataset(sp)
  d_ub <- d[-sample(seq_len(nrow(d)), floor(0.05 * nrow(d))), ]

  for (fm in c(TRUE, FALSE)) {
    expect_message(
      suppressWarnings(att_gt(yname = "Y", xformla = ~X, data = d_ub, tname = "period",
        idname = "id", gname = "G", allow_unbalanced_panel = TRUE,
        faster_mode = fm, bstrap = FALSE)),
      "You have an unbalanced panel"
    )
  }
})
