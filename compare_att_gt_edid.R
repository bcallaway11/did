###############################################################################
# compare_att_gt_edid.R
#
# Simulate a staggered adoption panel and compare estimates from:
#   (1) att_gt()  -- Callaway & Sant'Anna (2021)
#   (2) edid()    -- Chen, Sant'Anna & Xie (2025), efficient DiD
#
# True ATT = 1 in all post-treatment periods (homogeneous, no dynamics).
###############################################################################

devtools::load_all(quiet = TRUE)   # loads att_gt(), edid(), and all helpers

set.seed(2025)

# ── 1. Simulate balanced staggered panel ──────────────────────────────────────
sp <- reset.sim(time.periods = 6, n = 2000)
df <- build_sim_dataset(sp_list = sp, panel = TRUE)

# build_sim_dataset returns G = 0 for never-treated (att_gt convention).
# edid() now accepts G=0 directly and converts 0 -> Inf internally,
# so we no longer need a separate G_edid column.

cat("─── Dataset summary ────────────────────────────────────────────────────\n")
cat("Periods:", sort(unique(df$period)), "\n")
cat("Cohorts (att_gt coding, 0 = never):", sort(unique(df$G)), "\n")
cat("Units:", length(unique(df$id)), " | Obs:", nrow(df), "\n\n")


# ── 2. Fit att_gt() (Callaway-Sant'Anna, no covariates) ───────────────────────
cs <- att_gt(
  yname         = "Y",
  tname         = "period",
  idname        = "id",
  gname         = "G",
  data          = df,
  control_group = "nevertreated",
  est_method    = "reg",           # outcome-regression, closest to no-cov DiD
  bstrap        = FALSE,
  cband         = FALSE,
  print_details = FALSE
)

cs_overall    <- aggte(cs, type = "simple")
cs_eventstudy <- aggte(cs, type = "dynamic")
cs_group      <- aggte(cs, type = "group")


# ── 3. Fit edid() (Chen-Sant'Anna-Xie, no covariates, PT-All) ─────────────────
# Note: gname = "G" directly -- edid() auto-converts G=0 to Inf internally.
ed <- edid(
  data          = df,
  yname         = "Y",
  idname        = "id",
  tname         = "period",
  gname         = "G",
  pt_assumption = "all",
  control_group = "nevertreated",
  alp           = 0.05
)


# ── 4. Compare cell-level ATT(g, t) ──────────────────────────────────────────
cat("─── Cell-level ATT(g,t) comparison ─────────────────────────────────────\n")

# CS cell-level estimates
cs_gt <- data.frame(
  group  = cs$group,
  time   = cs$t,
  att_cs = round(cs$att, 4),
  se_cs  = round(cs$se, 4),
  is_pre = cs$t < cs$group
)

# EDiD cell-level estimates (post-treatment only, drop skipped cells)
ed_gt <- ed$att_gt[!is.na(ed$att_gt$att), c("group", "time", "att", "se", "is_pre")]
names(ed_gt)[3:4] <- c("att_ed", "se_ed")
ed_gt$att_ed <- round(ed_gt$att_ed, 4)
ed_gt$se_ed  <- round(ed_gt$se_ed,  4)

# Merge on (group, time)
cmp_gt <- merge(cs_gt, ed_gt, by = c("group", "time", "is_pre"), all = TRUE)
cmp_gt <- cmp_gt[order(cmp_gt$is_pre, cmp_gt$group, cmp_gt$time), ]

# Highlight post-treatment cells
post_gt <- subset(cmp_gt, !is_pre)
pre_gt  <- subset(cmp_gt,  is_pre)

cat("\nPost-treatment cells (ATT should \u2248 1):\n")
print(post_gt, row.names = FALSE)

cat("\nPre-treatment placebo cells (ATT should \u2248 0):\n")
print(pre_gt, row.names = FALSE)

cat("\nMean ATT diff (post): CS - EDiD =",
    round(mean(post_gt$att_cs - post_gt$att_ed, na.rm = TRUE), 4), "\n")


# ── 5. Compare overall ATT ────────────────────────────────────────────────────
cat("\n─── Overall ATT ─────────────────────────────────────────────────────────\n")
cat(sprintf("  %-20s  ATT = %6.4f   SE = %6.4f   95%% CI = [%6.4f, %6.4f]\n",
            "att_gt (CS-simple):",
            cs_overall$overall.att,
            cs_overall$overall.se,
            cs_overall$overall.att - qnorm(0.975) * cs_overall$overall.se,
            cs_overall$overall.att + qnorm(0.975) * cs_overall$overall.se))
cat(sprintf("  %-20s  ATT = %6.4f   SE = %6.4f   95%% CI = [%6.4f, %6.4f]\n",
            "edid (PT-Post):",
            ed$overall$att,
            ed$overall$se,
            ed$overall$ci_lower,
            ed$overall$ci_upper))
cat("  True ATT = 1\n")


# ── 6. Compare event-study ────────────────────────────────────────────────────
cat("\n─── Event-study (relative-time ATTs) ───────────────────────────────────\n")

# CS event-study
cs_es <- data.frame(
  rel_time = cs_eventstudy$egt,
  att_cs   = round(cs_eventstudy$att.egt, 4),
  se_cs    = round(cs_eventstudy$se.egt,  4)
)

# EDiD event-study
ed_es_names <- names(ed$event_study)
ed_es <- data.frame(
  rel_time = as.numeric(sub("^e", "", ed_es_names)),
  att_ed   = round(vapply(ed$event_study, `[[`, numeric(1), "att"),  4),
  se_ed    = round(vapply(ed$event_study, `[[`, numeric(1), "se"),   4)
)

cmp_es <- merge(cs_es, ed_es, by = "rel_time", all = TRUE)
cmp_es <- cmp_es[order(cmp_es$rel_time), ]
print(cmp_es, row.names = FALSE)


# ── 7. Compare group-level ATTs ───────────────────────────────────────────────
cat("\n─── Group-level ATTs ────────────────────────────────────────────────────\n")

cs_grp <- data.frame(
  group  = cs_group$egt,
  att_cs = round(cs_group$att.egt, 4),
  se_cs  = round(cs_group$se.egt,  4)
)

ed_grp_names <- names(ed$group)
ed_grp <- data.frame(
  group  = as.numeric(sub("^g", "", ed_grp_names)),
  att_ed = round(vapply(ed$group, `[[`, numeric(1), "att"), 4),
  se_ed  = round(vapply(ed$group, `[[`, numeric(1), "se"),  4)
)

cmp_grp <- merge(cs_grp, ed_grp, by = "group", all = TRUE)
cmp_grp <- cmp_grp[order(cmp_grp$group), ]
print(cmp_grp, row.names = FALSE)


# ── 8. SE comparison: who is more efficient? ─────────────────────────────────
cat("\n─── Efficiency comparison (post-treatment cells only) ──────────────────\n")
post_with_both <- subset(post_gt, !is.na(att_cs) & !is.na(att_ed))
post_with_both$se_ratio <- round(post_with_both$se_cs / post_with_both$se_ed, 3)
cat("SE ratio = se(CS) / se(EDiD).  Ratio > 1 means EDiD is more efficient.\n")
print(post_with_both[, c("group", "time", "se_cs", "se_ed", "se_ratio")],
      row.names = FALSE)
cat(sprintf("\nMedian SE ratio (CS / EDiD): %.3f\n",
            median(post_with_both$se_ratio, na.rm = TRUE)))


# ── 9. Quick plot: event-study side by side ───────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  es_plot <- rbind(
    transform(cmp_es[!is.na(cmp_es$att_cs), ],
              att      = att_cs,
              se       = se_cs,
              Estimator = "att_gt (CS)")[, c("rel_time", "att", "se", "Estimator")],
    transform(cmp_es[!is.na(cmp_es$att_ed), ],
              att      = att_ed,
              se       = se_ed,
              Estimator = "edid (CSX)"  )[, c("rel_time", "att", "se", "Estimator")]
  )
  es_plot$ci_lo <- es_plot$att - qnorm(0.975) * es_plot$se
  es_plot$ci_hi <- es_plot$att + qnorm(0.975) * es_plot$se

  # Dodge slightly so CIs don't overlap
  es_plot$rel_time_jit <- ifelse(es_plot$Estimator == "att_gt (CS)",
                                  es_plot$rel_time - 0.1,
                                  es_plot$rel_time + 0.1)

  p <- ggplot(es_plot, aes(x = rel_time_jit, y = att,
                            color = Estimator, fill = Estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black", alpha = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(aes(x = rel_time_jit)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15) +
    scale_x_continuous(breaks = sort(unique(es_plot$rel_time)),
                       labels = sort(unique(es_plot$rel_time))) +
    labs(
      title    = "Event-study: att_gt (CS) vs edid (CSX)",
      subtitle = "Dotted line = true ATT = 1",
      x        = "Relative time (periods since first treatment)",
      y        = "ATT estimate",
      color    = NULL, fill = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")

  print(p)
  cat("\nEvent-study plot printed.\n")
}

cat("\n─── Done ────────────────────────────────────────────────────────────────\n")


###############################################################################
# Part II: R edid() vs Python diff-diff on the canonical diff-diff dataset
#
# Dataset: generate_staggered_data(n_units=300, n_periods=10, ATT=2.0, seed=42)
#   Cohorts (first_treat): 0=never-treated, 3, 5, 7
#   Periods: 0–9 (0-indexed; period 0 = period_1, excluded from DiD pairs)
#   True ATT = 2.0
#
# Python benchmark CSVs are generated by:
#   python3 benchmark/generate_diffdiff_benchmark.py
#
# If the CSVs are missing this section is skipped with a warning.
###############################################################################

rm(list = ls())

cat("\n\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat(" PART II: R edid() vs Python diff-diff (canonical diff-diff dataset)  \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ── Locate benchmark CSVs ─────────────────────────────────────────────────────
bench_dir  <- file.path(getwd(), "benchmark", "data")
panel_path <- file.path(bench_dir, "diffdiff_panel.csv")

if (!file.exists(panel_path)) {
  warning(
    "Benchmark CSVs not found at '", bench_dir, "'.\n",
    "Run:  python3 benchmark/generate_diffdiff_benchmark.py\n",
    "Then re-run this script to see the diff-diff comparison."
  )
} else {

  # ── 10. Load data and fit R edid() ─────────────────────────────────────────
  cat("─── 10. Load diff-diff panel and fit R edid() ─────────────────────────\n")
  dd_panel <- read.csv(panel_path)

  cat("Dataset summary:\n")
  cat("  Periods  :", sort(unique(dd_panel$period)), "\n")
  cat("  Cohorts  :", sort(unique(dd_panel$first_treat)),
      " (0 = never-treated)\n")
  cat("  Units    :", length(unique(dd_panel$unit)),
      " | Obs:", nrow(dd_panel), "\n\n")

  # first_treat == 0 => never-treated; edid() auto-converts 0 -> Inf
  ed_dd <- edid(
    data          = dd_panel,
    yname         = "outcome",
    idname        = "unit",
    tname         = "period",
    gname         = "first_treat",
    pt_assumption = "all",
    control_group = "nevertreated",
    alp           = 0.05
  )

  # ── 11. Load Python results ────────────────────────────────────────────────
  py_gt  <- read.csv(file.path(bench_dir, "py_edid_all_gt.csv"))
  py_ov  <- read.csv(file.path(bench_dir, "py_edid_all_overall.csv"))
  py_es  <- read.csv(file.path(bench_dir, "py_edid_all_es.csv"))
  py_grp <- read.csv(file.path(bench_dir, "py_edid_all_grp.csv"))

  # ── 12. Cell-level ATT(g,t) comparison ─────────────────────────────────────
  cat("─── 11. Cell-level ATT(g,t): R edid() vs Python diff-diff ────────────\n")

  # Python GT columns: group, time, effect, se, t_stat, p_value, conf_int_lower, conf_int_upper
  py_gt_clean <- data.frame(
    group  = py_gt$group,
    time   = py_gt$time,
    att_py = round(py_gt$effect, 4),
    se_py  = round(py_gt$se,     4)
  )

  # R edid GT: group, time, att, se, is_pre
  r_gt_clean <- ed_dd$att_gt[!is.na(ed_dd$att_gt$att),
                              c("group", "time", "att", "se", "is_pre")]
  names(r_gt_clean)[3:4] <- c("att_r", "se_r")
  r_gt_clean$att_r <- round(r_gt_clean$att_r, 4)
  r_gt_clean$se_r  <- round(r_gt_clean$se_r,  4)

  cmp_dd <- merge(py_gt_clean, r_gt_clean, by = c("group", "time"), all = TRUE)
  cmp_dd <- cmp_dd[order(cmp_dd$is_pre, cmp_dd$group, cmp_dd$time), ]

  # Absolute difference
  cmp_dd$diff_att <- round(abs(cmp_dd$att_py - cmp_dd$att_r), 5)
  cmp_dd$diff_se  <- round(abs(cmp_dd$se_py  - cmp_dd$se_r),  5)

  cat("\nPost-treatment cells (ATT should \u2248 2.0):\n")
  print(subset(cmp_dd, !is_pre)[,
        c("group", "time", "att_py", "att_r", "diff_att",
          "se_py", "se_r", "diff_se")],
        row.names = FALSE)

  cat("\nPre-treatment placebo cells (ATT should \u2248 0):\n")
  print(subset(cmp_dd, is_pre)[,
        c("group", "time", "att_py", "att_r", "diff_att",
          "se_py", "se_r", "diff_se")],
        row.names = FALSE)

  post_dd <- subset(cmp_dd, !is_pre)
  cat(sprintf("\nMax |diff| in ATT (post): %.5f\n",
              max(post_dd$diff_att, na.rm = TRUE)))
  cat(sprintf("Max |diff| in SE  (post): %.5f\n",
              max(post_dd$diff_se,  na.rm = TRUE)))

  # ── 13. Overall ATT comparison ─────────────────────────────────────────────
  cat("\n─── 12. Overall ATT comparison ─────────────────────────────────────────\n")
  cat(sprintf("  %-25s  ATT = %7.4f   SE = %7.4f   95%% CI = [%7.4f, %7.4f]\n",
              "Python diff-diff (EDiD):",
              py_ov$att, py_ov$se, py_ov$ci_lower, py_ov$ci_upper))
  cat(sprintf("  %-25s  ATT = %7.4f   SE = %7.4f   95%% CI = [%7.4f, %7.4f]\n",
              "R edid() (PT-All):",
              ed_dd$overall$att,
              ed_dd$overall$se,
              ed_dd$overall$ci_lower,
              ed_dd$overall$ci_upper))
  cat("  True ATT = 2.0\n")
  cat(sprintf("  |diff| ATT : %.5f\n",
              abs(py_ov$att - ed_dd$overall$att)))
  cat(sprintf("  |diff| SE  : %.5f\n",
              abs(py_ov$se  - ed_dd$overall$se)))

  # ── 14. Event-study comparison ─────────────────────────────────────────────
  cat("\n─── 13. Event-study comparison ─────────────────────────────────────────\n")

  # Python event study columns: relative_period, effect, se, ...
  py_es_df <- data.frame(
    rel_time = py_es$relative_period,
    att_py   = round(py_es$effect, 4),
    se_py    = round(py_es$se,     4)
  )

  # R event study: named list, each element has $e, $att, $se
  r_es_df <- do.call(rbind, lapply(ed_dd$event_study, function(x) {
    data.frame(rel_time = x$e, att_r = round(x$att, 4), se_r = round(x$se, 4))
  }))
  r_es_df <- r_es_df[order(r_es_df$rel_time), ]

  cmp_es_dd <- merge(py_es_df, r_es_df, by = "rel_time", all = TRUE)
  cmp_es_dd <- cmp_es_dd[order(cmp_es_dd$rel_time), ]
  cmp_es_dd$diff_att <- round(abs(cmp_es_dd$att_py - cmp_es_dd$att_r), 5)
  cmp_es_dd$diff_se  <- round(abs(cmp_es_dd$se_py  - cmp_es_dd$se_r),  5)

  print(cmp_es_dd, row.names = FALSE)
  cat(sprintf("\nMax |diff| in event-study ATT: %.5f\n",
              max(cmp_es_dd$diff_att, na.rm = TRUE)))
  cat(sprintf("Max |diff| in event-study SE : %.5f\n",
              max(cmp_es_dd$diff_se,  na.rm = TRUE)))

  # ── 15. Group-level comparison ──────────────────────────────────────────────
  cat("\n─── 14. Group-level ATT comparison ─────────────────────────────────────\n")

  py_grp_df <- data.frame(
    group  = py_grp$group,
    att_py = round(py_grp$effect, 4),
    se_py  = round(py_grp$se,     4)
  )

  r_grp_df <- do.call(rbind, lapply(ed_dd$group, function(x) {
    data.frame(group = x$group, att_r = round(x$att, 4), se_r = round(x$se, 4))
  }))

  cmp_grp_dd <- merge(py_grp_df, r_grp_df, by = "group", all = TRUE)
  cmp_grp_dd$diff_att <- round(abs(cmp_grp_dd$att_py - cmp_grp_dd$att_r), 5)
  cmp_grp_dd$diff_se  <- round(abs(cmp_grp_dd$se_py  - cmp_grp_dd$se_r),  5)
  print(cmp_grp_dd, row.names = FALSE)

  cat("\n─── Summary ─────────────────────────────────────────────────────────────\n")
  cat("Differences < 1e-4 indicate numerical equivalence.\n")
  cat("Differences > 0.01 on estimates or SE warrant investigation.\n")

}  # end if (file.exists(panel_path))

cat("\n─── Done (Part II) ──────────────────────────────────────────────────────\n")
