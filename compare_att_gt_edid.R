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
  pt_assumption = "post",
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
