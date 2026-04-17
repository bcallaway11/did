#!/usr/bin/env python3
"""
generate_diffdiff_benchmark.py

Generate the canonical diff-diff benchmark dataset and run EfficientDiD
(Chen-Sant'Anna-Xie 2025) via the diff-diff Python library.

Outputs to benchmark/data/:
  - diffdiff_panel.csv      -- the raw panel data (unit, period, outcome, first_treat)
  - py_edid_all_gt.csv      -- ATT(g,t) group-time estimates  [PT-All]
  - py_edid_all_overall.csv -- overall ATT                    [PT-All]
  - py_edid_all_es.csv      -- event-study ATTs               [PT-All]
  - py_edid_all_grp.csv     -- group-level ATTs               [PT-All]

Usage:
  python3 benchmark/generate_diffdiff_benchmark.py

Requires: diff-diff (pip install diff-diff)
"""

import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd
from diff_diff import generate_staggered_data, EfficientDiD

# ── Output directory ──────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
out_dir = os.path.join(script_dir, "data")
os.makedirs(out_dir, exist_ok=True)

# ── 1. Generate canonical dataset ────────────────────────────────────────────
print("Generating staggered panel: n_units=300, n_periods=10, ATT=2.0, seed=42 ...")
data = generate_staggered_data(
    n_units=300, n_periods=10, treatment_effect=2.0,
    dynamic_effects=False, seed=42
)

print(f"  Shape      : {data.shape}")
print(f"  Cohorts    : {sorted(data['first_treat'].unique())}")
print(f"  Periods    : {sorted(data['period'].unique())}")
print(f"  Never-treat: first_treat == 0  ({(data['first_treat']==0).sum()//10} units)")

# Save the 4 columns that edid() needs
panel_cols = ["unit", "period", "outcome", "first_treat"]
panel_path = os.path.join(out_dir, "diffdiff_panel.csv")
data[panel_cols].to_csv(panel_path, index=False)
print(f"  Saved panel -> {panel_path}\n")

# ── 2. EDiD PT-All: group-time ────────────────────────────────────────────────
print("Running EfficientDiD(pt_assumption='all', aggregate='group_time') ...")
edid_gt = EfficientDiD(pt_assumption="all").fit(
    data, outcome="outcome", unit="unit", time="period",
    first_treat="first_treat", aggregate="all"
)
gt_df = edid_gt.to_dataframe(level="group_time")
gt_path = os.path.join(out_dir, "py_edid_all_gt.csv")
gt_df.to_csv(gt_path, index=False)
print(f"  {len(gt_df)} (g,t) cells  -> {gt_path}")

# ── 3. EDiD PT-All: overall ATT ──────────────────────────────────────────────
overall_df = pd.DataFrame({
    "att"      : [edid_gt.overall_att],
    "se"       : [edid_gt.overall_se],
    "ci_lower" : [edid_gt.overall_conf_int[0]],
    "ci_upper" : [edid_gt.overall_conf_int[1]],
})
ov_path = os.path.join(out_dir, "py_edid_all_overall.csv")
overall_df.to_csv(ov_path, index=False)
print(f"  Overall ATT = {edid_gt.overall_att:.4f}  SE = {edid_gt.overall_se:.4f}  -> {ov_path}")

# ── 4. EDiD PT-All: event-study ──────────────────────────────────────────────
print("Running EfficientDiD(pt_assumption='all', aggregate='event_study') ...")
edid_es = EfficientDiD(pt_assumption="all").fit(
    data, outcome="outcome", unit="unit", time="period",
    first_treat="first_treat", aggregate="event_study"
)
es_df = edid_es.to_dataframe(level="event_study")
es_path = os.path.join(out_dir, "py_edid_all_es.csv")
es_df.to_csv(es_path, index=False)
print(f"  {len(es_df)} relative-time points  -> {es_path}")

# ── 5. EDiD PT-All: group-level ──────────────────────────────────────────────
print("Running EfficientDiD(pt_assumption='all', aggregate='group') ...")
edid_grp = EfficientDiD(pt_assumption="all").fit(
    data, outcome="outcome", unit="unit", time="period",
    first_treat="first_treat", aggregate="group"
)
grp_df = edid_grp.to_dataframe(level="group")
grp_path = os.path.join(out_dir, "py_edid_all_grp.csv")
grp_df.to_csv(grp_path, index=False)
print(f"  {len(grp_df)} groups  -> {grp_path}")

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n── Python EDiD results (summary) ──────────────────────────────────────")
print(f"  Overall ATT = {edid_gt.overall_att:.4f}   SE = {edid_gt.overall_se:.4f}   "
      f"CI = [{edid_gt.overall_conf_int[0]:.4f}, {edid_gt.overall_conf_int[1]:.4f}]")
print(f"  True ATT   = 2.0")
print()
print("── Group-time ATT(g,t) ──────────────────────────────────────────────────")
print(gt_df.to_string(index=False))
print()
print("── Event-study ─────────────────────────────────────────────────────────")
print(es_df.to_string(index=False))
print()
print("── Group-level ─────────────────────────────────────────────────────────")
print(grp_df.to_string(index=False))
print()
print("Done.")
