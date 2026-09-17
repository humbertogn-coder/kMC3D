#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kp_focused_regimes_cfg.py
═════════════════════════
Answers Dr. Balbuena's questions 2 and 3 for the property (or properties)
that the sensitivity ranking flagged as most responsive to k_p:

  Q2  What happens when k_p is BELOW the reference (56.25)?
  Q3  What happens for INTERMEDIATE k_p between 56.25 and 1000?

Two kinds of output:

  (1) TREND across the whole sweep for the target property, with the three
      regimes shaded:  below-reference | reference→comparison | above.
      A monotonic / logistic fit is drawn so the intermediate behaviour is
      explicit (is it linear, saturating, or does it have a knee?).

  (2) PER-CYCLE curves for a few REPRESENTATIVE k_p in each regime, so you
      can see how the property evolves over cycling, not just its endpoint.
      Representatives (nearest available folder): the reference (56.25), two
      below it, three intermediate, and the comparison (1000).

The target property is selected automatically from kp_sensitivity_ranking.csv
(top high-sensitivity property whose per-cycle CSV exists). You can override
it with --property "<group:column>" and --per-cycle-col "<col in per-cycle CSV>".

Per-cycle sources understood out of the box (property group -> CSV, col):
  CE/Q      -> ce_capacity_per_cycle.csv        (CE, capacity_C)
  BuriedLi  -> li_buried_per_cycle.csv           (li_buried_kinetic,
                                                  n_li_buried_dynamic)
  SEI       -> sei_growth_per_cycle.csv          (SEI_thickness_A)
  Rxn       -> reaction_tracking_per_cycle.csv   (any Rxn*/*_xyz column)
  Anode     -> anode_consumption_per_cycle.csv   (total_sei_anode, ...)

Outputs to CFG['out_dir']:
  fig_<tag>_trend_regimes_vs_kp.png
  fig_<tag>_percycle_by_regime.png
  kp_focused_<tag>_summary.csv

Run:
  python kp_focused_regimes_cfg.py --config config_plating
  python kp_focused_regimes_cfg.py --config config_plating \
         --property "SEI:SEI_thickness_A_final" \
         --per-cycle-col SEI_thickness_A
"""

import argparse
import warnings
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import kmc_io as io

warnings.filterwarnings("ignore")

# group -> (per-cycle csv, default per-cycle column, default ML-feature col)
GROUP_MAP = {
    "CE/Q":     ("ce_capacity_per_cycle.csv",      "CE",
                 "CE_mean"),
    "BuriedLi": ("li_buried_per_cycle.csv",        "li_buried_kinetic",
                 "li_buried_kinetic_final"),
    "SEI":      ("sei_growth_per_cycle.csv",       "SEI_thickness_A",
                 "SEI_thickness_A_final"),
    "Rxn":      ("reaction_tracking_per_cycle.csv", "RxnStripping",
                 "RxnStripping"),
    "Anode":    ("anode_consumption_per_cycle.csv", "total_sei_anode",
                 "total_sei_anode_final"),
    "Barrier":  (None, None, "ea_asymmetry_ratio"),
}


def pick_target(cfg, out, override):
    """Return (group, ml_col, per_cycle_csv, per_cycle_col)."""
    if override:
        group = override.split(":")[0]
        ml_col = override
        pc_csv, pc_col_default, _ = GROUP_MAP.get(group, (None, None, None))
        pc_col = _derive_pc_col(out, pc_csv, ml_col, pc_col_default)
        return group, ml_col, pc_csv, pc_col
    rank_p = out / "kp_sensitivity_ranking.csv"
    if not rank_p.is_file():
        raise FileNotFoundError(
            "kp_sensitivity_ranking.csv not found. Run "
            "kp_sensitivity_ranking_cfg.py first, or pass --property.")
    rank = pd.read_csv(rank_p)
    # prefer high-sensitivity properties whose group has a per-cycle CSV
    cand = rank[rank.get("high_sensitivity", False) == True] \
        if "high_sensitivity" in rank.columns else rank
    for _, r in pd.concat([cand, rank]).drop_duplicates("property").iterrows():
        group = str(r["property"]).split(":")[0]
        pc_csv, pc_col, _ = GROUP_MAP.get(group, (None, None, None))
        if pc_csv and (out / pc_csv).is_file():
            pc_col = _derive_pc_col(out, pc_csv, r["property"], pc_col)
            return group, r["property"], pc_csv, pc_col
    # fallback: top row regardless of per-cycle availability
    r = rank.iloc[0]
    group = str(r["property"]).split(":")[0]
    pc_csv, pc_col, _ = GROUP_MAP.get(group, (None, None, None))
    pc_col = _derive_pc_col(out, pc_csv, r["property"], pc_col)
    return group, r["property"], pc_csv, pc_col


def _derive_pc_col(out, pc_csv, ml_col, default_pc_col):
    """If the target property's base name (e.g. 'RxnLiSurface' out of
    'Rxn:RxnLiSurface', or 'li_buried_kinetic' out of the '_final' feature)
    appears as a column in the per-cycle CSV, use it. Otherwise fall back
    to the group default."""
    if pc_csv is None:
        return default_pc_col
    p = out / pc_csv
    if not p.is_file():
        return default_pc_col
    cols = set(pd.read_csv(p, nrows=1).columns)
    base = ml_col.split(":", 1)[-1]
    # try the raw base, then strip common aggregation suffixes
    candidates = [base]
    for suf in ("_final", "_mean", "_max", "_at_db_cycle", "_A_final",
                "_A_mean"):
        if base.endswith(suf):
            candidates.append(base[: -len(suf)])
    candidates.append(base.replace("_final", "").replace("_mean", ""))
    for c in candidates:
        if c in cols:
            return c
    return default_pc_col


def logistic_log10(x, L, k, x0, b):
    return L / (1.0 + np.exp(-k * (np.log10(x) - x0))) + b


def nearest_cases(cfg, targets):
    """For each target k_p value, return the (case_key, actual_k0) whose
    k0 is closest."""
    cases = cfg["cases"]
    out = []
    for t in targets:
        key, k0 = min(cases, key=lambda ck: abs(ck[1] - t))
        out.append((key, k0))
    # dedupe preserving order
    seen = set(); uniq = []
    for key, k0 in out:
        if key not in seen:
            seen.add(key); uniq.append((key, k0))
    return uniq


def trend_figure(cfg, out, ml_col, tag):
    """Property (from ML-feature CSV) vs k_p across the full sweep, with
    the three regimes shaded and a logistic fit."""
    group = ml_col.split(":")[0]
    fname = {
        "CE/Q": "ce_capacity_ML_features.csv",
        "BuriedLi": "li_buried_ML_features.csv",
        "Rxn": "reaction_ML_features.csv",
        "Anode": "anode_ML_features.csv",
        "Barrier": "activation_barrier_ML_features.csv",
        "SEI": "sei_morphology_ML_features.csv",
    }.get(group)
    raw_col = ml_col.split(":", 1)[1]
    p = out / fname
    if not p.is_file():
        print(f"  trend: {fname} missing; skip"); return None
    df = pd.read_csv(p)
    if raw_col not in df.columns or "k0" not in df.columns:
        print(f"  trend: {raw_col} not in {fname}; skip"); return None
    d = df[["k0", raw_col]].dropna().sort_values("k0")
    x = d["k0"].values.astype(float); y = d[raw_col].values.astype(float)

    plt.rcParams.update(cfg["plot_style"])
    sym = cfg.get("rate_symbol_math", "$k_p$")
    kref = cfg["kp_reference"]; khigh = cfg["ref_high_k0"]
    fig, ax = plt.subplots(figsize=(9.2, 5.4))

    xmin, xmax = x.min() * 0.8, x.max() * 1.25
    ax.axvspan(xmin, kref, color="#AED6F1", alpha=0.25)
    ax.axvspan(kref, khigh, color="#ABEBC6", alpha=0.28)
    ax.axvspan(khigh, xmax, color="#F5B7B1", alpha=0.22)
    ax.text(np.sqrt(xmin * kref), ax.get_ylim()[1], "below ref",
            ha="center", va="bottom", fontsize=8.5, color="#2471A3")
    ax.text(np.sqrt(kref * khigh), ax.get_ylim()[1], "intermediate",
            ha="center", va="bottom", fontsize=8.5, color="#1E8449")

    ax.plot(x, y, "o", color="#1B2631", markersize=6,
            markeredgecolor="white", markeredgewidth=0.6, zorder=5,
            label="Per-$k_\\mathrm{p}$ value")
    # logistic fit if it behaves
    try:
        p0 = (np.nanmax(y) - np.nanmin(y), 2.0,
              np.log10(np.median(x)), np.nanmin(y))
        popt, _ = curve_fit(logistic_log10, x, y, p0=p0, maxfev=20000)
        xf = np.logspace(np.log10(xmin), np.log10(xmax), 300)
        yf = logistic_log10(xf, *popt)
        yhat = logistic_log10(x, *popt)
        ss = 1 - np.sum((y - yhat)**2) / np.sum((y - y.mean())**2)
        ax.plot(xf, yf, "-", color="#8E44AD", lw=2,
                label=f"Logistic fit (R²={ss:.3f})")
    except Exception:
        pass
    ax.axvline(kref, color="#27AE60", ls="--", lw=1.4)
    ax.axvline(khigh, color="#C0392B", ls=":", lw=1.4)
    ax.set_xscale("log"); ax.set_xlim(xmin, xmax)
    ax.set_xlabel(f"{sym}  ({cfg['rate_units_note']})")
    ax.set_ylabel(raw_col)
    ax.set_title(f"{raw_col} across the {sym} sweep", fontweight="bold")
    ax.grid(alpha=0.3, which="both", lw=0.5)
    ax.legend(loc="best")
    plt.tight_layout()
    op = out / f"fig_{tag}_trend_regimes_vs_kp.png"
    fig.savefig(op, bbox_inches="tight", facecolor="white"); plt.close()
    print(f"  saved {op.name}")
    return d


def percycle_figure(cfg, out, pc_csv, pc_col, tag):
    """Per-cycle curves for representative k_p in each regime."""
    if pc_csv is None:
        print("  per-cycle: property has no per-cycle CSV; skip"); return
    p = out / pc_csv
    if not p.is_file():
        print(f"  per-cycle: {pc_csv} missing; skip"); return
    df = pd.read_csv(p)
    if pc_col not in df.columns:
        print(f"  per-cycle: {pc_col} not in {pc_csv}; skip"); return

    kref = cfg["kp_reference"]; khigh = cfg["ref_high_k0"]
    below = sorted(cfg["kp_below_reference"])
    inter = sorted(cfg["kp_intermediate"])
    # representatives: 2 below, reference, 3 spread across intermediate, comparison
    reps_below = nearest_cases(cfg, [below[0], below[len(below)//2]])
    reps_inter = nearest_cases(cfg, [inter[len(inter)//4],
                                     inter[len(inter)//2],
                                     inter[3*len(inter)//4]])
    reps_ref  = nearest_cases(cfg, [kref])
    reps_high = nearest_cases(cfg, [khigh])

    plt.rcParams.update(cfg["plot_style"])
    sym = cfg.get("rate_symbol_math", "$k_p$")
    fig, ax = plt.subplots(figsize=(9.4, 5.6))

    # color families: blues for below, greens intermediate, ref black, high red
    blue = plt.cm.Blues(np.linspace(0.55, 0.9, len(reps_below)))
    green = plt.cm.Greens(np.linspace(0.5, 0.9, len(reps_inter)))

    def plot_set(reps, colors, ls="-"):
        for (key, k0), col in zip(reps, colors):
            g = io.select_by_k0(df, k0)
            if g.empty:
                continue
            ax.plot(g["cycle"], g[pc_col], ls, color=col, lw=1.8,
                    markersize=3, label=f"{sym}={k0:g}")

    plot_set(reps_below, blue)
    for (key, k0) in reps_ref:
        g = io.select_by_k0(df, k0)
        if not g.empty:
            ax.plot(g["cycle"], g[pc_col], "-", color="#111111", lw=2.6,
                    label=f"{sym}={k0:g} (reference)")
    plot_set(reps_inter, green)
    for (key, k0) in reps_high:
        g = io.select_by_k0(df, k0)
        if not g.empty:
            ax.plot(g["cycle"], g[pc_col], "-", color="#C0392B", lw=2.6,
                    label=f"{sym}={k0:g} (comparison)")

    ax.set_xlabel("Cycle number")
    ax.set_ylabel(pc_col)
    ax.set_title(f"{pc_col} vs cycle across {sym} regimes",
                 fontweight="bold")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(loc="best", fontsize=8.5, ncol=2)
    plt.tight_layout()
    op = out / f"fig_{tag}_percycle_by_regime.png"
    fig.savefig(op, bbox_inches="tight", facecolor="white"); plt.close()
    print(f"  saved {op.name}")


def regime_summary(cfg, trend_df, ml_col, tag, out):
    """Numeric summary: mean of the property in each regime + % change."""
    if trend_df is None or trend_df.empty:
        return
    raw_col = ml_col.split(":", 1)[1]
    kref = cfg["kp_reference"]; khigh = cfg["ref_high_k0"]
    x = trend_df["k0"].values; y = trend_df[raw_col].values
    below = y[x < kref]
    inter = y[(x >= kref) & (x <= khigh)]
    above = y[x > khigh]
    def stat(a):
        return (np.nan, np.nan) if len(a) == 0 else (np.mean(a), np.std(a))
    b_m, b_s = stat(below); i_m, i_s = stat(inter); a_m, a_s = stat(above)
    ref_val = float(np.interp(kref, x, y))
    rows = [{
        "property": raw_col,
        "reference_kp": kref, "reference_value": ref_val,
        "mean_below_ref": b_m, "std_below_ref": b_s, "n_below": int(len(below)),
        "mean_intermediate": i_m, "std_intermediate": i_s, "n_inter": int(len(inter)),
        "mean_above_comparison": a_m, "std_above": a_s, "n_above": int(len(above)),
        "pct_ref_to_intermediate_end":
            (float((y[x <= khigh][-1] - ref_val) / abs(ref_val) * 100)
             if abs(ref_val) > 1e-12 and (x <= khigh).any() else np.nan),
    }]
    op = out / f"kp_focused_{tag}_summary.csv"
    pd.DataFrame(rows).to_csv(op, index=False)
    print(f"  saved {op.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--property", default=None,
                    help='Override target, e.g. "SEI:SEI_thickness_A_final"')
    ap.add_argument("--per-cycle-col", default=None,
                    help="Override per-cycle column for the curves figure")
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]; out.mkdir(parents=True, exist_ok=True)

    group, ml_col, pc_csv, pc_col = pick_target(cfg, out, args.property)
    if args.per_cycle_col:
        pc_col = args.per_cycle_col
    tag = re.sub(r"[^A-Za-z0-9]+", "_", ml_col).strip("_")
    print(f"Focused regime analysis — target property: {ml_col}")
    print(f"  per-cycle: {pc_csv} column '{pc_col}'")

    trend_df = trend_figure(cfg, out, ml_col, tag)
    percycle_figure(cfg, out, pc_csv, pc_col, tag)
    regime_summary(cfg, trend_df, ml_col, tag, out)
    print("Done.")


if __name__ == "__main__":
    main()
