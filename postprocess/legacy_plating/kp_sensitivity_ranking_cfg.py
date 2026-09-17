#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kp_sensitivity_ranking_cfg.py
═════════════════════════════
STEP 0 for Dr. Balbuena's request: "First check which properties have the
largest changes and then do an analysis for that property."

This script does NOT run any new heavy analysis. It reads the per-k_p
ML-feature CSVs that the existing pipeline already produced across the
FULL 40-case sweep, and ranks every scalar property by how much it changes
across the k_p range. The output tells you which one or two properties are
worth analyzing in detail for the low / intermediate regimes, so you do
not have to repeat every analysis for every k_p.

Sensitivity metrics reported per property (each computed over the set of
per-k_p values, one value per case folder):
  - span_frac   : (max - min) / |mean|            (dimensionless dynamic range)
  - spearman_kp : Spearman rho vs k_p             (monotonicity, sign = direction)
  - rel_lo_hi   : (val@high_kp - val@low_kp)/|val@low_kp|  (end-to-end % change)
  - cv          : std / |mean|                    (overall scatter)

A property is flagged "high sensitivity" when it is BOTH large-swing
(span_frac in the top tertile) AND monotonic (|spearman_kp| >= 0.6). Those
are the properties for which the low/intermediate k_p analysis is most
informative.

Inputs (whichever exist in CFG['out_dir']):
  ce_capacity_ML_features.csv
  li_buried_ML_features.csv
  reaction_ML_features.csv
  anode_ML_features.csv
  activation_barrier_ML_features.csv
  sei_morphology_ML_features.csv

Outputs to CFG['out_dir']:
  kp_sensitivity_ranking.csv           full table, sorted by |span_frac|
  fig_kp_sensitivity_ranking.png       horizontal bar chart of top movers
  fig_kp_top_properties_vs_kp.png      the top-N properties, each vs k_p

Run:  python kp_sensitivity_ranking_cfg.py --config config_plating
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

import kmc_io as io

warnings.filterwarnings("ignore")

# ML-feature CSVs to scan, and a friendly label prefix so identical column
# names in different files do not collide.
FEATURE_FILES = [
    ("ce_capacity_ML_features.csv",        "CE/Q"),
    ("li_buried_ML_features.csv",          "BuriedLi"),
    ("reaction_ML_features.csv",           "Rxn"),
    ("anode_ML_features.csv",              "Anode"),
    ("activation_barrier_ML_features.csv", "Barrier"),
    ("sei_morphology_ML_features.csv",     "SEI"),
]

# Columns that are identifiers / not physical observables.
SKIP_COLS = {"case_key", "k0", "db_cycle"}


def load_merged(out_dir):
    """Merge every available per-k_p ML feature table on k0, prefixing
    each property with its source group so names stay unique."""
    merged = None
    for fname, group in FEATURE_FILES:
        p = out_dir / fname
        if not p.is_file():
            print(f"  (skip, not found: {fname})")
            continue
        df = pd.read_csv(p)
        if "k0" not in df.columns:
            print(f"  (skip, no k0 column: {fname})")
            continue
        keep = [c for c in df.columns if c not in SKIP_COLS]
        df = df[["k0"] + keep].copy()
        # collapse any duplicate k0 rows (shouldn't happen) by mean
        df = df.groupby("k0", as_index=False).mean(numeric_only=True)
        ren = {c: f"{group}:{c}" for c in keep}
        df = df.rename(columns=ren)
        merged = df if merged is None else merged.merge(df, on="k0", how="outer")
    return merged


def sensitivity_table(merged):
    kp = merged["k0"].astype(float).values
    order = np.argsort(kp)
    kp_sorted = kp[order]
    rows = []
    for col in merged.columns:
        if col == "k0":
            continue
        y = merged[col].astype(float).values[order]
        m = np.isfinite(y) & np.isfinite(kp_sorted)
        if m.sum() < 4:
            continue
        yv = y[m]; kv = kp_sorted[m]
        mean = np.mean(yv)
        denom = abs(mean) if abs(mean) > 1e-12 else np.nan
        span_frac = (np.max(yv) - np.min(yv)) / denom
        cv = np.std(yv) / denom
        try:
            rho, _ = spearmanr(kv, yv)
        except Exception:
            rho = np.nan
        lo, hi = yv[0], yv[-1]
        rel_lo_hi = (hi - lo) / abs(lo) if abs(lo) > 1e-12 else np.nan
        rows.append({
            "property": col,
            "n_points": int(m.sum()),
            "min": float(np.min(yv)), "max": float(np.max(yv)),
            "mean": float(mean),
            "span_frac": float(span_frac),
            "cv": float(cv),
            "spearman_kp": float(rho) if rho is not None else np.nan,
            "rel_lo_hi": float(rel_lo_hi),
        })
    tab = pd.DataFrame(rows)
    if tab.empty:
        return tab
    # High-sensitivity flag: big swing AND monotonic
    span_hi = tab["span_frac"].abs().quantile(2 / 3)
    tab["high_sensitivity"] = (
        (tab["span_frac"].abs() >= span_hi)
        & (tab["spearman_kp"].abs() >= 0.6)
    )
    return tab.sort_values("span_frac", key=lambda s: s.abs(),
                           ascending=False).reset_index(drop=True)


def fig_ranking(tab, cfg, out_path, top_n=18):
    plt.rcParams.update(cfg["plot_style"])
    d = tab.head(top_n).iloc[::-1]
    colors = ["#C0392B" if hs else "#5D6D7E"
              for hs in d["high_sensitivity"]]
    fig, ax = plt.subplots(figsize=(9.2, 0.42 * len(d) + 1.5))
    ax.barh(d["property"], d["span_frac"].abs(), color=colors)
    ax.set_xlabel("Dynamic range across the sweep  |max − min| / |mean|")
    sym = cfg.get("rate_symbol_math", "$k_p$")
    ax.set_title(f"Which properties respond most to {sym}?",
                 fontweight="bold")
    ax.grid(alpha=0.3, axis="x", lw=0.5)
    # legend proxy
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color="#C0392B", label="High sensitivity (big + monotonic)"),
        Patch(color="#5D6D7E", label="Other")],
        loc="lower right", fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close()


def fig_top_vs_kp(merged, tab, cfg, out_path, top_n=6):
    plt.rcParams.update(cfg["plot_style"])
    props = tab.head(top_n)["property"].tolist()
    n = len(props)
    if n == 0:
        return
    ncol = 2
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.0 * nrow))
    axes = np.atleast_1d(axes).ravel()
    kp = merged["k0"].astype(float).values
    order = np.argsort(kp)
    sym = cfg.get("rate_symbol_math", "$k_p$")
    kref = cfg.get("kp_reference", None)
    for ax, prop in zip(axes, props):
        y = merged[prop].astype(float).values
        ax.plot(kp[order], y[order], "o-", color="#2980B9",
                markersize=5, lw=1.6)
        if kref is not None:
            ax.axvline(kref, color="#27AE60", ls="--", lw=1.3, alpha=0.8)
            ax.text(kref, ax.get_ylim()[1], f" {sym}={kref:g}",
                    color="#27AE60", fontsize=8, va="top")
        ax.axvline(cfg["ref_high_k0"], color="#C0392B", ls=":", lw=1.3,
                   alpha=0.8)
        ax.set_xscale("log")
        ax.set_xlabel(sym)
        ax.set_ylabel(prop, fontsize=9)
        ax.grid(alpha=0.3, which="both", lw=0.4)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(f"Top responding properties vs {sym} "
                 "(green = reference, red = comparison)",
                 fontweight="bold", y=1.0)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]
    out.mkdir(parents=True, exist_ok=True)

    print(f"k_p sensitivity ranking — {cfg['sweep_name']}")
    merged = load_merged(out)
    if merged is None or merged.empty:
        print("No ML-feature CSVs found. Run the pipeline first "
              "(reaction_tracking, ce_capacity, li_buried, anode, "
              "porosity, activation_barrier).")
        return
    print(f"Merged {merged.shape[1]-1} properties over "
          f"{merged['k0'].nunique()} k_p values.")

    tab = sensitivity_table(merged)
    if tab.empty:
        print("Not enough finite data to rank."); return
    tab.to_csv(out / "kp_sensitivity_ranking.csv", index=False)
    print(f"Saved kp_sensitivity_ranking.csv ({len(tab)} properties)")

    print("\nTop responders (by dynamic range):")
    for _, r in tab.head(10).iterrows():
        flag = "  *** HIGH" if r["high_sensitivity"] else ""
        print(f"  {r['property']:<34s} span={r['span_frac']:+.2f}  "
              f"rho={r['spearman_kp']:+.2f}  loHi={r['rel_lo_hi']:+.1%}{flag}")

    fig_ranking(tab, cfg, out / "fig_kp_sensitivity_ranking.png")
    fig_top_vs_kp(merged, tab, cfg, out / "fig_kp_top_properties_vs_kp.png")
    print("\nSaved fig_kp_sensitivity_ranking.png and "
          "fig_kp_top_properties_vs_kp.png")
    print("\n>>> Use the top HIGH-sensitivity property as the target for "
          "kp_focused_regimes_cfg.py (Q2: below 56.25, Q3: intermediate).")


if __name__ == "__main__":
    main()
