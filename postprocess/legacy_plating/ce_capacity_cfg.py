#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ce_capacity_cfg.py
══════════════════
Per-cycle Coulombic Efficiency (CE) and Capacity (Q) for the sweep.

  CE  = 100 × re / ce  per cycle
  Q   = re_events × e   (per cycle, in Coulombs)

Outputs:
  ce_capacity_per_cycle.csv
  ce_capacity_ML_features.csv
  fig_CE_high_vs_low_k0.png
  fig_capacity_high_vs_low_k0.png
"""
import argparse, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import kmc_io as io
warnings.filterwarnings("ignore")

E_CHARGE = 1.602e-19


def per_cycle_ce_capacity(df_log):
    end_half = df_log.groupby("halfcycle").last().reset_index()
    rows = []
    for cyc, gcyc in end_half.groupby("cycle"):
        if gcyc.empty: continue
        cum_end = gcyc.iloc[-1]
        prev_mask = end_half["cycle"] == (cyc - 1)
        if prev_mask.any():
            cum_prev = end_half[prev_mask].iloc[-1]
        else:
            cum_prev = pd.Series({c: 0 for c in
                ["RxnPlating","RxnPlatingSEI","RxnStripping","RxnFSI",
                 "RxnSFO","RxnSOL","RxnF5D","RxnSOL2"]})
        ce_events = sum(cum_end[c] - cum_prev[c] for c in
            ["RxnPlating","RxnPlatingSEI","RxnFSI","RxnSFO",
             "RxnSOL","RxnF5D","RxnSOL2"])
        re_events = cum_end["RxnStripping"] - cum_prev["RxnStripping"]
        ce_pct = 100.0 * re_events / ce_events if ce_events > 0 else np.nan
        capacity_C = re_events * E_CHARGE if re_events > 0 else np.nan
        rows.append({
            "cycle": int(cyc),
            "ce_events": int(ce_events), "re_events": int(re_events),
            "CE": float(ce_pct) if pd.notna(ce_pct) else np.nan,
            "capacity_C": float(capacity_C) if pd.notna(capacity_C) else np.nan,
        })
    return pd.DataFrame(rows)


def process_case(case_key, k0, cfg):
    log = io.resolve_log(cfg["base_dir"] / case_key, cfg)
    if log is None:
        print(f"  [{case_key}] missing; skip", flush=True); return None
    df_log = io.parse_data2excel(log, full=True)
    per_cyc = per_cycle_ce_capacity(df_log)
    per_cyc.insert(0, "case_key", case_key); per_cyc.insert(1, "k0", k0)
    per_cyc = per_cyc[per_cyc["cycle"] <= cfg["max_cycle"]]
    print(f"  [{case_key}] {len(per_cyc)} cyc, CE_mean={per_cyc['CE'].mean():.1f}%",
          flush=True)
    return per_cyc


def fig_pair(df_all, ycol, ylabel, title, out_path, cfg, ylim=None):
    """Lowest / reference / highest k_p on the same axes, read from
    cfg["comparison_cases"] = [(key, k0, color, label), ...]."""
    plt.rcParams.update(cfg["plot_style"])
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    sym = cfg["rate_symbol_plain"]
    for key, k0, color, label in cfg["comparison_cases"]:
        g = df_all[df_all["case_key"] == key].sort_values("cycle")
        if g.empty: continue
        ax.plot(g["cycle"], g[ycol], "o-", lw=2.6 if label == "reference" else 2,
                color=color, markersize=4, alpha=0.9,
                label=f"{sym} = {k0:g} ({label})")
    ax.set_xlabel("Cycle number"); ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold"); ax.grid(alpha=0.3)
    ax.legend(loc="best")
    if ylim: ax.set_ylim(*ylim)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor="white"); plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]; out.mkdir(parents=True, exist_ok=True)

    print(f"CE & Capacity — {cfg['sweep_name']}")
    parts = []
    for case_key, k0 in cfg["cases"]:
        print(f"\n── case {case_key} (k0={k0}) ──", flush=True)
        df = process_case(case_key, k0, cfg)
        if df is not None: parts.append(df)
    if not parts: print("No data."); return

    df_all = pd.concat(parts, ignore_index=True)
    df_all.to_csv(out / "ce_capacity_per_cycle.csv", index=False)
    ml = (df_all.groupby(["case_key","k0"])
          .agg(CE_mean=("CE","mean"), CE_final=("CE","last"),
               capacity_C_mean=("capacity_C","mean"),
               capacity_C_final=("capacity_C","last"))
          .reset_index())
    ml.to_csv(out / "ce_capacity_ML_features.csv", index=False)
    print(f"\nSaved CSVs ({len(df_all)} per-cycle rows, {len(ml)} per-k0 rows)")

    fig_pair(df_all, "CE", "Coulombic efficiency (%)",
             f"Coulombic efficiency — lowest / reference / highest {cfg['rate_symbol_plain']}",
             out / "fig_CE_high_vs_low_k0.png", cfg, ylim=(0, 110))
    fig_pair(df_all, "capacity_C", "Capacity Q (C)",
             f"Reversible capacity — lowest / reference / highest {cfg['rate_symbol_plain']}",
             out / "fig_capacity_high_vs_low_k0.png", cfg)
    print("Saved figures.")


if __name__ == "__main__":
    main()
