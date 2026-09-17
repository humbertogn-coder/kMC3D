#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
activation_barrier_cfg.py
═════════════════════════
Computes the mean activation barrier (Ea) per kMC event, separated by
charge / discharge half-cycle, for every case in the config. Produces
the Ea-asymmetry analysis (charge/discharge ratio) vs k0.

Phase is determined by voltage: charge when voltage is the higher value
(typically 4.4), discharge when lower (typically 2.8). The split point is
the midpoint between the two most common voltages in the file.

Outputs to CFG['out_dir']:
  - activation_barrier_ML_features.csv   one row per k0
  - fig_activation_barrier_vs_k0.png     mean Ea charge & discharge vs k0
  - fig_activation_asymmetry_vs_k0.png   charge/discharge ratio vs k0

Run:  python activation_barrier_cfg.py --config config_plating
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import kmc_io as io

warnings.filterwarnings("ignore")


def split_voltage(df_log):
    """Return the voltage threshold separating charge from discharge."""
    vals = df_log["voltage"].value_counts().index.tolist()
    vals = sorted(vals)
    if len(vals) >= 2:
        return 0.5 * (vals[0] + vals[-1])
    return df_log["voltage"].median()


def process_case(case_key, k0, cfg):
    base = cfg["base_dir"] / case_key
    log  = io.resolve_log(base, cfg)
    if log is None:
        print(f"  [{case_key}] missing log; skip", flush=True)
        return None

    df_log = io.parse_data2excel(log)
    # Keep only events with a meaningful Ea (diffusion / surface events
    # carry the barrier; reaction events have Ea≈0 placeholder)
    df = df_log[(df_log["Ea"] > 0) & (df_log["cycle"] <= cfg["max_cycle"])]
    if df.empty:
        print(f"  [{case_key}] no Ea>0 events; skip", flush=True)
        return None

    thr = split_voltage(df_log)
    charge = df[df["voltage"] >= thr]
    disch  = df[df["voltage"] <  thr]

    ea_charge = float(charge["Ea"].mean()) if len(charge) else np.nan
    ea_disch  = float(disch["Ea"].mean())  if len(disch)  else np.nan
    ratio = ea_charge / ea_disch if (ea_disch and ea_disch > 0) else np.nan

    print(f"  [{case_key}] Ea_charge={ea_charge:.4f} "
          f"Ea_disch={ea_disch:.4f} ratio={ratio:.3f}", flush=True)
    return {
        "case_key": case_key, "k0": k0,
        "ea_charge_mean": ea_charge, "ea_disch_mean": ea_disch,
        "ea_asymmetry_ratio": ratio,
        "n_charge_events": int(len(charge)),
        "n_disch_events": int(len(disch)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]; out.mkdir(parents=True, exist_ok=True)

    print(f"Activation barrier — {cfg['sweep_name']}")
    rows = []
    for case_key, k0 in cfg["cases"]:
        print(f"\n── case {case_key} (k0={k0}) ──", flush=True)
        r = process_case(case_key, k0, cfg)
        if r is not None: rows.append(r)

    if not rows:
        print("No data."); return
    df_ml = pd.DataFrame(rows).sort_values("k0").reset_index(drop=True)
    df_ml.to_csv(out / "activation_barrier_ML_features.csv", index=False)
    print(f"\nSaved activation_barrier_ML_features.csv ({len(df_ml)} rows)")

    # ── Figure: Ea charge & discharge vs k0 ──
    plt.rcParams.update(cfg["plot_style"])
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    d = df_ml.sort_values("k0")
    ax.plot(d["k0"], d["ea_charge_mean"], "o-", color="#C0392B", lw=2,
            markersize=7, label="Charge")
    ax.plot(d["k0"], d["ea_disch_mean"], "s-", color="#2980B9", lw=2,
            markersize=7, label="Discharge")
    ax.set_xscale("log")
    ax.set_xlabel(cfg["rate_symbol_math"] + " (plating rate)")
    ax.set_ylabel("Mean $E_a$ per kMC event (eV)")
    ax.set_title(f"Activation barrier vs {cfg['rate_symbol_math']}", fontweight="bold")
    ax.grid(alpha=0.3, which="both", lw=0.5)
    ax.legend(loc="best")
    plt.tight_layout()
    fig.savefig(out / "fig_activation_barrier_vs_k0.png",
                bbox_inches="tight", facecolor="white")
    plt.close()

    # ── Figure: asymmetry ratio vs k0 ──
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.plot(d["k0"], d["ea_asymmetry_ratio"], "o-", color="#8E44AD", lw=2,
            markersize=7, label="Charge / discharge ratio")
    ax.axhline(1.0, ls="--", color="#444", lw=1.2, alpha=0.7,
               label=r"$E_{charge}=E_{discharge}$")
    ax.set_xscale("log")
    ax.set_xlabel(cfg["rate_symbol_math"] + " (plating rate)")
    ax.set_ylabel(r"$E_{charge}\,/\,E_{discharge}$")
    ax.set_title(f"Charge / discharge barrier asymmetry vs {cfg['rate_symbol_math']}",
                 fontweight="bold")
    ax.grid(alpha=0.3, which="both", lw=0.5)
    ax.legend(loc="best")
    plt.tight_layout()
    fig.savefig(out / "fig_activation_asymmetry_vs_k0.png",
                bbox_inches="tight", facecolor="white")
    plt.close()
    print("Saved figures.")


if __name__ == "__main__":
    main()
