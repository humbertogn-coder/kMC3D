#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_comparison_figures_cfg.py
══════════════════════════════
Final comparison figures for the meeting. Format matches the stripping
figures: plain title on top (no bold serif), sans-serif default font,
"k0 = X" legend labels, thick lines. Selection is by k0 value (robust
against case_key string/int round-trip).

Reproduces ALL reaction / consumption / morphology figures shown in the
presentation (no ML/SHAP/BO, no 3D):
  - CE and Capacity (high vs low k0)
  - SEI thickness vs cycle (high vs low) + vs k0
  - SEI porosity profile (high vs low)
  - Li buried dynamic + kinetic (vs cycle and vs k0)
  - Net Li balance vs k0 (Figure 4 style)
  - Activation barrier + asymmetry vs k0 (Figure 10)
  - Anode consumption: total, SFO, F5D, LiF (Figures 6-9 style)
  - Reaction products: O2-, F-, N3-(=a6), S2- per cycle (high vs low)
  - Reaction event rates: salt decomposition, solvent decomposition

Run:  python make_comparison_figures_cfg.py --config config_plating
"""
import argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import kmc_io as io
warnings.filterwarnings("ignore")

# Symbol used in all axis labels / legends (set from config in main()).
# PBB asked for k_p instead of k0 to avoid confusing the studies.
RATE_SYM = "k0"

# Stripping-style format: default sans-serif, plain title, no bold
STYLE = {
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "axes.linewidth": 1.0,
    "legend.frameon": True,
    "legend.framealpha": 0.9,
    "savefig.dpi": 200,
    "figure.facecolor": "white",
}
LW = 2.2


def safe_read(path):
    try: return pd.read_csv(path)
    except Exception: return None


def sel(df, k0):
    """Select rows for a k0 value robustly, sorted by cycle."""
    if df is None or df.empty or "k0" not in df.columns:
        return pd.DataFrame()
    return io.select_by_k0(df, k0)


def multi_case_fig(dfs_cols, cases, xlabel, ylabel, title, out_path,
                   ylim=None, marker=None):
    """Generic multi-case line figure. `cases` is the list of
    (key, k0, color, label) tuples from cfg["comparison_cases"], so the
    same figure shows lowest / reference / highest k_p together."""
    plt.rcParams.update(STYLE)
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    df, ycol = dfs_cols
    xcol = "cycle"
    for key, k0, color, label in cases:
        g = sel(df, k0)
        if g.empty or ycol not in g.columns: continue
        lw = LW + 0.6 if label == "reference" else LW
        ax.plot(g[xcol], g[ycol], "-" if marker is None else marker+"-",
                lw=lw, color=color, label=f"{RATE_SYM} = {k0:g} ({label})",
                markersize=4 if marker else 0)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
    if ylim: ax.set_ylim(*ylim)
    ax.grid(alpha=0.3); ax.legend(loc="best")
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight"); plt.close()
    print(f"  saved {out_path.name}")


def vs_k0_fig(df, ycol, ylabel, title, out_path, logx=True, fit=False):
    """Single quantity vs k0 across all cases (final/aggregated)."""
    if df is None or df.empty or ycol not in df.columns:
        return
    plt.rcParams.update(STYLE)
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    d = df.sort_values("k0")
    ax.plot(d["k0"], d[ycol], "o-", lw=LW, color="#1f5f8b", markersize=6)
    if logx: ax.set_xscale("log")
    ax.set_xlabel(f"{RATE_SYM} (plating rate)"); ax.set_ylabel(ylabel)
    ax.set_title(title); ax.grid(alpha=0.3, which="both", lw=0.5)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight"); plt.close()
    print(f"  saved {out_path.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]
    sym = cfg.get("rate_symbol_plain", "k0")   # PBB: use k_p in labels
    if not out.is_dir():
        print(f"Output dir not found: {out}"); return

    CASES = cfg["comparison_cases"]   # (key, k0, color, label) x3
    global RATE_SYM
    RATE_SYM = sym

    # Read CSVs
    libur   = safe_read(out / "li_buried_per_cycle.csv")
    libur_ml= safe_read(out / "li_buried_ML_features.csv")
    react   = safe_read(out / "reaction_tracking_per_cycle.csv")
    anode   = safe_read(out / "anode_consumption_per_cycle.csv")
    barr_ml = safe_read(out / "activation_barrier_ML_features.csv")
    cecap   = safe_read(out / "ce_capacity_per_cycle.csv")
    seig    = safe_read(out / "sei_growth_per_cycle.csv")
    profile = safe_read(out / "porosity_profiles.csv")

    print("Generating comparison figures (stripping format)...")

    # ── CE & Capacity ──────────────────────────────────────────────────
    multi_case_fig((cecap, "CE"), CASES,
                 "Cycle", "Coulombic Efficiency (%)",
                 "Comparative Coulombic Efficiency vs Cycle",
                 out / "cmp_CE_vs_cycle.png", ylim=(0, 110))
    multi_case_fig((cecap, "capacity_C"), CASES,
                 "Cycle", "Capacity (C)",
                 "Comparative Capacity vs Cycle",
                 out / "cmp_capacity_vs_cycle.png")

    # ── SEI thickness ──────────────────────────────────────────────────
    multi_case_fig((seig, "SEI_thickness_A"), CASES,
                 "Cycle", "SEI Thickness (Å)",
                 "Comparative SEI Thickness vs Cycle",
                 out / "cmp_SEI_thickness_vs_cycle.png")
    if seig is not None and not seig.empty:
        fin = (seig.sort_values("cycle").groupby("k0").last().reset_index())
        vs_k0_fig(fin, "SEI_thickness_A", "SEI thickness, final cycle (Å)",
                  f"SEI thickness vs {sym}", out / "cmp_SEI_thickness_vs_k0.png")

    # ── SEI porosity profile ───────────────────────────────────────────
    if profile is not None and not profile.empty:
        plt.rcParams.update(STYLE)
        fig, ax = plt.subplots(figsize=(8.8, 5.2))
        for key, k0, color, label in CASES:
            g = sel(profile, k0)
            if g.empty: continue
            cyc = g["cycle"].max()
            gc = g[g["cycle"] == cyc].sort_values("z_rel_A")
            # z_rel_A is already the distance from the anode surface
            ax.plot(gc["z_rel_A"], gc["porosity_corrected"], "-",
                    lw=LW + 0.6 if label == "reference" else LW,
                    color=color, label=f"{RATE_SYM} = {k0:g} ({label})")
        ax.set_xlabel("Distance from anode surface (Å)")
        ax.set_ylabel("Porosity"); ax.set_title("SEI porosity")
        ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend(loc="best")
        plt.tight_layout()
        fig.savefig(out / "cmp_SEI_porosity_profile.png", bbox_inches="tight")
        plt.close(); print("  saved cmp_SEI_porosity_profile.png")

    # ── Li buried (dynamic + kinetic) ──────────────────────────────────
    multi_case_fig((libur, "n_li_buried_dynamic"), CASES,
                 "Cycle", "Dynamic Li buried (voxels)",
                 "Dynamic Li Buried (absolute) vs Cycle",
                 out / "cmp_Li_buried_dynamic_vs_cycle.png", marker="o")
    multi_case_fig((libur, "frac_buried_dynamic"), CASES,
                 "Cycle", "Buried Li fraction",
                 "Dynamic Li Buried (fraction) vs Cycle",
                 out / "cmp_Li_buried_fraction_vs_cycle.png", marker="o")
    multi_case_fig((libur, "li_buried_kinetic"), CASES,
                 "Cycle", "Kinetic Li buried",
                 "Kinetic Li Buried vs Cycle",
                 out / "cmp_Li_buried_kinetic_vs_cycle.png", marker="o")
    if libur_ml is not None and not libur_ml.empty:
        vs_k0_fig(libur_ml, "n_li_buried_dynamic_at_db_cycle",
                  f"Dynamic Li buried, cycle {cfg['li_buried_cycle_for_db']}",
                  "Dynamic Li Buried vs k0 (Figure 5 style)",
                  out / "cmp_Li_buried_dynamic_vs_k0.png")
        # Net Li balance (Figure 4 style): kinetic mean over first cycles
        vs_k0_fig(libur_ml, "li_buried_kinetic_mean",
                  "Mean net Li per cycle",
                  "Net Li balance vs k0 (Figure 4 style)",
                  out / "cmp_net_Li_balance_vs_k0.png")

    # ── Activation barriers (Figure 10) ────────────────────────────────
    if barr_ml is not None and not barr_ml.empty:
        plt.rcParams.update(STYLE)
        d = barr_ml.sort_values("k0")
        fig, ax = plt.subplots(figsize=(8.8, 5.2))
        ax.plot(d["k0"], d["ea_charge_mean"], "o-", lw=LW, color="#d62728",
                markersize=6, label="Charge")
        ax.plot(d["k0"], d["ea_disch_mean"], "s-", lw=LW, color="#1f77b4",
                markersize=6, label="Discharge")
        ax.set_xscale("log"); ax.set_xlabel(f"{sym} (plating rate)")
        ax.set_ylabel("Mean Ea per event (eV)")
        ax.set_title(f"Activation barrier vs {sym}")
        ax.grid(alpha=0.3, which="both", lw=0.5); ax.legend(loc="best")
        plt.tight_layout()
        fig.savefig(out / "cmp_activation_barrier_vs_k0.png", bbox_inches="tight")
        plt.close(); print("  saved cmp_activation_barrier_vs_k0.png")
        vs_k0_fig(d, "ea_asymmetry_ratio", "Ea charge / Ea discharge",
                  f"Charge/discharge asymmetry vs {sym}",
                  out / "cmp_activation_asymmetry_vs_k0.png")

    # ── Anode consumption (Figures 6-9) ────────────────────────────────
    for ycol, ylabel, title, fname in [
        ("total_sei_anode", "Total SEI atoms in anode",
         "Anode consumption: total SEI (Figure 6)",
         "cmp_anode_total.png"),
        ("SFO_anode", "SFO species in anode",
         "SFO anode consumption (Figure 7)", "cmp_anode_SFO.png"),
        ("F5D_anode", "F5D species in anode",
         "F5D anode consumption (Figure 8)", "cmp_anode_F5D.png"),
        ("F_anode", "F (LiF precursor) in anode",
         "LiF anode consumption (Figure 9)", "cmp_anode_LiF.png"),
    ]:
        multi_case_fig((anode, ycol), CASES,
                     "Cycle", ylabel, title, out / fname, marker="o")

    # ── Reaction products per cycle (high vs low) ──────────────────────
    for ycol, ylabel, title, fname in [
        ("O_xyz", "O²⁻ count", "O²⁻ accumulation vs cycle",
         "cmp_rxn_O2.png"),
        ("F_xyz", "F⁻ count", "F⁻ accumulation vs cycle",
         "cmp_rxn_F.png"),
        ("N_xyz", "N³⁻ count (= a-6 events)",
         "N³⁻ accumulation vs cycle (Li₃N precursor)", "cmp_rxn_N.png"),
        ("S_xyz", "S²⁻ count", "S²⁻ accumulation vs cycle (Li₂S precursor)",
         "cmp_rxn_S.png"),
        ("SFO_xyz", "SFO count", "SFO accumulation vs cycle",
         "cmp_rxn_SFO.png"),
        ("F5D_xyz", "F5D count", "F5D accumulation vs cycle",
         "cmp_rxn_F5D.png"),
    ]:
        multi_case_fig((react, ycol), CASES,
                     "Cycle", ylabel, title, out / fname, marker="o")

    # ── Reaction event totals (cumulative paths) ───────────────────────
    for ycol, ylabel, title, fname in [
        ("RxnSFO", "Cumulative salt-decomp events (a-2…a-8)",
         "Salt decomposition events vs cycle", "cmp_rxn_RxnSFO.png"),
        ("RxnFSI", "Cumulative FSI reduction (a-1)",
         "FSI reduction events vs cycle", "cmp_rxn_RxnFSI.png"),
        ("RxnF5D", "Cumulative solvent-decomp events (b-2,b-3)",
         "Solvent decomposition events vs cycle", "cmp_rxn_RxnF5D.png"),
        ("RxnSOL", "Cumulative F5DEE reduction (b-1)",
         "F5DEE reduction events vs cycle", "cmp_rxn_RxnSOL.png"),
        ("RxnSOL2", "Cumulative SOL2 events (c-1)",
         "SOL2 events vs cycle", "cmp_rxn_RxnSOL2.png"),
        ("RxnPlating", "Cumulative plating events",
         "Plating events vs cycle", "cmp_rxn_RxnPlating.png"),
        ("RxnStripping", "Cumulative stripping events",
         "Stripping events vs cycle", "cmp_rxn_RxnStripping.png"),
        ("RxnLiSurface", "Cumulative Li surface diffusion events",
         "Li surface diffusion events vs cycle", "cmp_rxn_RxnLiSurface.png"),
        ("RxnPlatingSEI", "Cumulative plating-through-SEI events",
         "Plating through SEI events vs cycle", "cmp_rxn_RxnPlatingSEI.png"),
    ]:
        multi_case_fig((react, ycol), CASES,
                     "Cycle", ylabel, title, out / fname, marker="o")

    print(f"\nAll comparison figures in {out}")


if __name__ == "__main__":
    main()
