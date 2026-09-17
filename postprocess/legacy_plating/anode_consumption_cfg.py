#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
anode_consumption_cfg.py
════════════════════════
Counts SEI species (F, O, N, S, F5D, SFO) inside the original anode
region as a function of cycle, for every case in the config.

The anode region is defined from the calibration frame (kmc-coords-0):
  anode spans [z_center - half, z_center + half],
  where z_center = median(Li z) in calib and half = anode_z_width/2.

Outputs to CFG['out_dir']:
  - anode_consumption_per_cycle.csv   one row per (k0, cycle)
  - anode_ML_features.csv             one row per k0 (final + max)

Run:  python anode_consumption_cfg.py --config config_plating
"""

import argparse
import warnings
import numpy as np
import pandas as pd

import kmc_io as io

warnings.filterwarnings("ignore")


def anode_band_from_calib(calib_path, half):
    fr = io.parse_xyz_last_frame(calib_path)
    co = io.minimal_image(fr["coords"], fr["cell"])
    sp = fr["species"]
    z = co[sp == "Li", 2]
    if z.size == 0:
        raise RuntimeError(f"No Li in calib {calib_path}")
    c = float(np.median(z))
    return c - half, c + half


def count_in_band(path, z_lo, z_hi, species_list):
    fr = io.parse_xyz_last_frame(path)
    if fr is None:
        return {s: 0 for s in species_list}
    co = io.minimal_image(fr["coords"], fr["cell"])
    sp = fr["species"]
    mask = (co[:, 2] >= z_lo) & (co[:, 2] <= z_hi)
    spin = sp[mask]
    return {s: int(np.sum(spin == s)) for s in species_list}


def process_case(case_key, k0, cfg):
    base = cfg["base_dir"] / case_key
    traj = base / cfg["traj_subpath"]
    log  = io.resolve_log(base, cfg)
    half = cfg["anode_z_width"] / 2.0
    calib = traj / "kmc-coords-0.xyz"
    if not calib.is_file() or log is None:
        print(f"  [{case_key}] missing; skip", flush=True)
        return None

    z_lo, z_hi = anode_band_from_calib(calib, half)
    df_log = io.parse_data2excel(log)
    step_to_cycle = io.build_step_to_cycle(df_log)
    cyc2xyz = io.map_cycles_to_xyz(io.list_xyz_files(traj), step_to_cycle)

    cycles = list(range(cfg["cycle_step_reactions"],
                        cfg["max_cycle"] + 1, cfg["cycle_step_reactions"]))
    rows = []
    for cyc in cycles:
        if cyc not in cyc2xyz: continue
        counts = count_in_band(cyc2xyz[cyc][1], z_lo, z_hi, cfg["species_xyz"])
        inorg = sum(counts[s] for s in cfg["inorg_species"])
        amorph = sum(counts[s] for s in cfg["amorph_species"])
        total = inorg + amorph
        rows.append({
            "case_key": case_key, "k0": k0, "cycle": cyc,
            **{f"{s}_anode": counts[s] for s in cfg["species_xyz"] if s != "Li"},
            "inorg_anode": inorg, "amorph_anode": amorph,
            "total_sei_anode": total,
        })
    if not rows: return None
    df = pd.DataFrame(rows)
    print(f"  [{case_key}] final total SEI in anode = "
          f"{df['total_sei_anode'].iloc[-1]}", flush=True)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]; out.mkdir(parents=True, exist_ok=True)

    print(f"Anode consumption — {cfg['sweep_name']}")
    parts = []
    for case_key, k0 in cfg["cases"]:
        print(f"\n── case {case_key} (k0={k0}) ──", flush=True)
        df = process_case(case_key, k0, cfg)
        if df is not None: parts.append(df)

    if not parts:
        print("No data."); return
    df_all = pd.concat(parts, ignore_index=True)
    df_all.to_csv(out / "anode_consumption_per_cycle.csv", index=False)

    ml = (df_all.sort_values("cycle").groupby(["case_key", "k0"])
          .agg(total_sei_anode_final=("total_sei_anode", "last"),
               total_sei_anode_max=("total_sei_anode", "max"),
               inorg_anode_final=("inorg_anode", "last"),
               amorph_anode_final=("amorph_anode", "last"),
               SFO_anode_final=("SFO_anode", "last"),
               F5D_anode_final=("F5D_anode", "last"),
               F_anode_final=("F_anode", "last"))
          .reset_index())
    ml.to_csv(out / "anode_ML_features.csv", index=False)
    print(f"\nSaved anode_consumption_per_cycle.csv ({len(df_all)} rows)")
    print(f"Saved anode_ML_features.csv ({len(ml)} rows)")


if __name__ == "__main__":
    main()
