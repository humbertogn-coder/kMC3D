#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reaction_tracking_cfg.py
════════════════════════
Reaction-path tracking + mass-balance reconstruction for every case in
the config. Combines event counts from data2excel.txt with species
counts in the .xyz frames to infer per-path event counts.

EXACT inferences (conservation only):
    events_a6 = N_xyz                 (N only from a-6)
    events_a8 = S_xyz - N_xyz         (S from a-6 + a-8)
    events_a4+a5+a7 = O_xyz - events_a8
    events_a2+a3 = RxnSFO - events_a6 - events_a8 - events_a4+a5+a7

Outputs to CFG['out_dir']:
  - reaction_tracking_per_cycle.csv   one row per (k0, cycle)
  - reaction_ML_features.csv          one row per k0 (final cumulative)

Run:  python reaction_tracking_cfg.py --config config_plating
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import kmc_io as io

warnings.filterwarnings("ignore")


def process_case(case_key, k0, cfg):
    base = cfg["base_dir"] / case_key
    traj = base / cfg["traj_subpath"]
    log  = io.resolve_log(base, cfg)
    if not traj.is_dir() or log is None:
        print(f"  [{case_key}] missing; skip", flush=True)
        return None

    df_log = io.parse_data2excel(log, full=True)
    step_to_cycle = io.build_step_to_cycle(df_log)
    cyc2xyz = io.map_cycles_to_xyz(io.list_xyz_files(traj), step_to_cycle)
    per_cyc = io.per_cycle_last(
        df_log, ["RxnPlating", "RxnStripping", "RxnLiSurface",
                 "RxnPlatingSEI", "RxnFSI", "RxnSFO", "RxnSOL",
                 "RxnF5D", "RxnSOL2"])

    cycles = list(range(cfg["cycle_step_reactions"],
                        cfg["max_cycle"] + 1, cfg["cycle_step_reactions"]))
    rows = []
    for cyc in cycles:
        if cyc not in cyc2xyz: continue
        prow = per_cyc[per_cyc["cycle"] == cyc]
        if prow.empty: continue
        prow = prow.iloc[0]
        counts = io.count_species_last_frame(cyc2xyz[cyc][1],
                                              cfg["species_xyz"])
        N, S, O, F = counts["N"], counts["S"], counts["O"], counts["F"]
        RxnSFO = int(prow["RxnSFO"])
        a6 = N
        a8 = max(S - N, 0)
        a4_5_7 = max(O - a8, 0)
        a2_3 = max(RxnSFO - a6 - a8 - a4_5_7, 0)
        rows.append({
            "case_key": case_key, "k0": k0, "cycle": cyc,
            "F_xyz": F, "O_xyz": O, "N_xyz": N, "S_xyz": S,
            "F5D_xyz": counts["F5D"], "SFO_xyz": counts["SFO"],
            "RxnPlating": int(prow["RxnPlating"]),
            "RxnStripping": int(prow["RxnStripping"]),
            "RxnLiSurface": int(prow["RxnLiSurface"]),
            "RxnPlatingSEI": int(prow["RxnPlatingSEI"]),
            "RxnFSI": int(prow["RxnFSI"]), "RxnSFO": RxnSFO,
            "RxnSOL": int(prow["RxnSOL"]), "RxnF5D": int(prow["RxnF5D"]),
            "RxnSOL2": int(prow["RxnSOL2"]),
            "events_a6": a6, "events_a8": a8,
            "events_a4_a5_a7": a4_5_7, "events_a2_a3": a2_3,
        })
        print(f"  [{case_key}] cyc {cyc:3d}: N={N} S={S} O={O} "
              f"a6={a6} a8={a8}", flush=True)

    if not rows: return None
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]; out.mkdir(parents=True, exist_ok=True)

    print(f"Reaction tracking — {cfg['sweep_name']}")
    parts = []
    for case_key, k0 in cfg["cases"]:
        print(f"\n── case {case_key} (k0={k0}) ──", flush=True)
        df = process_case(case_key, k0, cfg)
        if df is not None: parts.append(df)

    if not parts:
        print("No data."); return
    df_all = pd.concat(parts, ignore_index=True)
    df_all.to_csv(out / "reaction_tracking_per_cycle.csv", index=False)

    # ML features: final cumulative per k0
    ml = (df_all.sort_values("cycle").groupby(["case_key", "k0"]).last()
          .reset_index()[["case_key", "k0", "RxnPlating", "RxnStripping",
                           "RxnLiSurface", "RxnPlatingSEI",
                           "RxnFSI", "RxnSFO", "RxnSOL",
                           "RxnF5D", "RxnSOL2", "events_a6", "events_a8",
                           "events_a4_a5_a7", "events_a2_a3",
                           "F_xyz", "O_xyz", "N_xyz", "S_xyz",
                           "F5D_xyz", "SFO_xyz"]])
    ml.to_csv(out / "reaction_ML_features.csv", index=False)
    print(f"\nSaved reaction_tracking_per_cycle.csv ({len(df_all)} rows)")
    print(f"Saved reaction_ML_features.csv ({len(ml)} rows)")


if __name__ == "__main__":
    main()
