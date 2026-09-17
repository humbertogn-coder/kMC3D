#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kmc_io.py
═════════
Shared I/O helpers for the kMC analysis pipeline. Imported by every
analysis script so parsers / cycle-mapping live in ONE place.

Data2Excel.txt column layout (verified from calculate_reactions.cpp):
  col 0:  seed
  col 1:  seed counter
  col 2:  step_kmc            (= xyz file index)
  col 3:  half_cycle          (cycle = (half_cycle+1)//2)
  col 4:  event_idx
  col 5:  timeCurrentV
  col 6:  currentTime
  col 7:  voltage             (4.4 = charge, 2.8 = discharge)
  col 8:  reaction type
  col 9:  reactant species
  col 10: Ea
  col 11: expValue
  col 12: reactRate
  col 13: time
  col 14: RxnPlating          (cumulative)
  col 15: RxnStripping
  col 16: RxnLiSurface
  col 17: RxnFSI
  col 18: RxnSFO
  col 19: RxnSOL
  col 20: RxnF5D
  col 21: RxnSOL2
  col 22: RxnPlatingSEI
  col 23: LiMetal
  col 24: LiIon
  col 25: LiMetalSEI
  col 26: LiIonSEI
  col 27: nO
  col 28: nF
"""

import importlib
import re
from pathlib import Path

import numpy as np
import pandas as pd

# ── Column indices ──────────────────────────────────────────────────────
C_STEP, C_HALFCYCLE, C_VOLTAGE = 2, 3, 7
C_RXN_TYPE, C_REACTANT, C_EA   = 8, 9, 10
C_RxnPlating, C_RxnStripping   = 14, 15
C_RxnLiSurface, C_RxnFSI       = 16, 17
C_RxnSFO, C_RxnSOL             = 18, 19
C_RxnF5D, C_RxnSOL2            = 20, 21
C_RxnPlatingSEI                = 22
C_LiMetal, C_LiIon             = 23, 24
C_LiMetalSEI, C_LiIonSEI       = 25, 26


def load_config(config_module):
    """Import a config_*.py module by name and return its CFG dict."""
    mod = importlib.import_module(config_module)
    return mod.CFG


def resolve_log(case_dir, cfg):
    """
    Return the Path of the Data2Excel log inside `case_dir`, trying all
    filename capitalizations in cfg['log_name_candidates'] (falls back to
    cfg['log_name'], then a case-insensitive glob). Returns None if no
    log file exists. This makes the pipeline robust to the uppercase /
    lowercase mismatch between the C++ default ("Data2Excel.txt") and
    renamed plating outputs ("data2excel.txt").
    """
    case_dir = Path(case_dir)
    candidates = cfg.get("log_name_candidates") or [cfg["log_name"]]
    for name in candidates:
        p = case_dir / name
        if p.is_file():
            return p
    # last resort: case-insensitive search
    for p in case_dir.iterdir() if case_dir.is_dir() else []:
        if p.is_file() and p.name.lower() == "data2excel.txt":
            return p
    return None


# ════════════════════════════════════════════════════════════════════════
# Data2Excel parsing
# ════════════════════════════════════════════════════════════════════════

def parse_data2excel(path, full=False):
    """
    Parse Data2Excel / data2excel.txt.
    If full=False (default) keeps a compact set of columns useful for most
    analyses. If full=True keeps all the cumulative Rxn* and Li* counters.
    Returns a DataFrame with at least: step, halfcycle, cycle, voltage,
    rxn_type, reactant, Ea.
    """
    # printStepfscreen2() in calculate_reactions.cpp writes 29 fields per
    # row. full=True needs up to col 26 (LiIonSEI); compact mode needs up
    # to col 22 (RxnPlatingSEI check). Rows shorter than the requirement
    # (truncated last line after a killed job) are skipped explicitly.
    min_cols = 27 if full else 23
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < min_cols:
                continue
            try:
                rec = {
                    "step":      int(parts[C_STEP]),
                    "halfcycle": int(parts[C_HALFCYCLE]),
                    "voltage":   float(parts[C_VOLTAGE]),
                    "rxn_type":  parts[C_RXN_TYPE],
                    "reactant":  parts[C_REACTANT],
                    "Ea":        float(parts[C_EA]),
                }
                if full:
                    rec.update({
                        "RxnPlating":    int(parts[C_RxnPlating]),
                        "RxnStripping":  int(parts[C_RxnStripping]),
                        "RxnLiSurface":  int(parts[C_RxnLiSurface]),
                        "RxnFSI":        int(parts[C_RxnFSI]),
                        "RxnSFO":        int(parts[C_RxnSFO]),
                        "RxnSOL":        int(parts[C_RxnSOL]),
                        "RxnF5D":        int(parts[C_RxnF5D]),
                        "RxnSOL2":       int(parts[C_RxnSOL2]),
                        "RxnPlatingSEI": int(parts[C_RxnPlatingSEI]),
                        "LiMetal":       int(parts[C_LiMetal]),
                        "LiIon":         int(parts[C_LiIon]),
                        "LiMetalSEI":    int(parts[C_LiMetalSEI]),
                        "LiIonSEI":      int(parts[C_LiIonSEI]),
                    })
                rows.append(rec)
            except (ValueError, IndexError):
                continue
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError(f"No valid rows parsed from {path}")
    df["cycle"] = (df["halfcycle"] + 1) // 2
    return df


def build_step_to_cycle(df_log):
    """{step_kmc: cycle} for mapping xyz file indices to electrochemical
    cycle numbers."""
    return dict(zip(df_log["step"], df_log["cycle"]))


# ════════════════════════════════════════════════════════════════════════
# XYZ parsing
# ════════════════════════════════════════════════════════════════════════

def parse_xyz_last_frame(path):
    """Read all frames of an extxyz file, return the LAST one as dict
    with 'species' (np array), 'coords' (Nx3), 'cell' (3x3)."""
    last = None
    with open(path, "r") as fh:
        while True:
            n_line = fh.readline()
            if not n_line:
                break
            n_line = n_line.strip()
            if not n_line:
                continue
            try:
                n = int(n_line)
            except ValueError:
                continue
            comment = fh.readline()
            m = re.search(r'Lattice="([^"]+)"', comment)
            cell = (np.array(list(map(float, m.group(1).split()))).reshape(3, 3)
                    if m else np.eye(3))
            sp, xyz = [], []
            for _ in range(n):
                p = fh.readline().split()
                if len(p) >= 4:
                    sp.append(p[0])
                    xyz.append([float(p[1]), float(p[2]), float(p[3])])
            last = {"species": np.array(sp),
                    "coords":  np.array(xyz, dtype=float),
                    "cell":    cell}
    return last


def count_species_last_frame(path, species_list):
    """Return {species: count} for the LAST frame of an xyz file."""
    counts = {s: 0 for s in species_list}
    fr = parse_xyz_last_frame(path)
    if fr is None:
        return counts
    for s in fr["species"]:
        if s in counts:
            counts[s] += 1
    return counts


def minimal_image(coords, cell):
    frac = coords @ np.linalg.inv(cell.T)
    frac -= np.floor(frac)
    return frac @ cell


def xyz_index(p):
    m = re.search(r"kmc-coords-(\d+)\.xyz$", p.name)
    return int(m.group(1)) if m else None


def list_xyz_files(folder):
    """Return [(step_index, Path), ...] sorted by step index."""
    files = [(xyz_index(p), p) for p in folder.glob("kmc-coords-*.xyz")
             if xyz_index(p) is not None]
    files.sort(key=lambda x: x[0])
    return files


def map_cycles_to_xyz(xyz_pairs, step_to_cycle):
    """
    Assign each xyz file to its electrochemical cycle, keeping the LAST
    (highest-step) xyz per cycle. Returns {cycle: (step, Path)}.
    """
    cycle_to_xyz = {}
    for step, p in xyz_pairs:
        if step == 0:
            continue
        c = step_to_cycle.get(step)
        if c is None:
            candidates = [s for s in step_to_cycle if s <= step]
            if not candidates:
                continue
            c = step_to_cycle[max(candidates)]
        prev = cycle_to_xyz.get(c)
        if prev is None or step > prev[0]:
            cycle_to_xyz[c] = (step, p)
    return cycle_to_xyz


# ════════════════════════════════════════════════════════════════════════
# Per-cycle aggregation of cumulative counters
# ════════════════════════════════════════════════════════════════════════

def per_cycle_last(df_log_full, cols):
    """For each cycle, take the LAST cumulative value of each column in
    `cols`. Returns DataFrame indexed by cycle (reset)."""
    last = df_log_full.groupby("cycle")[cols].last()
    return last.reset_index()


def select_by_k0(df, k0_value, tol=1e-6):
    """Robustly select rows for a given k0 value (avoids string/int
    case_key mismatches after CSV round-trip). Returns sorted-by-cycle."""
    sub = df[np.isclose(df["k0"].astype(float), float(k0_value), atol=tol)]
    if "cycle" in sub.columns:
        sub = sub.sort_values("cycle")
    return sub
