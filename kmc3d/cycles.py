"""
cycles.py
=========
Per-cycle electrochemical metrics from cycle_stats.csv (one file per seed),
with conservation checks and multi-seed statistics. This is the first stage
of the post-processing pipeline: everything downstream (figures, parameter
sweeps, ML / SHAP) consumes the tables produced here.

Definitions (half-cycle i is a charge at BeginV for even i, a discharge at
EndV for odd i; a cell built from S8 starts in the charged state):

    cycle n            = discharge half 2n+1 followed by charge half 2n+2
    Q_dis[n]           = Li+ consumed at the cathode in half 2n+1 (d_li_consumed)
    Q_ch[n]            = Li+ released by the cathode in half 2n+2 (d_li_released)
    CE_cathode[n]      = Q_dis[n] / Q_ch[n]        (charge out / charge in of
                         the cathode; the cathode-side coulombic efficiency)
    capacity_mAh_g[n]  = Q_dis[n] / (16 * n_S8_initial) * 1672
                         (one electron per Li+; 16 e- per S8 = 1672 mAh/g_S)
    utilization[n]     = Q_dis[n] / (16 * n_S8_initial)

The engine's own CE columns (CE_cycle, CE_mod, CE_shuttle) are passed through.
Conservation: li_total, s_total must be constant; deposit_unplaced must be 0.

Cathode composition is reported with the S8-unit labels mapped to
conventional formulas (species.S8_UNIT_MAP): cat_Li2S4 = 2 * cat_Li4S8, etc.
"""

from __future__ import annotations
import csv
import glob
import os
from typing import Dict, List, Optional

import numpy as np

from .species import S8_UNIT_MAP, li_s_per_site

THEORETICAL_MAH_G = 1672.0     # S -> Li2S, 16 e- per S8
E_PER_S8 = 16


def load_cycle_stats(path: str) -> Dict[str, np.ndarray]:
    """cycle_stats.csv -> {column: float array}; 'phase' kept as strings."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return {}
    out: Dict[str, np.ndarray] = {}
    for k in rows[0]:
        if k == "phase":
            out[k] = np.array([r[k] for r in rows])
            continue
        vals = []
        for r in rows:
            v = r.get(k, "")
            try:
                vals.append(float(v) if v != "" else np.nan)
            except ValueError:
                vals.append(np.nan)
        out[k] = np.array(vals)
    return out


def conservation_report(t: Dict[str, np.ndarray]) -> Dict[str, object]:
    """Constant li_total / s_total, zero deposit_unplaced. Returns a dict with
    booleans and the offending ranges; never raises."""
    rep: Dict[str, object] = {}
    for key in ("li_total", "s_total"):
        if key in t:
            col = t[key][~np.isnan(t[key])]
            rep[f"{key}_constant"] = bool(col.size and np.all(col == col[0]))
            rep[f"{key}_range"] = (float(col.min()), float(col.max())) if col.size else None
    if "deposit_unplaced" in t:
        col = t["deposit_unplaced"][~np.isnan(t["deposit_unplaced"])]
        rep["deposit_unplaced_zero"] = bool(col.size == 0 or col.max() == 0)
    rep["ok"] = all(v for k, v in rep.items() if k.endswith("_constant") or k.endswith("_zero"))
    return rep


def initial_s8_sites(t: Dict[str, np.ndarray]) -> int:
    """S8 sites at the first half-cycle close (cathode not yet cycled)."""
    if "cat_S8" in t and not np.isnan(t["cat_S8"][0]):
        return int(t["cat_S8"][0])
    if "n_cathode" in t:
        return int(t["n_cathode"][0])
    return 0


def per_cycle(t: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Cycle-resolved metrics for one seed (see module docstring)."""
    n_half = t["half_index"].size
    n_s8 = initial_s8_sites(t)
    q_theo = E_PER_S8 * n_s8
    dis = t.get("d_li_consumed", np.zeros(n_half))
    rel = t.get("d_li_released", np.zeros(n_half))
    cyc, qd, qc, ce, cap, util = [], [], [], [], [], []
    passthrough = {k: [] for k in ("CE_cycle", "CE_mod", "CE_shuttle") if k in t}
    extra = {k: [] for k in t if k.startswith(("cat_", "n_dis_", "sei_", "cei_"))
             or k in ("n_Li", "n_SEI", "n_cathode", "n_CEI", "li_pool", "li_bulk",
                      "n_ETH", "n_cathode_blocked", "s_cei", "li_cei",
                      "n_Li_dead", "sei_thickness", "surface_roughness")}
    n = 0
    for i in range(1, n_half, 2):          # discharge halves
        j = i + 1                           # following charge half
        Qd = float(dis[i]) if not np.isnan(dis[i]) else 0.0
        Qc = float(rel[j]) if j < n_half and not np.isnan(rel[j]) else np.nan
        cyc.append(n); qd.append(Qd); qc.append(Qc)
        ce.append(Qd / Qc if (Qc and not np.isnan(Qc) and Qc > 0) else np.nan)
        cap.append(Qd / q_theo * THEORETICAL_MAH_G if q_theo else np.nan)
        util.append(Qd / q_theo if q_theo else np.nan)
        for k in passthrough:
            v = t[k][i]
            passthrough[k].append(v if not np.isnan(v) else (t[k][j] if j < n_half else np.nan))
        for k in extra:                     # state at the end of the discharge
            extra[k].append(t[k][i])
        n += 1
    out = {"cycle": np.array(cyc), "Q_dis": np.array(qd), "Q_ch": np.array(qc),
           "CE_cathode": np.array(ce), "capacity_mAh_g": np.array(cap),
           "utilization": np.array(util), "n_S8_initial": np.full(len(cyc), n_s8)}
    for k, v in passthrough.items():
        out[k] = np.array(v, dtype=float)
    for k, v in extra.items():
        out[k] = np.array(v, dtype=float)
    # conventional cathode composition (formula units), from S8-unit labels
    for lab, (name, units, _li, _s) in S8_UNIT_MAP.items():
        col = f"cat_{lab}"
        if col in out:
            out[f"cathode_{name}_units"] = np.nan_to_num(out[col]) * units
    return out


def load_case(case_dir: str, pattern: str = "runs/seed_*/cycle_stats.csv"):
    """All seeds of a case -> {seed: per_cycle table}, plus conservation."""
    paths = sorted(glob.glob(os.path.join(case_dir, pattern)))
    seeds, reports = {}, {}
    for p in paths:
        seed = os.path.basename(os.path.dirname(p))
        t = load_cycle_stats(p)
        if not t:
            continue
        reports[seed] = conservation_report(t)
        seeds[seed] = per_cycle(t)
    return seeds, reports


def aggregate_seeds(seeds: Dict[str, Dict[str, np.ndarray]]) -> Dict[str, np.ndarray]:
    """Mean and std over seeds by cycle index (common prefix), n_seeds column."""
    if not seeds:
        return {}
    n_cyc = min(v["cycle"].size for v in seeds.values())
    keys = [k for k in next(iter(seeds.values())) if all(k in v for v in seeds.values())]
    out: Dict[str, np.ndarray] = {"cycle": np.arange(n_cyc), "n_seeds": np.full(n_cyc, len(seeds))}
    for k in keys:
        if k == "cycle":
            continue
        stack = np.vstack([v[k][:n_cyc] for v in seeds.values()])
        with np.errstate(all="ignore"):
            out[f"{k}_mean"] = np.nanmean(stack, axis=0)
            out[f"{k}_std"] = np.nanstd(stack, axis=0)
    return out


def write_table(table: Dict[str, np.ndarray], path: str) -> str:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    keys = list(table)
    n = len(table[keys[0]])
    with open(path, "w") as fh:
        fh.write(",".join(keys) + "\n")
        for i in range(n):
            row = []
            for k in keys:
                v = table[k][i]
                row.append("" if (isinstance(v, float) and np.isnan(v)) else f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v))
            fh.write(",".join(row) + "\n")
    return path


def summarize_case(case_dir: str, out_dir: Optional[str] = None) -> Dict[str, object]:
    """Load, check, aggregate and write per_cycle_<seed>.csv, per_cycle_mean.csv
    and conservation.txt under <case>/analysis (or out_dir)."""
    out_dir = out_dir or os.path.join(case_dir, "analysis")
    seeds, reports = load_case(case_dir)
    if not seeds:
        raise FileNotFoundError(f"no cycle_stats.csv under {case_dir}/runs/seed_*")
    for s, tbl in seeds.items():
        write_table(tbl, os.path.join(out_dir, f"per_cycle_{s}.csv"))
    agg = aggregate_seeds(seeds)
    write_table(agg, os.path.join(out_dir, "per_cycle_mean.csv"))
    with open(os.path.join(out_dir, "conservation.txt"), "w") as fh:
        for s, r in reports.items():
            fh.write(f"{s}: {'OK' if r['ok'] else 'VIOLATION'}  {r}\n")
    return {"seeds": seeds, "reports": reports, "aggregate": agg, "out_dir": out_dir}
