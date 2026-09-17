#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
li_buried_sweep_cfg.py
══════════════════════
Computes Li buried for every case in the config, using TWO criteria:

  (A) DYNAMIC TOPOLOGICAL  (MAINBODY, local_density_thr from config):
      multi-frame criterion — Li above original anode, present in all
      frames of the window, AND (not in main connected body OR low local
      density). Reported per cycle (cycle_step) and at the DB cycle.

  (B) KINETIC  (independent consistency check, no topology):
      Li_buried_kinetic(t) = (RxnPlating + RxnPlatingSEI - RxnStripping)
                             - (LiMetal(t) - LiMetal(0))
      i.e. cumulative net plating minus the change in metallic Li count.

Outputs to CFG['out_dir']:
  - li_buried_per_cycle.csv      one row per (k0, cycle): both criteria
  - li_buried_ML_features.csv    one row per k0: aggregated, for the DB
  - fig_li_buried_dynamic_vs_k0.png
  - fig_li_buried_kinetic_vs_k0.png

Run:  python li_buried_sweep_cfg.py --config config_plating
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation, convolve, label
from scipy.optimize import curve_fit

import kmc_io as io

warnings.filterwarnings("ignore")

# ── Voxelization constants (match v6/v7/v8) ─────────────────────────────
LI_SPECIES   = "Li"
IGNORE       = {"ETH", "SOL", "FSI"}
VOXEL        = 0.8
SIGMA_FRAC   = 0.45
CUTOFF_FRAC  = 2.5
GRID_PAD     = 2.0
DENSITY_THR  = 0.12
WINDOW_N     = 10
ALWAYS_TOL   = 1
MAIN_MIN_VOX = 50
LDR_VOX      = 3
RADIUS       = {"Li": 1.73}
DEF_R        = 2.0


def anode_top(coords, species, half):
    z = coords[species == LI_SPECIES, 2]
    return float(np.median(z)) + half if z.size else 0.0


def _kernel(r):
    sig = max(1e-6, SIGMA_FRAC * r); rc = CUTOFF_FRAC * sig
    h = int(np.ceil(rc / VOXEL)); ax = np.arange(-h, h + 1) * VOXEL
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    K = np.exp(-(R**2) / (2 * sig**2)); K[R > rc] = 0
    return K.astype(np.float32)


def _stamp(field, K, idx):
    hx, hy, hz = [s // 2 for s in K.shape]
    ix, iy, iz = idx
    x0, x1 = ix-hx, ix+hx+1; y0, y1 = iy-hy, iy+hy+1; z0, z1 = iz-hz, iz+hz+1
    fx0, fx1 = max(0,x0), min(field.shape[0],x1)
    fy0, fy1 = max(0,y0), min(field.shape[1],y1)
    fz0, fz1 = max(0,z0), min(field.shape[2],z1)
    if fx0>=fx1 or fy0>=fy1 or fz0>=fz1: return
    field[fx0:fx1, fy0:fy1, fz0:fz1] += K[fx0-x0:K.shape[0]-(x1-fx1),
                                           fy0-y0:K.shape[1]-(y1-fy1),
                                           fz0-z0:K.shape[2]-(z1-fz1)]


def build_li_density(coords, species, mn, dims):
    f = np.zeros(dims, np.float32); cache = {}
    for p, s in zip(coords, species):
        if s != LI_SPECIES: continue
        r = RADIUS.get(s, DEF_R)
        if r not in cache: cache[r] = _kernel(r)
        idx = tuple(np.floor((p - mn) / VOXEL).astype(int))
        if all(0 <= i < d for i, d in zip(idx, dims)):
            _stamp(f, cache[r], idx)
    return f


def conn_struct():
    s = np.zeros((3,3,3), np.int32)
    s[1,1,0]=s[1,1,2]=s[1,0,1]=s[1,2,1]=s[0,1,1]=s[2,1,1]=s[1,1,1]=1
    return s


def main_body(sol):
    if not sol.any(): return np.zeros_like(sol, bool)
    lbl, n = label(sol, structure=conn_struct())
    if n == 0: return np.zeros_like(sol, bool)
    cnt = np.bincount(lbl.ravel()); cnt[0] = 0
    m = int(np.argmax(cnt))
    return (lbl == m) if cnt[m] >= MAIN_MIN_VOX else np.zeros_like(sol, bool)


def low_local_density(sol, thr):
    if not sol.any(): return np.zeros_like(sol, bool)
    r = LDR_VOX; ax = np.arange(-r, r+1)
    X,Y,Z = np.meshgrid(ax,ax,ax,indexing="ij")
    sph = (X**2+Y**2+Z**2) <= r*r; nk = int(sph.sum())
    cnt = convolve(sol.astype(np.int32), sph.astype(np.int32),
                   mode="constant", cval=0)
    return sol & ((cnt.astype(np.float32)/max(nk,1)) < thr)


def compute_dynamic_buried(frame_files, half, thr):
    """MAINBODY dynamic buried for one cycle window. Returns n_buried,
    n_li_some_filtered, frac."""
    data = []
    for p in frame_files:
        fr = io.parse_xyz_last_frame(p)
        co = io.minimal_image(fr["coords"], fr["cell"])
        sp = fr["species"]
        at = anode_top(co, sp, half)
        mask = ~np.isin(sp, list(IGNORE)) & (co[:, 2] >= at)
        c = co[mask].copy(); c[:, 2] -= at
        data.append((c, sp[mask], at))
    if len(data) < 2: return None

    xy_min = np.array([np.inf, np.inf]); xy_max = -xy_min.copy()
    z_min, z_max = np.inf, -np.inf
    for c, _, _ in data:
        if c.shape[0] == 0: continue
        xy_min = np.minimum(xy_min, c[:, :2].min(0))
        xy_max = np.maximum(xy_max, c[:, :2].max(0))
        z_min = min(z_min, c[:, 2].min()); z_max = max(z_max, c[:, 2].max())
    if not np.isfinite(z_min): return None
    mn = np.array([xy_min[0], xy_min[1], z_min]) - GRID_PAD
    mx = np.array([xy_max[0], xy_max[1], z_max]) + GRID_PAD
    dims = np.maximum(np.ceil((mx - mn) / VOXEL).astype(int) + 1, 3)

    struct = conn_struct(); masks = []
    for c, sp, _ in data:
        f = build_li_density(c, sp, mn, dims)
        raw = f >= DENSITY_THR
        masks.append(binary_dilation(raw, structure=struct,
                                     iterations=ALWAYS_TOL) if ALWAYS_TOL else raw)

    li_some = np.zeros(dims, bool); li_always = np.ones(dims, bool)
    for m in masks: li_some |= m; li_always &= m

    last = masks[-1]
    iso = (last & ~main_body(last)) | low_local_density(last, thr)
    buried = li_always & iso

    return {
        "n_li_buried": int(buried.sum()),
        "n_li_some":   int(li_some.sum()),
        "frac_buried": float(buried.sum() / max(li_some.sum(), 1)),
    }


def window_files(cycle_to_xyz, cycle, n=WINDOW_N):
    """Collect up to n xyz files ending at `cycle` (extend backwards)."""
    collected = []; c = int(cycle)
    while c > 0 and len(collected) < n:
        if c in cycle_to_xyz:
            collected.insert(0, cycle_to_xyz[c][1])
        c -= 1
    return collected[-n:]


# ════════════════════════════════════════════════════════════════════════
def process_case(case_key, k0, cfg):
    base = cfg["base_dir"] / case_key
    traj = base / cfg["traj_subpath"]
    log  = io.resolve_log(base, cfg)
    half = cfg["anode_z_width"] / 2.0

    if not traj.is_dir() or log is None:
        print(f"  [{case_key}] missing traj or log; skip", flush=True)
        return None, None

    df_log = io.parse_data2excel(log, full=True)
    step_to_cycle = io.build_step_to_cycle(df_log)
    xyz_pairs = io.list_xyz_files(traj)
    cyc2xyz = io.map_cycles_to_xyz(xyz_pairs, step_to_cycle)

    # Per-cycle cumulative counters for the KINETIC criterion
    per_cyc = io.per_cycle_last(
        df_log, ["RxnPlating", "RxnPlatingSEI", "RxnStripping", "LiMetal"])
    li_metal_0 = int(df_log["LiMetal"].iloc[0])

    target_cycles = list(range(cfg["cycle_step"],
                               cfg["max_cycle"] + 1, cfg["cycle_step"]))
    rows = []
    for cyc in target_cycles:
        # ── kinetic ──
        prow = per_cyc[per_cyc["cycle"] == cyc]
        if prow.empty:
            kinetic = np.nan
        else:
            prow = prow.iloc[0]
            kinetic = ((prow["RxnPlating"] + prow["RxnPlatingSEI"]
                        - prow["RxnStripping"])
                       - (prow["LiMetal"] - li_metal_0))
        # ── dynamic topological ──
        wf = window_files(cyc2xyz, cyc)
        dyn = compute_dynamic_buried(wf, half, cfg["local_density_thr"]) \
            if len(wf) >= 2 else None

        rows.append({
            "case_key":            case_key,
            "k0":                  k0,
            "cycle":               cyc,
            "n_li_buried_dynamic": (dyn["n_li_buried"] if dyn else np.nan),
            "frac_buried_dynamic": (dyn["frac_buried"] if dyn else np.nan),
            "li_buried_kinetic":   kinetic,
        })
        bd = dyn["n_li_buried"] if dyn else -1
        print(f"  [{case_key}] cyc {cyc:3d}: dyn={bd}  kin={kinetic:.0f}",
              flush=True)

    df = pd.DataFrame(rows)

    # ── Suppress spurious single-cycle drops in the DYNAMIC metric ──────
    # The topological criterion occasionally relabels the connected main
    # body between frames, causing a one-cycle collapse (e.g. 15460 -> 6709
    # -> 16194). These are not physical. We add a rolling-median-smoothed
    # column (window 3) and use it for the aggregates, while keeping the
    # raw column for transparency / QA.
    if "n_li_buried_dynamic" in df.columns and len(df) >= 3:
        df["n_li_buried_dynamic_smooth"] = (
            df["n_li_buried_dynamic"]
            .rolling(window=3, center=True, min_periods=1).median())
    else:
        df["n_li_buried_dynamic_smooth"] = df.get("n_li_buried_dynamic",
                                                  np.nan)

    # Aggregate for ML DB (one row per k0). Dynamic aggregates use the
    # smoothed column; kinetic is already monotone so it is used as-is.
    db_cyc = cfg["li_buried_cycle_for_db"]
    at_db = df[df["cycle"] == db_cyc]
    dyn_db = (float(at_db["n_li_buried_dynamic_smooth"].iloc[0])
              if len(at_db) else np.nan)
    dyn_db_raw = (float(at_db["n_li_buried_dynamic"].iloc[0])
                  if len(at_db) else np.nan)
    kin_db = float(at_db["li_buried_kinetic"].iloc[0]) if len(at_db) else np.nan
    agg = {
        "case_key": case_key, "k0": k0,
        "n_li_buried_dynamic_mean": float(df["n_li_buried_dynamic_smooth"].mean()),
        "n_li_buried_dynamic_max":  float(df["n_li_buried_dynamic_smooth"].max()),
        "n_li_buried_dynamic_at_db_cycle": dyn_db,
        "n_li_buried_dynamic_at_db_cycle_raw": dyn_db_raw,
        "li_buried_kinetic_mean":   float(df["li_buried_kinetic"].mean()),
        "li_buried_kinetic_final":  float(df["li_buried_kinetic"].iloc[-1]),
        "li_buried_kinetic_at_db_cycle": kin_db,
        "db_cycle": db_cyc,
    }
    return df, agg


# ════════════════════════════════════════════════════════════════════════
def fig_vs_k0(df_ml, ycol, ylabel, title, out_path, cfg):
    plt.rcParams.update(cfg["plot_style"])
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    d = df_ml.sort_values("k0")
    x = d["k0"].values; y = d[ycol].values
    ax.plot(x, y, "o", color="#2980B9", markersize=7,
            markeredgecolor="white", markeredgewidth=0.8, label="Data")
    # exponential/logistic fit if enough points
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() >= 4:
        def f(k, a, b, c): return a * (1 - np.exp(-b * k)) + c
        try:
            p, _ = curve_fit(f, x[mask], y[mask],
                             p0=(np.nanmax(y), 1/np.median(x), np.nanmin(y)),
                             maxfev=10000)
            xf = np.logspace(np.log10(x.min()*0.7), np.log10(x.max()*1.4), 300)
            yf = f(xf, *p)
            ss = 1 - np.sum((y[mask]-f(x[mask],*p))**2)/np.sum((y[mask]-y[mask].mean())**2)
            ax.plot(xf, yf, "-", color="#27AE60", lw=2,
                    label=f"Exp. fit (R²={ss:.3f})")
        except Exception:
            pass
    ax.set_xscale("log")
    ax.set_xlabel(cfg["rate_symbol_math"] + " (plating rate)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold")
    ax.grid(alpha=0.3, which="both", lw=0.5)
    ax.legend(loc="best")
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--only-case", default=None,
                    help="Process ONE case folder and write partial CSVs to "
                         "out_dir/li_buried_parts/. Use for SLURM array jobs.")
    ap.add_argument("--merge", action="store_true",
                    help="Merge all partial CSVs in out_dir/li_buried_parts/ "
                         "into the final li_buried_*.csv + figures.")
    args = ap.parse_args()
    cfg = io.load_config(args.config)
    out = cfg["out_dir"]; out.mkdir(parents=True, exist_ok=True)
    parts_dir = out / "li_buried_parts"; parts_dir.mkdir(exist_ok=True)

    # ── ARRAY MODE: one case, write partial files, exit ─────────────────
    if args.only_case is not None:
        match = [(k, v) for k, v in cfg["cases"] if k == args.only_case]
        if not match:
            print(f"case_key {args.only_case} not in config; nothing to do")
            return
        case_key, k0 = match[0]
        print(f"Li buried (array mode) — case {case_key} (k0={k0})")
        df, agg = process_case(case_key, k0, cfg)
        if df is None:
            print(f"  case {case_key}: no data"); return
        df.to_csv(parts_dir / f"percycle_case_{case_key}.csv", index=False)
        pd.DataFrame([agg]).to_csv(
            parts_dir / f"agg_case_{case_key}.csv", index=False)
        print(f"  wrote partial CSVs for case {case_key}")
        return

    # ── MERGE MODE: gather partials, write finals + figures ─────────────
    if args.merge:
        pc = sorted(parts_dir.glob("percycle_case_*.csv"))
        ag = sorted(parts_dir.glob("agg_case_*.csv"))
        if not pc:
            print(f"No partial CSVs in {parts_dir}; run array jobs first.")
            return
        df_all = pd.concat([pd.read_csv(p) for p in pc], ignore_index=True)
        df_ml = (pd.concat([pd.read_csv(p) for p in ag], ignore_index=True)
                 .sort_values("k0").reset_index(drop=True))
        df_all.to_csv(out / "li_buried_per_cycle.csv", index=False)
        df_ml.to_csv(out / "li_buried_ML_features.csv", index=False)
        print(f"Merged {len(pc)} cases -> li_buried_per_cycle.csv "
              f"({len(df_all)} rows), li_buried_ML_features.csv "
              f"({len(df_ml)} rows)")
        fig_vs_k0(df_ml, "n_li_buried_dynamic_at_db_cycle",
                  f"Dynamic Li buried at cycle {cfg['li_buried_cycle_for_db']} (voxels)",
                  "Dynamic (topological) Li buried vs k_p",
                  out / "fig_li_buried_dynamic_vs_k0.png", cfg)
        fig_vs_k0(df_ml, "li_buried_kinetic_final",
                  "Kinetic Li buried (final cycle)",
                  "Kinetic Li buried vs k_p",
                  out / "fig_li_buried_kinetic_vs_k0.png", cfg)
        print("Saved figures.")
        return

    # ── SERIAL MODE (original): all cases in one process ────────────────
    print(f"Li buried sweep — {cfg['sweep_name']}")
    all_rows, all_agg = [], []
    for case_key, k0 in cfg["cases"]:
        print(f"\n── case {case_key} (k0={k0}) ──", flush=True)
        df, agg = process_case(case_key, k0, cfg)
        if df is not None:
            all_rows.append(df); all_agg.append(agg)

    if not all_rows:
        print("No data processed."); return

    df_all = pd.concat(all_rows, ignore_index=True)
    df_all.to_csv(out / "li_buried_per_cycle.csv", index=False)
    df_ml = pd.DataFrame(all_agg).sort_values("k0").reset_index(drop=True)
    df_ml.to_csv(out / "li_buried_ML_features.csv", index=False)
    print(f"\nSaved li_buried_per_cycle.csv ({len(df_all)} rows)")
    print(f"Saved li_buried_ML_features.csv ({len(df_ml)} rows)")

    fig_vs_k0(df_ml, "n_li_buried_dynamic_at_db_cycle",
              f"Dynamic Li buried at cycle {cfg['li_buried_cycle_for_db']} (voxels)",
              "Dynamic (topological) Li buried vs k_p",
              out / "fig_li_buried_dynamic_vs_k0.png", cfg)
    fig_vs_k0(df_ml, "li_buried_kinetic_final",
              "Kinetic Li buried (final cycle)",
              "Kinetic Li buried vs k_p",
              out / "fig_li_buried_kinetic_vs_k0.png", cfg)
    print("Saved figures.")


if __name__ == "__main__":
    main()
