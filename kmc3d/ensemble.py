"""
ensemble.py - multi-seed ensembles and their aggregation.

kMC trajectories are stochastic realisations: physical conclusions live in the
ensemble. This module (a) runs N seeds of the same inputs sequentially, and
(b) aggregates their cycle_stats.csv into mean +/- std tables.

On GRACE prefer the SLURM job-array template (goKMC_grace_array.slrm): each
array task runs one seed in parallel; then aggregate with:

    python -m kmc3d.ensemble --aggregate "runs/seed_*/cycle_stats.csv" --out agg

Local sequential use:

    python -m kmc3d.ensemble --dir base_case --out runs --nseeds 5 --base-seed 8597
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
from typing import Dict, List

import numpy as np


# ------------------------------------------------------------------ running
def run_seeds(input_dir: str, out_root: str, seeds: List[int]) -> List[str]:
    """Run one simulation per seed (sequential). Returns run directories."""
    from . import config, lattice
    from .mechanism import Mechanism
    from .engine import Engine
    from .output import Output

    dirs = []
    for sd in seeds:
        p = config.Parameters.read(os.path.join(input_dir, "PARAMETERS.in"))
        g = config.Geometry.read(os.path.join(input_dir, "GEOMETRY.in"))
        d = config.DecompositionKinetics.read(
            os.path.join(input_dir, "DECOMPOSITION.in"), p.numberDecompositionRxns)
        m = config.MobilityKinetics.read(
            os.path.join(input_dir, "MOBILITY.in"), p.nMobilityReactions)
        mech = Mechanism.read(os.path.join(input_dir, "MECHANISM.in"))
        p.seed = int(sd)
        run_dir = os.path.join(out_root, f"seed_{sd}")
        eng = Engine(lattice.generate(g), p, g, d, m, mech, Output(root=run_dir))
        eng.run()
        dirs.append(run_dir)
    return dirs


# -------------------------------------------------------------- aggregation
def load_cycle_stats(path: str) -> Dict[str, np.ndarray]:
    with open(path) as fh:
        rd = csv.DictReader(fh)
        rows = list(rd)
    if not rows:
        return {}
    out: Dict[str, np.ndarray] = {}
    for k in rows[0]:
        vals = []
        for r in rows:
            v = r.get(k, "")
            try:
                vals.append(float(v) if v != "" else np.nan)
            except ValueError:
                vals.append(np.nan)     # non-numeric (phase) -> skip numerically
        out[k] = np.array(vals)
    out["phase"] = np.array([r.get("phase", "") for r in rows])
    return out


def aggregate(pattern: str, out_dir: str) -> str:
    """Aggregate many cycle_stats.csv (glob pattern) into mean/std by half_index.

    Seeds may end at different half counts; aggregation uses the common prefix.
    Writes <out_dir>/cycle_stats_mean.csv and _std.csv; returns the mean path.
    """
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files match {pattern}")
    tables = [load_cycle_stats(p) for p in paths]
    tables = [t for t in tables if t]
    n_half = min(t["half_index"].size for t in tables)
    keys = [k for k in tables[0]
            if k != "phase" and all(k in t for t in tables)]

    os.makedirs(out_dir, exist_ok=True)
    mean_path = os.path.join(out_dir, "cycle_stats_mean.csv")
    std_path = os.path.join(out_dir, "cycle_stats_std.csv")
    stack = {k: np.vstack([t[k][:n_half] for t in tables]) for k in keys}
    with open(mean_path, "w") as fm, open(std_path, "w") as fs:
        fm.write(",".join(keys) + ",n_seeds\n")
        fs.write(",".join(keys) + ",n_seeds\n")
        for i in range(n_half):
            mvals, svals = [], []
            for k in keys:
                col = stack[k][:, i]
                ok = ~np.isnan(col)
                mvals.append(f"{np.nanmean(col):.6g}" if ok.any() else "")
                svals.append(f"{np.nanstd(col):.6g}" if ok.any() else "")
            fm.write(",".join(mvals) + f",{len(tables)}\n")
            fs.write(",".join(svals) + f",{len(tables)}\n")
    return mean_path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default="", help="input directory with *.in files")
    ap.add_argument("--out", default="runs", help="output root")
    ap.add_argument("--nseeds", type=int, default=0)
    ap.add_argument("--base-seed", type=int, default=8597)
    ap.add_argument("--seeds", type=int, nargs="*", default=None)
    ap.add_argument("--aggregate", default="",
                    help="glob of cycle_stats.csv files to aggregate (skips running)")
    args = ap.parse_args()

    if args.aggregate:
        p = aggregate(args.aggregate, args.out)
        print(f"aggregated -> {p}")
        return
    seeds = args.seeds or [args.base_seed + i for i in range(max(args.nseeds, 1))]
    dirs = run_seeds(args.dir, args.out, seeds)
    try:
        p = aggregate(os.path.join(args.out, "seed_*", "cycle_stats.csv"), args.out)
        print(f"ran {len(dirs)} seeds; aggregated -> {p}")
    except FileNotFoundError:
        print(f"ran {len(dirs)} seeds (no cycle_stats found to aggregate)")


if __name__ == "__main__":
    main()
