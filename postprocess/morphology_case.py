#!/usr/bin/env python3
"""
morphology_case.py - lattice-topological morphology descriptors for every
sampled frame of a case (buried Li, its coating, film extent, roughness,
porosity). Writes <case>/analysis/morphology.csv.

Usage:
    python postprocess/morphology_case.py cases/fullcell_cei [--frames all]
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmc3d.morphology import case_morphology   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--frames", default="first_last", choices=["first_last", "all"])
    a = ap.parse_args()
    rows = case_morphology(a.case, frames=a.frames)
    out_dir = os.path.join(a.case, "analysis"); os.makedirs(out_dir, exist_ok=True)
    keys = ["seed", "step"] + sorted({k for r in rows for k in r} - {"seed", "step"})
    path = os.path.join(out_dir, "morphology.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    print(f"{len(rows)} frames -> {path}")


if __name__ == "__main__":
    main()
