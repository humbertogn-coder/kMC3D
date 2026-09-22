#!/usr/bin/env python3
"""
compare_halfcell.py - one figure with several Li-metal half-cell runs from
their Data2Excel.txt (works for the original C++ engine output as well).

Usage:
    python postprocess/compare_halfcell.py out.png "original=../Original_kMC/Data2Excel.txt" \\
        "anode_physical=cases/anode_physical/runs/seed_8597/Data2Excel.txt" \\
        "A_SEI=100=cases/anode_physical_A100/runs/seed_8597/Data2Excel.txt"
Each positional argument after the output path is label=path (the label may
contain '=' as long as the path is the part after the LAST '=').
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmc3d.compare_runs import load_run, compare_figure   # noqa: E402


def main():
    if len(sys.argv) < 3:
        print(__doc__); return 1
    out = sys.argv[1]
    runs = {}
    for arg in sys.argv[2:]:
        label, path = arg.rsplit("=", 1)
        runs[label] = load_run(path)
    print(compare_figure(runs, out, title="Li-metal half cell: engine comparison"))
    for k, r in runs.items():
        print(f"{k}: halves {r['half'].size}, plating {r['plating'][-1]}, stripping {r['stripping'][-1]}, "
              f"decomposition {r['decomposition'][-1]}, Li {r['li_metal'][0]} -> {r['li_metal'][-1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
