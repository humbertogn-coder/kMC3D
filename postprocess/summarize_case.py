#!/usr/bin/env python3
"""
summarize_case.py - per-cycle tables, conservation check and overview figure
for one case directory holding runs/seed_*/cycle_stats.csv.

Usage (from the repository root, kmc3d environment active):
    python postprocess/summarize_case.py cases/fullcell_cei --title "series 2"

Writes <case>/analysis/{per_cycle_<seed>.csv, per_cycle_mean.csv,
conservation.txt, overview.png, overview.svg}. Exit code 1 if any seed
violates Li or S conservation (the figure is still written).
"""
import argparse
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmc3d.cycles import summarize_case          # noqa: E402
from kmc3d.figures import case_overview          # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", help="case directory, e.g. cases/fullcell_cei")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", default=None, help="output dir (default <case>/analysis)")
    a = ap.parse_args()
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    res = summarize_case(a.case, a.out)
    title = a.title or f"{os.path.basename(a.case.rstrip('/'))}, {len(res['seeds'])} seeds"
    png = case_overview(res["aggregate"], os.path.join(res["out_dir"], "overview.png"), title=title)
    bad = [s for s, r in res["reports"].items() if not r["ok"]]
    for s, r in res["reports"].items():
        print(f"{s}: {'OK' if r['ok'] else 'CONSERVATION VIOLATION'}  li_total {r.get('li_total_range')}  s_total {r.get('s_total_range')}")
    print(f"tables -> {res['out_dir']}  figure -> {png}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
