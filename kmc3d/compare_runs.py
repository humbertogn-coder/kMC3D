"""
compare_runs.py
===============
Side-by-side comparison of Li-metal half-cell runs from their Data2Excel.txt
(the only output format shared by the original C++ engine and kmc3d), so the
legacy engine (antecedent), the physically based engine and its accelerated
variant can be put on one figure.

Per half-cycle: Li metal count, cumulative plating / stripping / electrolyte
decomposition events, and per cycle CE_cycle = stripping / plating events.
Data2Excel columns (0-based, both engines): 3 half-cycle, 14 RxnPlating,
15 RxnStripping, 17 FSI, 18 SFO, 19 SOL, 20 F5D, 21 SOL2, 22 PlatingSEI,
23 LiMetal.

    from kmc3d.compare_runs import load_run, compare_figure
    runs = {"original (C++)": load_run("Original_kMC/Data2Excel.txt"),
            "anode_physical": load_run("cases/anode_physical/runs/seed_8597/Data2Excel.txt")}
    compare_figure(runs, "analysis/halfcell_compare.png")
"""

from __future__ import annotations
import os
from typing import Dict

import numpy as np

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"


def load_run(path: str) -> Dict[str, np.ndarray]:
    """Last row of every half-cycle -> arrays over half-cycle index."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 29:
                continue
            try:
                rows.append([int(p[3]), int(p[14]), int(p[15]), int(p[17]) + int(p[18])
                             + int(p[19]) + int(p[20]) + int(p[21]), int(p[22]), int(p[23])])
            except ValueError:
                continue
    a = np.array(rows)
    halves = np.unique(a[:, 0])
    last = np.array([a[a[:, 0] == h][-1] for h in halves])
    out = {"half": last[:, 0], "plating": last[:, 1] + last[:, 4], "stripping": last[:, 2],
           "decomposition": last[:, 3], "li_metal": last[:, 5]}
    # per cycle CE: stripping in discharge half (odd) / plating in the preceding charge half
    dpl = np.diff(np.concatenate([[0], out["plating"]]))
    dst = np.diff(np.concatenate([[0], out["stripping"]]))
    ce, cyc = [], []
    for i in range(1, len(halves), 2):
        cyc.append(i // 2)
        ce.append(dst[i] / dpl[i - 1] if dpl[i - 1] > 0 else np.nan)
    out["cycle"] = np.array(cyc); out["CE_cycle"] = np.array(ce)
    out["side_fraction"] = np.where(out["plating"] > 0, out["decomposition"] / np.maximum(out["plating"], 1), np.nan)
    return out


def _style(ax, title, ylabel, xlabel="half-cycle"):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=3)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8); ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=10, color=INK, pad=8)
    ax.set_ylabel(ylabel, fontsize=8, color=INK2); ax.set_xlabel(xlabel, fontsize=8, color=INK2)


def compare_figure(runs: Dict[str, Dict[str, np.ndarray]], out_path: str, title: str = "") -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(10, 7), facecolor=SURFACE)
    fig.subplots_adjust(hspace=0.5, wspace=0.3, left=0.08, right=0.97, top=0.9, bottom=0.08)
    if title:
        fig.suptitle(title, x=0.08, ha="left", fontsize=12, color=INK)
    panels = [("li_metal", "Li metal sites", "sites", axs[0, 0]),
              ("decomposition", "Cumulative electrolyte decomposition events", "events", axs[0, 1]),
              ("side_fraction", "Decomposition events per plating event", "fraction", axs[1, 0])]
    for key, ttl, yl, ax in panels:
        for i, (name, r) in enumerate(runs.items()):
            ax.plot(r["half"], r[key], color=SERIES[i], linewidth=2, label=name)
        _style(ax, ttl, yl)
        if key == "side_fraction":
            ax.set_yscale("log"); ax.axhline(2e-3, color=GRID, linewidth=1)
            ax.annotate("CE 0.998 target", (ax.get_xlim()[0], 2e-3), xytext=(4, 3),
                        textcoords="offset points", fontsize=8, color=INK2)
    ax = axs[1, 1]
    for i, (name, r) in enumerate(runs.items()):
        ax.plot(r["cycle"], r["CE_cycle"], color=SERIES[i], linewidth=2, label=name)
    _style(ax, "CE_cycle = stripping / plating events", "fraction", xlabel="cycle")
    ax.axhline(1.0, color=GRID, linewidth=1); ax.set_ylim(0, 1.6)
    axs[0, 0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=200); fig.savefig(os.path.splitext(out_path)[0] + ".svg")
    plt.close(fig)
    return out_path
