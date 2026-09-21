"""
figures.py
==========
Standard report figures from the per-cycle tables of cycles.py (matplotlib,
PNG + SVG). One axis per panel, fixed categorical hue order, thin marks,
recessive grid, text in neutral ink: the same conventions every figure of the
project follows so reports look like one system.

    from kmc3d.cycles import summarize_case
    from kmc3d.figures import case_overview
    res = summarize_case("cases/fullcell_cei")
    case_overview(res["aggregate"], "cases/fullcell_cei/analysis/overview.png",
                  title="fullcell_cei, 5 seeds")
"""

from __future__ import annotations
import os
from typing import Dict, Optional

import numpy as np

# fixed categorical order (validated reference palette, light surface)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"

# S8-unit states in lithiation order (S8 first) with conventional names
CATHODE_STATES = [("cat_S8", "S8"), ("cat_Li2S8", "Li2S8"), ("cat_Li4S8", "Li2S4 (Li4S8)"),
                  ("cat_Li8S8", "Li2S2 (Li8S8)"), ("cat_Li16S8", "Li2S (Li16S8)")]


def _style(ax, title, ylabel, xlabel="cycle"):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=3)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=10, color=INK, pad=8)
    ax.set_ylabel(ylabel, fontsize=8, color=INK2)
    ax.set_xlabel(xlabel, fontsize=8, color=INK2)


def _band(ax, x, mean, std, color, label, direct=True):
    ok = ~np.isnan(mean)
    ax.plot(x[ok], mean[ok], color=color, linewidth=2, label=label, solid_capstyle="round")
    if std is not None:
        s = np.nan_to_num(std)
        ax.fill_between(x[ok], (mean - s)[ok], (mean + s)[ok], color=color, alpha=0.15, linewidth=0)
    if direct and ok.any():   # direct label at the last point
        ax.annotate(label, (x[ok][-1], mean[ok][-1]), xytext=(4, 0), textcoords="offset points",
                    fontsize=8, color=INK2, va="center")


def case_overview(agg: Dict[str, np.ndarray], out_path: str,
                  title: Optional[str] = None, n_s8: Optional[int] = None) -> str:
    """Four panels: capacity, coulombic efficiencies, cathode composition,
    sulfur inventory. `agg` is cycles.aggregate_seeds output."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = agg["cycle"].astype(float)
    n_s8 = n_s8 or int(np.nan_to_num(agg.get("n_S8_initial_mean", [0])[0]))
    fig, axs = plt.subplots(2, 2, figsize=(10, 7), facecolor=SURFACE)
    fig.subplots_adjust(hspace=0.6, wspace=0.3, left=0.08, right=0.97, top=0.9, bottom=0.08)
    if title:
        fig.suptitle(title, x=0.08, ha="left", fontsize=12, color=INK)

    # (a) discharge capacity
    ax = axs[0, 0]
    _band(ax, x, agg["capacity_mAh_g_mean"], agg.get("capacity_mAh_g_std"), SERIES[0], "discharge")
    _style(ax, "Discharge capacity (theoretical 1672 mAh/g_S)", "mAh / g_S")
    ax.set_ylim(bottom=0)

    # (b) coulombic efficiencies
    ax = axs[0, 1]
    for i, (k, lab) in enumerate((("CE_cathode", "CE cathode (Q_dis / Q_ch)"),
                                  ("CE_shuttle", "CE shuttle (engine)"))):
        if f"{k}_mean" in agg:
            _band(ax, x, agg[f"{k}_mean"], agg.get(f"{k}_std"), SERIES[i], lab, direct=False)
    _style(ax, "Coulombic efficiency", "fraction")
    ax.set_ylim(0, 1.5)
    ax.axhline(1.0, color=GRID, linewidth=1)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower right")

    # (c) cathode composition (stacked, formula units per S8 site basis)
    ax = axs[1, 0]
    bottom = np.zeros_like(x)
    for i, (k, lab) in enumerate(CATHODE_STATES):
        if f"{k}_mean" in agg:
            y = np.nan_to_num(agg[f"{k}_mean"])
            ax.fill_between(x, bottom, bottom + y, color=SERIES[i], linewidth=0, label=lab, step=None)
            ax.plot(x, bottom + y, color=SURFACE, linewidth=0.8)   # 2px surface gap
            bottom = bottom + y
    _style(ax, "Cathode lattice sites by lithiation state", "S8-unit sites")
    if n_s8:
        ax.axhline(n_s8, color=GRID, linewidth=1)
        ax.annotate(f"initial {n_s8}", (x[0], n_s8), xytext=(4, -10), textcoords="offset points",
                    fontsize=8, color=INK2, va="center")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower left",
              bbox_to_anchor=(0.0, 1.0), ncol=3, borderaxespad=0.0)
    ax.set_title(ax.get_title(), pad=26)

    # (d) sulfur inventory as fractions of initial S
    ax = axs[1, 1]
    s_tot = 8.0 * n_s8 if n_s8 else np.nan
    parts = []
    cath = np.zeros_like(x)
    for k, _ in CATHODE_STATES:
        if f"{k}_mean" in agg:
            cath += 8 * np.nan_to_num(agg[f"{k}_mean"])
    parts.append(("cathode lattice", cath))
    dis = np.zeros_like(x)
    for k, ns in (("n_dis_Li2S8_d", 8), ("n_dis_Li2S6_d", 6), ("n_dis_Li2S4_d", 4)):
        if f"{k}_mean" in agg:
            dis += ns * np.nan_to_num(agg[f"{k}_mean"])
    parts.append(("dissolved", dis))
    if "sei_Li2S2_an_mean" in agg:
        parts.append(("anode deposit Li2S2", 2 * np.nan_to_num(agg["sei_Li2S2_an_mean"])))
    if "s_cei_mean" in agg:
        parts.append(("CEI film", np.nan_to_num(agg["s_cei_mean"])))
    bottom = np.zeros_like(x)
    for i, (lab, y) in enumerate(parts):
        f = y / s_tot if s_tot else y
        ax.fill_between(x, bottom, bottom + f, color=SERIES[i], linewidth=0, label=lab)
        ax.plot(x, bottom + f, color=SURFACE, linewidth=0.8)
        bottom = bottom + f
    _style(ax, "Sulfur inventory (fraction of initial S; total must be 1)", "fraction")
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0, color=GRID, linewidth=1)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="lower left",
              bbox_to_anchor=(0.0, 1.0), ncol=2, borderaxespad=0.0)
    ax.set_title(ax.get_title(), pad=26)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=200)
    root, _ = os.path.splitext(out_path)
    fig.savefig(root + ".svg")
    plt.close(fig)
    return out_path
