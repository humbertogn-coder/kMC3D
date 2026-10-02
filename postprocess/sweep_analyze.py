#!/usr/bin/env python3
"""
sweep_analyze.py - analyse a one-parameter sweep built by sweep_make.py:
per-run metrics from the half-cycle ledger (cycle_stats.csv) and from the
first/last trajectory frames, mean and spread over seeds for every member,
a sensitivity table with significance tests, per-cycle curves, and a set of
slide-ready figures with the swept parameter in physical units.

Usage (repository root, kmc3d environment active):

    python postprocess/sweep_analyze.py cases/sweep_plating_j0 [--skip 10] [--no-xyz]

Outputs under <sweep>/analysis/:
    per_run.csv        one row per (member, seed): every metric
    summary.csv        mean, std, n over seeds per member
    sensitivity.csv    slope vs factor, % change per +50 %, p-values, ranking
    curves_<m>.csv     per-cycle mean over seeds, one column per member
    porosity_profiles.csv
    figures/*.png      F01 .. F11 (16:9-friendly, one per slide)
    README.md          what each file and figure contains, with the numbers

Metric definitions (half cell, one row of the ledger per half-cycle):
    plated / stripped      Plating + PlatingSEI events in the charge half,
                           LiStripping events in the discharge half
    side events            FSI + SFO + SOL + SOL2 + F5D events (both halves)
    side_per_1000_plated   side events per 1000 plating events over the window
    sei_rate               slope of n_SEI vs cycle over the window (sites/cycle)
    li_lost_per_1000       increase of n_Li_ion (Li bound to SEI, the proxy
                           for dead Li) per 1000 Li plated over the window
    loss_mod_permille      (1 - CE_mod) x 1000, CE_mod = stripped / (plated +
                           decomposition electrons)
    j_plate / j_strip      e x events / (dt_half x A) with A = both faces of
                           the slab (nx a x ny a each), in mA/cm2
    xyz descriptors        kmc3d.morphology on the last frame: topological
                           buried Li, islands, film extent, film porosity
The window is cycles >= --skip (default 10) up to the shortest run.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys
import warnings
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmc3d.cycles import load_cycle_stats, recover_cycle_stats   # noqa: E402

E_CHARGE = 1.602176634e-19
SIDE = ("rxn_FSI", "rxn_SFO", "rxn_SOL", "rxn_SOL2", "rxn_F5D")
INORG = ("sei_F", "sei_O", "sei_N", "sei_S")
BC_FRACTION = 2.0 / 16.0      # body-centre sites per cube: what dense Li metal fills
ORG = ("sei_SFO", "sei_F5D")

LABELS = {"j_plate_mA_cm2": "plating current density", "j_strip_mA_cm2": "stripping current density",
          "dt_charge_s": "time per charge half", "side_per_1000_plated": "side events / 1000 plated",
          "FSI_per_1000_plated": "FSI reduction / 1000", "SOL_per_1000_plated": "F5DEE reduction / 1000",
          "SFO_per_1000_plated": "SFO step / 1000", "SOL2_per_1000_plated": "F5DEE on Li2O / 1000",
          "sei_rate": "SEI growth (sites/cycle)", "sei_final": "SEI sites, last cycle",
          "sei_thickness_mean_A": "SEI z extent", "sei_inorganic_fraction": "SEI inorganic fraction",
          "li_lost_per_1000_plated": "Li bound to SEI / 1000 plated", "li_ion_final": "Li bound to SEI, last cycle",
          "li_dead_final": "dead Li (disconnected)", "roughness_mean_A": "Li front roughness",
          "surface_diff_per_cycle": "Li surface diffusion / cycle", "ce_event": "CE (electron budget)",
          "loss_event_permille": "charge lost (permille)", "solvent_consumed": "solvent consumed",
          "salt_consumed": "salt consumed", "xyz_li_buried": "buried Li (islands)", "xyz_islands": "Li islands",
          "xyz_film_extent_A": "film extent", "xyz_film_porosity_rel": "film porosity",
          "cycle_ce_below_0.9": "cycle where CE_cycle < 0.9", "cycle_failure": "cycle of anode failure",
          "stripped_last10_over_first10": "stripping retention (last/first 10)",
          "li_ion_fraction_final": "fraction of anode Li bound to SEI"}

# figure conventions (same system as kmc3d.figures)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"


# ----------------------------------------------------------------- helpers
def read_geometry(path: str) -> Dict[str, float]:
    g: Dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                t = s.split()
                try:
                    g[t[0]] = float(t[1])
                except (ValueError, IndexError):
                    g[t[0]] = t[1] if len(t) > 1 else ""
    return g


def _col(t, k, n):
    v = t.get(k)
    if v is None:
        return np.zeros(n)
    return np.nan_to_num(v.astype(float))


def _slope(x, y):
    ok = ~np.isnan(y)
    if ok.sum() < 3:
        return np.nan
    return float(np.polyfit(x[ok], y[ok], 1)[0])


def run_metrics(run_dir: str, area_cm2: float, skip: int, max_cycle: Optional[int]) -> Optional[Dict[str, float]]:
    recover_cycle_stats(run_dir)
    p = os.path.join(run_dir, "cycle_stats.csv")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    t = load_cycle_stats(p)
    if not t:
        return None
    n = t["half_index"].size
    phase = t["phase"]
    ch = [i for i in range(n - 1) if str(phase[i]).startswith("ch") and str(phase[i + 1]).startswith("dis")]
    if not ch:
        return None
    time = _col(t, "time", n)
    dt = np.diff(np.concatenate([[0.0], time]))
    plated_h = _col(t, "rxn_Plating", n) + _col(t, "rxn_PlatingSEI", n)
    strip_h = _col(t, "rxn_LiStripping", n)
    side_h = sum(_col(t, k, n) for k in SIDE)
    per_side = {k: _col(t, k, n) for k in SIDE}
    surf_h = _col(t, "rxn_LiSurface", n)
    m: Dict[str, np.ndarray] = {}
    idx_c = np.array(ch); idx_d = idx_c + 1
    m["cycle"] = np.arange(1, idx_c.size + 1)
    m["plated"] = plated_h[idx_c]
    m["stripped"] = strip_h[idx_d]
    m["side"] = side_h[idx_c] + side_h[idx_d]
    for k in SIDE:
        m[k] = per_side[k][idx_c] + per_side[k][idx_d]
    m["surface_diff"] = surf_h[idx_c] + surf_h[idx_d]
    m["dt_charge"] = dt[idx_c]
    m["dt_discharge"] = dt[idx_d]
    with np.errstate(divide="ignore", invalid="ignore"):
        m["j_plate"] = np.where(m["dt_charge"] > 0, E_CHARGE * m["plated"] / (m["dt_charge"] * area_cm2) * 1e3, np.nan)
        m["j_strip"] = np.where(m["dt_discharge"] > 0, E_CHARGE * m["stripped"] / (m["dt_discharge"] * area_cm2) * 1e3, np.nan)
    for k in ("n_SEI", "n_Li_ion", "n_Li_dead", "n_Li_interface", "n_Li", "sei_thickness",
              "surface_roughness", "li_front_z", "n_SOL", "n_FSI"):
        m[k] = _col(t, k, n)[idx_d]           # state at the end of the cycle
    for k in INORG + ORG:
        if k in t:
            m[k] = _col(t, k, n)[idx_d]
    for k in ("CE_cycle", "CE_mod"):
        v = t.get(k)
        m[k] = v[idx_d].astype(float) if v is not None else np.full(idx_c.size, np.nan)
    m["cum_side"] = np.cumsum(m["side"])

    nc = m["cycle"].size
    if max_cycle:
        nc = min(nc, max_cycle)
    w = (m["cycle"] >= skip) & (m["cycle"] <= nc)
    if w.sum() < 3:
        w = m["cycle"] <= nc
    out: Dict[str, float] = {"n_cycles": float(nc), "window_from": float(m["cycle"][w][0]),
                             "window_to": float(m["cycle"][w][-1])}
    plated_w = m["plated"][w].sum()
    out["plated_per_cycle"] = float(m["plated"][w].mean())
    out["stripped_per_cycle"] = float(m["stripped"][w].mean())
    out["side_per_cycle"] = float(m["side"][w].mean())
    out["side_per_1000_plated"] = float(m["side"][w].sum() / plated_w * 1e3) if plated_w else np.nan
    for k in SIDE:
        out[f"{k[4:]}_per_1000_plated"] = float(m[k][w].sum() / plated_w * 1e3) if plated_w else np.nan
    out["surface_diff_per_cycle"] = float(m["surface_diff"][w].mean())
    out["dt_charge_s"] = float(np.nanmean(m["dt_charge"][w]))
    out["dt_discharge_s"] = float(np.nanmean(m["dt_discharge"][w]))
    out["j_plate_mA_cm2"] = float(np.nanmean(m["j_plate"][w]))
    out["j_strip_mA_cm2"] = float(np.nanmean(m["j_strip"][w]))
    out["time_total_s"] = float(time[idx_d[nc - 1]])
    out["sei_final"] = float(m["n_SEI"][nc - 1])
    out["sei_rate"] = _slope(m["cycle"][w].astype(float), m["n_SEI"][w])
    out["sei_thickness_final_A"] = float(m["sei_thickness"][nc - 1])
    out["sei_thickness_mean_A"] = float(m["sei_thickness"][w].mean())
    tot_in = sum(m[k][nc - 1] for k in INORG if k in m)
    tot_org = sum(m[k][nc - 1] for k in ORG if k in m)
    out["sei_inorganic_final"] = float(tot_in)
    out["sei_organic_final"] = float(tot_org)
    out["sei_inorganic_fraction"] = float(tot_in / (tot_in + tot_org)) if (tot_in + tot_org) else np.nan
    for k in INORG + ORG:
        if k in m:
            out[f"{k}_final"] = float(m[k][nc - 1])
    out["li_ion_final"] = float(m["n_Li_ion"][nc - 1])
    out["li_ion_rate"] = _slope(m["cycle"][w].astype(float), m["n_Li_ion"][w])
    d_ion = m["n_Li_ion"][w][-1] - m["n_Li_ion"][w][0]
    out["li_lost_per_1000_plated"] = float(d_ion / plated_w * 1e3) if plated_w else np.nan
    out["li_dead_final"] = float(m["n_Li_dead"][nc - 1])
    out["li_interface_final"] = float(m["n_Li_interface"][nc - 1])
    out["roughness_mean_A"] = float(m["surface_roughness"][w].mean())
    out["roughness_final_A"] = float(m["surface_roughness"][nc - 1])
    out["li_front_mean_A"] = float(m["li_front_z"][w].mean())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out["ce_cycle_mean"] = float(np.nanmean(m["CE_cycle"][w]))
        out["ce_mod_mean"] = float(np.nanmean(m["CE_mod"][w]))
        out["loss_mod_permille"] = float(np.nanmean((1.0 - m["CE_mod"][w]) * 1e3))
    # CE implied by the electron budget: every side event costs one electron
    # that plates no Li. With half-cycles closed by an event budget, CE_cycle
    # and CE_mod are fixed by the budget (stripped = budget), so this is the
    # CE that carries the physics of the sweep.
    tot_side = m["side"][w].sum()
    out["ce_event"] = float(plated_w / (plated_w + tot_side)) if (plated_w + tot_side) else np.nan
    out["loss_event_permille"] = float(1e3 * tot_side / (plated_w + tot_side)) if (plated_w + tot_side) else np.nan
    out["solvent_consumed"] = float(_col(t, "n_SOL", n)[0] - m["n_SOL"][nc - 1]) if "n_SOL" in t else np.nan
    out["salt_consumed"] = float(_col(t, "n_FSI", n)[0] - m["n_FSI"][nc - 1]) if "n_FSI" in t else np.nan
    # aging: the discharge half strips less than the charge half plated once the
    # strippable Li runs short (end_half_when_idle closes the half early)
    with np.errstate(divide="ignore", invalid="ignore"):
        ce_c = np.where(m["plated"] > 0, m["stripped"] / m["plated"], np.nan)
    m["CE_cycle_calc"] = ce_c
    out["ce_cycle_calc_mean"] = float(np.nanmean(ce_c[w]))
    out["stripped_last10_over_first10"] = float(m["stripped"][max(nc - 10, 0):nc].mean()
                                                / max(m["stripped"][:min(10, nc)].mean(), 1e-9))
    for thr in (0.9, 0.5):
        below = np.where(np.convolve((ce_c < thr).astype(float), np.ones(3), "valid") >= 3)[0]
        out[f"cycle_ce_below_{thr}"] = float(m["cycle"][below[0]]) if below.size else np.nan
    out["li_ion_fraction_final"] = float(m["n_Li_ion"][nc - 1] / max(m["n_Li"][nc - 1], 1))
    lt = t.get("li_total")
    out["li_conserved"] = float(np.nanmax(lt) == np.nanmin(lt)) if lt is not None and np.isfinite(lt).any() else np.nan
    log = os.path.join(run_dir, "kmc_info_log.txt")
    out["stalled"] = 0.0
    if os.path.exists(log):
        with open(log, errors="replace") as fh:
            out["stalled"] = float("stalled" in fh.read())
    # failure cycle: sustained CE_cycle < 0.5, or the run stalling (no strippable
    # Li left, the anode surface fully bound to SEI) at its last cycle
    fc = out["cycle_ce_below_0.5"]
    if not np.isfinite(fc) and out["stalled"] == 1.0:
        fc = float(nc)
    out["cycle_failure"] = fc
    return {"scalars": out, "curves": {k: v[:nc] for k, v in m.items()}}


def xyz_metrics(run_dir: str, geometry_file: str) -> Dict[str, object]:
    """Descriptors of the last frame (and the porosity profile) via kmc3d.morphology."""
    from kmc3d.morphology import frame_descriptors, lattice_site_z, porosity_profile
    from kmc3d.postprocess import load_frame
    files = sorted(glob.glob(os.path.join(run_dir, "trajectory", "kmc-coords-*.xyz")),
                   key=lambda p: int(os.path.basename(p)[11:-4]))
    if not files:
        return {}
    zs = lattice_site_z(geometry_file)
    fr = load_frame(files[-1])
    d = frame_descriptors(fr, zs)
    out = {"xyz_step": d["step"], "xyz_li_buried": d["n_li_buried"], "xyz_islands": d["n_islands"],
           "xyz_frac_buried": d["frac_buried"],
           "xyz_film_extent_A": 0.5 * (d["film_extent_top"] + d["film_extent_bottom"]),
           "xyz_film_porosity": d.get("film_porosity_mean", np.nan),
           "xyz_roughness_A": 0.5 * (d["roughness_top"] + d["roughness_bottom"]),
           "xyz_n_sei": d["n_sei"]}
    for k, v in d.items():
        if k.startswith("buried_coat_"):
            out[f"xyz_{k}"] = v
    prof = porosity_profile(fr, zs)
    # porosity_profile counts every lattice site; a dense Li slab fills only
    # the BC sites (2 of 16 per cube), so rescale to "0 = dense Li metal"
    with np.errstate(invalid="ignore"):
        prof["porosity_rel"] = np.clip(1.0 - (1.0 - prof["porosity"]) / BC_FRACTION, 0.0, 1.0)
    out["xyz_film_porosity_rel"] = float(np.clip(1.0 - (1.0 - out["xyz_film_porosity"]) / BC_FRACTION, 0.0, 1.0)) \
        if np.isfinite(out["xyz_film_porosity"]) else np.nan
    return {"scalars": out, "profile": prof}


# ------------------------------------------------------------- statistics
def sensitivity(per_run: List[Dict[str, float]], metrics: List[str], ref_factor: float = 1.0) -> List[Dict[str, object]]:
    from scipy import stats
    rows = []
    x_all = np.array([r["factor"] for r in per_run], float)
    fmin, fmax = x_all.min(), x_all.max()
    for mname in metrics:
        y_all = np.array([r.get(mname, np.nan) for r in per_run], float)
        ok = np.isfinite(y_all)
        if ok.sum() < 4 or np.unique(x_all[ok]).size < 3:
            continue
        x, y = x_all[ok], y_all[ok]
        lr = stats.linregress(x, y)
        ref = y[np.isclose(x, ref_factor)]
        ref_mean = float(ref.mean()) if ref.size else float(lr.intercept + lr.slope * ref_factor)
        lo, hi = y[np.isclose(x, fmin)], y[np.isclose(x, fmax)]
        if lo.size > 1 and hi.size > 1:
            tt = stats.ttest_ind(lo, hi, equal_var=False)
            p_t = float(tt.pvalue)
        else:
            p_t = np.nan
        sp = stats.spearmanr(x, y)
        # pooled within-member noise, to compare with the effect size
        groups = [y[np.isclose(x, f)] for f in np.unique(x)]
        noise = float(np.sqrt(np.mean([g.var(ddof=1) for g in groups if g.size > 1]))) if any(g.size > 1 for g in groups) else np.nan
        # power-law exponent on the member means (side/plating is expected ~ factor^-1)
        fs = np.unique(x)
        means = np.array([y[np.isclose(x, f)].mean() for f in fs])
        okp = means > 0
        if okp.sum() >= 3:
            pl = stats.linregress(np.log(fs[okp]), np.log(means[okp]))
            exponent, exp_err = float(pl.slope), float(pl.stderr)
        else:
            exponent, exp_err = np.nan, np.nan
        pct50 = 100.0 * lr.slope * 0.5 / ref_mean if ref_mean else np.nan
        span = 100.0 * lr.slope * (fmax - fmin) / ref_mean if ref_mean else np.nan
        rows.append({"metric": mname, "ref_mean": ref_mean, "slope_per_unit_factor": float(lr.slope),
                     "pct_change_per_plus50": pct50, "pct_change_over_sweep": span,
                     "power_exponent": exponent, "power_exponent_se": exp_err,
                     "r": float(lr.rvalue), "p_linear": float(lr.pvalue), "spearman_rho": float(sp.correlation),
                     "p_spearman": float(sp.pvalue), "p_welch_low_vs_high": p_t,
                     "seed_noise_sd": noise, "effect_to_noise": (abs(lr.slope) * (fmax - fmin) / noise) if noise else np.nan,
                     "significant_p05": bool(lr.pvalue < 0.05)})
    rows.sort(key=lambda r: -abs(r["pct_change_per_plus50"]) if np.isfinite(r["pct_change_per_plus50"]) else 0)
    return rows


# ---------------------------------------------------------------- figures
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13,
                         "legend.fontsize": 11, "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "figure.facecolor": "white", "savefig.dpi": 200, "figure.dpi": 100})
    return plt


def _style(ax, title, ylabel, xlabel):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, length=3)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", color=INK, pad=10)
    ax.set_ylabel(ylabel, color=INK2)
    ax.set_xlabel(xlabel, color=INK2)


def _phys_axis(ax, manifest):
    """Secondary x axis in physical units when the swept parameter has them."""
    ref = manifest.get("reference_physical", {})
    if "j0_mA_cm2" in ref:
        j0 = ref["j0_mA_cm2"]
        sec = ax.secondary_xaxis("top", functions=(lambda f: f * j0, lambda j: j / j0))
        sec.set_xlabel("j0 (mA/cm2)", color=INK2)
        sec.tick_params(colors=INK2, length=3)
        for s in sec.spines.values():
            s.set_color(GRID)
        return sec
    return None


def _xlabel(manifest):
    p = manifest["param"]
    if p == "bv_k0_site":
        return "factor on the plating exchange rate k0 (1.0 = reference)"
    return f"factor on {p} (1.0 = reference)"


def _errorbar(ax, summ, metric, color=SERIES[0], label=None, scale=1.0):
    f = np.array([s["factor"] for s in summ])
    y = np.array([s.get(f"{metric}_mean", np.nan) for s in summ]) * scale
    e = np.array([s.get(f"{metric}_std", np.nan) for s in summ]) * scale
    ax.errorbar(f, y, yerr=e, fmt="o-", color=color, ecolor=color, elinewidth=1.2, capsize=4,
                markersize=7, linewidth=2, label=label)
    ref = np.isclose(f, 1.0)
    if ref.any():
        ax.plot(f[ref], y[ref], "o", color=INK, markersize=9, zorder=5)


def _seq_colors(n):
    import matplotlib
    cmap = matplotlib.colormaps["Blues"]
    return [cmap(0.35 + 0.6 * i / max(n - 1, 1)) for i in range(n)]


def make_figures(out_dir, manifest, summ, per_run, curves, sens, profiles):
    plt = _mpl()
    fig_dir = os.path.join(out_dir, "figures"); os.makedirs(fig_dir, exist_ok=True)
    xl = _xlabel(manifest)
    members = [s["name"] for s in summ]
    factors = [s["factor"] for s in summ]
    cols = _seq_colors(len(members))
    made = []

    def save(fig, name):
        p = os.path.join(fig_dir, name)
        fig.tight_layout()
        fig.savefig(p, bbox_inches="tight"); plt.close(fig); made.append(name)

    # F01 design: swept value in physical units and the resulting current density
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    _errorbar(ax, summ, "j_plate_mA_cm2", SERIES[0], "plating (charge half)")
    _errorbar(ax, summ, "j_strip_mA_cm2", SERIES[1], "stripping (discharge half)")
    _style(ax, "Current density of the box at each setting", "current density (mA/cm2)", xl)
    _phys_axis(ax, manifest); ax.legend(frameon=False)
    save(fig, "F01_current_density.png")

    # F02 side-reaction events per 1000 plating events, by channel (stacked)
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    chans = [k[4:] for k in SIDE]
    bottom = np.zeros(len(summ))
    for i, c in enumerate(chans):
        v = np.array([s.get(f"{c}_per_1000_plated_mean", 0.0) for s in summ])
        v = np.nan_to_num(v)
        if v.sum() == 0:
            continue
        ax.bar(range(len(summ)), v, bottom=bottom, color=SERIES[i % len(SERIES)], label=c, width=0.7,
               edgecolor="white", linewidth=1.5)
        bottom += v
    tot_e = np.array([s.get("side_per_1000_plated_std", 0.0) for s in summ])
    ax.errorbar(range(len(summ)), bottom, yerr=tot_e, fmt="none", ecolor=INK, capsize=4)
    # expectation if side reactions simply compete with plating: ratio ~ 1/factor
    ref = [s for s in summ if np.isclose(s["factor"], 1.0)]
    if ref and np.isfinite(ref[0].get("side_per_1000_plated_mean", np.nan)):
        exp_v = ref[0]["side_per_1000_plated_mean"] / np.array(factors)
        ax.plot(range(len(summ)), exp_v, "D--", color=INK, markersize=6, linewidth=1.2,
                label="expected (reference / factor)")
    ax.set_xticks(range(len(summ))); ax.set_xticklabels([f"{f:g}" for f in factors])
    _style(ax, "Side reactions per 1000 plating events", "events per 1000 plating events", xl)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    save(fig, "F02_side_reactions.png")

    # F03 SEI growth curves + rate
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0), gridspec_kw={"width_ratios": [1.5, 1]})
    for m, f, c in zip(members, factors, cols):
        cv = curves.get(m, {}).get("n_SEI")
        if cv is None:
            continue
        a1.plot(cv["cycle"], cv["mean"], color=(INK if np.isclose(f, 1.0) else c),
                linewidth=(2.6 if np.isclose(f, 1.0) else 1.8), label=f"{f:g}")
    _style(a1, "SEI sites vs cycle (mean over seeds)", "SEI sites in the box", "cycle")
    a1.legend(frameon=False, title="factor", ncol=3, fontsize=9)
    _errorbar(a2, summ, "sei_rate", SERIES[0])
    _style(a2, "SEI growth rate", "SEI sites per cycle", xl.split(" (")[0])
    _phys_axis(a2, manifest)
    save(fig, "F03_sei_growth.png")

    # F04 SEI composition (inorganic vs organic), final
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    keys = [k for k in INORG + ORG if any(f"{k}_final_mean" in s for s in summ)]
    bottom = np.zeros(len(summ))
    for i, k in enumerate(keys):
        v = np.nan_to_num(np.array([s.get(f"{k}_final_mean", 0.0) for s in summ]))
        if v.sum() == 0:
            continue
        ax.bar(range(len(summ)), v, bottom=bottom, color=SERIES[i % len(SERIES)], label=k[4:], width=0.7,
               edgecolor="white", linewidth=1.5)
        bottom += v
    ax.set_xticks(range(len(summ))); ax.set_xticklabels([f"{f:g}" for f in factors])
    _style(ax, "SEI composition at the last cycle", "SEI sites (F, O, N, S inorganic; SFO, F5D organic)", xl)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    save(fig, "F04_sei_composition.png")

    # F05 SEI thickness and roughness
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0))
    _errorbar(a1, summ, "sei_thickness_mean_A", SERIES[0])
    _style(a1, "SEI z extent (window mean)", "angstrom", xl.split(" (")[0]); _phys_axis(a1, manifest)
    _errorbar(a2, summ, "roughness_mean_A", SERIES[1])
    _style(a2, "Li front roughness (window mean)", "std of column height (angstrom)", xl.split(" (")[0]); _phys_axis(a2, manifest)
    save(fig, "F05_thickness_roughness.png")

    # F06 Li inventory: Li bound to SEI vs cycle, and loss per 1000 plated
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0), gridspec_kw={"width_ratios": [1.5, 1]})
    for m, f, c in zip(members, factors, cols):
        cv = curves.get(m, {}).get("n_Li_ion")
        if cv is None:
            continue
        a1.plot(cv["cycle"], cv["mean"], color=(INK if np.isclose(f, 1.0) else c),
                linewidth=(2.6 if np.isclose(f, 1.0) else 1.8), label=f"{f:g}")
    _style(a1, "Li bound to SEI (ionised, cannot strip)", "Li sites", "cycle")
    a1.legend(frameon=False, title="factor", ncol=3, fontsize=9)
    _errorbar(a2, summ, "li_lost_per_1000_plated", SERIES[0])
    _style(a2, "Li lost per 1000 Li plated", "Li bound per 1000 plated", xl.split(" (")[0]); _phys_axis(a2, manifest)
    save(fig, "F06_li_inventory.png")

    # F07 coulombic efficiency implied by the electron budget
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0))
    _errorbar(a1, summ, "ce_event", SERIES[0], None, 100.0)
    _style(a1, "CE = plating / (plating + side reactions)", "CE (%)", xl.split(" (")[0]); _phys_axis(a1, manifest)
    _errorbar(a2, summ, "loss_event_permille", SERIES[1], "measured")
    if ref and np.isfinite(ref[0].get("loss_event_permille_mean", np.nan)):
        fx = np.linspace(min(factors), max(factors), 50)
        a2.plot(fx, ref[0]["loss_event_permille_mean"] / fx, "--", color=INK, linewidth=1.2,
                label="expected (reference / factor)")
        a2.legend(frameon=False)
    _style(a2, "Charge lost to side reactions", "(1 - CE) x 1000", xl.split(" (")[0]); _phys_axis(a2, manifest)
    save(fig, "F07_coulombic_efficiency.png")

    # F08 morphology from the last frame
    if any("xyz_li_buried_mean" in s for s in summ):
        fig, axs = plt.subplots(1, 3, figsize=(13.3, 4.6))
        _errorbar(axs[0], summ, "xyz_li_buried", SERIES[0])
        _style(axs[0], "Buried Li (isolated islands)", "Li sites", "factor"); _phys_axis(axs[0], manifest)
        _errorbar(axs[1], summ, "xyz_film_extent_A", SERIES[1])
        _style(axs[1], "Film extent above the slab", "angstrom", "factor"); _phys_axis(axs[1], manifest)
        _errorbar(axs[2], summ, "xyz_film_porosity_rel", SERIES[2])
        _style(axs[2], "Film porosity (0 = dense Li)", "porosity", "factor"); _phys_axis(axs[2], manifest)
        save(fig, "F08_morphology_last_frame.png")

    # F09 porosity profiles, low / reference / high
    if profiles:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        pick = [members[0]] + ([m for m, f in zip(members, factors) if np.isclose(f, 1.0)]) + [members[-1]]
        lab = {members[0]: f"lowest ({factors[0]:g})", members[-1]: f"highest ({factors[-1]:g})"}
        for m in pick:
            if m not in profiles:
                continue
            z, por = profiles[m]
            c = INK if lab.get(m) is None else (SERIES[0] if m == members[0] else SERIES[1])
            ax.plot(z, por, color=c, linewidth=2.2, label=lab.get(m, "reference (1.0)"))
        _style(ax, "Porosity profile across the slab (last frame, mean over seeds)",
               "porosity (0 = dense Li metal, 1 = electrolyte)", "z from the anode centre (angstrom)")
        ax.set_xlim(-45, 45); ax.legend(frameon=False, loc="lower left")
        save(fig, "F09_porosity_profile.png")

    # F10 time per half-cycle (what the clock says)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    _errorbar(ax, summ, "dt_charge_s", SERIES[0], "charge half")
    _errorbar(ax, summ, "dt_discharge_s", SERIES[1], "discharge half")
    _style(ax, "kMC time per half-cycle", "seconds", xl); _phys_axis(ax, manifest)
    ax.legend(frameon=False)
    save(fig, "F10_time_per_half.png")

    # F11 sensitivity ranking
    if sens:
        rows = [r for r in sens if np.isfinite(r["pct_change_per_plus50"])][:14]
        fig, ax = plt.subplots(figsize=(9, 6))
        names = [LABELS.get(r["metric"], r["metric"]) for r in rows][::-1]
        vals = [r["pct_change_per_plus50"] for r in rows][::-1]
        sig = [r["significant_p05"] for r in rows][::-1]
        colors = [SERIES[0] if s else GRID for s in sig]
        ax.barh(names, vals, color=colors, edgecolor="white", height=0.7)
        vmax = max(1.0, float(np.nanmax(np.abs(vals))))
        ax.set_xlim(-1.35 * vmax, 1.35 * vmax)
        for i, (v, r) in enumerate(zip(vals, rows[::-1])):
            ax.annotate(f"p={r['p_linear']:.2g}", (v, i), xytext=(5 if v >= 0 else -5, 0),
                        textcoords="offset points", va="center",
                        ha="left" if v >= 0 else "right", fontsize=9, color=INK2)
        ax.axvline(0, color=INK2, linewidth=1)
        _style(ax, "Sensitivity: % change per +50 % in the parameter (blue = p < 0.05)", "", "% change per +50 %")
        ax.grid(True, axis="x", color=GRID); ax.grid(False, axis="y")
        save(fig, "F11_sensitivity.png")
    # F12 aging: CE per cycle (stripped / plated) and cycles to failure
    any_fail = any(np.isfinite(r.get("cycle_failure", np.nan)) or np.isfinite(r.get("cycle_ce_below_0.9", np.nan))
                   for r in per_run)
    if any(("CE_cycle_calc" in curves.get(m, {})) for m in members):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0), gridspec_kw={"width_ratios": [1.5, 1]})
        for m, f, c in zip(members, factors, cols):
            cv = curves.get(m, {}).get("CE_cycle_calc")
            if cv is None:
                continue
            a1.plot(cv["cycle"], cv["mean"], color=(INK if np.isclose(f, 1.0) else c),
                    linewidth=(2.6 if np.isclose(f, 1.0) else 1.8), label=f"{f:g}")
        _style(a1, "CE per cycle = stripped / plated (mean over seeds)", "CE_cycle", "cycle")
        a1.legend(frameon=False, title="factor", ncol=3, fontsize=9)
        if any_fail:
            _errorbar(a2, summ, "cycle_ce_below_0.9", SERIES[0], "CE_cycle < 0.9")
            _errorbar(a2, summ, "cycle_failure", SERIES[1], "failure (CE < 0.5 or stall)")
            _style(a2, "Cycles to failure (failed runs only)", "cycle", xl.split(" (")[0]); _phys_axis(a2, manifest)
            a2.legend(frameon=False)
        else:
            _errorbar(a2, summ, "li_ion_fraction_final", SERIES[0])
            _style(a2, "Anode Li bound to SEI at the end", "fraction of n_Li", xl.split(" (")[0]); _phys_axis(a2, manifest)
        save(fig, "F12_aging.png")
    return made


# ------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep", help="sweep directory with sweep.json and f_*/runs/seed_*")
    ap.add_argument("--skip", type=int, default=10, help="first cycle of the analysis window")
    ap.add_argument("--max-cycle", type=int, default=None)
    ap.add_argument("--no-xyz", action="store_true", help="skip trajectory-frame descriptors")
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args(argv)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    with open(os.path.join(a.sweep, "sweep.json")) as fh:
        manifest = json.load(fh)
    out_dir = os.path.join(a.sweep, "analysis"); os.makedirs(out_dir, exist_ok=True)
    members = sorted(manifest["members"], key=lambda m: m["factor"])

    per_run: List[Dict[str, float]] = []
    curves: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    profiles: Dict[str, tuple] = {}
    curve_keys = ("n_SEI", "n_Li_ion", "n_Li_dead", "sei_thickness", "surface_roughness", "CE_mod", "cum_side", "j_plate",
                  "CE_cycle_calc", "stripped")
    for m in members:
        mdir = os.path.join(a.sweep, m["name"])
        geo = read_geometry(os.path.join(mdir, "GEOMETRY.in"))
        a0 = float(geo.get("latticeConstant", 4.0))
        area = 2.0 * (geo["nx"] * a0) * (geo["ny"] * a0) * 1e-16     # both faces, cm2
        runs = sorted(glob.glob(os.path.join(mdir, "runs", "seed_*")))
        seed_curves: Dict[str, List[np.ndarray]] = {k: [] for k in curve_keys}
        prof_acc = []
        for r in runs:
            res = run_metrics(r, area, a.skip, a.max_cycle)
            if res is None:
                print(f"  {m['name']} {os.path.basename(r)}: no ledger, skipped")
                continue
            row = {"member": m["name"], "factor": m["factor"], "value": m["value"],
                   "seed": os.path.basename(r)}
            row.update(m.get("physical", {}) if isinstance(m.get("physical"), dict) else {})
            row.pop("note", None)
            row.update(res["scalars"])
            if not a.no_xyz:
                try:
                    xm = xyz_metrics(r, os.path.join(mdir, "GEOMETRY.in"))
                    if xm:
                        row.update(xm["scalars"])
                        prof_acc.append(xm["profile"])
                except Exception as exc:          # noqa: BLE001
                    print(f"  {m['name']} {os.path.basename(r)}: xyz descriptors failed ({exc})")
            per_run.append(row)
            for k in curve_keys:
                if k in res["curves"]:
                    seed_curves[k].append(res["curves"][k])
        curves[m["name"]] = {}
        for k, lst in seed_curves.items():
            if not lst:
                continue
            n = min(len(v) for v in lst)
            st = np.vstack([v[:n] for v in lst]).astype(float)
            curves[m["name"]][k] = {"cycle": np.arange(1, n + 1), "mean": np.nanmean(st, 0), "std": np.nanstd(st, 0)}
        if prof_acc:
            z = prof_acc[0]["z_rel"]
            por = np.nanmean(np.vstack([p["porosity_rel"] for p in prof_acc]), 0)
            profiles[m["name"]] = (z, por)
        print(f"{m['name']}: {sum(1 for r in per_run if r['member'] == m['name'])} runs")

    if not per_run:
        sys.exit("no runs found")

    # per_run.csv
    keys = []
    for r in per_run:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(out_dir, "per_run.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in per_run:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) and np.isfinite(v) else v) for k, v in r.items()})

    # summary.csv: mean / std / n per member
    num_keys = [k for k in keys if k not in ("member", "seed") and any(isinstance(r.get(k), (int, float, np.floating)) for r in per_run)]
    summ = []
    for m in members:
        rows = [r for r in per_run if r["member"] == m["name"]]
        if not rows:
            continue
        s = {"name": m["name"], "factor": m["factor"], "value": m["value"], "n_seeds": len(rows)}
        for k in num_keys:
            v = np.array([r.get(k, np.nan) for r in rows], float)
            s[f"{k}_mean"] = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
            s[f"{k}_std"] = float(np.nanstd(v, ddof=1)) if np.isfinite(v).sum() > 1 else np.nan
        summ.append(s)
    skeys = list(summ[0])
    with open(os.path.join(out_dir, "summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=skeys); w.writeheader()
        for s in summ:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) and np.isfinite(v) else v) for k, v in s.items()})

    # sensitivity
    metrics = ["j_plate_mA_cm2", "j_strip_mA_cm2", "dt_charge_s", "side_per_1000_plated",
               "FSI_per_1000_plated", "SOL_per_1000_plated", "SFO_per_1000_plated", "SOL2_per_1000_plated",
               "sei_rate", "sei_final", "sei_thickness_mean_A", "sei_inorganic_fraction",
               "li_lost_per_1000_plated", "li_ion_final", "li_dead_final", "roughness_mean_A",
               "surface_diff_per_cycle", "ce_event", "loss_event_permille",
               "cycle_ce_below_0.9", "cycle_failure", "stripped_last10_over_first10", "li_ion_fraction_final",
               "xyz_li_buried", "xyz_islands", "xyz_film_extent_A", "xyz_film_porosity_rel"]
    metrics = [k for k in metrics if k in num_keys]
    sens = sensitivity(per_run, metrics)
    if sens:
        with open(os.path.join(out_dir, "sensitivity.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sens[0])); w.writeheader()
            for r in sens:
                w.writerow({k: (f"{v:.4g}" if isinstance(v, float) and np.isfinite(v) else v) for k, v in r.items()})

    # curves
    for k in curve_keys:
        names = [m["name"] for m in members if k in curves.get(m["name"], {})]
        if not names:
            continue
        n = min(curves[nm][k]["cycle"].size for nm in names)
        with open(os.path.join(out_dir, f"curves_{k}.csv"), "w", newline="") as fh:
            fh.write("cycle," + ",".join(f"{nm}_mean,{nm}_std" for nm in names) + "\n")
            for i in range(n):
                fh.write(f"{i + 1}," + ",".join(f"{curves[nm][k]['mean'][i]:.6g},{curves[nm][k]['std'][i]:.6g}" for nm in names) + "\n")
    if profiles:
        names = list(profiles)
        with open(os.path.join(out_dir, "porosity_profiles.csv"), "w", newline="") as fh:
            fh.write("z_rel_A," + ",".join(names) + "\n")
            z = profiles[names[0]][0]
            for i in range(z.size):
                fh.write(f"{z[i]:.3f}," + ",".join(f"{profiles[nm][1][i]:.4f}" for nm in names) + "\n")

    made = []
    if not a.no_figures:
        try:
            made = make_figures(out_dir, manifest, summ, per_run, curves, sens, profiles)
        except ImportError as exc:
            print(f"figures skipped: {exc}")

    # README with the headline numbers
    with open(os.path.join(out_dir, "README.md"), "w") as fh:
        fh.write(f"# Sweep analysis: {manifest['param']} around {manifest['reference_value']:g}\n\n")
        fh.write(f"Base case `{manifest['base']}`; window cycles {int(per_run[0]['window_from'])} to "
                 f"{int(min(r['window_to'] for r in per_run))}; {len(per_run)} runs.\n\n")
        rp = manifest.get("reference_physical", {})
        if "j0_mA_cm2" in rp:
            fh.write(f"Reference: k0 = {rp['k0_site_per_s']:g} s^-1 per site = j0 {rp['j0_mA_cm2']:.1f} mA/cm2; "
                     f"plating rate per growth site {rp['k_plate_per_site_per_s']:.0f} s^-1 at eta_charge.\n\n")
        fh.write("## Summary (mean +- sd over seeds)\n\n| factor | value | n | j_plate mA/cm2 | side/1000 | SEI/cycle | SEI z (A) | Li lost/1000 | roughness (A) | 1-CE (permille) |\n|---|---|---|---|---|---|---|---|---|---|\n")
        for s in summ:
            def c(k, d=2):
                return f"{s.get(k + '_mean', np.nan):.{d}f} +- {s.get(k + '_std', np.nan):.{d}f}"
            fh.write(f"| {s['factor']:g} | {s['value']:.4g} | {s['n_seeds']} | {c('j_plate_mA_cm2')} | {c('side_per_1000_plated')} | "
                     f"{c('sei_rate', 3)} | {c('sei_thickness_mean_A', 1)} | {c('li_lost_per_1000_plated')} | {c('roughness_mean_A')} | {c('loss_event_permille')} |\n")
        fh.write("\n## Sensitivity (linear fit over all runs; % change per +50 % of the parameter; "
                 "power-law exponent from the member means, metric ~ factor^n)\n\n"
                 "| metric | reference | % per +50 % | exponent n | p (linear) | p (Welch low vs high) | effect / seed noise |\n|---|---|---|---|---|---|---|\n")
        for r in sens:
            fh.write(f"| {r['metric']} | {r['ref_mean']:.4g} | {r['pct_change_per_plus50']:+.1f} | "
                     f"{r['power_exponent']:+.2f} +- {r['power_exponent_se']:.2f} | {r['p_linear']:.3g} | "
                     f"{r['p_welch_low_vs_high']:.3g} | {r['effect_to_noise']:.2f} |\n")
        fh.write("\n## Figures\n\n" + "\n".join(f"- figures/{m}" for m in made) + "\n")
        bad = [r for r in per_run if r.get("li_conserved") == 0.0 or r.get("stalled") == 1.0]
        fh.write(f"\nConservation or stall flags: {len(bad)} run(s)" + (": " + ", ".join(f"{r['member']}/{r['seed']}" for r in bad) if bad else "") + "\n")
    print(f"tables and figures -> {out_dir}")


if __name__ == "__main__":
    main()
