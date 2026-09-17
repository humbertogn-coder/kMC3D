#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
config_plating.py
═════════════════
Central configuration for the PLATING k0 sweep analysis pipeline.

This is the SINGLE SOURCE OF TRUTH for paths, k0 mapping, reference
cases, and analysis parameters. Every analysis script imports from here,
so to analyze a NEW parameter sweep in the future you only edit THIS file
(or make a copy like config_stripping.py / config_newparam.py) and
re-run the master orchestrator.

USAGE in any analysis script:
    from config_plating import CFG
    base = CFG["base_dir"]
    for case_key, k0 in CFG["cases"]:
        ...
"""

from pathlib import Path

CFG = {
    # ── Identity of this sweep ──────────────────────────────────────────
    "sweep_name": "plating",
    "varied_event": "Plating",   # which kMC event the coefficient controls
    "fixed_note":  "stripping held at its reference case",

    # ── Paths on HPRC Grace ─────────────────────────────────────────────
    "base_dir":  Path("/scratch/user/humbertogn/Plating_kMC"),
    "out_dir":   Path("/scratch/user/humbertogn/Plating_kMC/Analysis_Results"),
    "log_dir":   Path("/scratch/user/humbertogn/Plating_kMC/logs"),
    "venv_path": Path("/scratch/user/humbertogn/ase_env"),

    # Where the trajectory xyz files live inside each case folder
    "traj_subpath": "output/trajectory",
    # Log file name. The C++ engine writes "Data2Excel.txt" (verified in
    # declarations.cpp: ofstream fscreen2("Data2Excel.txt")). Some plating
    # folders use lowercase "data2excel.txt", so ALL candidates are listed
    # and kmc_io.resolve_log() picks whichever exists in each case folder.
    "log_name": "Data2Excel.txt",
    "log_name_candidates": ["Data2Excel.txt", "data2excel.txt",
                             "DATA2EXCEL.TXT"],

    # ── Module load line for Slurm (matches ase_env) ────────────────────
    "module_load": "GCCcore/13.2.0 Python/3.11.5",

    # ── Case mapping: (folder_name, k_p_value) ──────────────────────────
    # VERY NARROW plating sweep, all clustered around the reference 1e-3.
    # Folders are NOT in value order (folder 1 is the reference; 2-5 step
    # DOWN, 6-10 step UP). Every script sorts by k_p, so this is fine.
    # Sorted low->high: 5,4,3,2,1,6,7,8,9,10.
    "cases": [
        ("1",  1.0e-3),    # <- REFERENCE
        ("2",  0.9e-3),
        ("3",  0.8e-3),
        ("4",  0.7e-3),
        ("5",  0.6e-3),    # lowest
        ("6",  1.1e-3),
        ("7",  1.2e-3),
        ("8",  1.3e-3),
        ("9",  1.4e-3),
        ("10", 1.5e-3),    # highest
    ],

    # ── THREE comparison cases shown together in every comparison figure ─
    # (key, k_p, colour, label). lowest = blue, reference = black,
    # highest = red.
    "comparison_cases": [
        ("5",  0.6e-3, "#1f77b4", "lowest"),
        ("1",  1.0e-3, "#111111", "reference"),
        ("10", 1.5e-3, "#d62728", "highest"),
    ],

    # ── Legacy two-point pair (kept for scripts that still use it) ───────
    "ref_low_key":  "5",        # k_p = 0.6e-3 (lowest)
    "ref_low_k0":   0.6e-3,
    "ref_high_key": "10",       # k_p = 1.5e-3 (highest)
    "ref_high_k0":  1.5e-3,
    "ref_mid_key":  "1",        # k_p = 1.0e-3 (reference)
    "ref_mid_k0":   1.0e-3,

    # ── Nomenclature (PBB request: use k_p, not k0) ─────────────────────
    "rate_symbol_plain": "k_p",
    "rate_symbol_math":  r"$k_{\mathrm{p}}$",
    "rate_units_note":   "plating mobility prefactor (model units, from MOBILITY.in)",

    # ── Regimes for PBB's questions ─────────────────────────────────────
    "kp_reference": 1.0e-3,
    # (Q2) k_p BELOW the reference
    "kp_below_reference": [0.6e-3, 0.7e-3, 0.8e-3, 0.9e-3],
    # (Q3) k_p ABOVE the reference
    "kp_intermediate": [1.1e-3, 1.2e-3, 1.3e-3, 1.4e-3, 1.5e-3],

    # ── Cycle sampling ──────────────────────────────────────────────────
    "max_cycle":  140,
    "cycle_step": 10,     # for Li buried sweep (heavy); 10 cycles
    "cycle_step_reactions": 2,   # finer for reaction/species curves

    # ── Li buried parameters (MAINBODY criterion) ───────────────────────
    "local_density_thr": 0.5,   # user-validated optimum
    "li_buried_cycle_for_db": 110,   # cycle at which dynamic Li buried is
                                      # reported for the ML database
                                      # (matches "cycle 110" in slide Fig 5)

    # ── Dynamic Li buried fitted expression (from STRIPPING sweep) ──────
    # Logistic-style fit used to FILL the dynamic Li buried column in the
    # ML database where direct values are unavailable.
    #   n_li_buried(k0) = L / (1 + exp(-k*(log10(k0) - x0))) + b
    # NOTE: these coefficients come from the stripping analysis. When the
    # stripping DB is provided we will re-fit / confirm them. Placeholders
    # below are flagged and must be confirmed before trusting DB values.
    "dyn_li_buried_fit": {
        "form":  "logistic_log10",   # see fill_dynamic_li_buried() in ML script
        "L":     None,   # amplitude         ← CONFIRM from stripping fit
        "k":     None,   # steepness         ← CONFIRM
        "x0":    None,   # log10(k0) midpoint← CONFIRM
        "b":     None,   # baseline          ← CONFIRM
        "_status": "PLACEHOLDER — confirm coefficients from stripping fit",
    },

    # ── Anode region (for anode-consumption analysis) ───────────────────
    "anode_z_width": 20.2,   # Å (anodeThickness from PARAMETERS.in)

    # ── Species conventions ─────────────────────────────────────────────
    "species_xyz":   ["F", "O", "N", "S", "F5D", "SFO", "Li"],
    "inorg_species": ["F", "O", "N", "S"],
    "amorph_species": ["F5D", "SFO"],

    # ── Plot conventions ────────────────────────────────────────────────
    # Labels NOT on lines (use legend), per user request.
    "plot_style": {
        "font.family":         "serif",
        "font.serif":          ["DejaVu Serif"],
        "font.size":           11,
        "axes.linewidth":      1.1,
        "axes.labelsize":      12,
        "axes.titlesize":      13,
        "xtick.direction":     "in",
        "ytick.direction":     "in",
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.frameon":      True,
        "legend.framealpha":   0.92,
        "savefig.dpi":         300,
    },
    # Color map for the 10 plating cases (sequential, low→high k0)
    "cmap_name": "viridis",
    # Highlight colors for reference cases
    "ref_low_color":  "#2980B9",   # blue  (reference k0=56.25)
    "ref_high_color": "#C0392B",   # red   (highest k0=10000)
}


def get_case_k0(case_key):
    for k, v in CFG["cases"]:
        if k == case_key:
            return v
    raise KeyError(f"case_key {case_key} not in CFG['cases']")


def list_case_keys():
    return [k for k, _ in CFG["cases"]]


def list_k0_values():
    return [v for _, v in CFG["cases"]]


if __name__ == "__main__":
    print(f"Sweep: {CFG['sweep_name']}")
    print(f"Base dir: {CFG['base_dir']}")
    print(f"Cases ({len(CFG['cases'])}):")
    for k, v in CFG["cases"]:
        tag = ""
        if k == CFG["ref_mid_key"]:
            tag = "  <- REFERENCE"
        elif k == CFG["ref_low_key"]:
            tag = "  <- LOWEST"
        elif k == CFG["ref_high_key"]:
            tag = "  <- HIGHEST"
        print(f"  folder {k:>3s}:  k_p = {v:>10.2e}{tag}")
