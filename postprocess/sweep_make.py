#!/usr/bin/env python3
"""
sweep_make.py - build a one-parameter sweep (a family of case directories)
from a base case, with the physical meaning of the swept parameter recorded
in a manifest, and write the launch script for GRACE.

Usage (repository root, kmc3d environment active):

    # plating exchange rate, +-50 % around the physical reference (anode_physical)
    python postprocess/sweep_make.py --base cases/anode_physical \
        --out cases/sweep_plating_j0 --param bv_k0_site \
        --factors 0.5,0.625,0.75,0.875,1.0,1.125,1.25,1.375,1.5 \
        --set maxCycles=280 --seeds 10

    # the same idea on the legacy network (DECOMPOSITION.in k0 of one reaction)
    python postprocess/sweep_make.py --base validation/anode_small \
        --out cases/sweep_plating_legacy --param decomp:Plating --factors 0.5,1,1.5

--param names either a PARAMETERS.in keyword (bv_k0_site, eta_charge, ...)
or decomp:<REACTION> for the k0 column of that reaction in DECOMPOSITION.in.
Every factor gets its own directory <out>/f_<factor>/ holding the five *.in
files of the base case with only the swept value (and any --set overrides)
changed, so goKMC_prod_array.slrm runs each one unchanged. <out>/sweep.json
records the design; sweep_analyze.py reads it.

Nothing in the engine is touched: a sweep is only a set of input decks.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys

E_CHARGE = 1.602176634e-19      # C
K_B_EV = 8.617333262e-5         # eV/K
A_SITE_CM2 = 4.0e-16            # (2.0 A)^2, the surface area of one lattice site
IN_FILES = ("PARAMETERS.in", "GEOMETRY.in", "DECOMPOSITION.in", "MOBILITY.in", "MECHANISM.in")


def read_param(path: str, key: str) -> str | None:
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            if tok[0] == key:
                return " ".join(tok[1:])
    return None


def set_param(path: str, key: str, value: str) -> bool:
    """Replace the value of `key` in a PARAMETERS.in-style file (append if absent)."""
    with open(path) as fh:
        lines = fh.readlines()
    done = False
    for i, line in enumerate(lines):
        s = line.strip()
        if s and not s.startswith("#") and s.split()[0] == key:
            lines[i] = f"{key} {value}\n"
            done = True
    if not done:
        lines.append(f"{key} {value}\n")
    with open(path, "w") as fh:
        fh.writelines(lines)
    return done


def decomp_k0(path: str, reaction: str) -> float:
    """k0 of the first channel of `reaction` in DECOMPOSITION.in."""
    with open(path) as fh:
        lines = [l for l in fh]
    for i, line in enumerate(lines):
        s = line.strip()
        if s and not s.startswith("#") and s.split()[0] == reaction:
            for j in range(i + 1, len(lines)):
                t = lines[j].strip()
                if t and not t.startswith("#"):
                    return float(t.split()[1])
    raise KeyError(f"{reaction} not found in {path}")


def scale_decomp_k0(path: str, reaction: str, factor: float) -> float:
    with open(path) as fh:
        lines = [l for l in fh]
    for i, line in enumerate(lines):
        s = line.strip()
        if s and not s.startswith("#") and s.split()[0] == reaction:
            for j in range(i + 1, len(lines)):
                t = lines[j].strip()
                if t and not t.startswith("#"):
                    tok = t.split()
                    new = float(tok[1]) * factor
                    tok[1] = f"{new:.6g}"
                    lines[j] = " ".join(tok) + "\n"
                    with open(path, "w") as fh:
                        fh.writelines(lines)
                    return new
    raise KeyError(f"{reaction} not found in {path}")


def physical_mapping(param: str, value: float, params_path: str) -> dict:
    """Physical units of the swept value, when the parameter has them."""
    T = float(read_param(params_path, "temperature") or 298.0)
    kT = K_B_EV * T
    if param == "bv_k0_site":
        alpha = float(read_param(params_path, "bv_alpha") or 0.5)
        eta_c = float(read_param(params_path, "eta_charge") or -0.03)
        eta_d = float(read_param(params_path, "eta_discharge") or 0.03)
        j0 = value * E_CHARGE / A_SITE_CM2 * 1e3          # mA/cm2
        return {"k0_site_per_s": value, "j0_mA_cm2": j0,
                "k_plate_per_site_per_s": value * __import__("math").exp(-alpha * eta_c / kT),
                "k_strip_per_site_per_s": value * __import__("math").exp(alpha * eta_d / kT),
                "note": "j0 = k0_site e / A_site with A_site = (2.0 A)^2; plating "
                        "rate per growth site at eta_charge, stripping at eta_discharge"}
    if param.startswith("decomp:"):
        return {"k0_model_units": value,
                "note": "legacy network: k0 multiplies globalA and sigma, model units"}
    return {"value": value}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="base case directory (five *.in files)")
    ap.add_argument("--out", required=True, help="sweep directory to create")
    ap.add_argument("--param", required=True, help="PARAMETERS.in key or decomp:<REACTION>")
    ap.add_argument("--factors", required=True, help="comma-separated multipliers of the base value")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="extra PARAMETERS.in overrides applied to every member (repeatable)")
    ap.add_argument("--seeds", type=int, default=5, help="seeds per member (array size in launch.sh)")
    ap.add_argument("--base-seed", type=int, default=8597)
    ap.add_argument("--time", default="04:00:00", help="SLURM --time per seed")
    ap.add_argument("--force", action="store_true", help="overwrite an existing sweep directory")
    a = ap.parse_args(argv)

    base = a.base.rstrip("/\\")
    for f in IN_FILES:
        if not os.path.exists(os.path.join(base, f)):
            sys.exit(f"missing {f} in {base}")
    if os.path.exists(a.out):
        if not a.force:
            sys.exit(f"{a.out} exists (use --force to rebuild the input decks; runs/ are kept)")
    os.makedirs(a.out, exist_ok=True)

    factors = [float(x) for x in a.factors.split(",")]
    params_path = os.path.join(base, "PARAMETERS.in")
    if a.param.startswith("decomp:"):
        ref = decomp_k0(os.path.join(base, "DECOMPOSITION.in"), a.param.split(":", 1)[1])
    else:
        v = read_param(params_path, a.param)
        if v is None:
            sys.exit(f"{a.param} not found in {params_path}")
        ref = float(v.split()[0])

    members = []
    for f in factors:
        name = f"f_{f:.3f}"
        d = os.path.join(a.out, name)
        os.makedirs(d, exist_ok=True)
        for fn in IN_FILES:
            shutil.copy(os.path.join(base, fn), os.path.join(d, fn))
        if a.param.startswith("decomp:"):
            new = scale_decomp_k0(os.path.join(d, "DECOMPOSITION.in"), a.param.split(":", 1)[1], f)
        else:
            new = ref * f
            set_param(os.path.join(d, "PARAMETERS.in"), a.param, f"{new:.6g}")
        for kv in a.set:
            k, v = kv.split("=", 1)
            set_param(os.path.join(d, "PARAMETERS.in"), k, v)
        members.append({"name": name, "factor": f, "value": new,
                        "physical": physical_mapping(a.param, new, os.path.join(d, "PARAMETERS.in"))})
        print(f"{name}: {a.param} = {new:.6g}")

    manifest = {"base": base, "param": a.param, "reference_value": ref,
                "reference_physical": physical_mapping(a.param, ref, params_path),
                "overrides": a.set, "seeds": a.seeds, "base_seed": a.base_seed,
                "members": members}
    with open(os.path.join(a.out, "sweep.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    # launch script: one SLURM array per member, then the analysis job after all
    sh = os.path.join(a.out, "launch.sh")
    out_rel = a.out.rstrip("/\\")
    with open(sh, "w", newline="\n") as fh:
        fh.write("#!/bin/bash\n# generated by sweep_make.py; run from the repository root on GRACE\n")
        fh.write("set -e\nmkdir -p slurm_logs\nJOBS=\"\"\n")
        for m in members:
            fh.write(f"J=$(sbatch --parsable --job-name=sw_{m['factor']:.3f} --time={a.time} "
                     f"--array=0-{a.seeds - 1} --export=ALL,CASE={out_rel}/{m['name']},BASE_SEED={a.base_seed} "
                     f"slurm/goKMC_prod_array.slrm)\n")
            fh.write("echo \"submitted $J for %s\"\nJOBS=\"$JOBS:$J\"\n" % m['name'])
        fh.write(f"sbatch --job-name=sw_analyze --dependency=afterany${{JOBS}} "
                 f"--export=ALL,SWEEP={out_rel} slurm/sweep_analyze.slrm\n")
        fh.write("echo \"analysis job queued after${JOBS}\"\n")
    os.chmod(sh, 0o755)
    print(f"manifest -> {os.path.join(a.out, 'sweep.json')}\nlaunch   -> {sh}")


if __name__ == "__main__":
    main()
