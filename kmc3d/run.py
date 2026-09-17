"""
run.py
======
Command-line driver.  Reads the five input files, builds (or loads) the
superlattice, runs the kMC engine, and writes the trajectory + Data2Excel.txt.

Usage
-----
    python -m kmc3d.run [--dir INPUT_DIR] [--out OUTPUT_DIR]
                        [--params PARAMETERS.in] [--geom GEOMETRY.in]
                        [--decomp DECOMPOSITION.in] [--mob MOBILITY.in]
                        [--mech MECHANISM.in]

By default every file is looked up in --dir (default: current directory).
This is what goKMC_grace.slrm invokes on GRACE.
"""

from __future__ import annotations
import argparse
import os
import sys
import time

from . import config
from . import lattice as L
from .mechanism import Mechanism
from .engine import Engine
from .output import Output


def _resolve(d, name, override):
    return override if override else os.path.join(d, name)


def build_engine(input_dir=".", out_dir="output",
                 params_file=None, geom_file=None, decomp_file=None,
                 mob_file=None, mech_file=None):
    """Assemble an Engine from the input files (also usable from notebooks)."""
    p_path = _resolve(input_dir, "PARAMETERS.in", params_file)
    g_path = _resolve(input_dir, "GEOMETRY.in", geom_file)
    d_path = _resolve(input_dir, "DECOMPOSITION.in", decomp_file)
    m_path = _resolve(input_dir, "MOBILITY.in", mob_file)
    x_path = _resolve(input_dir, "MECHANISM.in", mech_file)

    params = config.Parameters.read(p_path)
    geom = config.Geometry.read(g_path)
    decomp = config.DecompositionKinetics.read(d_path, params.numberDecompositionRxns)

    # number of mobility reactions is auto-detected from PARAMETERS or counted
    n_mob = params.nMobilityReactions
    if n_mob <= 0:
        n_mob = _count_mobility_blocks(m_path)
    mob = config.MobilityKinetics.read(m_path, n_mob)
    mech = Mechanism.read(x_path)

    # ---- lattice: generate programmatically or read an existing POSCAR -----
    if geom.poscar_in:
        poscar = geom.poscar_in if os.path.isabs(geom.poscar_in) \
            else os.path.join(input_dir, geom.poscar_in)
        lat = L.read_poscar(poscar)
    else:
        lat = L.generate(geom)
        if geom.write_poscar:
            L.write_poscar(lat, os.path.join(out_dir, "POSCAR_generated.vasp"))

    out = Output(root=out_dir)
    out.log(f"sites={lat.n}  box={lat.box.diagonal().tolist()}  "
            f"mobility_reactions={n_mob}")
    eng = Engine(lat, params, geom, decomp, mob, mech, out)
    return eng, out


def _count_mobility_blocks(path):
    """Count leading mobility-reaction blocks before the INTERACTIONS label."""
    toks = config._tokens(path)
    i, n = 0, 0
    while i < len(toks):
        if toks[i].upper().startswith("INTERACTION"):
            break
        # <rxn> <count> then count*5 tokens
        cnt = int(toks[i + 1])
        i += 2 + cnt * 5
        n += 1
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generalized 3D kMC (kmc3d)")
    ap.add_argument("--dir", default=".", help="directory holding the *.in files")
    ap.add_argument("--out", default="output", help="output directory")
    ap.add_argument("--params"); ap.add_argument("--geom")
    ap.add_argument("--decomp"); ap.add_argument("--mob"); ap.add_argument("--mech")
    ap.add_argument("--restart", default="", help="path to checkpoint.npz to resume from")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    eng, out = build_engine(args.dir, args.out, args.params, args.geom,
                            args.decomp, args.mob, args.mech)
    t0 = time.time()
    try:
        eng.run(restart=args.restart)
    finally:
        out.log(f"wall_time_s={time.time() - t0:.2f}")
        out.close()
    print(f"Done. steps={eng.currentStep} cycles={eng.cycleNumber} "
          f"-> {args.out}/  ({time.time() - t0:.1f}s)")


if __name__ == "__main__":
    sys.exit(main())
