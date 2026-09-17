"""
kmc3d
=====
Generalized 3D kinetic Monte Carlo (kMC) for electrochemical cells: SEI
formation and Li electrodeposition (Li-metal half cell, Perez-Beltran et al.
2024) extended to a Li-S full cell with polysulfide shuttle, shared Li+ pool
and cathode passivation. Five input files: PARAMETERS.in, GEOMETRY.in,
DECOMPOSITION.in, MOBILITY.in, MECHANISM.in.

Public surface
--------------
    config       PARAMETERS.in / GEOMETRY.in / DECOMPOSITION.in / MOBILITY.in parsers
    stats        per-half-cycle ledger (cycle_stats.csv, CE metrics, conservation)
    mechanism    MECHANISM.in parser (reaction products -> data driven)
    lattice      OVITO-free superlattice generation + POSCAR I/O
    engine       the kMC engine (BKL; reproduces trajectory + Data2Excel.txt)
    output       writers for the exact original output formats
    postprocess  hooks for z-warp / density-field / morphology / ML pipelines
    run          CLI entry point (used by goKMC_grace.slrm)
"""

from . import config, mechanism, lattice, engine, output, postprocess  # noqa: F401

__all__ = ["config", "mechanism", "lattice", "engine", "output", "postprocess"]
__version__ = "1.1.0"
