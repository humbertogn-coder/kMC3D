"""
structio.py - read kMC snapshot POSCARs and export structures for the
characterisation tools (Zeo++, RASPA3).

A kMC snapshot (written by Engine._write_snapshot) is a VASP/POSCAR file
containing only the frozen solid skeleton (SEI + metal), electrolyte excluded.
Both Zeo++ and RASPA can consume a CIF, so CIF is the single export format.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List
import numpy as np


@dataclass
class Structure:
    box: np.ndarray            # 3x3 cell (Angstrom, rows = lattice vectors)
    species: np.ndarray        # (N,) element/species symbols
    cart: np.ndarray           # (N,3) cartesian coordinates (Angstrom)

    @property
    def frac(self) -> np.ndarray:
        return self.cart @ np.linalg.inv(self.box)

    @property
    def lengths(self):
        a, b, c = (np.linalg.norm(self.box[i]) for i in range(3))
        return a, b, c

    @property
    def angles(self):
        # alpha (b,c), beta (a,c), gamma (a,b) in degrees
        def ang(u, v):
            cs = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
            return float(np.degrees(np.arccos(np.clip(cs, -1, 1))))
        a, b, c = self.box
        return ang(b, c), ang(a, c), ang(a, b)


def read_poscar(path: str) -> Structure:
    """Read a (Cartesian or Direct) POSCAR/VASP snapshot."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    scale = float(lines[1].split()[0])
    box = np.array([[float(x) for x in lines[i].split()[:3]]
                    for i in (2, 3, 4)]) * scale
    syms = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    mode = lines[7].strip().lower()
    cartesian = mode.startswith("c") or mode.startswith("k")
    species: List[str] = []
    for s, n in zip(syms, counts):
        species += [s] * n
    n_atoms = sum(counts)
    coords = np.array([[float(x) for x in lines[8 + i].split()[:3]]
                       for i in range(n_atoms)])
    cart = coords if cartesian else coords @ box
    return Structure(box=box, species=np.array(species), cart=cart)


def write_cif(struct: Structure, path: str, name: str = "kmc_snapshot") -> str:
    """Write a P1 CIF (explicit atoms, no symmetry) readable by Zeo++ and RASPA."""
    a, b, c = struct.lengths
    al, be, ga = struct.angles
    frac = struct.frac % 1.0
    with open(path, "w") as fh:
        fh.write(f"data_{name}\n")
        fh.write(f"_cell_length_a    {a:.6f}\n")
        fh.write(f"_cell_length_b    {b:.6f}\n")
        fh.write(f"_cell_length_c    {c:.6f}\n")
        fh.write(f"_cell_angle_alpha {al:.4f}\n")
        fh.write(f"_cell_angle_beta  {be:.4f}\n")
        fh.write(f"_cell_angle_gamma {ga:.4f}\n")
        fh.write("_symmetry_space_group_name_H-M 'P 1'\n")
        fh.write("_symmetry_Int_Tables_number 1\n")
        fh.write("loop_\n_symmetry_equiv_pos_as_xyz\n  'x, y, z'\n")
        fh.write("loop_\n_atom_site_label\n_atom_site_type_symbol\n"
                 "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
        counter = {}
        for sym, (fx, fy, fz) in zip(struct.species, frac):
            counter[sym] = counter.get(sym, 0) + 1
            label = f"{sym}{counter[sym]}"
            fh.write(f"  {label} {sym} {fx:.6f} {fy:.6f} {fz:.6f}\n")
    return path


def write_cssr(struct: Structure, path: str) -> str:
    """Write a Zeo++ .cssr (alternative to CIF; some Zeo++ builds prefer it)."""
    a, b, c = struct.lengths
    al, be, ga = struct.angles
    frac = struct.frac % 1.0
    with open(path, "w") as fh:
        fh.write(f"{a:.4f} {b:.4f} {c:.4f}\n")
        fh.write(f"{al:.2f} {be:.2f} {ga:.2f} SPGR =  1 P 1\n")
        fh.write(f"{len(struct.species)}   0\n")
        fh.write("kmc snapshot\n")
        for i, (sym, (fx, fy, fz)) in enumerate(zip(struct.species, frac), 1):
            fh.write(f"{i} {sym} {fx:.5f} {fy:.5f} {fz:.5f} "
                     "0 0 0 0 0 0 0 0 0.00\n")
    return path
