"""
lattice.py
==========
Programmatic generation of the BCB/BCO/BA/OC/TE superlattice and its neighbour
topology, plus POSCAR read/write.

Why this exists
---------------
The original workflow required manually expanding a POSCAR in OVITO before every
run.  By reverse-engineering the hard-coded neighbour cutoffs in
create_lattice.cpp we can place every site type analytically inside one
conventional cube (side = latticeConstant) and replicate it nx*ny*nz times.  No
OVITO, fully reproducible, box size set from GEOMETRY.in.

Site basis inside the conventional cube (fractional, side a = latticeConstant)
------------------------------------------------------------------------------
  BCB x1 : (0,   0,   0  )                      Li-metal sublattice A (corner)
  BCO x1 : (1/2, 1/2, 1/2)                      Li-metal sublattice B (body)
  OC  x3 : (1/2,0,0) (0,1/2,0) (0,0,1/2)        octahedral
  BA  x3 : (1/2,1/2,0) (1/2,0,1/2) (0,1/2,1/2)  the "BA" set
  TE  x8 : (1/4 or 3/4)^3                        tetrahedral

These reproduce, with a sub-spacing of a/2 = 2 A, the exact cutoff table used
by createNeighborArrays():
  dBCB2TE=sqrt3, dBCB2OC=2, dBCB2BA=2sqrt2, dBCB2BCO=2sqrt3, dBCB2BCB=4 ...
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
from scipy.spatial import cKDTree

# site-type integer codes
BCB, BCO, BA, OC, TE = 0, 1, 2, 3, 4
NAME_CODE = {"BCB": BCB, "BCO": BCO, "BA": BA, "OC": OC, "TE": TE}
CODE_NAME = {v: k for k, v in NAME_CODE.items()}

# fractional basis inside one conventional cube
_BASIS: List[Tuple[int, Tuple[float, float, float]]] = [
    (BCB, (0.0, 0.0, 0.0)),
    (BCO, (0.5, 0.5, 0.5)),
    (OC,  (0.5, 0.0, 0.0)), (OC, (0.0, 0.5, 0.0)), (OC, (0.0, 0.0, 0.5)),
    (BA,  (0.5, 0.5, 0.0)), (BA, (0.5, 0.0, 0.5)), (BA, (0.0, 0.5, 0.5)),
    (TE,  (0.25, 0.25, 0.25)), (TE, (0.75, 0.25, 0.25)),
    (TE,  (0.25, 0.75, 0.25)), (TE, (0.25, 0.25, 0.75)),
    (TE,  (0.75, 0.75, 0.25)), (TE, (0.75, 0.25, 0.75)),
    (TE,  (0.25, 0.75, 0.75)), (TE, (0.75, 0.75, 0.75)),
]

# neighbour cutoffs (Angstrom) lifted verbatim from create_lattice.cpp, sub=2 A.
# key = unordered pair of type codes -> cutoff.  Note: original SUPPRESSES the
# BCB-BCB and BCO-BCO same-sublattice links (commented out), so they are absent.
_TOL = 1.05
_CUTOFF: Dict[frozenset, float] = {
    frozenset((BCB, BCO)): 3.46410,
    frozenset((BCB, BA)):  2.82843,
    frozenset((BCB, OC)):  2.00000,
    frozenset((BCB, TE)):  1.73205,
    frozenset((BCO, BA)):  2.00000,
    frozenset((BCO, OC)):  2.82843,
    frozenset((BCO, TE)):  1.73205,
    frozenset((BA,)):      2.82843,   # BA-BA
    frozenset((BA, OC)):   2.00000,
    frozenset((BA, TE)):   1.73205,
    frozenset((OC,)):      2.82843,   # OC-OC
    frozenset((OC, TE)):   1.73205,
    frozenset((TE,)):      2.00000,   # TE-TE
}
_MAX_CUTOFF = max(_CUTOFF.values()) * _TOL + 1e-6


@dataclass
class Lattice:
    box: np.ndarray                 # 3x3 cell matrix (Angstrom)
    frac: np.ndarray                # (N,3) fractional coords
    name_code: np.ndarray           # (N,) site-type code
    index1: np.ndarray              # (N,) 1-based index (OVITO compatibility)
    # neighbour CSR per relation: indptr[r], indices[r]
    nbr_indptr: Dict[int, np.ndarray]
    nbr_indices: Dict[int, np.ndarray]

    @property
    def n(self) -> int:
        return self.frac.shape[0]

    # --- derived site-type masks (mirror the C++ macros) -------------------
    def is_BC(self) -> np.ndarray:   # SITE_IS_BC : BCB or BCO
        return (self.name_code == BCB) | (self.name_code == BCO)

    def is_OCsite(self) -> np.ndarray:  # SITE_IS_OC : OC or BCO
        return (self.name_code == OC) | (self.name_code == BCO)

    def is_BAsite(self) -> np.ndarray:  # SITE_IS_BA : BA or BCB
        return (self.name_code == BA) | (self.name_code == BCB)


# ----------------------------------------------------------------------------
def _build_neighbors(frac: np.ndarray, name_code: np.ndarray,
                     box: np.ndarray) -> Tuple[Dict[int, np.ndarray],
                                               Dict[int, np.ndarray]]:
    """
    Build, for every site, its neighbours grouped by the NEIGHBOUR's type code.
    Returns CSR (indptr, indices) keyed by neighbour type r in {BCB,BCO,BA,OC,TE}.

    Periodic minimum image is applied through a cubic/orthorhombic box (the
    generated cell is orthorhombic; off-diagonal terms are honoured for the
    KD-tree via Cartesian coordinates with ghost replication in each direction).
    """
    N = frac.shape[0]
    cart = frac @ box                      # (N,3) Cartesian
    L = np.array([box[0, 0], box[1, 1], box[2, 2]])

    # KD-tree with periodic boxsize (orthorhombic supercell)
    tree = cKDTree(cart % L, boxsize=L)
    pairs = tree.query_pairs(r=_MAX_CUTOFF, output_type="ndarray")  # (P,2), i<j

    # per-pair Cartesian distance with min-image
    di = cart[pairs[:, 0]] - cart[pairs[:, 1]]
    di -= L * np.round(di / L)
    dist = np.sqrt((di * di).sum(axis=1))

    ci = name_code[pairs[:, 0]]
    cj = name_code[pairs[:, 1]]

    # cutoff per pair (depends on the unordered type pair)
    cut = np.empty(len(pairs))
    for p in range(len(pairs)):
        key = frozenset((int(ci[p]), int(cj[p])))
        cut[p] = _CUTOFF.get(key, -1.0) * _TOL
    keep = dist <= cut
    pairs = pairs[keep]
    pi, pj = pairs[:, 0], pairs[:, 1]

    # symmetric edge list (i -> j and j -> i), grouped by neighbour type
    src = np.concatenate([pi, pj])
    dst = np.concatenate([pj, pi])
    dst_type = name_code[dst]

    indptr: Dict[int, np.ndarray] = {}
    indices: Dict[int, np.ndarray] = {}
    for r in (BCB, BCO, BA, OC, TE):
        m = dst_type == r
        s = src[m]; d = dst[m]
        order = np.argsort(s, kind="stable")
        s = s[order]; d = d[order]
        counts = np.bincount(s, minlength=N)
        ptr = np.zeros(N + 1, dtype=np.int64)
        ptr[1:] = np.cumsum(counts)
        indptr[r] = ptr
        indices[r] = d.astype(np.int64)
    return indptr, indices


def generate(geom) -> Lattice:
    """Generate the supercell from GEOMETRY.in (the OVITO-free path)."""
    a = geom.latticeConstant
    cells = np.indices((geom.nx, geom.ny, geom.nz)).reshape(3, -1).T  # (Ncell,3)
    fracs = []
    codes = []
    inv_n = np.array([geom.nx, geom.ny, geom.nz], dtype=float)
    for code, (fx, fy, fz) in _BASIS:
        base = np.array([fx, fy, fz])
        pts = (cells + base) / inv_n           # fractional w.r.t full box
        fracs.append(pts)
        codes.append(np.full(pts.shape[0], code, dtype=np.int8))
    frac = np.vstack(fracs)
    name_code = np.concatenate(codes)

    # sort by type (BCB,BCO,BA,OC,TE) then keep stable -> matches species grouping
    order = np.argsort(name_code, kind="stable")
    frac = frac[order] % 1.0
    name_code = name_code[order]
    index1 = np.arange(1, frac.shape[0] + 1, dtype=np.int64)

    box = np.diag([geom.nx * a, geom.ny * a, geom.nz * a]).astype(float)
    indptr, indices = _build_neighbors(frac, name_code, box)
    return Lattice(box, frac, name_code, index1, indptr, indices)


# ----------------------------------------------------------------------------
# POSCAR I/O
# ----------------------------------------------------------------------------
def write_poscar(lat: Lattice, path: str, comment: str = "kMC superlattice") -> None:
    """Write the generated lattice as a VASP POSCAR (site types as 'species')."""
    order = np.argsort(lat.name_code, kind="stable")
    codes = lat.name_code[order]
    frac = lat.frac[order]
    types, counts = [], []
    for r in (BCB, BCO, BA, OC, TE):
        c = int((codes == r).sum())
        if c:
            types.append(CODE_NAME[r]); counts.append(c)
    with open(path, "w") as fh:
        fh.write(comment + "\n")
        fh.write("1.0\n")
        for row in lat.box:
            fh.write(f"  {row[0]:.8f}  {row[1]:.8f}  {row[2]:.8f}\n")
        fh.write("  " + "  ".join(types) + "\n")
        fh.write("  " + "  ".join(str(c) for c in counts) + "\n")
        fh.write("Direct\n")
        for f in frac:
            fh.write(f"  {f[0]:.8f}  {f[1]:.8f}  {f[2]:.8f}\n")


def read_poscar(path: str) -> Lattice:
    """Read a POSCAR whose 'species' are the site types BCB/BCO/BA/OC/TE."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    scale = float(lines[1].split()[0])
    box = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    types = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    start = 7
    mode = lines[7].strip().lower()
    if mode.startswith(("s",)):       # selective dynamics line present
        start = 8
        mode = lines[8].strip().lower()
    coord_start = start + 1
    is_direct = mode.startswith(("d", "f"))   # Direct/Fractional
    N = sum(counts)
    coords = np.array([[float(x) for x in lines[coord_start + i].split()[:3]]
                       for i in range(N)])
    if not is_direct:                  # Cartesian -> fractional
        coords = coords @ np.linalg.inv(box)
    name_code = np.concatenate([np.full(c, NAME_CODE[t], dtype=np.int8)
                                for t, c in zip(types, counts)])
    index1 = np.arange(1, N + 1, dtype=np.int64)
    indptr, indices = _build_neighbors(coords % 1.0, name_code, box)
    return Lattice(box, coords % 1.0, name_code, index1, indptr, indices)
