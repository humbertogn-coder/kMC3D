"""
postprocess.py
==============
Adapter layer so the kMC output plugs straight into the user's existing
analysis modules (z-warp, gaussian density field, morphology-parameter
extraction, ML/SHAP models).

Nothing here changes the simulation; these are *readers* and light
transforms over the standard outputs:

    output/trajectory/kmc-coords-<step>.xyz     (extended XYZ)
    Data2Excel.txt                              (reaction-count log)

Design
------
Each external module typically wants one of two things:
  (a) a per-frame point cloud  -> `load_frame` / `iter_frames` give you
      (species, positions[N,3], charge, site_name) arrays.
  (b) the scalar time series   -> `load_data2excel` returns a dict of
      columns (numpy arrays), keyed by the DATA_COLUMNS names.

The z-warp, density-field and morphology helpers below are deliberately
thin and dependency-light (numpy only).  They are meant as *connection
points*: if you already have your own implementation, call `iter_frames`
and feed your function; if not, the defaults give a reasonable baseline.
"""

from __future__ import annotations
import glob
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterator, List, Tuple

import numpy as np

from .output import DATA_COLUMNS


# ---------------------------------------------------------------------------
# Trajectory readers
# ---------------------------------------------------------------------------
@dataclass
class Frame:
    step: int
    box: np.ndarray                 # 3x3
    species: np.ndarray             # (N,) str
    pos: np.ndarray                 # (N,3) cartesian
    index: np.ndarray               # (N,) 1-based site index
    charge: np.ndarray              # (N,) int
    site: np.ndarray                # (N,) str  (BCB/BCO/BA/OC/TE)


_LAT_RE = re.compile(r'Lattice="([^"]+)"')
_STEP_RE = re.compile(r"kmc-coords-(\d+)\.xyz$")


def load_frame(path: str) -> Frame:
    """Parse one extended-XYZ frame written by output.Output.write_xyz."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    m = _LAT_RE.search(lines[1])
    box = (np.array([float(x) for x in m.group(1).split()]).reshape(3, 3)
           if m else np.eye(3))
    species, pos, index, charge, site = [], [], [], [], []
    for ln in lines[2:2 + n]:
        f = ln.split()
        species.append(f[0])
        pos.append([float(f[1]), float(f[2]), float(f[3])])
        index.append(int(f[4]) if len(f) > 4 else 0)
        charge.append(int(f[5]) if len(f) > 5 else 0)
        site.append(f[6] if len(f) > 6 else "")
    sm = _STEP_RE.search(os.path.basename(path))
    step = int(sm.group(1)) if sm else -1
    return Frame(step, box, np.array(species), np.array(pos, float),
                 np.array(index, int), np.array(charge, int), np.array(site))


def iter_frames(traj_dir: str = "output/trajectory") -> Iterator[Frame]:
    """Yield frames in step order."""
    paths = glob.glob(os.path.join(traj_dir, "kmc-coords-*.xyz"))
    paths.sort(key=lambda p: int(_STEP_RE.search(os.path.basename(p)).group(1)))
    for p in paths:
        yield load_frame(p)


# ---------------------------------------------------------------------------
# Reaction-count log reader
# ---------------------------------------------------------------------------
def load_data2excel(path: str = "Data2Excel.txt") -> Dict[str, np.ndarray]:
    """Return {column_name: array}.  Tolerates the optional '# ' header row."""
    rows: List[List[str]] = []
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            rows.append(s.split())
    if not rows:
        return {c: np.array([]) for c in DATA_COLUMNS}
    ncol = len(rows[0])
    cols = DATA_COLUMNS[:ncol]
    out: Dict[str, np.ndarray] = {}
    for j, c in enumerate(cols):
        raw = [r[j] for r in rows if len(r) > j]
        try:
            out[c] = np.array([float(x) for x in raw])
        except ValueError:
            out[c] = np.array(raw)            # e.g. the textual 'type' column
    return out


# ---------------------------------------------------------------------------
# z-warp
# ---------------------------------------------------------------------------
def z_warp(frame: Frame, ref_species: Tuple[str, ...] = ("Li",),
           nbins: int = 200) -> Dict[str, np.ndarray]:
    """
    Flatten the electrode/SEI interface by remapping z onto the local mean
    interface height, so morphology is measured relative to the growing
    front rather than the absolute box frame.

    Returns dict with z (original), z_warped, and the per-(x,y) reference
    surface used.  Swap in your own z-warp by consuming `frame.pos` directly.
    """
    sel = np.isin(frame.species, list(ref_species))
    p = frame.pos[sel]
    if p.shape[0] == 0:
        return dict(z=frame.pos[:, 2], z_warped=frame.pos[:, 2],
                    surface=np.array([]))
    # coarse (x,y) grid; reference height = max z of ref atoms per cell
    lx, ly = frame.box[0, 0], frame.box[1, 1]
    gx = np.clip((p[:, 0] / lx * nbins).astype(int), 0, nbins - 1)
    gy = np.clip((p[:, 1] / ly * nbins).astype(int), 0, nbins - 1)
    surf = np.full((nbins, nbins), np.nan)
    for cx, cy, cz in zip(gx, gy, p[:, 2]):
        if np.isnan(surf[cx, cy]) or cz > surf[cx, cy]:
            surf[cx, cy] = cz
    mean_h = np.nanmean(surf)
    surf = np.where(np.isnan(surf), mean_h, surf)
    # warp every atom by its column reference
    ax = np.clip((frame.pos[:, 0] / lx * nbins).astype(int), 0, nbins - 1)
    ay = np.clip((frame.pos[:, 1] / ly * nbins).astype(int), 0, nbins - 1)
    z_ref = surf[ax, ay]
    return dict(z=frame.pos[:, 2], z_warped=frame.pos[:, 2] - z_ref, surface=surf)


# ---------------------------------------------------------------------------
# Gaussian density field
# ---------------------------------------------------------------------------
def gaussian_density_field(frame: Frame, grid=(48, 48, 96), sigma: float = 1.5,
                           species: Tuple[str, ...] | None = None) -> np.ndarray:
    """
    Deposit selected atoms onto a 3D grid and convolve with a Gaussian to get
    a smooth density field (useful for porosity / dendrite descriptors and as
    a CNN input).  numpy-only separable Gaussian; pass species=None for all.
    """
    sel = (np.ones(frame.species.shape, bool) if species is None
           else np.isin(frame.species, list(species)))
    p = frame.pos[sel]
    nx, ny, nz = grid
    L = frame.box.diagonal()
    field = np.zeros(grid, float)
    if p.shape[0]:
        ix = np.clip((p[:, 0] / L[0] * nx).astype(int), 0, nx - 1)
        iy = np.clip((p[:, 1] / L[1] * ny).astype(int), 0, ny - 1)
        iz = np.clip((p[:, 2] / L[2] * nz).astype(int), 0, nz - 1)
        np.add.at(field, (ix, iy, iz), 1.0)
    # separable Gaussian blur in grid units
    s = max(sigma / (L[2] / nz), 1e-3)
    rad = int(3 * s) + 1
    k = np.exp(-0.5 * (np.arange(-rad, rad + 1) / s) ** 2)
    k /= k.sum()
    for ax in range(3):
        field = np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"),
                                    ax, field)
    return field


# ---------------------------------------------------------------------------
# Morphology parameters
# ---------------------------------------------------------------------------
def morphology_parameters(frame: Frame,
                          sei_exclude=None) -> Dict[str, float]:
    """
    Scalar descriptors per frame, ready to be stacked into a feature table for
    ML/SHAP.  SEI = the anode film (species.is_anode_film: anode SEI products
    plus *_an deposits); cathode material and CEI film are NOT anode SEI.
    Passing `sei_exclude` (a tuple of labels) restores the legacy rule
    "everything not excluded is SEI".  Extend freely; the keys become your
    ML feature names.  Cathode / CEI counts are reported separately.
    """
    from .species import classify_array
    sp = frame.species
    is_li = sp == "Li"
    cls = classify_array(sp)
    if sei_exclude is not None:
        is_sei = ~np.isin(sp, list(sei_exclude))
    else:
        is_sei = (cls == "sei") | (cls == "deposit")
    is_cath = cls == "cathode"
    is_cei = cls == "cei"
    z = frame.pos[:, 2]
    Lz = frame.box[2, 2]
    out = {
        "n_atoms": float(sp.size),
        "n_Li": float(is_li.sum()),
        "n_Li_ion": float((is_li & (frame.charge != 0)).sum()),
        "n_SEI": float(is_sei.sum()),
        "sei_fraction": float(is_sei.sum() / max(sp.size, 1)),
        "li_mean_z": float(z[is_li].mean()) if is_li.any() else 0.0,
        "li_front_z": float(z[is_li].max()) if is_li.any() else 0.0,
        "sei_mean_z": float(z[is_sei].mean()) if is_sei.any() else 0.0,
        "sei_thickness": (float(z[is_sei].max() - z[is_sei].min())
                          if is_sei.any() else 0.0),
        "roughness_std_z": float(z[is_li].std()) if is_li.any() else 0.0,
        "coverage_frac": float(np.unique(
            (frame.pos[is_sei, :2] / frame.box.diagonal()[:2] * 32)
            .astype(int), axis=0).shape[0] / (32 * 32)) if is_sei.any() else 0.0,
        "box_z": float(Lz),
        # full-cell extras (0 for anode-only runs)
        "n_cathode": float(is_cath.sum()),
        "n_CEI": float(is_cei.sum()),
        "n_deposit": float((cls == "deposit").sum()),
    }
    return out


def morphology_table(traj_dir: str = "output/trajectory") -> Dict[str, np.ndarray]:
    """Build a {feature: array_over_frames} table for the whole trajectory."""
    rows, steps = [], []
    for fr in iter_frames(traj_dir):
        rows.append(morphology_parameters(fr)); steps.append(fr.step)
    if not rows:
        return {}
    keys = rows[0].keys()
    table = {k: np.array([r[k] for r in rows]) for k in keys}
    table["step"] = np.array(steps)
    return table


# ---------------------------------------------------------------------------
# ML feature matrix
# ---------------------------------------------------------------------------
def ml_feature_matrix(traj_dir: str = "output/trajectory",
                      data2excel: str = "Data2Excel.txt"
                      ) -> Tuple[np.ndarray, List[str]]:
    """
    Join per-frame morphology descriptors with the nearest reaction-count row
    to produce (X, feature_names) for a model / SHAP analysis.
    """
    morph = morphology_table(traj_dir)
    if not morph:
        return np.empty((0, 0)), []
    feat_names = [k for k in morph if k != "step"]
    X = np.column_stack([morph[k] for k in feat_names])
    if os.path.exists(data2excel):
        d = load_data2excel(data2excel)
        if "currentStep" in d and d["currentStep"].size:
            cs = d["currentStep"]
            add = ["LiMetal", "LiIon", "nO", "nF", "RxnPlating", "RxnStripping"]
            extra = []
            for name in add:
                if name in d:
                    col = np.interp(morph["step"], cs, d[name])
                    extra.append(col); feat_names.append(name)
            if extra:
                X = np.column_stack([X] + extra)
    return X, feat_names


# ---------------------------------------------------------------------------
# Unified ML table: cycle ledger + pore/adsorption descriptors
# ---------------------------------------------------------------------------
def run_feature_table(out_dir: str,
                      characterization: bool = False,
                      **char_kw) -> Dict[str, np.ndarray]:
    """One-call ML table for a finished run directory.

    Loads <out_dir>/cycle_stats.csv (per-half-cycle targets and descriptors,
    written natively by the engine) and, optionally, joins Zeo++/RASPA
    descriptors from <out_dir>/snapshots on the 'cycle' column.

    Returns {column: np.ndarray}. Feed directly into pandas / sklearn / SHAP:
        import pandas as pd
        df = pd.DataFrame(run_feature_table("output", characterization=True))
    """
    from .ensemble import load_cycle_stats
    path = os.path.join(out_dir, "cycle_stats.csv")
    tbl = load_cycle_stats(path)
    if not tbl:
        return {}
    if characterization:
        from .characterize import characterize_snapshots
        snap_dir = os.path.join(out_dir, "snapshots")
        if os.path.isdir(snap_dir):
            char = characterize_snapshots(snap_dir, **char_kw)
            cyc_char = np.array(char.get("cycle", []), dtype=float)
            if cyc_char.size:
                base_cyc = tbl["cycle"]
                for k, vals in char.items():
                    if k == "cycle":
                        continue
                    col = np.full(base_cyc.size, np.nan)
                    for cc, vv in zip(cyc_char, vals):
                        j = np.where(base_cyc == cc)[0]
                        if j.size:
                            col[j[0]] = vv
                    tbl[f"char_{k}"] = col
    return tbl
