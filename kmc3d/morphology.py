"""
morphology.py
=============
Lattice-topological morphology descriptors of the Li anode from one .xyz
frame (extended XYZ written by output.write_xyz). Complements the Gaussian
density route of postprocess.py with descriptors that need no voxelization:
on a fixed lattice, connectivity is exact and cheap.

Descriptors (all relative to the anode, so full cells with a cathode wrapping
the z faces are handled):

    li_main_body        Li sites connected (Li-Li network, first neighbours)
                        to the largest Li component
    li_buried_topo      Li sites NOT in the main body: isolated Li islands,
                        the topological definition of dead / buried Li
                        (Perez-Beltran 2024 and the group's 2026 draft use
                        the same idea on a voxel field; here on the lattice)
    li_buried_coating   composition of the first shell around buried Li
                        (fraction SEI, deposit, CEI, electrolyte, cathode)
    porosity_profile    per z-slab (relative to the anode centre) fraction of
                        lattice sites that are NOT solid (empty or electrolyte)
    film_thickness      z extent of the anode film (SEI + deposits) above and
                        below the initial anode band
    surface_roughness   std of the Li front height on an xy grid, per side

Solid = metal, sei, deposit, cei, cathode (species.classify).

Neighbour cutoff: 3.5 A joins BC-BC (3.46 A body-diagonal in the 4 A cube),
BC-BA/OC (2.0 to 2.83 A) and BA-OC shells while excluding 4.0 A second
neighbours. Periodic in x, y, z (the cell is periodic in all three).
"""

from __future__ import annotations
import glob
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .postprocess import Frame, load_frame
from .species import classify_array

NBR_CUTOFF = 3.5


def _tree(frame: Frame):
    L = frame.box.diagonal()
    pos = np.mod(frame.pos, L)
    return cKDTree(pos, boxsize=L), pos


def li_components(frame: Frame) -> Tuple[np.ndarray, np.ndarray, int]:
    """Connected components of the Li-Li first-neighbour graph.
    Returns (labels over ALL atoms, -1 for non-Li; size of each label;
    label of the main body)."""
    is_li = frame.species == "Li"
    idx = np.where(is_li)[0]
    labels = np.full(frame.species.size, -1, dtype=int)
    if idx.size == 0:
        return labels, np.array([], int), -1
    tree, pos = _tree(frame)
    sub = cKDTree(pos[idx], boxsize=frame.box.diagonal())
    pairs = sub.query_pairs(NBR_CUTOFF, output_type="ndarray")
    n = idx.size
    if pairs.size:
        a = coo_matrix((np.ones(pairs.shape[0]), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
        ncomp, lab = connected_components(a, directed=False)
    else:
        ncomp, lab = n, np.arange(n)
    labels[idx] = lab
    sizes = np.bincount(lab, minlength=ncomp)
    return labels, sizes, int(np.argmax(sizes))


def buried_li(frame: Frame) -> Dict[str, object]:
    labels, sizes, main = li_components(frame)
    is_li = labels >= 0
    if not is_li.any():
        return {"n_Li": 0, "n_li_main": 0, "n_li_buried": 0, "frac_buried": 0.0,
                "n_islands": 0, "coating": {}}
    in_main = labels == main
    buried = is_li & ~in_main
    out = {"n_Li": int(is_li.sum()), "n_li_main": int(in_main.sum()),
           "n_li_buried": int(buried.sum()),
           "frac_buried": float(buried.sum() / is_li.sum()),
           "n_islands": int((sizes > 0).sum() - 1)}
    # first-shell coating of buried Li (non-Li neighbours by class)
    coating: Dict[str, int] = {}
    if buried.any():
        tree, pos = _tree(frame)
        cls = classify_array(frame.species)
        nb = tree.query_ball_point(pos[buried], NBR_CUTOFF)
        flat = np.unique(np.concatenate([np.asarray(x, int) for x in nb])) if nb else np.array([], int)
        flat = flat[frame.species[flat] != "Li"]
        for c in np.unique(cls[flat]):
            coating[str(c)] = int((cls[flat] == c).sum())
        tot = max(sum(coating.values()), 1)
        coating = {k: v / tot for k, v in coating.items()}
    out["coating"] = coating
    return out


def anode_center_z(frame: Frame) -> float:
    """z of the anode centre: median z of the main Li body (robust to
    dendrites and to Li far from the slab)."""
    labels, sizes, main = li_components(frame)
    z = frame.pos[labels == main, 2]
    return float(np.median(z)) if z.size else 0.5 * frame.box[2, 2]


def porosity_profile(frame: Frame, all_sites_z: np.ndarray, slab: float = 4.0,
                     z0: Optional[float] = None) -> Dict[str, np.ndarray]:
    """Fraction of lattice sites per z-slab that are not solid. all_sites_z:
    z of every lattice site (from lattice.generate on the case GEOMETRY.in),
    needed because empty sites are not written to the .xyz. z is measured
    from the anode centre (z0) with periodic wrap."""
    Lz = frame.box[2, 2]
    z0 = anode_center_z(frame) if z0 is None else z0
    cls = classify_array(frame.species)
    solid = np.isin(cls, ["metal", "sei", "deposit", "cei", "cathode"])
    rel_all = (all_sites_z - z0 + 0.5 * Lz) % Lz - 0.5 * Lz
    rel_sol = (frame.pos[solid, 2] - z0 + 0.5 * Lz) % Lz - 0.5 * Lz
    edges = np.arange(-0.5 * Lz, 0.5 * Lz + slab, slab)
    n_all, _ = np.histogram(rel_all, edges)
    n_sol, _ = np.histogram(rel_sol, edges)
    with np.errstate(divide="ignore", invalid="ignore"):
        por = np.where(n_all > 0, 1.0 - n_sol / np.maximum(n_all, 1), np.nan)
    return {"z_rel": 0.5 * (edges[:-1] + edges[1:]), "porosity": por,
            "n_sites": n_all, "n_solid": n_sol}


def film_and_roughness(frame: Frame, grid_xy: int = 10) -> Dict[str, float]:
    """Anode film extent and Li front roughness on both sides of the slab."""
    Lz = frame.box[2, 2]
    z0 = anode_center_z(frame)
    cls = classify_array(frame.species)
    rel = (frame.pos[:, 2] - z0 + 0.5 * Lz) % Lz - 0.5 * Lz
    film = (cls == "sei") | (cls == "deposit")
    out: Dict[str, float] = {}
    for side, sgn in (("top", 1), ("bottom", -1)):
        f = film & (sgn * rel > 0)
        # 90th percentile of |z|: robust to stray film sites near the cathode
        out[f"film_extent_{side}"] = float(np.percentile(np.abs(rel[f]), 90)) if f.any() else 0.0
        li = (cls == "metal") & (sgn * rel > 0)
        if li.any():
            L = frame.box.diagonal()
            gx = np.clip((frame.pos[li, 0] / L[0] * grid_xy).astype(int), 0, grid_xy - 1)
            gy = np.clip((frame.pos[li, 1] / L[1] * grid_xy).astype(int), 0, grid_xy - 1)
            key = gx * grid_xy + gy
            h = np.full(grid_xy * grid_xy, -np.inf)
            np.maximum.at(h, key, sgn * rel[li])
            h = h[np.isfinite(h)]
            out[f"li_front_{side}"] = float(h.mean()) if h.size else 0.0
            out[f"roughness_{side}"] = float(h.std()) if h.size else 0.0
        else:
            out[f"li_front_{side}"] = 0.0; out[f"roughness_{side}"] = 0.0
    out["anode_center_z"] = z0
    return out


def frame_descriptors(frame: Frame, all_sites_z: Optional[np.ndarray] = None) -> Dict[str, float]:
    """Flat dict of scalar descriptors for one frame (ML-ready)."""
    d: Dict[str, float] = {"step": float(frame.step)}
    b = buried_li(frame)
    for k in ("n_Li", "n_li_main", "n_li_buried", "frac_buried", "n_islands"):
        d[k] = float(b[k])
    for k, v in b["coating"].items():
        d[f"buried_coat_{k}"] = float(v)
    d.update(film_and_roughness(frame))
    cls = classify_array(frame.species)
    for c in ("sei", "deposit", "cei", "cathode", "electrolyte"):
        d[f"n_{c}"] = float((cls == c).sum())
    if all_sites_z is not None:
        p = porosity_profile(frame, all_sites_z)
        Lz = frame.box[2, 2]
        film_top = d["film_extent_top"]; film_bot = d["film_extent_bottom"]
        sel = ((p["z_rel"] > 0) & (p["z_rel"] <= max(film_top, 4.0))) | \
              ((p["z_rel"] < 0) & (p["z_rel"] >= -max(film_bot, 4.0)))
        por = p["porosity"][sel]
        d["film_porosity_mean"] = float(np.nanmean(por)) if por.size else np.nan
        d["total_solid_fraction"] = float(p["n_solid"].sum() / max(p["n_sites"].sum(), 1))
    return d


def lattice_site_z(geometry_file: str) -> np.ndarray:
    """z of every lattice site for the case geometry (empty sites included)."""
    from . import config, lattice as L
    g = config.Geometry.read(geometry_file)
    lat = L.generate(g)
    return (lat.frac @ lat.box)[:, 2]


def case_morphology(case_dir: str, frames: str = "first_last") -> List[Dict[str, float]]:
    """Descriptors for every seed of a case. frames: 'first_last' (the sampled
    frames brought from GRACE) or 'all'."""
    zs = lattice_site_z(os.path.join(case_dir, "GEOMETRY.in"))
    rows = []
    for seed_dir in sorted(glob.glob(os.path.join(case_dir, "runs", "seed_*"))):
        files = sorted(glob.glob(os.path.join(seed_dir, "trajectory", "kmc-coords-*.xyz")),
                       key=lambda p: int(os.path.basename(p)[11:-4]))
        if frames == "first_last" and len(files) > 2:
            files = [files[0], files[-1]]
        for f in files:
            d = frame_descriptors(load_frame(f), zs)
            d["seed"] = os.path.basename(seed_dir)
            rows.append(d)
    return rows
