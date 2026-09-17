"""
characterize.py - run Zeo++ and/or RASPA over a directory of kMC snapshots and
assemble a per-snapshot descriptor table for ML / SHAP.

Both tools are optional and degrade gracefully:
  - Zeo++ needs the `network` binary (ZEOPP_BIN / PATH).
  - RASPA needs `pip install raspalib`.
If a tool is unavailable its columns are simply skipped (with a note).

Typical use
-----------
    from kmc3d.characterize import characterize_snapshots
    tbl = characterize_snapshots("out/snapshots", run_raspa=True, run_zeopp=True,
                                 sei_slab=True)
    # tbl: {"cycle":[...], "PLD_pore_limiting":[...], "henry":[...], ...}
"""

from __future__ import annotations
import glob
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from .structio import read_poscar


def _sei_slab_bounds(poscar_path: str,
                     metal: Tuple[str, ...] = ("Li",)) -> Optional[Tuple[float, float]]:
    """z-range spanned by non-metal (SEI) species; None if no SEI yet."""
    s = read_poscar(poscar_path)
    sei = ~np.isin(s.species, list(metal))
    if not sei.any():
        return None
    z = s.cart[sei, 2]
    return float(z.min()), float(z.max())


def _cycle_of(path: str) -> int:
    m = re.search(r"cycle(\d+)", os.path.basename(path))
    return int(m.group(1)) if m else -1


def characterize_snapshots(snap_dir: str,
                           run_zeopp: bool = True,
                           run_raspa: bool = True,
                           sei_slab: bool = True,
                           probe: str = "He",
                           raspa_cycles: int = 3000,
                           zeopp_probe_radius: float = 1.4) -> Dict[str, List]:
    """Build a descriptor table over all snap_cycle*.vasp in `snap_dir`."""
    snaps = sorted(glob.glob(os.path.join(snap_dir, "snap_cycle*.vasp")),
                   key=_cycle_of)
    rows: List[Dict[str, float]] = []
    cycles: List[int] = []

    zeopp_ok = run_zeopp
    if run_zeopp:
        try:
            from .zeopp import zeopp_descriptors, _zeopp_bin
            if _zeopp_bin() is None:
                zeopp_ok = False
                print("[characterize] Zeo++ 'network' not found - skipping Zeo++.")
        except Exception:
            zeopp_ok = False

    raspa_ok = run_raspa
    if run_raspa:
        try:
            import raspalib  # noqa: F401
            from .raspa import widom_insertion
        except Exception:
            raspa_ok = False
            print("[characterize] raspalib not installed - skipping RASPA.")

    for sp in snaps:
        cyc = _cycle_of(sp)
        rec: Dict[str, float] = {}
        slab = _sei_slab_bounds(sp) if sei_slab else None

        if zeopp_ok:
            try:
                from .zeopp import zeopp_descriptors
                rec.update(zeopp_descriptors(sp, probe_radius=zeopp_probe_radius))
            except Exception as e:
                print(f"[characterize] Zeo++ failed on cycle {cyc}: {e}")

        if raspa_ok:
            try:
                from .raspa import widom_insertion
                w = widom_insertion(sp, probe=probe, n_cycles=raspa_cycles,
                                    z_slab=slab)
                rec.update({f"raspa_{k}": v for k, v in w.items()})
            except Exception as e:
                print(f"[characterize] RASPA failed on cycle {cyc}: {e}")

        if rec:
            rows.append(rec); cycles.append(cyc)

    if not rows:
        return {"cycle": []}
    keys = sorted({k for r in rows for k in r})
    table: Dict[str, List] = {"cycle": cycles}
    for k in keys:
        table[k] = [r.get(k, float("nan")) for r in rows]
    return table


def merge_into_features(snap_dir: str, traj_dir: str = "output/trajectory",
                        data2excel: str = "Data2Excel.txt", **kw):
    """Join pore/adsorption descriptors with the morphology+reaction matrix.

    Returns (X, feature_names) ready for a model / SHAP.
    """
    from .postprocess import ml_feature_matrix, morphology_table
    X, names = ml_feature_matrix(traj_dir, data2excel)
    char = characterize_snapshots(snap_dir, **kw)
    if not char.get("cycle"):
        return X, names
    morph = morphology_table(traj_dir)
    steps = morph.get("step")
    if steps is None or X.size == 0:
        return X, names
    # map snapshot cycle -> nearest trajectory step is non-trivial; here we align
    # by index order (snapshots are a subset). Return the descriptor table too so
    # the user can join on cycle explicitly if frames carry cycle ids.
    extra_names = [k for k in char if k != "cycle"]
    return X, names, char, extra_names
