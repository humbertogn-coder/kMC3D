"""
zeopp.py - geometric pore analysis of kMC snapshots with Zeo++.

Zeo++ ( http://www.zeoplusplus.org ) is a standalone C++ binary (`network`).
It is NOT pip-installable; install separately and point ZEOPP_BIN at the
`network` executable (or put it on PATH).

What it gives you (all purely geometric, hard-sphere):
  Di  = largest included sphere      (LCD, largest cavity diameter)
  Df  = largest free sphere          (PLD, pore-limiting diameter / bottleneck)
  Dif = largest included sphere along the free path
  plus accessible surface area and accessible volume for a chosen probe.

These complement the Gaussian+warp descriptors with an interpretable
bottleneck/connectivity measure that text-free morphology cannot give.
"""

from __future__ import annotations
import os
import shutil
import subprocess
import tempfile
from typing import Dict, Optional

from .structio import read_poscar, write_cssr
from .ff_data import SPECIES_PROPS


def _zeopp_bin() -> Optional[str]:
    return os.environ.get("ZEOPP_BIN") or shutil.which("network")


def write_radii_file(path: str) -> str:
    """Write a Zeo++ radii file mapping every kMC species label to its radius."""
    with open(path, "w") as fh:
        for label, p in SPECIES_PROPS.items():
            fh.write(f"{label} {p.radius:.3f}\n")
    return path


def _parse_res(path: str) -> Dict[str, float]:
    # .res line:  <name>   Di   Df   Dif
    with open(path) as fh:
        toks = fh.read().split()
    vals = [t for t in toks if _is_float(t)]
    if len(vals) >= 3:
        di, df, dif = (float(vals[0]), float(vals[1]), float(vals[2]))
        return {"LCD_largest_cavity": di, "PLD_pore_limiting": df,
                "Dif_free_path_sphere": dif}
    return {}


def _parse_sa(path: str) -> Dict[str, float]:
    out = {}
    with open(path) as fh:
        txt = fh.read()
    for key, tag in (("ASA_A2", "ASA_A^2:"), ("ASA_m2_cm3", "ASA_m^2/cm^3:"),
                     ("ASA_m2_g", "ASA_m^2/g:"), ("NASA_A2", "NASA_A^2:")):
        v = _after(txt, tag)
        if v is not None:
            out[key] = v
    return out


def _parse_vol(path: str) -> Dict[str, float]:
    out = {}
    with open(path) as fh:
        txt = fh.read()
    for key, tag in (("AV_A3", "AV_A^3:"), ("AV_cm3_g", "AV_cm^3/g:"),
                     ("AV_volume_fraction", "AV_Volume_fraction:"),
                     ("NAV_A3", "NAV_A^3:")):
        v = _after(txt, tag)
        if v is not None:
            out[key] = v
    return out


def zeopp_descriptors(poscar_path: str,
                      probe_radius: float = 1.4,
                      chan_radius: float = 1.4,
                      n_samples: int = 2000,
                      keep_files: bool = False) -> Dict[str, float]:
    """Run Zeo++ -res/-sa/-vol on one snapshot. Returns a descriptor dict.

    probe_radius : probe sphere for surface area / volume (e.g. 1.4 = N2-ish).
    chan_radius  : channel probe for accessibility.
    Raises RuntimeError if the `network` binary is not found.
    """
    binpath = _zeopp_bin()
    if binpath is None:
        raise RuntimeError(
            "Zeo++ 'network' binary not found. Install Zeo++ and set ZEOPP_BIN "
            "or add it to PATH. (http://www.zeoplusplus.org)")

    struct = read_poscar(poscar_path)
    work = tempfile.mkdtemp(prefix="zeopp_")
    try:
        cssr = write_cssr(struct, os.path.join(work, "snap.cssr"))
        rad = write_radii_file(os.path.join(work, "kmc.rad"))
        base = os.path.join(work, "snap")

        def run(args, out_ext):
            subprocess.run([binpath, "-r", rad, *args, cssr],
                           cwd=work, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return base + out_ext

        desc: Dict[str, float] = {}
        desc.update(_parse_res(run(["-res"], ".res")))
        desc.update(_parse_sa(run(
            ["-sa", str(chan_radius), str(probe_radius), str(n_samples)], ".sa")))
        desc.update(_parse_vol(run(
            ["-vol", str(chan_radius), str(probe_radius), str(n_samples)], ".vol")))
        return desc
    finally:
        if not keep_files:
            shutil.rmtree(work, ignore_errors=True)


# ------------------------------------------------------------------ helpers
def _is_float(t: str) -> bool:
    try:
        float(t); return True
    except ValueError:
        return False


def _after(txt: str, tag: str):
    i = txt.find(tag)
    if i < 0:
        return None
    rest = txt[i + len(tag):].split()
    return float(rest[0]) if rest and _is_float(rest[0]) else None
