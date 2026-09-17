"""
ff_data.py - per-species characterisation properties for the kMC pseudo-species.

The kMC uses coarse pseudo-species (F, O, S, N, Li, plus composite fragments
SFO and F5D). Zeo++ needs a hard-sphere RADIUS per label; RASPA needs a
PseudoAtom (mass, charge) + VDW (epsilon, sigma) per label.

Composite kMC fragments are mapped to a representative real atom:
    SFO  -> S   (S-F-O fragment, dominant heavy atom S)
    F5D  -> F   (fluorinated decomposition fragment, F-like)
These proxies are APPROXIMATIONS; refine with your own DFT/empirical data.

UFF Lennard-Jones parameters: Rappe et al., J. Am. Chem. Soc. 1992, 114, 10024.
  epsilon in Kelvin, sigma in Angstrom.
Zeo++ radii here are vdW-like effective radii (Angstrom); tune for your system.
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass
class SpeciesProp:
    element: str        # real-atom proxy used for FF/radii lookup
    mass: float         # amu
    charge: float       # partial charge (e); placeholder, refine with CHELPG/DFT
    radius: float       # Zeo++ hard-sphere radius (Angstrom)
    eps_K: float        # UFF epsilon (Kelvin)
    sigma_A: float      # UFF sigma (Angstrom)


# label -> properties.  Keys are the kMC species symbols as written in snapshots.
SPECIES_PROPS = {
    "Li":  SpeciesProp("Li",  6.941,  0.0, 1.82,  12.58, 2.184),
    "F":   SpeciesProp("F",  18.998, -0.5, 1.47,  36.48, 2.997),
    "O":   SpeciesProp("O",  15.999, -0.5, 1.52,  30.19, 3.118),
    "S":   SpeciesProp("S",  32.06,  -0.3, 1.80, 137.90, 3.595),
    "N":   SpeciesProp("N",  14.007, -0.4, 1.55,  38.29, 3.261),
    # composite fragments (proxy element noted)
    "SFO": SpeciesProp("S",  60.07,   0.0, 1.80, 137.90, 3.595),
    "F5D": SpeciesProp("F",  50.0,    0.0, 1.47,  36.48, 2.997),
}

# Common probe molecules (single-site) for Widom / pore accessibility.
#   radius for Zeo++ accessibility; UFF eps/sigma + mass for RASPA Widom.
PROBES = {
    "He":  dict(radius=1.40, mass=4.0026,  eps_K=10.9,  sigma_A=2.640, charge=0.0),
    "Li":  dict(radius=0.76, mass=6.941,   eps_K=12.58, sigma_A=2.184, charge=1.0),
    "N2":  dict(radius=1.82, mass=28.013,  eps_K=34.7,  sigma_A=3.298, charge=0.0),
    # united-atom DME-like solvent probe (very coarse, for solvent uptake trend)
    "DME": dict(radius=2.40, mass=90.12,   eps_K=200.0, sigma_A=4.50,  charge=0.0),
}


def props_for(label: str) -> SpeciesProp:
    if label in SPECIES_PROPS:
        return SPECIES_PROPS[label]
    raise KeyError(f"No characterisation properties for species '{label}'. "
                   f"Add it to ff_data.SPECIES_PROPS.")
