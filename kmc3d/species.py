"""
species.py
==========
Single source of truth for classifying kMC species labels in post-processing.

Every analysis module (morphology, Zeo++, RASPA, ledger readers) must use these
helpers instead of hard-coded lists, so that new labels (S8-unit cathode
states, anode deposits, CEI films) are handled consistently.

Classes
-------
    metal        Li (metallic or ionic Li on lattice sites)
    electrolyte  ETH, SOL, FSI            (mobile, excluded from skeletons)
    cathode      S8, Li2S8, Li4S8, Li8S8, Li16S8 (S8-unit model) and the legacy
                 Li2S6, Li2S4, Li2S2, Li2S lattice labels
    deposit      *_an  (insoluble polysulfide deposits on the anode, SEI-like)
    cei          CEI_* (cathode electrolyte interphase film)
    sei          everything else (anode SEI: F, O, N, S, SFO, F5D, ...)
Dissolved reservoir species (*_d) never appear on the lattice.

S8-unit lattice labels map to conventional formulas for reporting:
    Li4S8 -> 2 Li2S4, Li8S8 -> 4 Li2S2, Li16S8 -> 8 Li2S.
"""

from __future__ import annotations
import re
from typing import Dict, Iterable, Tuple

import numpy as np

ELECTROLYTE = ("ETH", "SOL", "FSI")
METAL = ("Li",)

_CATHODE_RE = re.compile(r"^(S8|Li\d*S\d*)$")

# S8-unit lattice label -> (conventional formula, formula units per site,
#                           Li per site, S per site)
S8_UNIT_MAP: Dict[str, Tuple[str, int, int, int]] = {
    "S8":     ("S8",    1,  0, 8),
    "Li2S8":  ("Li2S8", 1,  2, 8),
    "Li4S8":  ("Li2S4", 2,  4, 8),
    "Li8S8":  ("Li2S2", 4,  8, 8),
    "Li16S8": ("Li2S",  8, 16, 8),
}
# legacy single-site cascade labels (cases/fullcell_shuttle)
LEGACY_CATHODE = ("Li2S6", "Li2S4", "Li2S2", "Li2S")


def classify(label: str) -> str:
    if label in METAL:
        return "metal"
    if label in ELECTROLYTE:
        return "electrolyte"
    if label.endswith("_an"):
        return "deposit"
    if label.startswith("CEI"):
        return "cei"
    if label in S8_UNIT_MAP or label in LEGACY_CATHODE or _CATHODE_RE.match(label):
        return "cathode"
    return "sei"


def classify_array(labels: Iterable[str]) -> np.ndarray:
    labels = np.asarray(list(labels), dtype=str)
    uniq = np.unique(labels)
    lut = {u: classify(str(u)) for u in uniq}
    return np.array([lut[str(l)] for l in labels], dtype=object)


def is_anode_film(labels: Iterable[str]) -> np.ndarray:
    """SEI + anode deposits: the film that grows on the Li anode."""
    cls = classify_array(labels)
    return (cls == "sei") | (cls == "deposit")


def li_s_per_site(label: str) -> Tuple[int, int]:
    """(Li, S) atoms represented by one lattice site with this label."""
    if label in S8_UNIT_MAP:
        _, _, li, s = S8_UNIT_MAP[label]
        return li, s
    m = re.match(r"^Li(\d*)S(\d*)", label)
    if m:
        return (int(m.group(1)) if m.group(1) else 1,
                int(m.group(2)) if m.group(2) else 1)
    if label == "S8":
        return 0, 8
    return 0, 0


def conventional_name(label: str) -> str:
    """Report-friendly name: Li4S8 -> Li2S4, Li2S2_an -> Li2S2 (anode)."""
    if label in S8_UNIT_MAP:
        return S8_UNIT_MAP[label][0]
    if label.endswith("_an"):
        return label[:-3] + " (anode deposit)"
    return label
