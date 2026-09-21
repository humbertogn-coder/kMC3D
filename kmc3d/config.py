"""
config.py
=========
Reads the scalar / tabular input files of the generalized 3D kMC:

  PARAMETERS.in   - global scalars (seed, voltages, steps, T, molarities ...)
  GEOMETRY.in     - box replication (nx,ny,nz), lattice constant, electrode
                    regions (anode / separator / cathode) -> replaces the
                    manual OVITO + POSCAR-expansion step.
  DECOMPOSITION.in- Arrhenius KINETICS for each decomposition reaction type
                    (sigma,k0,Ea,alpha,E0), one or more rate channels per type.
  MOBILITY.in     - Arrhenius KINETICS for mobility reactions (Diffusion,
                    LiStripping, LiSurface) + pairwise interaction energies
                    (ICOHP-derived) used to build the environment-dependent
                    diffusion barrier.
  MECHANISM.in    - PRODUCTS (stoichiometry) of every reaction type.  This is
                    the file that was previously hard-coded inside
                    electrolyteDecom() in calculate_reactions.cpp.

The split DECOMPOSITION (kinetics) vs MECHANISM (products) mirrors the C++
exactly: a type may expose several FRM rate-candidates while firing a single
(possibly stochastic) product channel.

All parsers are dependency-free (stdlib only) so they run on GRACE without a
conda environment.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List
import os


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
def _strip_comments(line: str) -> str:
    for c in ("#", "!", ";"):
        i = line.find(c)
        if i >= 0:
            line = line[:i]
    return line.strip()


def _tokens(path: str) -> List[str]:
    """Flatten a whitespace-delimited file into a token stream (comments removed)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")
    toks: List[str] = []
    with open(path) as fh:
        for raw in fh:
            line = _strip_comments(raw)
            if line:
                toks.extend(line.split())
    return toks


# ----------------------------------------------------------------------------
# PARAMETERS.in
# ----------------------------------------------------------------------------
@dataclass
class Parameters:
    seed: int = 12345
    temperature: float = 298.0          # K
    saltMolar: float = 1.0              # mol/L
    solvMolar: float = 1.0              # mol/L
    scanIntervalType: str = "time"      # only "time" implemented (matches C++)
    potentialType: str = "step"         # only "step" implemented
    scanInterval: float = 1.0e-3        # half-cycle duration (time units)
    maxInterval: float = 1.0e6          # max steps per half cycle
    BeginV: float = -0.1                # plating / reduction potential
    EndV: float = 0.1                   # stripping / oxidation potential
    totalSteps: int = 100000
    maxCycles: int = 200
    printStepFreq: int = 0
    XYZprintFreq: int = 100             # write trajectory frame every N steps
    snapshotEveryCycles: int = 0        # dump a POSCAR snapshot every N half-cycles
                                        #   (0 = off) for Zeo++/RASPA post-processing
    checkpointEveryCycles: int = 0      # save a restartable checkpoint every N
                                        #   half-cycles (0 = off)
    globalA: float = 1.0                # global Arrhenius prefactor
    # number-of-entries are auto-derived from the files; kept for compatibility
    numberDecompositionRxns: int = 0
    nMobilityReactions: int = 0
    # ---- Li-S full cell: shared Li pool + well-mixed polysulfide reservoir
    #   li_pool_mode  fixed  -> legacy (plating/cathode Li are pure counters;
    #                           anode-only runs stay byte-identical)
    #                 shared -> one Li+ pool: plating/FSI/cathode reduction
    #                           debit it, stripping/cathode oxidation credit it;
    #                           a reaction needing n Li+ is ineligible if the
    #                           pool holds fewer, and plating scales with
    #                           pool/pool0.
    li_pool_mode: str = "fixed"
    li_pool_init: int = -1              # -1 -> saltMolar * electrolyte volume
    reservoir_sites: int = -1           # normalisation volume (in lattice
                                        #   sites) for dissolved-species
                                        #   concentrations; -1 -> n_ETH at t=0
    #   li_bulk_init  anode bulk-metal reservoir (Li foil below the lattice).
    #                 In shared mode the anode framework (wrap / updateLiMetal)
    #                 draws its Li from here instead of creating it, and
    #                 returns removed floating Li here. -1 = unlimited foil
    #                 (accounting only, li_bulk goes negative = net drawn);
    #                 N >= 0 = finite foil, framework growth stops at 0
    #                 (anode depletion; sets an effective N/P ratio).
    li_bulk_init: int = -1
    #   li_pool_max   upper bound of the shared Li+ pool (shared mode only).
    #                 -1 = unlimited (legacy). N > 0: stripping and cathodic
    #                 oxidation (RELEASE_LI) become ineligible once the pool
    #                 holds N Li+. Mean-field electroneutrality / salt
    #                 solubility limit: forces the anodic current to follow the
    #                 cathodic one (a real full cell has one current). Without
    #                 it the pool ran away (74 -> 4500 Li+ in fullcell_cei,
    #                 2026-09-20) and the anode aged ~8x faster per cycle.
    #                 Suggested: 2 * li_pool0 (supersaturation margin).
    #                 MODEL LIMITATION: a hard cap, not a Nernstian penalty.
    li_pool_max: int = -1
    #   stop_electrolyte_fraction  clean termination when the electrolyte sites
    #                 (ETH + SOL + FSI) fall below this fraction of their
    #                 initial number ("cell dried out"). 0 = off (legacy).
    #                 Mirrors the natural stall of the original half-cell
    #                 engine, which the full cell never reaches on its own.
    stop_electrolyte_fraction: float = 0.0

    # field name -> caster
    _CAST = {
        "seed": int, "temperature": float, "saltMolar": float, "solvMolar": float,
        "scanIntervalType": str, "potentialType": str, "scanInterval": float,
        "maxInterval": float, "BeginV": float, "EndV": float, "totalSteps": int,
        "maxCycles": int, "printStepFreq": int, "XYZprintFreq": int,
        "snapshotEveryCycles": int, "checkpointEveryCycles": int,
        "globalA": float, "numberDecompositionRxns": int, "nMobilityReactions": int,
        "li_pool_mode": str, "li_pool_init": int, "reservoir_sites": int,
        "li_bulk_init": int, "li_pool_max": int,
        "stop_electrolyte_fraction": float,
    }

    @classmethod
    def read(cls, path: str) -> "Parameters":
        toks = _tokens(path)
        p = cls()
        i = 0
        while i < len(toks):
            key = toks[i]
            if i + 1 >= len(toks):
                break                       # trailing key with no value -> keep default
            if key not in cls._CAST:
                # unknown token: skip one value defensively
                i += 2
                continue
            val = toks[i + 1]
            setattr(p, key, cls._CAST[key](val))
            i += 2
        return p


# ----------------------------------------------------------------------------
# GEOMETRY.in  (replaces manual POSCAR expansion in OVITO)
# ----------------------------------------------------------------------------
@dataclass
class Geometry:
    """
    Programmatic supercell.  The conventional cubic cell (side = latticeConstant)
    hosts the BCB/BCO/BA/OC/TE superlattice that can simultaneously accommodate
    BCC Li, rock-salt LiF and antifluorite Li2O.  The cell is replicated
    nx*ny*nz times - no OVITO needed.

    Electrode regions are specified as fractional z-windows so the same engine
    handles a half cell (anode only) or a full Li||S cell (anode + separator +
    sulfur cathode).
    """
    nx: int = 8
    ny: int = 8
    nz: int = 24
    latticeConstant: float = 4.0        # Angstrom (conventional cube side)

    # electrode geometry (all as FRACTION of total box height in z)
    anode_center: float = 0.5           # matches original (slab centred at 0.5)
    anode_thickness: float = 0.10       # original used anodeThickness (Angstrom);
                                        # here expressed as a z-fraction for generality
    anode_thickness_is_fraction: bool = True
    anode_thickness_A: float = 8.0      # used only if *_is_fraction == False
    anode_species: str = "Li"           # base anode metal (e.g. Li, or Zn)
    anode_vacancy_fraction: float = 0.0  # fraction of anode BC sites that are holes
    anode_vacancy_species: str = ""     # "" -> empty hole; else fill holes with this
                                        #   species (e.g. "Zn" for a Li/Zn mixed anode)

    # optional sulfur cathode (disabled by default; primary goal is the anode model)
    cathode_enabled: bool = False
    cathode_species: str = "S"          # simple cathode active material
    cathode_center: float = 0.90
    cathode_thickness: float = 0.10     # z-fraction
    cathode_charge: int = 0
    # static electronic/ionic accessibility mask (fraction of cathode sites
    # that can undergo conversion; 1.0 = all, legacy behaviour). Drawn once
    # at cell creation from the engine RNG only when < 1.0.
    cathode_access_fraction: float = 1.0
    # dynamic electronic passivation by insulating discharge products (opt-in,
    # point 3). A cathode-region site whose neighbour shell (same shells as the
    # nLi count: BCB, BCO, BA, OC) holds >= cathode_passivation_nmin sites of
    # any species in cathode_passivating_species is electronically blocked:
    # it cannot fire reactions that transfer Li+ (CONSUME_LI / RELEASE_LI).
    # Dissolution / precipitation (no electron transfer) are NOT blocked.
    # cathode_passivation_nmin 0 (default) disables the feature completely
    # (legacy code path, byte-identical output).
    #   cathode_passivating_species  comma-separated, e.g. "Li2S,Li2S2"
    #   cathode_passivation_nmin     integer threshold (0 = off)
    # MODEL LIMITATION: mean-field electronic isolation by neighbour count; no
    # explicit electron-conduction network (percolation) is solved, and the
    # threshold is a placeholder to be calibrated (Li2S is a wide-gap insulator,
    # Li2S2 is less so; both are lumped as "blocking" here).
    cathode_passivating_species: str = ""
    cathode_passivation_nmin: int = 0

    # POSCAR i/o
    write_poscar: bool = True           # dump the generated supercell as POSCAR
    poscar_in: str = ""                 # if set, read the lattice from this POSCAR
                                        # instead of generating it programmatically

    _CAST = {
        "nx": int, "ny": int, "nz": int, "latticeConstant": float,
        "anode_center": float, "anode_thickness": float,
        "anode_thickness_is_fraction": lambda s: s.lower() in ("1", "true", "yes"),
        "anode_thickness_A": float,
        "anode_species": str,
        "anode_vacancy_fraction": float,
        "anode_vacancy_species": str,
        "cathode_enabled": lambda s: s.lower() in ("1", "true", "yes"),
        "cathode_species": str, "cathode_center": float,
        "cathode_thickness": float, "cathode_charge": int,
        "cathode_access_fraction": float,
        "cathode_passivating_species": str,
        "cathode_passivation_nmin": int,
        "write_poscar": lambda s: s.lower() in ("1", "true", "yes"),
        "poscar_in": str,
    }

    @classmethod
    def read(cls, path: str) -> "Geometry":
        toks = _tokens(path)
        g = cls()
        i = 0
        while i < len(toks):
            key = toks[i]
            if i + 1 >= len(toks):
                break                       # trailing key with no value -> keep default
            if key not in cls._CAST:
                i += 2
                continue
            setattr(g, key, cls._CAST[key](toks[i + 1]))
            i += 2
        return g


# ----------------------------------------------------------------------------
# DECOMPOSITION.in   (kinetics only)
# ----------------------------------------------------------------------------
# Format (whitespace, free-form), repeated per type:
#     <TYPE> <nRateChannels>
#         <sigma> <k0> <Ea> <alpha> <E0>     # channel 0
#         ...                                # channel nRateChannels-1
@dataclass
class DecompositionKinetics:
    # type -> list of dict(sigma,k0,Ea,alpha,E0)
    channels: Dict[str, List[Dict[str, float]]] = field(default_factory=dict)

    def n_channels(self, rtype: str) -> int:
        return len(self.channels.get(rtype, []))

    @classmethod
    def read(cls, path: str, n_reactions: int = -1) -> "DecompositionKinetics":
        toks = _tokens(path)
        d = cls()
        i = 0
        read = 0
        while i + 1 < len(toks):
            if n_reactions >= 0 and read >= n_reactions:
                break                          # stop after the declared count
            rtype = toks[i]
            try:
                n = int(toks[i + 1])
            except ValueError:
                break                          # hit trailing notes/prose -> stop
            i += 2
            rows = []
            ok = True
            for _ in range(n):
                try:
                    sigma, k0, Ea, alpha, E0 = (float(toks[i + k]) for k in range(5))
                except (ValueError, IndexError):
                    ok = False
                    break
                i += 5
                rows.append(dict(sigma=sigma, k0=k0, Ea=Ea, alpha=alpha, E0=E0))
            if not ok:
                break
            d.channels[rtype] = rows
            read += 1
        return d


# ----------------------------------------------------------------------------
# MOBILITY.in   (kinetics + pair interactions)
# ----------------------------------------------------------------------------
# Format:
#     <RxnType> <nSpecies>
#         <species> <sigma> <k0> <alpha> <E0>
#         ...
#     ... (repeated nMobilityReactions times)
#     INTERACTIONS <nPairs>
#         <spcsA> <spcsB> <Eb>
#         ...
@dataclass
class MobilityKinetics:
    # rxn -> species -> dict(sigma,k0,alpha,E0)
    params: Dict[str, Dict[str, Dict[str, float]]] = field(default_factory=dict)
    # interaction[A][B] = Eb
    interaction: Dict[str, Dict[str, float]] = field(default_factory=dict)

    def get(self, rxn: str, species: str) -> Dict[str, float]:
        return self.params.get(rxn, {}).get(species, {})

    def Eb(self, a: str, b: str) -> float:
        return self.interaction.get(a, {}).get(b, 0.0)

    @classmethod
    def read(cls, path: str, n_mobility_reactions: int) -> "MobilityKinetics":
        toks = _tokens(path)
        m = cls()
        i = 0
        for _ in range(n_mobility_reactions):
            rxn = toks[i]; n = int(toks[i + 1]); i += 2
            m.params[rxn] = {}
            for _ in range(n):
                spcs = toks[i]
                sigma, k0, alpha, E0 = (float(toks[i + 1 + k]) for k in range(4))
                i += 5
                m.params[rxn][spcs] = dict(sigma=sigma, k0=k0, alpha=alpha, E0=E0)
        # interactions block
        if i < len(toks):
            _label = toks[i]; npairs = int(toks[i + 1]); i += 2
            for _ in range(npairs):
                a, b = toks[i], toks[i + 1]; eb = float(toks[i + 2]); i += 3
                m.interaction.setdefault(a, {})[b] = eb
        return m
