"""
raspa.py - thermodynamic characterisation of kMC snapshots with RASPA3 (raspalib).

This is a POST-PROCESSING / frozen-snapshot tool, NOT a coupled loop: the kMC
structure is held rigid and RASPA runs Widom test-particle insertions to get,
at infinite dilution:
    Henry coefficient        (adsorption affinity)
    heat of adsorption Qst   (binding strength)

Install with:  pip install raspalib

IMPORTANT practical note
------------------------
Full-cell snapshots are dominated by DENSE Li metal, which has essentially no
probe-accessible volume. Characterise a REPRESENTATIVE SEI SLAB instead (use
`z_slab=(zmin,zmax)` to keep only the decomposition layer between the metal and
the bulk electrolyte). On a dense block the Henry coefficient collapses to ~0,
which is itself a meaningful "non-porous" signal rather than an error.
"""

from __future__ import annotations
import math
from typing import Dict, Optional, Tuple

import numpy as np

from .structio import Structure, read_poscar
from .ff_data import SPECIES_PROPS, PROBES


def _slab(struct: Structure, z_slab: Optional[Tuple[float, float]]) -> Structure:
    if z_slab is None:
        return struct
    lo, hi = z_slab
    z = struct.cart[:, 2]
    keep = (z >= lo) & (z <= hi)
    box = struct.box.copy()
    box[2, 2] = max(hi - lo, 1.0)
    cart = struct.cart[keep].copy()
    cart[:, 2] -= lo
    return Structure(box=box, species=struct.species[keep], cart=cart)


def widom_insertion(poscar_path: str,
                    probe: str = "He",
                    temperature: float = 298.0,
                    n_cycles: int = 5000,
                    n_init: int = 500,
                    z_slab: Optional[Tuple[float, float]] = None) -> Dict[str, float]:
    """Run a RASPA3 Widom insertion on one frozen snapshot.

    Returns {henry_mol_kg_Pa, qst_kJ_mol, rosenbluth, accessible}. If the
    structure is non-porous to the probe, rosenbluth ~ 0 and henry -> 0 with
    accessible=0 (reported, not raised).
    """
    try:
        import raspalib
    except ImportError as e:
        raise RuntimeError("raspalib not installed. `pip install raspalib`.") from e

    struct = _slab(read_poscar(poscar_path), z_slab)
    uniq = sorted(set(struct.species.tolist()))
    pp = PROBES[probe]

    pseudo = [raspalib.PseudoAtom(name=s, framework_type=True,
                                  mass=SPECIES_PROPS[s].mass, charge=0.0,
                                  atomic_number=6) for s in uniq]
    pseudo.append(raspalib.PseudoAtom(name=probe, framework_type=False,
                                      mass=pp["mass"], charge=pp["charge"],
                                      atomic_number=2))
    params = [raspalib.VDWParameters(SPECIES_PROPS[s].eps_K,
                                     SPECIES_PROPS[s].sigma_A) for s in uniq]
    params.append(raspalib.VDWParameters(pp["eps_K"], pp["sigma_A"]))

    ff = raspalib.ForceField(
        pseudo_atoms=pseudo, parameters=params,
        mixing_rule=raspalib.ForceField.MixingRule.LORENTZ_BERTHELOT,
        cutoff_framework_vdw=12.0, cutoff_molecule_vdw=12.0,
        cutoff_coulomb=12.0, shifted=True, tail_corrections=False,
        use_charge=False)

    type_of = {s: i for i, s in enumerate(uniq)}
    a, b, c = struct.lengths
    frac = struct.frac % 1.0
    atoms = [raspalib.Atom(position=tuple(frac[i]), charge=0.0,
                           type=type_of[struct.species[i]])
             for i in range(struct.species.size)]
    framework = raspalib.Framework(
        force_field=ff, component_name="kmc_snap",
        simulation_box=raspalib.SimulationBox(a=a, b=b, c=c),
        space_group_hall_number=1, defined_atoms=atoms,
        number_of_unit_cells=[1, 1, 1])

    comp = raspalib.Component(
        force_field=ff, component_name=probe,
        critical_temperature=5.2, critical_pressure=2.3e5,
        acentric_factor=-0.39,
        defined_atoms=[raspalib.Atom(position=(0, 0, 0),
                                     charge=pp["charge"], type=len(uniq))],
        particle_probabilities=raspalib.MCMoveProbabilities(widom_probability=1.0),
        number_of_blocks=5)

    system = raspalib.System(
        force_field=ff, external_temperature=temperature, external_pressure=1e4,
        helium_void_fraction=1.0, framework_components=framework,
        components=[comp], initial_number_of_molecules=[0])

    mc = raspalib.MonteCarlo(
        number_of_cycles=n_cycles, number_of_initialization_cycles=n_init,
        number_of_equilibration_cycles=0, print_every=max(n_cycles, 1),
        systems=[system])
    mc.run()

    # results populate the system's COPY of the component, not our handle
    run_comp = mc.systems[0].components[0]
    run_sys = mc.systems[0]
    rosen, _err = run_comp.average_rosenbluth_weights.result()
    accessible = 1.0 if (rosen == rosen and rosen > 1e-9) else 0.0  # nan check
    henry = float(rosen) if accessible else 0.0
    qst = float("nan")
    if accessible:
        try:
            ea = run_sys.average_enthalpies_of_adsorption.result()
            if ea:
                qst = float(np.ravel(ea)[0])
        except Exception:
            pass
    return {"henry": henry, "qst_kJ_mol": qst,
            "rosenbluth": float(rosen) if rosen == rosen else 0.0,
            "accessible": accessible}
