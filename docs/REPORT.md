# kmc3d: technical report (state at milestone 4, 2026-09-29)

Author: Carlos Guerrero (Balbuena group, Texas A&M, Battery500).
Scope: the Python engine `kmc3d`, its relation to the group's original C++
kinetic Monte Carlo (kMC) code, the changes made between 2026-09-04 and
2026-09-29, the results obtained so far, the known limitations and the
applications that the engine now makes possible. Every number quoted here is
traceable to `docs/MODEL_NOTES.md` (results, numbered items 0000 to 0007),
`docs/CHANGELOG.md` (engine changes) and `docs/KINETICS_TABLE.md`
(parameters and references). This document is updated at every milestone.

## 1. What the engine is

`kmc3d` is a rejection-free (BKL / Gillespie) kinetic Monte Carlo model of a
lithium-sulfur cell on a fixed superlattice. Every lattice site holds one
"species": Li metal, an electrolyte unit (ether background ETH, solvent SOL,
salt anion FSI), a solid electrolyte interphase (SEI) fragment, a cathode
unit (S8 with 0 to 16 Li) or a cathode electrolyte interphase (CEI) unit.
Events are lattice moves and reactions with rates; the algorithm picks one
event per step with probability proportional to its rate and advances the
clock by an exponentially distributed waiting time. A cycling protocol
alternates charge and discharge half-cycles.

The model has two purposes that must be kept apart: (a) capacity, Coulombic
efficiency and cell life with rates taken from the literature, and (b) the
morphology of the SEI (thickness, porosity, buried Li), which needs an
explicit acceleration of the SEI-forming reactions because a nanometre box
moves in one "cycle" about 1/4000 of the charge of a real 1 mAh/cm2 cycle
(section 7).

## 2. The original engine (antecedent)

The group's C++ code models a Li-metal half cell: a Li slab, ether
electrolyte with randomly placed solvent and salt units, plating and
stripping at the slab surface, decomposition of salt and solvent at the Li
surface into SEI fragments, diffusion of the fragments (MOBILITY.in, group
DFT values), and a semi-infinite Li foil that refills the slab from below.
Li within a cut-off distance of an SEI fragment is flagged "ionised"
(bound in LiF, Li2O, ...) and can no longer be stripped. Half-cycles last a
fixed number of events and the run ends when the electrolyte or the free Li
surface is consumed, typically around 100 cycles.

Its plating and decomposition rates are placeholders of the same order
(1e-4 : 1e-4), which gives about one decomposition per plating event. That
is the SEI-morphology regime: an implicit acceleration of roughly 500 on the
side reactions relative to a cell with Coulombic efficiency 99.8 %. It has
no cathode, no mass balance on sulfur, and its voltage is a label with no
effect on any rate.

`kmc3d` reproduces the C++ code byte for byte on its own input files
(`validation/anode_small`) and every addition is behind an opt-in keyword,
so with the keywords absent the two engines are the same model.

## 3. What was added (all opt-in; keywords in parentheses)

Full cell. A cathode built from S8 units on one lattice layer, each unit
storing 8 S and 0 to 16 Li (S8, Li2S8, Li4S8 = 2 Li2S4, Li8S8 = 4 Li2S2,
Li16S8 = 8 Li2S). Reduction consumes Li+ from a shared pool, oxidation
returns it. Li and S are conserved exactly at every step (`li_total`,
`s_total` columns). (MECHANISM.in REGION cathode, CONVERT, CONSUME_LI,
RELEASE_LI; GEOMETRY.in cathode_*.)

Shared Li+ pool and finite Li inventory (`li_pool_mode shared`,
`li_bulk_init`). Plating, stripping and the cathode exchange Li+ through one
counter; the semi-infinite foil can be switched off or made finite.

Dissolved polysulfides and shuttle. A well-mixed reservoir of Li2S8_d,
Li2S6_d, Li2S4_d with dissolution, precipitation, disproportionation and
corrosion of Li0 at the anode surface (REGION anode_surface, DISSOLVE,
PRECIP, REQUIRE, RATE_SCALE, STRIP_LI0, DEPOSIT).

CEI. Reduction of FSI by long-chain polysulfides at the cathode surface into
an inert film species with its own S and Li accounting (REGION
cathode_surface, CEI, S_LOSS, LI_LOSS), and dynamic passivation of cathode
sites by Li2S / Li2S2 / film neighbours (cathode_passivating_species,
cathode_passivation_nmin).

Physically based kinetics. Butler-Volmer plating and stripping with the
exchange current density of Boyle 2020 and an overpotential that sets the
current (`anode_kinetics bv`); the cathode cascade evaluated at the cathode
potential with the standard potentials of Kumaresan 2008 and exchange
current densities of Marinescu 2016 (`cathode_kinetics bv`); a
galvanostatic mode in which the cathode potential is solved at every scan
so that one current flows through both electrodes (`cathode_kinetics
galvanostatic`, cathode_V_min / max), which yields a voltage profile per
half-cycle. Electrolyte-reduction rates anchored to the experimental
Coulombic efficiency of LiFSI / F5DEE (Yu 2022) by the Li lost per Li plated
(KINETICS_TABLE.md section 7b), with the DFT barrier of Tan 2024 as the
relative attenuation of the F5DEE-on-Li2O route.

Protocol. Start with a discharge (`first_half end`), idle cut-off
(`end_half_when_idle`), current-efficiency cut-off (`end_half_min_ce`, the
analogue of the end-of-charge voltage cut-off), structural brakes
(`li_pool_max`, `stop_electrolyte_fraction`, `stop_anode_inactive_halves`).

Robustness. Clean stall detection with a diagnostic line, atomic writes of
the per-half-cycle ledger, ledger recovery from the checkpoint side file,
SIGTERM handling and a SLURM script that waits for the process before the
job ends; restart from checkpoint.

Speed. Incremental neighbour counts, edge gathering through a CSR view,
one-pass electrolyte refresh, vectorised lattice construction: 7x faster
than the first Python version on the physical cases, byte-identical output;
a 640000-site box builds in 20 s and runs at about 1 s per event.

Post-processing (`kmc3d/` modules and `postprocess/` scripts). Species
classification; per-cycle capacity, Coulombic efficiency (four definitions:
CE_cycle, CE_mod with a configurable electron-channel list, CE_cathode,
CE_shuttle), sulfur inventory and multi-seed statistics; four-panel overview
figure; lattice morphology (connected Li components, buried Li and its
coating, film extent, roughness, porosity profile); Gaussian density fields
with per-species radii; comparison of runs and of the two engines; Zeo++ and
RASPA input writers.

## 4. How every change was validated

Two reference cases (`validation/anode_small`, a legacy half cell, and
`validation/fullcell_small`, a small full cell with every keyword active)
have stored md5 sums of Data2Excel.txt, cycle_stats.csv and every xyz
frame. Every engine change is accepted only if both cases reproduce the
sums byte for byte; the log of validations is `validation/README.md`
(13 entries). The physical cases are additionally checked on three bench
runs (half cell, accelerated half cell, small-box full cell) for the speed
and trajectory-output changes. Conservation of Li and S is checked on every
production run (`analysis/conservation.txt`).

## 5. Parameters and references (summary; full table in KINETICS_TABLE.md)

Every rate is `globalA sigma k0 exp(-(Ea - alpha (V - E0)) / kT)` with tags
[V] verified in a PDF in hand, [V-assumed] taken from a source that states it
as an assumption, [V-derived] derived from a verified quantity, [P]
placeholder, [G] group value never modified.

- Li plating / stripping: j0 = 29.8 mA/cm2 (LiFSI in DME, Boyle 2020,
  Table 1) -> 74 s^-1 per site; alpha 0.5; eta +-30 mV [V, P for eta].
- Electrolyte reduction (FSI, SFO, SOL): anchored to CE 99.74 to 99.90 %
  (Yu 2022, LiFSI / F5DEE) by Li lost per Li plated, with the candidate
  multiplicities of the box [V-derived]. F5DEE on Li2O: barrier 0.364 eV
  (Tan 2024) as relative attenuation [V-derived].
- Cathode: E0 2.39, 2.37, 2.24, 2.04, 2.01 V (Kumaresan 2008, Table II);
  k0 from iH,0 10 and iL,0 5 A/m2 (Marinescu 2016, Table 1); alpha 0.5
  [V-assumed]; cycling window 1.7 to 2.8 V [V-assumed].
- Shuttle: ks 0.53 h^-1 at 298 K (Mikhaylik 2004) [V]. Precipitation and
  dissolution constants: Marinescu 2016 / Zhang 2016 [V-assumed].
- Surface diffusion and pair interactions: MOBILITY.in [G].
- Lattice criteria without a literature value (passivation threshold,
  deposit rules, reservoir normalisation): [P], listed as such.

The Li2Sx steps are the weakest link: their standard potentials and exchange
currents come from cell-level fits in DOL / DME, not F5DEE, and the
conversion to per-site rates is an assumption. They are the natural target
for the group's DFT / AIMD work (section 9).

## 6. Results so far

Run 1 (2026-09-23). Physical half cell: side reactions 0.4 to 0.6 % of
platings, stable; but the semi-infinite foil filled the box in 4 of 5 seeds.
Full cell: capacity 560 to 610 mAh/g_S for 5 to 9 cycles, then the anode
buried by the F5DEE-on-Li2O reaction (a transition-state prefactor is not an
effective electrode rate). Fixes: finite anode, SOL2 anchoring, stall
detection, atomic ledger, SLURM wait.

Run 2 (2026-09-24). 550 to 600 mAh/g_S for 13 cycles, CE 1.00, then death by
ionised Li: the SEI took 15 to 40 % of the events at the end of charge
because a trickle of cathode release kept charge halves alive. Fix:
current-efficiency cut-off.

Run 3 (2026-09-27). Life 3x run 2 (50 % capacity at cycle 57); side events
still 2 to 3x the half-cell anchor because at a fixed cathode potential the
last oxidation step (Li2S8 -> S8) is slow and the charge is supply-limited.
Fix: galvanostatic cathode potential.

Run 4 (2026-09-29). Voltage profile 2.47 to 2.66 V (charge) and 2.06 to
2.31 V (discharge), pool populated, side events 0.9 % of platings; gradual
fade with a large seed spread. Measurement: 2.7 fragments per reduction and
2.1 Li ionised per fragment, i.e. 0.034 Li lost per Li plated against 0.002
at CE 99.8 %. Fix: side rates divided by 17 (anchoring by Li lost).

Run 5 (2026-09-29), milestone 4. 100 cycles at 550 to 600 mAh/g_S with no
fade in 5 seeds, CE 1.00, Li lost per Li plated 0.0024, SEI 0.8 sites per
cycle, 4 % of S in the CEI, 21 % of S dissolved (stable). Figure:
`cases/fullcell_physical/analysis/overview.png`.

Accelerated half cell (A_SEI 100, 2026-09-27). SEI of 300 to 800 sites,
film porosity 0.86 to 0.88, roughness 3 to 4 A, 1 to 3 % of Li buried in
7 to 27 islands (`cases/anode_physical_A100/analysis/morphology.csv`).

## 7. Comparison with the original engine

| aspect | original C++ | kmc3d |
|---|---|---|
| system | Li half cell | Li-S full cell (half cell still available) |
| Li inventory | semi-infinite foil | finite slab, finite or infinite foil, shared Li+ pool |
| cathode | none | S8-unit lattice, 5 lithiation states, exact S balance |
| polysulfides | none | dissolved reservoir, shuttle, precipitation, CEI |
| voltage | label, no effect | sets every electrode rate; galvanostatic mode gives a profile |
| rates | proportional placeholders (1e-4 : 1e-4) | literature values with references and tags |
| SEI per plating | ~1 (accelerated ~500x) | 0.002 Li lost per Li plated (CE anchor), or A_SEI x that |
| cycle life | ~100 kMC cycles until the box is consumed | > 100 cycles without fade at the CE anchor |
| conservation | not enforced | Li and S exact, checked every run |
| termination | electrolyte / surface consumed | same, plus stall diagnostic, cut-offs, checkpoints |
| speed | compiled, fast per event | Python, ~0.2 s per event (10x10x30), ~1 s (20x20x100) |
| output | Data2Excel, xyz | same plus cycle_stats ledger, POSCAR snapshots, analysis tables and figures |
| validation | none stored | byte-identical reference cases, 13 logged validations |

What the original does better: raw speed per event and simplicity. What it
cannot do: predict capacity, Coulombic efficiency or the effect of voltage,
keep a mass balance, or represent the cathode.

Scale. One kMC cycle in the 4 x 4 nm box moves about 5 atomic layers of Li;
a real 1 mAh/cm2 cycle moves about 20000. The acceleration factor A_SEI is
therefore the ratio of the real charge per cycle to the box charge per
cycle (about 4000 for 1 mAh/cm2): with it, one kMC cycle forms the SEI of a
fraction of a real cycle. A realistic SEI thickness per real cycle (about
10 nm at CE 99.8 %) needs a box tens of nanometres tall and represents 2 to
3 real cycles; the 8 x 8 x 40 nm case `fullcell_tall` (A_SEI 100, finite foil)
is that box.

## 8. Known limitations (state of the model, not bugs)

- Li next to an SEI fragment cannot be stripped and plating needs an
  electrolyte neighbour: there is no Li+ transport through the SEI. The rule
  is a proxy for Li bound in SEI compounds and for buried Li; with the side
  rates anchored by Li lost it no longer limits life within 100 cycles, but
  a through-SEI transfer with an attenuation exp(-beta d) remains the
  physically complete option (decision pending with the advisor).
- The voltage is constant within a half-cycle (galvanostatic mode moves the
  cathode potential between scans, not the anode overpotential).
- The cathode has no volume change and one lattice layer; about 30 % of its
  sites stay lithiated after a charge behind the passivation mask.
- The utilisation is capped by the event budget per half-cycle (35 % at 800
  events in the 10x10x30 box); raising the budget is safe with the cut-off.
- The CEI is one species, one layer, no dissolution.
- Every event is a lattice site: molecular sizes enter only through the
  post-processing (Gaussian density with per-species radii).
- The Li2Sx kinetics are borrowed from DOL / DME cell models.

## 9. Applications now within reach

- Capacity and Coulombic efficiency versus cycle for a given electrolyte,
  with the electrolyte entering through its decomposition mechanism: a
  solvent with a known mechanism is a change of DECOMPOSITION.in and
  MECHANISM.in, and the CE becomes an output.
- Voltage and current as variables: overpotential (C-rate) sweeps, cycling
  window sweeps, cathode potential profiles.
- SEI morphology in the accelerated mode: thickness, porosity, inorganic
  fraction, buried Li and its coating, with the Gaussian density giving the
  fragments a size; Zeo++ / RASPA on the POSCAR snapshots.
- Geometry: anode thickness, N/P ratio, cathode loading, imperfections.
- Parameter sweeps and ML / SHAP on the per-cycle tables.
- Parameters for continuum (DFN) models: exchange current versus SEI
  thickness, SEI yield per Ah, film porosity, precipitation and shuttle
  constants.

## 10. Roadmap

1. Run 6 of `fullcell_physical`: 500 cycles (maxCycles 1000).
2. `fullcell_tall`: the SEI-morphology movie (A_SEI 100, 8 x 8 x 40 nm,
   about 1 h per cycle, checkpoints every 10 half-cycles).
3. Advisor decisions: through-SEI Li transport; DFT / AIMD of the Li2Sx
   steps on the cathode host and of polysulfide solvation in F5DEE.
4. Full utilisation (larger event budget), cathode re-oxidation deficit,
   variable anode overpotential (CC-CV).
5. Sweeps, ML / SHAP, DFN parameter extraction, SPAN cathode (optional).
