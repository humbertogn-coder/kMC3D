# kmc3d model notes: Li-S full cell (engine state, keywords, placeholders, limits)

Last updated: 2026-09-17. Engine modules: config.py, mechanism.py, engine.py,
stats.py. Reference case: cases/fullcell_shuttle/.

## 1. What the engine does

Fixed-lattice 3D kinetic Monte Carlo (BKL) ported from the confidential C++
half-cell SEI engine (Perez-Beltran et al.) and generalized to a full Li-S cell:
Li anode at the centre of the box, S8 cathode wrapping both z faces (periodic
band). Reactions are tagged by REGION (anode / cathode / anode_surface) and
VOLTAGE (begin / end / any). The anode-only legacy path must stay byte-identical
after any change (see validation/README.md).

## 2. Input files (five)

PARAMETERS.in, GEOMETRY.in, DECOMPOSITION.in, MOBILITY.in, MECHANISM.in.
Golden rule: on any unexpected divergence, first verify the five .in files are
exactly the intended ones. Historically most "bugs" were misplaced parameters.

## 3. Opt-in keywords (absent keyword = legacy path)

PARAMETERS.in
    li_pool_mode     fixed | shared   (default fixed = legacy counters)
    li_pool_init     -1 -> saltMolar * electrolyte volume (74 Li+ in the 10x10x30 box)
    reservoir_sites  -1 -> number of ETH sites at t=0; normalises
                     c_x = N_dis[x] / reservoir_sites in RATE_SCALE
    li_bulk_init     -1 = unlimited Li foil (framework Li accounted, not capped);
                     N >= 0 models a finite Li excess (effective N/P ratio)
    li_pool_max      -1 = unlimited (legacy). N > 0: stripping and cathodic
                     oxidation are ineligible while the pool holds >= N Li+.
                     Mean-field electroneutrality: one current for both
                     electrodes. Suggested 2 * li_pool0. (structural brake 1)
    stop_electrolyte_fraction  0 = off. f > 0: the run stops cleanly when
                     ETH + SOL + FSI sites fall below f * initial ("cell dried
                     out"), the explicit form of the original engine's natural
                     stall. Suggested 0.05. (structural brake 2)
GEOMETRY.in
    cathode_access_fraction      1.0 = all cathode sites convert; < 1.0 draws a
                                 static contact mask at cell creation (point 4, static)
    cathode_passivating_species  comma-separated, e.g. Li2S,Li2S2 (point 3)
    cathode_passivation_nmin     0 = off; N > 0 blocks a cathode site whose
                                 neighbour shell holds >= N passivating sites.
                                 A cathode BC site has 8 BC neighbours.
    cathode_center 0.0 + anode_center 0.5: periodic bands (_zband helper)

MECHANISM.in keywords
    REGION anode_surface      live Li0 surface (same criterion as LiStripping)
    TRIGGER EMPTY             fires on empty BC sites of the region
    REQUIRE <spc_d> <n>       reservoir must hold >= n
    RATE_SCALE <spc_d>        rate *= N_dis[spc_d] / reservoir_sites
    DISSOLVE <spc_d>          empties the trigger site, reservoir += 1
    PRECIP <spc> [q]          fills the empty trigger site with <spc>
    RES <spc_d> <delta>       changes the reservoir
    STRIP_LI0 <n>             corrodes n surface Li0 (trigger + n-1)
    DEPOSIT <spc> <q> <n>     places n insoluble <spc> on BA sites of the anode
                              surface (SEI class, does not diffuse)
    CONVERT / CONSUME_LI / RELEASE_LI  cathode cascade with Li+ transfer
    REGION cathode_surface    electrolyte sites (ETH/SOL/FSI on OC/BA) with at
                              least one neighbour of cathode material (CEI zone)
    CEI <spc> [q]             converts the trigger electrolyte site into the
                              inert film species <spc> (excluded from the anode
                              SEI class; may be listed as passivating)
    S_LOSS <n>                n S atoms sequestered in the CEI (enters s_total)
    LI_LOSS <n>               n Li trapped in the CEI (li_cei column, reporting)

Dissolved species (reservoir only, never on the lattice): Li2S8_d, Li2S6_d,
Li2S4_d. Insoluble on the lattice: Li2S2, Li2S (cathode), Li2S2_an (anode deposit).

Stoichiometries (Li and S balanced; S balance closes exactly in tests):
    3 Li2S8 + 2 Li0 -> 4 Li2S6        2 Li2S6 + 2 Li0 -> 3 Li2S4
    Li2S4 + 2 Li0 -> 2 Li2S2(an)      irreversible: Li2Sx + 2 Li0 -> Li2S(x-2) + Li2S2(an)

## 4. Cathode passivation (point 3, partial)

A cathode-region site with >= cathode_passivation_nmin neighbours (BCB, BCO,
BA, OC shells, same edge set as the nLi count) occupied by a passivating
species is electronically blocked: it cannot fire reactions carrying
CONSUME_LI or RELEASE_LI. Dissolution and precipitation are not gated.
The mask is recomputed from the lattice state each step (no RNG, not stored
in the checkpoint). Column n_cathode_blocked appears in cycle_stats.csv only
when active. Still missing for point 3: CEI (electrolyte decomposition at the
cathode, REGION cathode).

## 5. Literature placeholders (DECOMPOSITION.in; sigma k0 Ea alpha E0)

Per-site rate = globalA * sigma * k0 * exp(-(Ea - alpha (V - E0))/kT); with
globalA 0.5 and Ea = alpha = 0 the rate is 0.5 * k0. DO NOT calibrate yet:
the engine must be complete and conservative first.

| Reaction          | k0   | Basis (PLACEHOLDER)                                              |
|-------------------|------|------------------------------------------------------------------|
| cat_diss_Li2S8/6  | 0.05 | S as Li2Sx > 10 M in ethers (Rauh 1977; Mikhaylik and Akridge, JES 151 (2004) A1969): fast high-order dissolution |
| cat_diss_Li2S4    | 0.02 | Li2S4 marginally soluble; Li2S2/Li2S have no dissolution channel |
| cat_prec_Li2Sx    | 5.0  | Precipitation constant k_p of the 0D model, Marinescu et al., JES 165 (2018) A6107. Calibrate the ratio k_diss/k_prec, not each value |
| sh_Li2S8/6        | 0.5  | Calibrate against the shuttle factor f_c = k_s q_H [S]/I (exp. 1.1 to 4) and target CE: 0.95 to 0.99 with LiNO3, 0.6 to 0.9 without |
| sh_Li2S4          | 0.2  | idem, low order less abundant in solution                        |
| p_irrev (weights) | 0.05 | Marinescu 2018: fraction of shuttled material lost irreversibly  |
| passivation nmin  | n/a  | not set in production yet; suggested first trial Li2S,Li2S2 with nmin 4 to 6 of 8 |

Observed orders of magnitude: k_sh 5/5/2 gives CE_shuttle ~0.33 and drains the
pool; 0.5/0.5/0.2 gives CE_shuttle ~0.44 to 0.49. With reservoir_sites 18211
re-precipitation never fires (the shuttle, with hundreds of Li surface sites,
beats 30 to 60 empty cathode sites 10:1); reservoir_sites 2000 or k_prec 50
activates it.

## 6. cycle_stats.csv columns (only when the feature is active)

    li_pool             Li+ in the pool at the end of the half-cycle
    li_bulk             foil reservoir
    li_total            CONSERVATION CHECK = li_pool + li_bulk + n_Li
                        + li_consumed - li_released + li_shuttled. Must be
                        constant over the whole run. Check it before reading any CE.
    d_li_plated_pool    Li0 placed through the pool in the half-cycle
    d_li_framework      Li added by wrap()/updateLiMetal() from the foil
    n_dis_<spc>         reservoir at the end of the half-cycle
    d_li_shuttled       Li0 corroded by the shuttle
    d_li_deposit        Li2S2_an sites deposited
    sei_Li2S2_an        anode deposit inventory
    n_cathode_blocked   electronically blocked cathode sites (passivation)
    s_total             SULFUR CONSERVATION CHECK = S on lattice cathode
                        species + S in anode deposits + S in the reservoir
                        + s_cei. Must be constant (= 8 * initial S8 sites).
    s_cei, n_CEI,       S sequestered in the CEI, CEI film sites, per-species
    cei_<spc>           CEI inventory (only when a CEI reaction exists)
    deposit_unplaced    DEPOSIT sites that found no eligible BA site (audit;
                        must stay 0, otherwise s_total drifts)

CE definitions (always present):

    CE_cycle   = stripped(discharge) / plated(charge)                  [legacy]
    CE_mod     = stripped(discharge) / [plated + FSI + SFO + SOL + SOL2 + F5D](charge)
    CE_shuttle = Li released by cathodic oxidation (charge) /
                 [that + Li0 corroded by the shuttle in charge and discharge]

Honest warning: CE_cycle and CE_mod can exceed 1 because stripping is a
per-site kinetic process, not capacity limited. They are count ratios, not
charge balances. CE_shuttle is an electron balance at the cathode and is the
one to calibrate first.

## 7. Output compatibility

Extended .xyz (same format, new symbols Li2S2_an and cathode species; dissolved
species are not sites and do not appear). POSCAR snapshots for Zeo++/RASPA over
the SEI z-slice (Li2S2_an is part of the skeleton). cycle_stats.csv keeps the
legacy columns; new ones are appended by name.

## 9. S8-unit lattice model of the cathode (cases/fullcell_cei, 2026-09-19)

Every cathode lattice site holds one S8 unit; the label is its lithiation
state. Post-processing must map labels to conventional formulas:

    lattice label   formula equivalent   Li per site   soluble
    S8              S8                    0            no
    Li2S8           Li2S8                 2            yes -> 1 Li2S8_d
    Li4S8           2 Li2S4               4            yes -> 2 Li2S4_d
    Li8S8           4 Li2S2               8            no  (passivating)
    Li16S8          8 Li2S               16            no  (passivating)

Cascade: 2 + 2 + 4 + 8 = 16 Li+ per S8 (theoretical 1672 mAh/g). Li2S6 is
not a lattice state (non-integer Li per S8); Li2S6_d exists in solution only,
from the shuttle (3 Li2S8 + 2 Li -> 4 Li2S6). Precipitation of Li2S6_d goes
through the disproportionation 2 Li2S6_d -> Li2S8(site) + Li2S4_d, exact in
Li and S. Known simplification: no direct Li2S6 dissolution from the cathode.
Rates remain placeholders (0.10/0.06/0.04/0.02 discharge, mirrored on charge).

## 8. Known limits (next increments)

00. ROOT CAUSE of the fast anode saturation (2026-09-20 analysis vs the
   original half-cell run): the two electrodes were not current-coupled.
   Stripping (k0 56 per surface Li0 site) outran cathodic reduction (k0 0.05
   per site) by ~1000x, the shared pool grew from 74 to ~4500 Li+, and the
   400-event half-cycles delivered ~8x more anode events per cycle than the
   original's 50-event half-cycles (cycle 20 here ~ cycle 140 there in
   wear). Brakes li_pool_max and stop_electrolyte_fraction address this
   structurally; the rate calibration remains pending.
0. First fullcell_cei production run (2026-09-20, placeholder rates): the
   shuttle consumes ~90 % of the S8 sites over 160 half-cycles, ~73 % of all
   sulfur ends as Li2S2_an on the anode, re-precipitation rarely fires, the
   surface FSI salt is exhausted (CEI self-limits) and the anode framework
   draws ~15x its initial Li from the foil (li_bulk about -15000; li_total is
   conserved). All of this is the calibration problem, not a bug: k_sh vs
   k_prec vs reservoir_sites must be set so that the cathode survives.

1. REQUIRE uses the max over channels (sh_Li2S8 demands 3 even though the
   irreversible channel uses 1): small bias toward the reversible channel at low N_dis.
2. CEI: only the FSI + long-chain polysulfide chemical path is modelled
   (cases/fullcell_cei). Solvent (F5DEE) path disabled for lack of evidence;
   no electrochemical oxidation of the electrolyte (outside the 1.7 to 2.8 V
   window of ethers). Rates are placeholders (docs/CEI_LITERATURE.md).
3. Accessibility only as a static mask; passivation is mean field by neighbour
   count, no electronic percolation.
4. SULFUR BALANCE. The legacy single-site cascade of cases/fullcell_shuttle
   (Li2S8 -> Li2S6 -> Li2S4 -> Li2S2 -> Li2S by relabelling one site) is NOT
   S-conservative: each step drops S atoms (s_total drifts). It is kept only as
   the milestone-1 legacy/validation case. Since 2026-09-19 the reference case
   cases/fullcell_cei uses the S8-UNIT LATTICE MODEL (section 9), which
   conserves S by construction; s_total must be constant there.
5. Engine speed ~0.5 s per event regardless of box size: the cost is in the
   per-candidate Python loops, not in the lattice. Optimization pending.
5. (fixed 2026-09-18) cycle_stats.csv is flushed at every checkpoint and on
   SIGTERM, so a SLURM time-limit kill no longer loses it.
