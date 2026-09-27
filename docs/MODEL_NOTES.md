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
    stop_anode_inactive_halves  0 = off. N > 0: stop after N consecutive
                     half-cycles with zero plating and zero stripping (anode
                     dead; the cathode alone exchanging pool Li+ is not cycling).
    anode_kinetics   legacy | bv. bv: Butler-Volmer plating/stripping with
                     bv_k0_site (74 s^-1), bv_alpha (0.5), eta_charge (-0.03 V),
                     eta_discharge (+0.03 V); one aggregated plating candidate
                     (k x eligible growth sites). docs/ANODE_KINETICS.md.
    cathode_kinetics legacy | bv. bv: REGION cathode rates evaluated at
                     cathode_V_charge / cathode_V_discharge with signed alpha.
    first_half       begin (legacy, charge first) | end (discharge first)
    end_half_when_idle  0 = off. N: the half-cycle ends after N consecutive
                     events without Li transfer (zero-current cut-off).
                     fullcell_physical uses 10 (50 in run 1 let a trickle of
                     cathode release keep dead charge halves alive).
    stall_attempts   200 (retries with electrolyte refresh when no event can fire)
    stall_reruns     10000 (rejected slow draws, dt > scanInterval/5, per step)
    stall_p_accept_min 1e-4 (a draw is hopeless when 1 - exp(-W scanInterval/5)
                     is below this). Reaching any limit ends the run cleanly
                     with a "Simulation stalled" diagnostic (W, pool, eligible
                     plating sites, electrolyte). Legacy values were 1000 /
                     100000 with a traceback; the reference cases never reach
                     them, so their output is unchanged.
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

MILESTONE 3 (2026-09-23): PHYSICALLY BASED FULL CELL CYCLES. Smoke test of
cases/fullcell_physical on a 6x6x24 box (12 half-cycles): with anode_kinetics
bv, cathode_kinetics bv, first_half end, end_half_when_idle 50 and a FINITE Li
inventory (li_bulk_init 0), the cell reaches a steady state: anode Li 403
(charged) <-> 185 (discharged), electrolyte constant, cathode fully converted
S8 <-> Li8S8/Li16S8 each half-cycle, ~205 Li+ moved per half-cycle both ways,
pool bounded, li_total and s_total exact, SEI 22 sites in 12 half-cycles
(~0.5 % of platings). Every rate traces to docs/KINETICS_TABLE.md. Two legacy
behaviours had to be switched off for this: the semi-infinite foil (which
refilled stripped vacancies and inflated the slab by the plated amount every
cycle) and the charge-first protocol (a cell assembled charged must discharge
first, otherwise the event budget of the first half goes to electrolyte
decomposition).

0003. RUN 2 OF THE PHYSICAL CASES (GRACE, 5 seeds each, 2026-09-24/27).
   All 15 runs closed cleanly (5 "Simulation terminated", 10 "Interrupted"
   at the time limit with ledger and checkpoint flushed); conservation exact.
   anode_physical (finite anode, 153 half-cycles, 1.4 h per seed): Li 1150
   <-> 1100, pool 27 <-> 77, 18 decompositions per 3832 platings (0.47 %,
   the CE anchor), 40 SEI sites after 76 cycles. Reference half cell of the
   physical line. Figure cases/anode_physical/analysis/
   halfcell_original_vs_physical_run2.png.
   anode_physical_A100 (A_SEI 100, 25 to 115 half-cycles in 6 h): 16 % side
   events per plating, SEI 300 to 800 sites, film porosity 0.86 to 0.88,
   roughness 3 to 4 A, 1 to 3 % of Li buried in 7 to 27 islands
   (analysis/morphology.csv). Slow: 4 to 16 s per counted step because SEI
   fragments diffuse (about 8 diffusion moves per reaction event, each a full
   framework update); see the profile in item 0004.
   fullcell_physical (10x10x30, 800 events per half, 75 to 89 half-cycles in
   12 h): 550 to 600 mAh/g_S (35 % utilisation, event-budget limited) with
   CE_cathode 0.99 to 1.01 for 13 cycles, then fade to ~10 mAh/g by cycle 37
   (figure analysis/overview.png). Death mechanism, seed 8597: the Li flagged
   as ionised by the legacy charge rule (a Li site within minR of an SEI
   species, chg = 1) grows 24 -> 1670 of 1714 Li sites while neutral surface
   Li (LiMetalS) falls 264 -> 13; stripping candidates require chg = 0, so
   the discharge halves end idle and the cathode stays at S8. The anode is
   NOT covered (electrolyte-facing Li rises 229 -> 870, the surface roughens):
   it is deactivated site by site, about one Li per SEI fragment. Sources of
   SEI: 5 to 15 electrolyte reductions per charge half = 1 to 2.5 % of
   platings, 2 to 5 times the anchored 0.47 % of the half cell, because the
   pool empties at the end of every charge (cathode fully oxidised) and a
   trickle of Li release keeps the half alive while only side reactions can
   fire (SFO second-step reductions reach 35 to 60 per half late in life).
   Each event yields ~4 fragment sites, each ionising ~1 Li: ~8 % of the
   transferred Li per cycle is deactivated against 0.2 to 0.5 % in the real
   cell. Two conclusions: (1) the protocol needs a real cut-off: end the
   charge half when the cathode has nothing left to oxidise and the discharge
   half when it has nothing left to reduce (no candidate of the RELEASE_LI /
   CONSUME_LI reactions), instead of an idle counter that a trickle resets;
   (2) MODEL DECISION for the advisor: Li next to SEI fragments cannot be
   stripped and plating needs an electrolyte neighbour, i.e. there is no Li+
   transport through the SEI; every fragment permanently removes ~1 Li from
   the active inventory. Options: through-SEI transfer with an attenuation
   exp(-beta d) (KINETICS_TABLE, tunnelling row, [P]), or anchoring the side
   rates on Li deactivated per cycle (fragments) instead of events (factor
   ~4 lower). The three fullcell_physical checkpoints can be resumed but the
   cells are dead; a run 3 needs (1) first.

0004. SPEED PROFILE (cloud bench, 2 cores, 6 half-cycles, cProfile).
   anode_physical: 0.41 s per counted step; A100: 0.077 s per step() call
   but ~9 step() calls per counted step (SEI fragment diffusion moves), so
   0.67 s per counted step on the bench and 4 to 16 s on GRACE as the SEI
   grows. Hot spots (share of wall time): _count_all_classes 30 to 40 %
   (neighbour class counts recomputed from scratch at every refresh, 15000 to
   32000 calls), _eact_all 25 to 40 % (activation energies of every Li site
   recomputed per diffusion scan), updateCharges 10 %, createEther 10 %,
   _nonsei_mask rebuilt 324000 times inside the wrap fill loop (5 %). A
   vectorisation pass (incremental neighbour counts, cached masks per step,
   Eact only for candidate sites) should give 3 to 5x with byte-identical
   output, to be validated on validation/ as every engine change.

0002. RUN 1 OF THE PHYSICAL CASES (GRACE, 5 seeds each, 2026-09-23/24).
   anode_physical (half cell, 152 half-cycles): side fraction 4 to 6 x 10^-3
   decompositions per plating, stable over 150 half-cycles (CE_cycle 1.00;
   original C++ engine: 1 to 2 decompositions per plating). Li metal sites
   grow 1150 -> 5900 in both engines with net plating ~0: this is the
   framework fill of the semi-infinite foil (legacy behaviour, kept for the
   half cell). One seed finished (1.5 h); four stalled at half-cycle 116-151
   while plating: the last frames show 6114 Li sites of 6000 BC sites (plus
   BA/OC), 3 ETH sites left. The semi-infinite foil had filled the whole box,
   so no electrolyte-facing site remained (this is the "cell consumed" end
   of the original engine). The stall then burned 4 h in the retry loop and
   the final ledger write failed (0-byte cycle_stats.csv), which motivated
   the stall detection, the atomic write and the recovery. Fix for the
   physical half cell: li_pool_mode shared + li_bulk_init 0 (finite anode,
   as in fullcell_physical); bench test: nLi 1150 <-> 1100, pool 27 <-> 77,
   li_total exact, steady. The legacy fixed mode stays for the antecedent.
   fullcell_physical (10x10x30, N/P 1.4, 800 events per half): capacity
   560 to 610 mAh/g_S (35 % utilisation, event-budget limited) with
   CE_cathode 0.99 to 1.00 for 5 to 9 cycles and exact conservation; then the
   anode is buried by SOL2 (F5DEE on Li2O): 53, 239, 700, 668 ... events per
   charge half once the first O sites exist, 2950 F/F5D sites by cycle 13,
   stripping stops, the discharge halves end idle, and the run stops at
   W = 0 ("no acceptable reaction", 3 seeds) or hangs (2 seeds). Cause: the
   kT/h prefactor made SOL2 4.5e6 s^-1 per site, 2e5 times the aggregated
   plating rate, and the idle cut-off (50) was kept alive by a trickle of
   cathode release. Fixes: SOL2 anchored to the SOL prefactor with Tan's
   barrier as relative attenuation, end_half_when_idle 10. To re-run.

0000. SERIES 2 RESULT (fullcell_cei, corrected engine, 5 seeds, 200 half-
   cycles, 2026-09-21): conservation exact; cathode survives (425/700 sites at
   cycle 100), dissolution/precipitation in equilibrium, CEI is the sulfur sink
   (35 % of S) after exhausting the surface FSI. BUT the anode dies by cycle
   ~30 (framework Li 1200 -> 5330 sites, electrolyte 17800 -> 1800, zero
   plating/stripping afterwards) and capacity plateaus at 22 mAh/g (1.3 %):
   the 400-event half-cycle can move at most ~200 Li while the cathode needs
   11200 (N/P 0.10), and the legacy anode kinetics (tuned for 50-event
   half-cycles) age the anode 8x per cycle. Levers: (1) geometry N/P ~ 1
   (cases/fullcell_cei_np1, series 3), (2) larger maxInterval (kills the anode
   faster, rejected alone), (3) physical anode kinetics (Butler-Volmer plating/
   stripping with experimental j0, alpha 0.5; DFT Ea for decomposition) as a
   new case anode_physical: the real calibration, in progress.
   Morphology of series 2 (kmc3d/morphology.py): the Li slab grows from 20 A
   to 104 A thick (5320 Li sites of 6000 BC sites), a single connected body
   with NO buried Li, and only 9 to 21 anode-SEI sites in 200 half-cycles
   (the original half cell formed ~1300 in 152). Anode decomposition (k0
   1e-4 to 1e-6 per site) is frozen out of the shared BKL by cathode events
   1000x faster: another face of the rate-scale mismatch. Legacy Li growth
   by the framework (foil refill + wrap) is what ends both the original and
   the full cell, in ~150 half-cycles of 50 events there, ~60 of 400 here.

000. FIXED 2026-09-21: cathode sites are no longer overwritten by the anode
   framework (Li fill / ETH fill). Earlier full-cell runs (fullcell_shuttle
   and the first fullcell_cei production) carry this artifact: their cathode
   loss and part of their li_bulk draw are not physical.

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
