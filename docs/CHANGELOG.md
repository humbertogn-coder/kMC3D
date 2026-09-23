# Changelog (engine changes, each validated byte-identical on validation/)

## 2026-09-23
- Cycling protocol keywords (opt-in): first_half end (start with a discharge)
  and end_half_when_idle N (zero-current cut-off: the half-cycle ends after N
  consecutive events without Li transfer). Legacy path byte-identical.
- fullcell_physical smoke test: steady cycling (see MODEL_NOTES milestone 3).
  Case set to li_bulk_init 0 and anode 40.4 A (N/P ~ 1.4).
- cathode_kinetics bv (opt-in): REGION cathode reactions evaluated at the
  cathode potential of the half-cycle (cathode_V_charge / cathode_V_discharge)
  through the existing Arrhenius-Butler-Volmer term; alpha is signed by the
  electron direction (negative reduction, positive oxidation). Legacy path
  byte-identical (both validation cases).
- cases/fullcell_physical: thin-cathode full cell with BV anode (j0 Boyle
  2020, CE-anchored side rates), BV cathode cascade (E0 Kumaresan 2008 Table
  II, k0 per site from Marinescu 2016 iH,0 / iL,0), shuttle 0.01 (Mikhaylik
  2004), CEI, brakes and stops; 800-event half-cycles.

## 2026-09-22 (post-processing)
- kmc3d/compare_runs.py + postprocess/compare_halfcell.py: side-by-side
  figure of half-cell runs from Data2Excel.txt (original C++ output included):
  Li metal, cumulative decomposition, side events per plating, CE_cycle.
- gaussian_density_field: periodic wrap and optional sigma per species
  (sigma_i = 0.45 r_i from ff_data radii; solids only by default), so size
  heterogeneity enters the density field at no pair-wise cost.
- cases/anode_physical_A100: anode_physical with A_SEI = 100 on FSI/SFO/SOL
  only (morphology variant; KINETICS_TABLE.md section 7).
- anode_kinetics bv (opt-in, PARAMETERS.in): Butler-Volmer plating and
  stripping sharing one per-site exchange rate (bv_k0_site = j0 A_site / e,
  bv_alpha, eta_charge, eta_discharge); stripping keeps the MOBILITY.in
  coordination dependence as a relative modulation. Legacy path byte-identical
  (both validation cases). New case cases/anode_physical (half cell, 10x10x30,
  152 half-cycles of 50 events like the original run): j0 29.8 mA/cm2 (Boyle
  2020), CE-anchored side-reaction rates (CE 0.998, Yu 2022), SOL2 with kT/h
  prefactor and Ea 0.364 eV (Tan 2024).
- bv plating candidate is aggregated (k_bv x number of growth sites) instead
  of salt-site based: in the legacy scheme a slow side reaction fired (with a
  100 to 600 s time jump) whenever no salt site touched Li0, so the event CE
  was set by salt geometry rather than kinetics. Legacy path untouched.
  First anode_physical result (22 half-cycles): plating 565, stripping 550,
  side reactions 0, Li inventory stable at ~1100: consistent with CE 0.998.
- docs/KINETICS_TABLE.md verified against six PDFs (Boyle 2020, Kumaresan
  2008, Marinescu 2016, Mikhaylik 2004, Zhang 2016, Yue 2018): no [M] entries
  left. Key numbers: j0(Li, LiFSI/DME) 29.8 mA/cm2 = 74 s^-1 per site; cathode
  i0 per Kumaresan Table II (assumed by the authors) and Marinescu Table 1;
  ks 0.53 h^-1 (Mikhaylik) = 0.009 s^-1 per Li0 site in the 10x10x30 box.
- kmc3d/morphology.py + postprocess/morphology_case.py: lattice-topological
  descriptors per frame (Li connected components, buried Li and its coating,
  film extent and roughness per side, porosity profile relative to the anode
  centre, total solid fraction). No voxelization needed; 0.6 s per frame.
  Output <case>/analysis/morphology.csv.
- docs/KINETICS_TABLE.md: every rate of the physically based engine with its
  formula, value at 298 K, source and verification tag ([V] verified, [M]
  from memory to verify, [P] pending, [G] group). Nothing [M]/[P] enters a
  production case. docs/ANODE_KINETICS.md: Butler-Volmer formulation.
- stop_anode_inactive_halves (opt-in): clean stop when the anode has had no
  plating or stripping for N consecutive half-cycles. Both validation cases
  byte-identical.
- cases/fullcell_cei_np1: thin cathode (one BC layer, 100 S8 sites, N/P ~ 0.7)
  + all three stops. Series 3 candidate.
- Series 2 of fullcell_cei analyzed (docs/MODEL_NOTES.md section 8).

## 2026-09-21 (post-processing)
- kmc3d/species.py: single classification of labels (metal, electrolyte,
  cathode, deposit, cei, sei) and S8-unit label mapping; used by
  characterize (_sei_slab_bounds = anode film only), postprocess
  (morphology_parameters + n_cathode/n_CEI/n_deposit), ff_data (properties
  for all full-cell labels + proxy fallback), zeopp radii and raspa.
- kmc3d/cycles.py: per-cycle Q_dis, Q_ch, CE_cathode, capacity (mAh/g_S),
  utilization, cathode composition in conventional formulas, conservation
  report, multi-seed mean/std. kmc3d/figures.py: 4-panel case overview.
  postprocess/summarize_case.py: CLI (exit 1 on conservation violation).
- cases/fullcell_cei: shuttle k0 lowered 10x (0.05/0.05/0.02, calibration
  step 1, CE_shuttle 0.63 to 0.80 in the small-box test) and brakes enabled
  (li_pool_max 148, stop_electrolyte_fraction 0.05). Series 2 to be run on
  GRACE with the corrected engine; series 1 results are artifact-laden.
- Cathode lattice sites (BC sites in the cathode z-band) excluded from the
  three anode-framework fills: updateLiMetal "fill vacuum", wrap and
  createEther. Before, a cathode site emptied by dissolution was refilled with
  foil Li (or ETH) on the next step: Li metal inside the cathode, dissolution
  irreversible, re-precipitation never fired, li_bulk inflated. Found with the
  shuttle-rate sweep (k_sh/10, /50, /100 all lost the cathode identically).
  anode_small byte-identical; fullcell_small reference regenerated (intended
  change). Shuttle sweep conclusion: rates were not the cause of cathode loss.

## 2026-09-20
- Structural brakes (opt-in, PARAMETERS.in): li_pool_max (pool cap: no
  stripping / RELEASE_LI when full; mean-field current coupling of the two
  electrodes) and stop_electrolyte_fraction (clean stop when the cell dries
  out). Column pool_full_blocks when the cap is active. Both validation
  cases byte-identical.
- Sulfur leak fixed: a shuttle reaction with DEPOSIT is a candidate only when
  enough eligible BA surface sites exist (the Li surface saturated with
  Li2S2_an is passivated toward the shuttle). Found with s_total in the
  fullcell_cei production run (seeds 8597/8599/8601, after half-cycle ~120):
  deposits that found no room silently dropped their S. Audit counter and
  column deposit_unplaced (must stay 0). anode_small byte-identical;
  fullcell_small identical except the new column.
- First production results of fullcell_cei (5 seeds, 140 to 180 half-cycles,
  killed at the 24 h limit; cycle_stats.csv survived thanks to the checkpoint
  flush). li_total constant. Placeholder rates let the shuttle consume ~90 %
  of the cathode (700 -> 42 to 59 S8 sites; ~73 % of S as Li2S2_an on the
  anode); surface FSI exhausted, CEI 120 to 186 sites. Calibration target
  defined: the cathode must not vanish.

## 2026-09-19
- cases/fullcell_cei switched to the S8-unit lattice model of the cathode
  (S8, Li2S8, Li4S8, Li8S8, Li16S8; dissolution and precipitation rewritten
  exact in Li and S; Li2S6_d precipitates by disproportionation). Input files
  only, no engine change. Acceptance: s_total constant with everything on.
- CEI reactions for LiFSI/F5DEE in the new case cases/fullcell_cei
  (fullcell_shuttle + FSI reduced by long-chain polysulfides + film
  passivation). Literature basis in docs/CEI_LITERATURE.md (Small Methods
  2026, Nat. Energy 2022, JACS 2024). New MECHANISM op LI_LOSS and ledger
  column li_cei. Both validation cases byte-identical.

## 2026-09-18
- CEI infrastructure: REGION cathode_surface, MECHANISM ops CEI and S_LOSS,
  ledger columns s_total (S conservation, whenever the reservoir is active),
  s_cei, n_CEI, cei_<spc>. No CEI reactions in any case yet (rates pending
  literature). anode_small byte-identical; fullcell_small identical except
  the added s_total column (reference regenerated).
- s_total exposed the non-conservative single-site cathode cascade (model
  decision pending, MODEL_NOTES section 8).
- cycle_stats.csv is rewritten at every checkpoint and run.py turns SIGTERM
  (SLURM time limit) into a clean interrupt, so the ledger and a final
  checkpoint survive a killed job. Restart from that checkpoint continues the
  ledger without gaps. Both validation cases byte-identical.

## 2026-09-17
- Cathode Li2S passivation restored (GEOMETRY.in: cathode_passivating_species,
  cathode_passivation_nmin). Gates only Li-transfer reactions. New column
  n_cathode_blocked. anode_small and fullcell_small byte-identical with the
  keywords absent.
- validation/ created: tiny reference cases + reference md5 + README.
- Repository layout: cases/, slurm/, docs/, examples/, postprocess/.
- Cleanup: __pycache__, SLURM logs, superseded test outputs removed.

## 2026-09-04 (milestone 1, before this repository)
- Shared Li+ pool (li_pool_mode shared), Li foil reservoir (li_bulk_init),
  li_total conservation column.
- Well-mixed polysulfide reservoir and shuttle keywords (REGION anode_surface,
  TRIGGER EMPTY, REQUIRE, RATE_SCALE, DISSOLVE, PRECIP, RES, STRIP_LI0, DEPOSIT).
- Static cathode accessibility mask (cathode_access_fraction).
- CE_mod and CE_shuttle metrics.
