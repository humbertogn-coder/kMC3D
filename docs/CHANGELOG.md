# Changelog (engine changes, each validated byte-identical on validation/)

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
