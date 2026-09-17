# validation/ - byte-identity reference cases for kmc3d engine changes

Rule: any change to the engine is accepted only if the cases below reproduce
the reference md5 sums (Data2Excel.txt, cycle_stats.csv and every .xyz frame)
when the new feature keyword is NOT present in the input files.

Cases (tiny, seconds to minutes, NOT physical):

* anode_small/     anode-only, legacy Li pool (li_pool_mode fixed), 6 anode
                   reactions, cathode disabled. This is the legacy code path.
* fullcell_small/  Li-S full cell with shared pool, reservoir and shuttle
                   (same mechanism as fullcell_shuttle/, 6x6x24 box, 30 events
                   per half-cycle, 6 half-cycles).

How to run (from kMC_Alpha/, with the kmc3d environment active):

    set PYTHONDONTWRITEBYTECODE=1
    python -m kmc3d.run --dir validation/anode_small    --out val_out/anode_small
    python -m kmc3d.run --dir validation/fullcell_small --out val_out/fullcell_small

Then compare, per case, the md5 of Data2Excel.txt, cycle_stats.csv and
trajectory/*.xyz against reference_md5/<case>.md5. All must match.

Reference sums were produced on 2026-09-17 with Python 3.10.12 / numpy 2.2.6.
A different numpy major version may change float formatting in Data2Excel.txt;
if only Data2Excel.txt differs, regenerate the reference with the unchanged
engine on the new platform before judging a code change.

Feature-on checks (not byte-identity, physical sanity):

* li_total column of cycle_stats.csv must be constant over the whole run.
* Sulfur balance: lattice Li2Sx sites + reservoir counts must close.
* cathode passivation: n_cathode_blocked column appears only when
  cathode_passivation_nmin > 0; dissolution/precipitation counts must not be
  suppressed by it (only Li-transfer reactions are gated).

Change log of validated engine changes:

* 2026-09-17  cathode Li2S passivation (GEOMETRY.in: cathode_passivating_species,
              cathode_passivation_nmin). Both cases byte-identical with the
              keywords absent; feature-on smoke test: li_total constant,
              n_cathode_blocked > 0, dissolution unaffected.
