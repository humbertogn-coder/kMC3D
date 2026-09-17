# kmc3d

Generalized 3D fixed-lattice kinetic Monte Carlo for electrochemical cells,
written in Python. Ported from a half-cell Li-metal SEI engine and extended to a
full Li-S cell (Li anode + S8 cathode, polysulfide shuttle, shared Li+ pool,
cathode passivation). Battery500 consortium, Balbuena group, Texas A&M.

## Layout

    kmc3d/          engine package (config, mechanism, lattice, engine, stats,
                    output, run) and post-processing modules (postprocess,
                    characterize, zeopp, raspa, structio, ensemble, ff_data)
    cases/          input decks (the five .in files per case). Outputs go to
                    cases/<case>/runs/ and are NOT versioned
    validation/     tiny byte-identity reference cases + reference md5 sums
    postprocess/    analysis scripts (CE, capacity, morphology, ML, SHAP)
    slurm/          GRACE job scripts (test and multi-seed production array)
    docs/           MODEL_NOTES.md (keywords, placeholders, limits), CHANGELOG.md
    examples/       minimal sample outputs (one cycle_stats.csv, two .xyz frames)

## Quick start (local)

    conda env create -f environment.yml
    conda activate kmc3d
    python -m kmc3d.run --dir cases/fullcell_shuttle --out cases/fullcell_shuttle/runs/seed_8597

## Validation (mandatory after any engine change)

See validation/README.md. Without the new keyword the engine must reproduce
the reference md5 sums byte for byte; with it, li_total must stay constant.

## GRACE (HPRC, Texas A&M)

    cd /scratch/user/humbertogn/kmc3d && git pull
    mkdir -p slurm_logs
    sbatch slurm/goKMC_test.slrm            # short test
    sbatch slurm/goKMC_prod_array.slrm      # 5 seeds, one core each

Results stay on /scratch. Bring back only cycle_stats.csv, Data2Excel.txt and
the first and last .xyz of each seed for analysis.

## Model notes

Rates for the cathode cascade, dissolution, precipitation and shuttle are
literature placeholders and are NOT calibrated yet (docs/MODEL_NOTES.md,
section 5). The confidential C++ reference engine is not part of this repository.
