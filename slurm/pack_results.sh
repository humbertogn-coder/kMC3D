#!/bin/bash
# pack_results.sh  -  bundle a LIGHT sample of a multi-seed run for local analysis.
# Usage (on GRACE, from the repository root):
#     bash slurm/pack_results.sh cases/fullcell_shuttle
#     bash slurm/pack_results.sh cases/fullcell_cei
# Produces <case>_sample.tar.gz in the repository root containing, per seed:
#     the five *.in actually used, cycle_stats.csv, kmc_info_log.txt,
#     Data2Excel.txt, checkpoint.npz.json (scalars only, small),
#     trajectory/: FIRST and LAST frame only, snapshots/: all POSCARs.
# Full trajectories stay on /scratch (hundreds of MB per seed).
set -e
CASE="${1:?usage: pack_results.sh <case dir>}"
NAME=$(basename "$CASE")
LIST=$(mktemp)
for d in "$CASE"/runs/seed_*; do
    [ -d "$d" ] || continue
    ls "$d"/*.in "$d"/cycle_stats.csv "$d"/kmc_info_log.txt "$d"/Data2Excel.txt \
       "$d"/checkpoint.npz.json 2>/dev/null >> "$LIST" || true
    ls "$d"/snapshots/*.vasp 2>/dev/null >> "$LIST" || true
    # first and last trajectory frame by numeric step index
    ls "$d"/trajectory/kmc-coords-*.xyz 2>/dev/null \
        | sed 's/.*kmc-coords-\([0-9]*\)\.xyz/\1 &/' | sort -n \
        | awk 'NR==1{print $2} {last=$2} END{if (NR>1) print last}' >> "$LIST"
done
tar -czf "${NAME}_sample.tar.gz" -T "$LIST"
rm -f "$LIST"
echo "wrote ${NAME}_sample.tar.gz ($(du -h "${NAME}_sample.tar.gz" | cut -f1))"
tar -tzf "${NAME}_sample.tar.gz" | wc -l | xargs echo "files:"
