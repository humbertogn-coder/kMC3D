#!/bin/bash
# pack_sweep.sh - bundle a LIGHT sample of a parameter sweep for local analysis.
# Usage (on GRACE, from the repository root):
#     bash slurm/pack_sweep.sh cases/sweep_plating_j0
# Produces <sweep name>_sample.tar.gz in the repository root containing
# sweep.json, launch.sh, analysis/ (if present) and, per member and seed, the
# same light sample as pack_results.sh (the *.in used, cycle_stats.csv,
# kmc_info_log.txt, Data2Excel.txt, checkpoint.npz.json, snapshots, and the
# first and last trajectory frame). Full trajectories stay on /scratch.
set -e
SWEEP="${1:?usage: pack_sweep.sh <sweep dir>}"
NAME=$(basename "$SWEEP")
LIST=$(mktemp)
ls "$SWEEP"/sweep.json "$SWEEP"/launch.sh 2>/dev/null >> "$LIST" || true
[ -d "$SWEEP/analysis" ] && find "$SWEEP/analysis" -type f >> "$LIST"
for m in "$SWEEP"/f_*; do
    [ -d "$m" ] || continue
    ls "$m"/*.in 2>/dev/null >> "$LIST" || true
    for d in "$m"/runs/seed_*; do
        [ -d "$d" ] || continue
        ls "$d"/*.in "$d"/cycle_stats.csv "$d"/kmc_info_log.txt "$d"/Data2Excel.txt \
           "$d"/checkpoint.npz.json 2>/dev/null >> "$LIST" || true
        ls "$d"/snapshots/*.vasp 2>/dev/null >> "$LIST" || true
        ls "$d"/trajectory/kmc-coords-*.xyz 2>/dev/null \
            | sed 's/.*kmc-coords-\([0-9]*\)\.xyz/\1 &/' | sort -n \
            | awk 'NR==1{print $2} {last=$2} END{if (NR>1) print last}' >> "$LIST"
    done
done
tar -czf "${NAME}_sample.tar.gz" -T "$LIST"
rm -f "$LIST"
echo "wrote ${NAME}_sample.tar.gz ($(du -h "${NAME}_sample.tar.gz" | cut -f1))"
tar -tzf "${NAME}_sample.tar.gz" | wc -l | xargs echo "files:"
