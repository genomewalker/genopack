#!/usr/bin/env bash
# Rebuild GTDB r226 multi-k GPK (4 parts, k=16,21,31) across 4 chaos nodes.
# kmer-sort and taxon-group are now defaults in genopack — no explicit flags needed.
# Partition TSVs must already exist (see build_gtdb_r226.sh steps 1-3).
set -euo pipefail

GENOPACK=/maps/projects/fernandezguerra/apps/repos/genopack/build/genopack
PARTS_DIR=/maps/projects/fernandezguerra/apps/repos/genopack/gtdb_5195k_dist_parts/work/parts
OUT_DIR=/maps/projects/caeg/scratch/kbd606/geodesic/gtdb/226/all_genomes_r226_mp
LOG_DIR="${OUT_DIR}/_logs"
NODES=(dandycomp01fl dandycomp02fl dandycomp03fl dandycomp06fl)
N_PARTS=4

mkdir -p "${LOG_DIR}"

LOG="${LOG_DIR}/master_build.log"
echo "=== Rebuild started $(date) on $(hostname -f) ===" | tee "${LOG}"

WORKER_PIDS=()
for i in $(seq 0 $(( N_PARTS - 1 ))); do
    PART="${PARTS_DIR}/part_${i}.tsv"
    OUT_GPK="${OUT_DIR}/part_${i}.gpk"
    NODE="${NODES[$i]}"
    LOG_W="${LOG_DIR}/part_${i}_build.log"

    # kmer-sort and taxon-group are defaults; explicit here for clarity
    BUILD_CMD="${GENOPACK} build"
    BUILD_CMD+=" -i ${PART}"
    BUILD_CMD+=" -o ${OUT_GPK}"
    BUILD_CMD+=" --sketch-kmers 16,21,31"
    BUILD_CMD+=" --mem-delta"
    BUILD_CMD+=" --threads 48"

    echo "  Part ${i}: SSH → ${NODE}" | tee -a "${LOG}"
    ssh -o StrictHostKeyChecking=no -o BatchMode=yes \
        "${NODE}" "nohup bash -c '${BUILD_CMD} > ${LOG_W} 2>&1' </dev/null &" &
    WORKER_PIDS+=($!)
done

echo "  All ${N_PARTS} builds launched. Polling every 60s..." | tee -a "${LOG}"

# Poll until all part GPKs appear (genopack writes manifest on completion)
DONE=0
while [[ "${DONE}" -lt "${N_PARTS}" ]]; do
    sleep 60
    DONE=0
    for i in $(seq 0 $(( N_PARTS - 1 ))); do
        # genopack build is complete when the .gpk directory contains a manifest
        if [[ -f "${OUT_DIR}/part_${i}.gpk/manifest.json" ]]; then
            (( DONE++ )) || true
        fi
    done
    echo "  $(date '+%Y-%m-%d %H:%M:%S') ${DONE}/${N_PARTS} parts complete" | tee -a "${LOG}"
done

for pid in "${WORKER_PIDS[@]}"; do
    wait "${pid}" 2>/dev/null || true
done

echo "=== Rebuild complete: ${DONE}/${N_PARTS} parts done ($(date)) ===" | tee -a "${LOG}"
