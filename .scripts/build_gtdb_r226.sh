#!/usr/bin/env bash
# Build GTDB r226 all-genomes GPK (5.2M genomes) — NFS manifest + direct pwrite mode.
#
# Architecture (no TCP, no firewall issues, no double-write):
#   - genopack coordinator starts on this node, watches a manifest dir on NFS
#   - Each worker builds its partition to a temp .gpk, pwrite()s sections
#     directly into the final output file at coordinator-allocated offsets,
#     then signals completion via an NFS manifest file
#   - Coordinator writes unified TOC once all workers are done
#   - Workers launch via SSH (no reverse tunnel needed — data flows via NFS)
#   - Sketches (SKCH) are always built (--no-sketch to disable)
#
# Usage:
#   Local (all workers on this node):
#     ./build_gtdb_r226.sh -o /path/to/output.gpk
#
#   Distributed (SSH to chaos nodes):
#     ./build_gtdb_r226.sh -o /path/to/output.gpk \
#       --nodes "dandycomp01fl dandycomp02fl dandycomp03fl dandycomp06fl"
#
#   With spill dir for large intermediate data:
#     ./build_gtdb_r226.sh -o /path/to/output.gpk \
#       --spill-dir /scratch/genopack_spill \
#       --nodes "dandycomp01fl dandycomp02fl dandycomp03fl dandycomp06fl"
set -euo pipefail

GENOPACK=/maps/projects/fernandezguerra/apps/repos/genopack/build/genopack
SLICES_DIR=/maps/projects/fernandezguerra/apps/repos/genopack/gtdb_5195k_dist_parts
N_PARTS=4
OUT_GPK=""
NODES=()
SPILL_DIR=""

usage() {
    echo "Usage: $0 -o <output.gpk> [-n <parts>]"
    echo "          [--nodes \"node1 node2 ...\"] [--spill-dir /path]"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o)          OUT_GPK="$2";                   shift 2 ;;
        -n)          N_PARTS="$2";                   shift 2 ;;
        --nodes)     read -ra NODES <<< "$2";        shift 2 ;;
        --spill-dir) SPILL_DIR="$2";                 shift 2 ;;
        *)           usage ;;
    esac
done
[[ -z "${OUT_GPK}" ]] && usage

WORK_DIR="${SLICES_DIR}/work"
mkdir -p "${WORK_DIR}"

LOG="${WORK_DIR}/build.log"
exec > >(tee -a "${LOG}") 2>&1
echo "=== Build started $(date) on $(hostname -f) ==="
echo "    output:  ${OUT_GPK}"
echo "    parts:   ${N_PARTS}"
[[ ${#NODES[@]} -gt 0 ]] && echo "    nodes:   ${NODES[*]}" || echo "    mode:    local"

# ── Step 1: combine slices ────────────────────────────────────────────────────
if [[ ! -f "${WORK_DIR}/combined.tsv" ]]; then
    echo "[1/4] Combining slices..."
    head -1 "${SLICES_DIR}/slice_0.tsv" > "${WORK_DIR}/combined.tsv"
    for i in 0 1 2; do
        tail -n +2 "${SLICES_DIR}/slice_${i}.tsv"
    done >> "${WORK_DIR}/combined.tsv"
    echo "      $(( $(wc -l < "${WORK_DIR}/combined.tsv") - 1 )) genomes combined."
else
    echo "[1/4] combined.tsv exists, skipping."
fi

# ── Step 2: normalize taxonomy ────────────────────────────────────────────────
if [[ ! -f "${WORK_DIR}/normalized.tsv" ]]; then
    echo "[2/4] Normalizing taxonomy..."
    "${GENOPACK}" taxonomy normalize \
        -i "${WORK_DIR}/combined.tsv" \
        -o "${WORK_DIR}/normalized.tsv"
else
    echo "[2/4] normalized.tsv exists, skipping."
fi

# ── Step 3: genus-balanced partition ─────────────────────────────────────────
if [[ ! -f "${WORK_DIR}/parts/part_0.tsv" ]]; then
    echo "[3/4] Partitioning into ${N_PARTS} genus-balanced parts..."
    mkdir -p "${WORK_DIR}/parts"
    "${GENOPACK}" taxonomy partition \
        -i "${WORK_DIR}/normalized.tsv" \
        -n "${N_PARTS}" \
        -o "${WORK_DIR}/parts" \
        -r g
else
    echo "[3/4] parts/part_0.tsv exists, skipping."
fi

# ── Step 4: NFS coordinator + parallel workers ────────────────────────────────
MANIFEST_DIR="${WORK_DIR}/manifest"
mkdir -p "${MANIFEST_DIR}"

# Clean up any leftover manifest files from a previous failed run.
rm -f "${MANIFEST_DIR}"/*.pending "${MANIFEST_DIR}"/*.alloc \
       "${MANIFEST_DIR}"/*.done   "${MANIFEST_DIR}"/*.tmp

echo "[4/4] Starting NFS coordinator (${N_PARTS} workers) → ${OUT_GPK}"
"${GENOPACK}" coordinator \
    -o "${OUT_GPK}" \
    --nfs-dir "${MANIFEST_DIR}" \
    --workers "${N_PARTS}" \
    >> "${WORK_DIR}/coordinator.log" 2>&1 &
COORD_PID=$!

# Give the coordinator a moment to create the output file before workers start.
sleep 3

WORKER_PIDS=()
for i in $(seq 0 $(( N_PARTS - 1 ))); do
    PART="${WORK_DIR}/parts/part_${i}.tsv"
    PART_TMP="${WORK_DIR}/parts/part_${i}_tmp"
    LOG_W="${WORK_DIR}/worker_${i}.log"
    COORD_ADDR="nfs:${MANIFEST_DIR}:${OUT_GPK}"

    BUILD_CMD="${GENOPACK} build"
    BUILD_CMD+=" -i ${PART}"
    BUILD_CMD+=" -o ${PART_TMP}"
    BUILD_CMD+=" --taxon-group --kmer-sort --mem-delta"
    BUILD_CMD+=" --coordinator ${COORD_ADDR}"
    [[ -n "${SPILL_DIR}" ]] && BUILD_CMD="GENOPACK_SPILL_DIR=${SPILL_DIR} ${BUILD_CMD}"

    WRAPPED="${BUILD_CMD} > ${LOG_W} 2>&1"

    if [[ ${#NODES[@]} -gt 0 ]]; then
        NODE="${NODES[$i % ${#NODES[@]}]}"
        echo "  Worker ${i}: SSH → ${NODE}"
        ssh -o "StrictHostKeyChecking=no" -o "BatchMode=yes" \
            "${NODE}" "bash -c '${WRAPPED}'" &
    else
        echo "  Worker ${i}: local"
        eval "${WRAPPED}" &
    fi
    WORKER_PIDS+=($!)
done

# ── Wait for all workers ──────────────────────────────────────────────────────
FAILED=0
for i in "${!WORKER_PIDS[@]}"; do
    if wait "${WORKER_PIDS[$i]}"; then
        echo "  Worker ${i}: done ✓"
    else
        echo "  ERROR: worker ${i} failed — see ${WORK_DIR}/worker_${i}.log"
        FAILED=$(( FAILED + 1 ))
    fi
done

if [[ "${FAILED}" -gt 0 ]]; then
    echo "=== ${FAILED}/${N_PARTS} workers failed. Fix and re-run. ==="
    kill "${COORD_PID}" 2>/dev/null || true
    exit 1
fi

# Wait for coordinator to finish writing TOC.
echo "  Waiting for coordinator to finalize..."
wait "${COORD_PID}" || { echo "=== coordinator failed — see ${WORK_DIR}/coordinator.log ==="; exit 1; }

echo "=== Build complete: ${OUT_GPK} ($(date)) ==="
