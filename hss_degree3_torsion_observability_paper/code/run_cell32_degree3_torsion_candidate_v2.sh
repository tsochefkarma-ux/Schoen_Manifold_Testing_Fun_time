#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

mkdir -p logs results
QMAX="${1:-24}"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="logs/32_degree3_A3_B3_hss_torsion_candidate_v2_${STAMP}.log"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Running corrected Cell 32 with q_max=${QMAX}"
echo "Log: ${LOG}"
sage code/32_degree3_A3_B3_hss_torsion_candidate_v2.sage "${QMAX}" 2>&1 | tee "$LOG"
