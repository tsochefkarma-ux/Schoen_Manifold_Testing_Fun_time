#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
LOG="logs/41_external_data_bridge_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to $LOG"
sage code/41_degree3_minimal_external_data_bridge.sage 2>&1 | tee "$LOG"
