#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
LOG="logs/36_selectivity_ledger_$(date +%Y%m%d_%H%M%S).log"
echo "Writing log to $LOG"
sage code/36_degree2_full_space_baseline_and_degree3_selectivity_ledger.sage 2>&1 | tee "$LOG"
