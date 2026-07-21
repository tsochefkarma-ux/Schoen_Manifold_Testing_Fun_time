#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs results

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

LOG="logs/33_evaluation_kernel_$(date +%Y%m%d_%H%M%S).log"
echo "Running Cell 33"
echo "Log: $LOG"

sage code/33_degree3_evaluation_kernel_and_sector_lift.sage 2>&1 | tee "$LOG"
