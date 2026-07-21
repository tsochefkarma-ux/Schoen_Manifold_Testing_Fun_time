#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
LOG="logs/37d_minimal_observable_lift_$(date +%Y%m%d_%H%M%S).log"
echo "Running Cell 37D"
echo "Log: $LOG"
sage code/37d_degree3_minimal_observable_lift.sage 2>&1 | tee "$LOG"
