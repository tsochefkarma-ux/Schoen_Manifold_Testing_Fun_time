#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs results

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

LOG="logs/40_degree3_reconstruction_atlas_$(date +%Y%m%d_%H%M%S).log"
echo "Running Cell 40; log: $LOG"
sage code/40_degree3_exact_reconstruction_atlas.sage 2>&1 | tee "$LOG"
