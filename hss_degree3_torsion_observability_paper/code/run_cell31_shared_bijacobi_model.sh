#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs results

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

LOG="logs/31_shared_bijacobi_model_$(date +%Y%m%d_%H%M%S).log"
echo "Running Cell 31"
echo "Log: $LOG"
sage code/31_shared_bijacobi_torsion_model.sage 2>&1 | tee "$LOG"
