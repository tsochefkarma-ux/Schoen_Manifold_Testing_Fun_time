#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
LOG="logs/37f_exchange_parity_$(date +%Y%m%d_%H%M%S).log"
echo "Running Cell 37F; log: $LOG"
sage code/37f_degree3_exchange_parity_decomposition.sage 2>&1 | tee "$LOG"
