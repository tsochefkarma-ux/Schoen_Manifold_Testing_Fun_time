#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
LOG="logs/35b_torsion_modularity_audit_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to $LOG"
sage code/35b_degree3_torsion_modularity_audit.sage 2>&1 | tee "$LOG"
