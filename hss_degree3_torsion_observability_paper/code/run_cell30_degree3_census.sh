#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs results

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

STAMP=$(date +%Y%m%d_%H%M%S)
LOG="logs/30_degree3_bootstrap_census_${STAMP}.log"

if [[ -f code/30_degree3_bootstrap_dimension_census.sage ]]; then
    SCRIPT="code/30_degree3_bootstrap_dimension_census.sage"
elif [[ -f 30_degree3_bootstrap_dimension_census.sage ]]; then
    SCRIPT="30_degree3_bootstrap_dimension_census.sage"
else
    echo "Cannot find Cell 30 Sage script." >&2
    exit 1
fi

echo "Running $SCRIPT"
echo "Log: $LOG"
sage "$SCRIPT" 2>&1 | tee "$LOG"
