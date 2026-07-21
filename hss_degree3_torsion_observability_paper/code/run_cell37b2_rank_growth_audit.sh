#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
PREC_BITS="${1:-320}"
LAURENT_PREC="${2:-64}"
EXP_MAX="${3:-32}"
TOL_DIGITS="${4:-60}"
LOG="logs/37b2_rank_growth_${PREC_BITS}b_U${EXP_MAX}_$(date +%Y%m%d_%H%M%S).log"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
sage code/37b2_degree3_response_rank_growth_audit.sage "$PREC_BITS" "$LAURENT_PREC" "$EXP_MAX" "$TOL_DIGITS" 2>&1 | tee "$LOG"
