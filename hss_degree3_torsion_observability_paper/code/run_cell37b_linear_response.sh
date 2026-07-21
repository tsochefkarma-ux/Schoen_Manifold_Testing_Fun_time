#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
PREC_BITS="${1:-256}"
LAURENT_PREC="${2:-32}"
STAMP="$(date +%Y%m%d_%H%M%S)"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
sage code/37b_degree3_linear_curvature_response.sage "$PREC_BITS" "$LAURENT_PREC" 2>&1 \
  | tee "logs/37b_linear_response_${STAMP}.log"
