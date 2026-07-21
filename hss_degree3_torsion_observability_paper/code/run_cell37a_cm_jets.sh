#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
PREC_BITS="${1:-256}"
THETA_MAX="${2:-6}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
sage code/37_degree3_cm_jet_extraction.sage "$PREC_BITS" "$THETA_MAX" 2>&1 \
  | tee "logs/37a_cm_jets_$(date +%Y%m%d_%H%M%S).log"
