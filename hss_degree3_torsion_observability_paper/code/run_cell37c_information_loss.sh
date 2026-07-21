#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
PREC="${1:-320}"
LPREC="${2:-64}"
UMAX="${3:-24}"
TOLDIG="${4:-60}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sage code/37c_degree3_information_loss_and_observable_stack.sage "$PREC" "$LPREC" "$UMAX" "$TOLDIG" 2>&1 | tee "logs/37c_information_loss_$(date +%Y%m%d_%H%M%S).log"
