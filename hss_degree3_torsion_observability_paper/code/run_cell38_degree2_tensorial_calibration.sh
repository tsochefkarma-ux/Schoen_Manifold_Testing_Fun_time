#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
PREC_BITS="${1:-320}"
LAURENT_PREC="${2:-64}"
EXP_MAX="${3:-24}"
TOL_DIGITS="${4:-50}"
Q_MAX="${5:-80}"
LOG="logs/38_degree2_tensorial_calibration_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to $LOG"
sage code/38_degree2_tensorial_parity_calibration.sage \
  "$PREC_BITS" "$LAURENT_PREC" "$EXP_MAX" "$TOL_DIGITS" "$Q_MAX" 2>&1 | tee "$LOG"
