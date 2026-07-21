#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs results

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SCRIPT="code/29_exact_polar_dependency_certificate.sage"
if [[ ! -f "$SCRIPT" ]]; then
  echo "Missing $SCRIPT" >&2
  exit 1
fi

sage "$SCRIPT" 2>&1 | tee logs/29_exact_polar_dependency_certificate.log

echo
echo "Complete. Summary:"
echo "  results/degree2_exact_polar_dependency_summary.txt"
