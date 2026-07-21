#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p logs results/28e_specialize_first

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

WORKER="code/28e_specialize_first_worker.sage"
AGGREGATOR="code/28f_aggregate_specialize_first.sage"
POINTS="${POINTS:-8}"

if [[ ! -f "$WORKER" ]]; then
  echo "Missing $WORKER" >&2
  exit 1
fi
if [[ ! -f "$AGGREGATOR" ]]; then
  echo "Missing $AGGREGATOR" >&2
  exit 1
fi

echo "Specialize-first inverse-selection run"
echo "Workers: $POINTS"
echo "No formal checkpoint will be loaded."

for ((i=0; i<POINTS; i++)); do
  out=$(printf 'results/28e_specialize_first/point_%02d.sobj' "$i")
  log=$(printf 'logs/28e_specialize_first_point_%02d.log' "$i")
  if [[ -f "$out" ]]; then
    echo "[skip] point $i already complete"
    continue
  fi
  echo "[run] point $i"
  sage "$WORKER" "$i" 2>&1 | tee "$log"
  echo "[done] point $i"
done

echo "[aggregate]"
sage "$AGGREGATOR" 2>&1 | tee logs/28f_aggregate_specialize_first.log

echo "Complete. Summary:"
echo "  results/degree2_inverse_selection_specialize_first_summary.txt"
