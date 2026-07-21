#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
sage code/39_degree3_torsion_recombination_constraint_census.sage 2>&1 \
  | tee "logs/39_torsion_constraint_census_$(date +%Y%m%d_%H%M%S).log"
