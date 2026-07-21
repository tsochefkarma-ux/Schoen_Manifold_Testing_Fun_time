#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
sage code/37e_degree3_parity_block_mixing_certificate_v2.sage 2>&1 | tee "logs/37e_parity_block_v2_$(date +%Y%m%d_%H%M%S).log"
