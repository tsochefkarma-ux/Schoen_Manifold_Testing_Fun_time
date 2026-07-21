#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
export PYTHONUNBUFFERED=1
sage code/34_degree3_certified_15d_quotient_v3.sage 2>&1 | tee "logs/34_certified_15d_quotient_v3.log"
