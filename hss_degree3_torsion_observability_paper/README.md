# Paper IV - Three-Torsion Observability and the Limits of Residue Transport

This folder contains the fourth paper in the HSS/Schoen spectral-curvature series.

## Main files

- `three_torsion_observability.pdf` - compiled paper
- `three_torsion_observability.tex` - LaTeX source
- `certificate_summary.txt` - concise result ledger
- `RUN_ORDER.md` - Sage reproduction order
- `code/` - final Sage scripts and launchers
- `results/` - compact certificate/result objects available in this bundle
- `supplement/` - research notes and output summaries

## Scope

Paper IV proves the exact degree-two zero-selectivity baseline, the all-orders HSS three-torsion specialization relation, the exact 15-dimensional observable quotient, and the exact finite coefficient reconstruction atlas. The degree-three geometric rank ladder and the 14+1 parity split are high-precision computational certificates and are marked accordingly.

The paper does **not** claim that the unknown mixed torsion sector has already been derived. The next stage is to certify the BKOS-to-HSS variable and multicover dictionary and fill the six requested mixed-sector coefficients.

## SageMath

The scripts were developed and run with SageMath on Linux/WSL. Later cells depend on result files created by earlier cells. Follow `RUN_ORDER.md`.

The large Part III formal checkpoint is not included because the Paper IV specialize-first workflow avoids it.
