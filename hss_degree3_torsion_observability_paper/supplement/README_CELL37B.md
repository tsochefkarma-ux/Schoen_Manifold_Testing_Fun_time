# Cell 37B — Degree-three linear curvature response

## Purpose

Cell 37B converts the fifteen symmetric CM-jet tensors from Cell 37A into the actual linear degree-three HSS curvature-defect columns.

For every observable coefficient and each parked sector it computes

\[
G_j(U)=\Delta S_{\mathrm{curv},3}(U)-\alpha_0\sigma_{3,j}(U)-\beta_0D_{3,j}(U).
\]

The calculation is performed directly as a truncated Laurent expansion at the frame divisor

\[
U=3X+1=0.
\]

It records the common Laurent support and the numerical rank of the fifteen shared ambiguity columns.

## Important scope restriction

The script also constructs a deliberately maximal support-matched convolution module from shifts of the two leading channels. This is only an **absorption diagnostic**. It is not asserted to be the physically admissible degree-three connection.

Cell 37C must derive the actual degree-three pole depth, allowed Laurent powers, universal drift and residue relations before a physical selectivity rank can be quoted.

## Required inputs

Place these files under `results/`:

- `degree3_cm_jets_cell37a.sobj`
- `degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj`

## Install

Copy:

- `37b_degree3_linear_curvature_response.sage` to `code/`
- `run_cell37b_linear_response.sh` to the workspace root

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell37b_linear_response.sh
bash run_cell37b_linear_response.sh 256 32
```

The second argument is the Laurent precision. Increase it to `40` if any reported support reaches the upper census edge.

## Outputs

- `results/degree3_linear_curvature_response_cell37b.sobj`
- `results/degree3_linear_curvature_response_cell37b.txt`
- `results/degree3_transport_selectivity_input_cell37b.sobj`

## Interpretation

The important outputs are:

- the observed defect support in `U`;
- `numerical rank of fifteen ambiguity responses`;
- the maximal support-matched selectivity diagnostic.

A zero maximal diagnostic is expected if shifts of the coprime leading channels span the whole observed Laurent row space. That would not yet decide physical selectivity. The physical number requires the geometrically derived module in Cell 37C.
