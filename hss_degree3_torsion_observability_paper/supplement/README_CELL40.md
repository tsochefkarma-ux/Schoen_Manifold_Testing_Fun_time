# Cell 40 — exact degree-three reconstruction atlas

## Purpose

Cell 39 proved that the ordinary torsion-summed sector has rank 9 in the
all-orders 15-dimensional observable quotient, while one additional character
sector is enough to reach rank 15. Cell 40 constructs the exact inverse maps.

It does **not** assume that the extra sector is already known. It answers the
linear-algebra question: which extra sector would suffice, which six of its
q-coefficients provide the missing information, and how are the fifteen shared
coefficients reconstructed exactly?

## Inputs

Place these existing Cell 32 and Cell 34 outputs in `results/`:

- `degree3_A3_B3_hss_torsion_candidate.sobj`
- `degree3_certified_15d_quotient_v3.sobj`

## Install

Copy:

- `40_degree3_exact_reconstruction_atlas.sage` to `code/`
- `run_cell40_reconstruction_atlas.sh` to the workspace root

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell40_reconstruction_atlas.sh
bash run_cell40_reconstruction_atlas.sh
```

## Outputs

- `results/degree3_exact_reconstruction_atlas_cell40.sobj`
- `results/degree3_exact_reconstruction_atlas_cell40.txt`

For every viable extra character sector, the `.sobj` contains:

- the earliest certified q cutoff reaching rank 15;
- nine independent torsion-summed coefficient rows;
- six independent rows from the extra sector;
- the exact 15 by 15 system matrix;
- its exact inverse;
- an exact inverse-verification flag.

The convention is

```text
y = system_matrix * c
c = reconstruction_inverse * y
```

where `c` is the 15-component observable ambiguity vector and `y` contains the
selected q-series coefficients.

## Interpretation

A successful run proves a data-sufficiency statement: the torsion sum plus one
independent additional character sector determines the complete observable
potential. It does not prove that transport supplies that additional sector.
