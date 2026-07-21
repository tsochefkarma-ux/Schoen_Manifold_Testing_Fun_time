# Cell 39 — degree-three torsion recombination constraint census

Cell 38 proves that exchange parity improves observability but supplies no
physical degree-two selectivity.  Cell 39 therefore inventories independent
enumerative information before any new matrix transport law is proposed.

The calculation works in the audited all-orders 15-dimensional observable
quotient and rebuilds exact character-evaluation matrices from the Cell-32
cyclotomic q-series.

It measures:

- the rank of the ordinary torsion-summed series, which is the trivial-character
  evaluation `(one,one)`;
- the maximum possible rank after adding the current scalar-transport ceiling
  of two;
- the smallest sets of character evaluations that determine all 15 observable
  coefficients;
- the ranks of the nine inverse-Fourier torsion-class projectors;
- the ranks of natural diagonal, off-diagonal, and total character sums.

No new closure law is assumed.  All ranks are exact over `Q(zeta_12)` and use
only the certified Cell-32 q range.

## Inputs

- `results/degree3_A3_B3_hss_torsion_candidate.sobj`
- `results/degree3_certified_15d_quotient_v3.sobj`

## Install

Copy:

- `39_degree3_torsion_recombination_constraint_census.sage` to
  `/home/will/research/code/`
- `run_cell39_torsion_constraint_census.sh` to `/home/will/research/`

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell39_torsion_constraint_census.sh
bash run_cell39_torsion_constraint_census.sh
```

## Outputs

- `results/degree3_torsion_recombination_constraint_census_cell39.sobj`
- `results/degree3_torsion_recombination_constraint_census_cell39.txt`
