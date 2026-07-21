# Cell 38 — degree-two tensorial / exchange-parity calibration

Cell 37F found an all-orders observable decomposition of the degree-three
three-torsion HSS response into a 14-dimensional exchange-even sector and a
one-dimensional exchange-odd quotient direction.  Cell 38 calibrates the same
idea at degree two before any matrix-valued closure is proposed.

The complete physical degree-two ambiguity is the three-dimensional symmetric
square of the one-sided basis `P2` and `B(q^2)`.  Cell 38 adds one deliberately
unphysical ordered-space control,

```text
F_- = P2(y) B(z^2) - B(y^2) P2(z),
```

and computes the linear degree-two Kähler, metric, curvature, spectral, scalar,
and odd/even block responses.

The calculation tests whether:

- the three physical columns are exchange-even;
- the added control is exchange-odd;
- scalar/eigenvalue observables discard the odd direction;
- the odd-to-even metric covector restores the fourth ordered-space direction;
- the exact degree-two scalar selectivity remains zero on the complete physical
  space, as proved in Cell 36.

A successful run means that parity and the extra matrix block improve
**observability**, but do not themselves create an enumerative constraint.
Any future degree-three matrix-valued closure must therefore justify additional
universal relations independently of the unknown q-series.

## Install

Copy:

- `38_degree2_tensorial_parity_calibration.sage` to `/home/will/research/code/`
- `run_cell38_degree2_tensorial_calibration.sh` to `/home/will/research/`

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell38_degree2_tensorial_calibration.sh
bash run_cell38_degree2_tensorial_calibration.sh 320 64 24 50 80
```

## Inputs

- `results/degree3_cm_jets_cell37a.sobj`
- `results/degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj`

## Outputs

- `results/degree2_tensorial_parity_calibration_cell38.sobj`
- `results/degree2_tensorial_parity_calibration_cell38.txt`
- `results/degree3_matrix_closure_acceptance_template_cell38.sobj`
