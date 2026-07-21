# Cell 30 — Degree-three bootstrap dimension census

This is the first degree-three cell. It is a fast counting certificate and does
not run the HSS geometry.

## What it computes

For genus zero and base degree three, the one-sided numerator has weight 16 and
E8 index 3. The script enumerates the fixed-weight part of the holomorphic
index-three module generated over level-one modular forms by:

- `A3`
- `B3`
- `A1*A2`
- `A1*B2`
- `A1^3`

It independently cross-checks the answer by raw monomial enumeration, builds
the symmetric E8 x E8 tensor-square basis, and reports the reconstruction rank
needed for the two parked torsion sectors `(omega,omega)` and
`(omega,omega^2)`.

## Installation

Copy:

```text
30_degree3_bootstrap_dimension_census.sage
```

to:

```text
/home/will/research/code/
```

Copy:

```text
run_cell30_degree3_census.sh
```

to:

```text
/home/will/research/
```

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell30_degree3_census.sh
bash run_cell30_degree3_census.sh
```

This should finish in seconds and use negligible memory.

## Outputs

```text
results/degree3_bootstrap_dimension_census.sobj
results/degree3_bootstrap_dimension_census.txt
logs/30_degree3_bootstrap_census_*.log
```

## Expected one-sided basis

```text
E6*A1*B2
E6^2*A3
E4*A1^3
E4*E6*B3
E4^2*A1*A2
E4^3*A3
```

The expected dimensions are:

```text
one-sided:                  6
ordered E8 x E8:           36
symmetric E8 x E8:         21
two one-sided sectors:     12 unknowns
two symmetric bi-sectors:  42 unknowns
```

The next cell will use this manifest to construct the actual index-three
Sakai generators and their HSS-specialized q-jets.
