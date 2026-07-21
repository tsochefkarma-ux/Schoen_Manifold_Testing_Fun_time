# Cell 32 — Degree-three A3/B3 HSS torsion candidate

This cell constructs the six weight-16, index-3 one-sided basis series along

\[
\tau=3t,\qquad z_r=t\gamma+\frac r3\gamma,\qquad r=0,1,2,
\]

with \(\gamma=(1,1,1,1,1,1,1,-1)\).  The tags `one`, `omega`, and `omega2`
refer to \(r=0,1,2\).

This is the canonical order-three shift along the HSS vector.  Its identification
with the geometric BKOS torsion labels is a hypothesis to be tested, not an input
treated as proved.

## Install

Copy:

- `32_degree3_A3_B3_hss_torsion_candidate.sage` to `/home/will/research/code/`
- `run_cell32_degree3_torsion_candidate.sh` to `/home/will/research/`

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell32_degree3_torsion_candidate.sh
bash run_cell32_degree3_torsion_candidate.sh 24
```

The final argument is the maximum sector q-power.  Start with 24.  A larger
cutoff such as 36 can be used later to check rank stabilization.

## Outputs

- `results/degree3_one_sided_torsion_basis.sobj`
- `results/degree3_A3_B3_hss_torsion_candidate.sobj`
- `results/degree3_A3_B3_hss_torsion_candidate.txt`

The cell verifies the five zero-elliptic-variable normalizations, constructs the
six basis series in exact cyclotomic arithmetic, and computes the individual and
combined ranks of the `(omega,omega)` and `(omega,omega2)` shared-sector maps.
