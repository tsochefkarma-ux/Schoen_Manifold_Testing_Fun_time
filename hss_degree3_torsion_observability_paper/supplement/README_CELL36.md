# Cell 36 — Full-space baseline and degree-three selectivity ledger

This cell closes the degree-two calibration and creates the exact matrix
interface for the degree-three inverse bootstrap.

## Inputs

Required:

- `results/degree3_torsion_modularity_audit_35b.sobj`

Recommended:

- `results/degree3_certified_15d_quotient_v3.sobj`
- `results/degree2_exact_polar_dependency_certificate.sobj`

## Install

Copy:

- `36_degree2_full_space_baseline_and_degree3_selectivity_ledger.sage`
  to `/home/will/research/code/`
- `run_cell36_selectivity_ledger.sh`
  to `/home/will/research/`

Run:

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell36_selectivity_ledger.sh
bash run_cell36_selectivity_ledger.sh
```

## Exact content

At degree two the cell enumerates the complete weight-10/index-2 module,
verifies that `P2` and `B(q^2)` are independent using an exact `q^0/q^1`
minor, and concludes that `FA,FB,FC` span the full symmetric ambiguity.
The exact 9-by-9 transport is surjective onto the Laurent support, so bare
selectivity is zero.

At degree three the cell requires the passed Cell-35b audit, constructs the
fifteen symmetric observable coordinates, and writes:

`results/degree3_transport_selectivity_input_template.sobj`

Future geometry workers must supply exact matrices `T_admissible` and
`G_ambiguity`. The selective rank is

```text
rank([T_admissible | G_ambiguity]) - rank(T_admissible).
```

A rank of 0 means total absorption; a rank of 15 permits a unique observable
reconstruction.
