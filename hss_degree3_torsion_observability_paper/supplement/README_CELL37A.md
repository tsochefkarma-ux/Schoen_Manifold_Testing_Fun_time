# Cell 37A — degree-three CM-jet extraction

This cell converts the audited canonical three-torsion q-series from Cell 32
into high-precision CM jets at `t=i`, using only the certified q-range.  It
prepares the five reduced one-sided jet families and the fifteen symmetric jet
tensors for `(omega,omega)` and `(omega,omega2)`.

## Required local result files

- `results/degree3_A3_B3_hss_torsion_candidate.sobj`
- `results/degree3_torsion_modularity_audit_35b.sobj`
- `results/degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj`

## Install

Copy `37_degree3_cm_jet_extraction.sage` to `code/` and
`run_cell37a_cm_jets.sh` to the workspace root.

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell37a_cm_jets.sh
bash run_cell37a_cm_jets.sh 256 6
```

## Outputs

- `results/degree3_cm_jets_cell37a.sobj`
- `results/degree3_cm_jets_cell37a.txt`

The output is a high-precision bridge, not the final exact formal-CM
certificate.  It is intended for the first degree-three selectivity census.
