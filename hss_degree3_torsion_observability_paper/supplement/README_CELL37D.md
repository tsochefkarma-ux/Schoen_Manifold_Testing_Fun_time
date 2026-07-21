# Cell 37D — minimal local-observable lift

Cell 37C found that the five curvature primitives retain rank 14 of the 15
observable degree-three ambiguity directions, while the final scalar defect has
rank 2. Cell 37D identifies the smallest local observable that restores the
missing fifteenth direction.

## Inputs

The following local result must exist:

```text
results/degree3_information_loss_cell37c.sobj
```

## Install

Copy:

```text
37d_degree3_minimal_observable_lift.sage
```

to `/home/will/research/code/`, and copy

```text
run_cell37d_minimal_observable_lift.sh
```

to `/home/will/research/`.

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell37d_minimal_observable_lift.sh
bash run_cell37d_minimal_observable_lift.sh
```

The cell tests frame-metric components first, then coordinate metric,
connection-vector, Kähler-jet, and raw-tensor components. It requires a lift to
be stable at tolerances `1e-30`, `1e-40`, and `1e-50`.

## Outputs

```text
results/degree3_minimal_observable_lift_cell37d.sobj
results/degree3_minimal_observable_lift_cell37d.txt
results/degree3_selected_observable_stack_cell37d.sobj
```

This cell is an information census only. A full-rank stack does not yet imply a
physical inverse theorem; the added observable still needs an independently
justified transport or covariance law.
