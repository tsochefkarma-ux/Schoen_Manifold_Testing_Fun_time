# Cell 37C — Degree-three information-loss ladder

This is a diagnostic stage, not a physical transport-module assumption.

It reconstructs the same 30 degree-three columns as Cell 37B2 and measures the
rank of the fifteen shared ambiguity directions after each geometric layer:

1. raw symmetric CM tensor;
2. Kähler-potential jet;
3. metric perturbation;
4. curvature connection vector;
5. eigenvalue channels;
6. spectral channels `(sigma,D)`;
7. curvature primitives;
8. scalar curvature;
9. the transport triplet `(Scurv,sigma,D)`;
10. the final leading-channel defect.

It also verifies the observed rank-one factorization separately in each parked
sector and checkpoints every completed column.

## Install

Copy:

- `37c_degree3_information_loss_and_observable_stack.sage` to `code/`
- `run_cell37c_information_loss.sh` to the workspace root

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell37c_information_loss.sh
bash run_cell37c_information_loss.sh 320 64 24 60
```

Use `tmux` for a persistent run. The output identifies the earliest geometric
operation at which rank drops below 15. That determines whether the next
observable family should use metric components, curvature primitives, or a
new CM/slice evaluation.
