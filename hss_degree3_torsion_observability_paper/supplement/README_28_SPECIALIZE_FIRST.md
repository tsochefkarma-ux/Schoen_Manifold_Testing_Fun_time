# Cell 28 specialize-first replacement

This version avoids the formal 467 MB checkpoint and never constructs the
multivariate rational connection map.

Each worker substitutes the seven CM parameters into a finite field first,
then rebuilds the original degree-two geometry over `GF(p)(X)`. It solves the
small Laurent system, writes a tiny witness, and exits. The aggregator combines
nonzero-minor witnesses into exact generic rank/consistency certificates.

## Install

Copy:

- `28e_specialize_first_worker.sage`
- `28f_aggregate_specialize_first.sage`

into `/home/will/research/code/`, and copy
`run_calibration_v4_specialize_first.sh` into `/home/will/research/`.

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_calibration_v4_specialize_first.sh
tmux new -s hss
bash run_calibration_v4_specialize_first.sh
```

Detach with `Ctrl-b`, release, then `d`.

## Restarting

Completed point files are skipped automatically. Each point is a fresh Sage
process, so memory is fully released between points.

## Output

```text
results/degree2_inverse_selection_specialize_first_summary.txt
results/degree2_inverse_selection_specialize_first.sobj
```
