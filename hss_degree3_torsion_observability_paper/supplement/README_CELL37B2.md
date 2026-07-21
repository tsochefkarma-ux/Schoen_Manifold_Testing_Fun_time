# Cell 37B2 — Degree-three response rank-growth audit

This replaces the premature move to a physical degree-three connection module.
Cell 37B reported support `[-4,12]` while its census stopped at `U^12`, so the
upper support boundary and the rank `2/15` must be tested in a wider, safely
truncated Laurent window.

## Install

Copy:

- `37b2_degree3_response_rank_growth_audit.sage` to `~/research/code/`
- `run_cell37b2_rank_growth_audit.sh` to `~/research/`

Run in `~/research` with the Sage environment active:

```bash
bash run_cell37b2_rank_growth_audit.sh 320 64 32 60
```

The parameters are precision bits, Laurent precision, maximum audited
exponent, and display-tolerance digits. The script requires at least twelve
unused Laurent orders above the audited maximum.

It checkpoints each of the 30 response columns. A repeat with identical
parameters resumes automatically.

## Interpretation

The script prints the response rank in each parked sector and the combined
rank as the upper Laurent cutoff grows, at three numerical tolerances. It also
reports whether nonzero coefficients remain at the upper boundary.

If the combined rank remains `2/15`, then any selectivity computed from this
single CM curvature observable is at most two. A unique fifteen-parameter
bootstrap would then require additional CM points, off-symmetric slices, or
independent curvature/connection observables.

Outputs:

- `results/degree3_response_rank_growth_cell37b2.sobj`
- `results/degree3_response_rank_growth_cell37b2.txt`
- `results/degree3_response_rank_growth_cell37b2_partial.sobj`
