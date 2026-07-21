# Cell 33 — Degree-three evaluation kernel and sector lift

Cell 32 found exact rank 15 out of 21 for the two parked candidate sectors
`(omega,omega)` and `(omega,omega2)`. Cell 33 computes the resulting exact
six-dimensional coefficient kernel and asks which additional character-sector
evaluations see those six invisible directions.

## Required input

The Cell-32 result must exist at:

```text
results/degree3_A3_B3_hss_torsion_candidate.sobj
```

## Install

Copy:

```text
33_degree3_evaluation_kernel_and_sector_lift.sage
```

to:

```text
/home/will/research/code/
```

Copy:

```text
run_cell33_evaluation_kernel.sh
```

to:

```text
/home/will/research/
```

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell33_evaluation_kernel.sh
bash run_cell33_evaluation_kernel.sh
```

## Outputs

```text
results/degree3_evaluation_kernel_and_sector_lift.sobj
results/degree3_evaluation_kernel_and_sector_lift.txt
```

## Interpretation

The script reports:

- exact one-sided ranks at `one`, `omega`, and `omega2`;
- exact ranks of all six symmetric character pairs;
- a six-vector exact basis for the parked-sector kernel;
- how many kernel directions each additional sector lifts;
- the smallest additional sector subsets attaining full rank 21, when full
  rank is available at this candidate specialization;
- q-cutoff profiles showing whether the lift has stabilized.

The calculation is exact over the cyclotomic coefficient field saved by Cell
32. It does not load the large formal-CM checkpoint.
