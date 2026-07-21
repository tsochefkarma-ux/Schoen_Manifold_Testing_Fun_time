# Cell 29: Exact polar-dependency certificate

This cell resolves the 13 subsets left open by the specialize-first Cell 28
summary. It does **not** load the 467 MB formal checkpoint and does not rebuild
the degree-two geometry.

It uses the existing files:

```text
results/28e_specialize_first/point_*.sobj
```

and the exact Laurent equations at orders `U^-4` and `U^-3`.

## Install

Copy:

```text
29_exact_polar_dependency_certificate.sage
```

into:

```text
/home/will/research/code/
```

Copy:

```text
run_cell29_exact_polar_certificate.sh
```

into:

```text
/home/will/research/
```

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell29_exact_polar_certificate.sh
bash run_cell29_exact_polar_certificate.sh
```

This should be a lightweight run. Its outputs are:

```text
results/degree2_exact_polar_dependency_summary.txt
results/degree2_exact_polar_dependency_certificate.sobj
logs/29_exact_polar_dependency_certificate.log
```

## Expected hierarchy

The structural model predicts:

- no cubic pole: consistent, dimension 2;
- no cubic or quadratic poles: consistent, dimension 1;
- completely pole-free transport: inconsistent;
- pure universal transport: inconsistent;
- unresolved polar subsets: 0.
