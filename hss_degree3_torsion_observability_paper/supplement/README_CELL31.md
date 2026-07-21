# Cell 31 — Shared symmetric bi-Jacobi torsion model

Cell 30 counted two possibilities: 12 independent one-sided coefficients or 42
independent symmetric bi-sector coefficients.  The general Schoen potential is,
however, one symmetric `E8 x E8` bi-quasi-Jacobi object.  Character sectors are
evaluations of that common object, so the two parked sectors share one symmetric
`6 x 6` coefficient matrix.

The correct shared ambiguity count is therefore

```text
6 * 7 / 2 = 21
```

rather than `2 * 21 = 42`.

## Install

Copy:

```text
31_shared_bijacobi_torsion_model.sage
```

to:

```text
/home/will/research/code/
```

Copy:

```text
run_cell31_shared_bijacobi_model.sh
```

to:

```text
/home/will/research/
```

Run:

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell31_shared_bijacobi_model.sh
bash run_cell31_shared_bijacobi_model.sh
```

## Outputs

```text
results/degree3_shared_bijacobi_model.sobj
results/degree3_shared_bijacobi_model.txt
results/degree3_one_sided_torsion_basis_template.sobj
```

The template defines the exact data schema required from Cell 32.  Cell 32 must
construct the six weight-16, index-3 one-sided basis series at the characters
`one`, `omega`, and `omega2`, including the new `A3` and `B3` specializations.

When the completed file

```text
results/degree3_one_sided_torsion_basis.sobj
```

exists, rerunning Cell 31 automatically builds the two parked sector matrices,
computes their exact combined rank over `QQ`, and reports how many of the 21
shared coefficients remain undetermined.
