# Cell 35b — Torsion modularity audit

This cell replaces the provisional group assumption in Cell 35 with an explicit audit.

## What it proves

It distinguishes the global Jacobi form from its torsion specialization:

\[
\Phi(\tau,0)=12096E_4\Delta\ne0,
\]

so the relation is **not** an identity in the full weight-16, index-3 E8 Jacobi-form space.

After restricting to the HSS line `z=s*gamma`, the scalar Jacobi index is

\[
3\,Q(\gamma)=3\cdot4=12.
\]

For `s=(tau+r)/3`, the intrinsic torsion correction is

\[
\exp(2\pi i(4/3)\tau).
\]

With `tau=3t`, this is exactly `q_t^4`. It is an exponential translation correction, not multiplication by Delta or eta, so the modular weight remains 16.

The code derives the common stabilizer of the three pairs `(1,r) mod 3` as `Gamma(3)` in the `tau` variable. It also verifies the conjugate group in the `t` variable as

\[
\Gamma_0(9)\cap\Gamma_1(3)=\Gamma_H(9,\langle4\rangle),
\]

with index 24.

The cusp width at infinity for `Gamma(3)` is 3, so the local parameter is

\[
\exp(2\pi i\tau/3)=q_t.
\]

The standard torsion-specialization theorem for holomorphic Jacobi forms supplies holomorphy at every cusp and permits a finite character. The code does not assume that character is trivial. It uses the conservative bound

\[
\left\lfloor16[\mathrm{SL}_2(\mathbb Z):\Gamma(3)]/12\right\rfloor=32.
\]

Only certified Cell-32 coefficients are used; the unsafe convolution buffer is explicitly excluded.

## Install

Copy:

```text
35b_degree3_torsion_modularity_audit.sage
```

to:

```text
/home/will/research/code/
```

Copy:

```text
run_cell35b_torsion_modularity_audit.sh
```

to:

```text
/home/will/research/
```

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell35b_torsion_modularity_audit.sh
bash run_cell35b_torsion_modularity_audit.sh
```

Or directly:

```bash
sage code/35b_degree3_torsion_modularity_audit.sage
```

## Input

```text
results/degree3_A3_B3_hss_torsion_candidate.sobj
```

The current certified cutoff `q^40` is sufficient.

## Outputs

```text
results/degree3_torsion_modularity_audit_35b.sobj
results/degree3_torsion_modularity_audit_35b.txt
```

## Expected final line

```text
CELL 35b AUDIT: PASS
```

## Mathematical dependency

The all-orders conclusion invokes the standard theorem that a holomorphic Jacobi form evaluated at a rational torsion section, with its intrinsic exponential correction, is a holomorphic modular form of the same weight on the torsion stabilizer, possibly with finite character. The code audits every problem-specific input to that theorem and the subsequent Sturm calculation.
