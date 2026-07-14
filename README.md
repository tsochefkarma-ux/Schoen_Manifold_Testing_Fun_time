# Square-Lattice Spectral–Curvature Closure in an HSS-Motivated Schoen Sector

This repository contains a three-part computational study of spectral–curvature closure in a symmetry-reduced Hessian model motivated by the one-sectional Hosono–Saito–Stienstra (HSS) sector of Schoen’s Calabi–Yau threefold.

The project follows a single mathematical arc:

1. an exact constant-coefficient closure appears at the square-lattice CM point;
2. the first degree-two correction obstructs that same constant-coefficient closure;
3. the obstruction is resolved by an exact finite Laurent connection over the classical volume divisor.

The main one-variable series is

$$
B(q)=9\prod_{n\ge1}(1-q^n)^{-4}
    =9q^{1/6}\eta(\tau)^{-4},
\qquad q=e^{2\pi i\tau}.
$$

All three papers study the symmetric line and the square-lattice point

$$
v=(X,1,1),\qquad \tau=i.
$$

---

## Papers

### Part I — Exact CM closure

**[Square-Lattice Spectral–Curvature Closure in an HSS-Motivated Leading Sector](square_lattice_spectral_curvature_note.pdf)**

The first paper proves that the leading even–odd primitive bisectional-curvature response closes exactly on two symmetry-adapted spectral-response channels:

$$
\mathcal S_1(X)=\alpha_0\bar{\sigma}(X)+\beta_0\bar{D}(X),
$$

with coefficients independent of the radial modulus $X$.

The closure is driven by the square-lattice fixed-point identities

$$
LE_2(i)=6,\qquad E_6(i)=0,\qquad L=2\pi,
$$

and is certified by exact finite-jet and rational-function calculations.

**Main message:** a special rank reduction occurs at the CM point $\tau=i$.

---

### Part II — Degree-two obstruction

**[Modular Fixed-Point Rank Reduction and Degree-Two Obstructions in an HSS-Motivated Schoen Sector](modular_fixed_point_degree2_obstruction.pdf)**

The second paper studies the rigidity of the leading closure and the first nontrivial degree-two correction.

It proves exact fixed-point ideal membership,

$$
N\in\langle LE_2-6,E_6\rangle,
$$

establishes local rigidity and transversality, and constructs the exact finite-jet obstruction map.

For the natural factorized degree-two family, a high-precision rank calculation gives

$$
\mathrm{rank}(A)=4,
\qquad
\mathrm{rank}([A\mid b])=5,
$$

showing that the leading closure cannot be extended by merely adding constant degree-two coefficients within that ansatz.

**Main message:** the rank-two response system survives, but constant correction coefficients are too restrictive.

> The final factorized degree-two obstruction in Part II is presented as stable high-precision evidence for the stated ansatz. It is not claimed as a theorem about the complete relative Schoen potential.

---

### Part III — Exact residue transport

**[Residue Transport at the Volume Divisor and Exact Degree-Two Closure in an HSS-Motivated Schoen Sector](exact_degree2_residue_transport.pdf)**

The third paper resolves the degree-two obstruction.

The intrinsic coordinate

$$
U=3X+1
$$

is simultaneously:

- the classical volume coordinate, since $V(X,1,1)=3U$;
- the degeneration divisor of the background Hessian metric;
- the collapse divisor of the symmetry-adapted frame.

The two response channels are finite Laurent polynomials in $U$, and the full factorized degree-two defect closes exactly once the degree-two coefficients are promoted to a finite Laurent connection:

$$
R_{\mathrm{quad}} + a r_A + b r_B + c r_C
= \alpha_1(U)\bar{\sigma}(U) + \beta_1(U)\bar{D}(U).
$$

The connection has the form

$$
\alpha_1(U)
=
\alpha_{-3}U^{-3}
+\alpha_{-2}U^{-2}
+\alpha_{-1}U^{-1}
+\alpha_{\mathrm{const}}
+\alpha_U U
+\alpha_{U^2}U^2,
$$

$$
\beta_1(U)
=
\beta_{-3}U^{-3}
+\beta_{-2}U^{-2}
+\beta_{\mathrm{const}}.
$$

The deepest residues satisfy

$$
\beta_{-3}=\rho_3\alpha_{-3},
\qquad
\rho_3=
\frac{2\bigl(L^2E_4(i)+108\bigr)}
{3\bigl(L^2E_4(i)+36\bigr)},
$$

and the universal positive-power drift is

$$
\alpha_{\mathrm{univ}}(U)
=
\frac{\alpha_0C_0^2}{243}
U\left(\frac3L-U\right),
\qquad
\beta_{\mathrm{univ}}(U)=0.
$$

The final exact certificate is

```text
EXACT DEGREE-TWO FAST CERTIFICATE
intercept Rquad: True (degree <= 13)
slope rA: True (degree <= 10)
slope rB: True (degree <= 10)
slope rC: True (degree <= 10)
```

**Main message:** the degree-two obstruction is the signature of residue transport over the unique volume/frame divisor, not the failure of the rank-two response system.

---

## Conceptual progression

| Part | Question | Result | Status |
|---|---|---|---|
| I | Does leading curvature close on the spectral channels at $\tau=i$? | Yes, with constant coefficients | Exact |
| II | Does the same constant-coefficient closure persist at degree two? | No, for the tested natural factorized ansatz | Exact structural results plus high-precision obstruction |
| III | What replaces constant degree-two coefficients? | A finite Laurent connection in $U=3X+1$ | Exact within the stated formal model |

The progression can be summarized as

$$
\text{CM fixed-point closure}
\longrightarrow
\text{constant-coefficient obstruction}
\longrightarrow
\text{exact residue transport}.
$$

---

## Reproducibility

The calculations were performed in **SageMath**, primarily using CoCalc.

Recommended environment:

- SageMath 10.x or a current CoCalc Sage kernel;
- substantial memory for the full formal-CM degree-two reconstruction;
- exact arithmetic enabled throughout the final certificate.

The Part III supplement contains several execution paths:

- a standalone frame and Laurent-matching certificate;
- a full formal-CM reconstruction of the degree-two geometry;
- a checkpointed CoCalc implementation;
- restartable one-column exact verifiers;
- numerical discovery scripts leading to the exact result.

The final rational identities are proved by clearing denominators, obtaining rigorous degree bounds, and checking exact vanishing at sufficiently many distinct rational values of $X$. No floating-point tolerance is used in the final certificate.

---

## Suggested reading order

1. Read Part I for the model, the two response channels, and the original square-lattice closure.
2. Read Part II for fixed-point rigidity, finite-jet deformation theory, and the degree-two obstruction.
3. Read Part III for the volume-divisor geometry, finite Laurent channels, residue recursion, and the exact degree-two theorem.

Readers mainly interested in the final result can begin with the introduction and theorem sections of Part III, then consult Parts I and II for the origin of the channels and the obstruction that motivated the Laurent connection.

---

## Scope

The exact theorems apply to the stated three-dimensional HSS-motivated Hessian model and, at degree two, to the factorized family generated by the components $F_A,F_B,F_C$ used in the papers.

The repository does **not** claim that this factorized family is the complete base-degree-two relative Schoen potential. Part III proves that the complete family considered here has an exact Laurent spectral–curvature closure over the formal square-lattice/level-two CM coefficient field.

---

## Repository layout

A typical layout is:

```text
.
├── square_lattice_spectral_curvature_note.pdf
├── modular_fixed_point_degree2_obstruction.pdf
├── exact_degree2_residue_transport.pdf
├── code/                  # SageMath scripts and certificates
├── supplements/           # LaTeX sources, logs, manifests, and checkpoints
└── README.md
```

Adjust the links above if the PDF filenames in the repository differ.

---

## Citation

These manuscripts are part of one continuous project. Until formal bibliographic records are available, cite the relevant paper by title and repository URL.

A repository-level citation may be written as:

```text
[Author], “Square-Lattice Spectral–Curvature Closure in an HSS-Motivated
Schoen Sector,” Parts I–III, GitHub repository, [year].
```

Individual BibTeX entries can be added once author, year, and archival identifiers are finalized.

---

## Status

This is an active research repository. The main Part III degree-two identity has an exact formal-CM certificate. Future work may include:

- identifying the full global geometric meaning of the Laurent connection;
- comparing the factorized family with the complete relative/contact Schoen potential;
- extending residue transport to higher base degree;
- relating the volume-divisor connection to broader BCOV, quasi-Jacobi, or variation-of-Hodge-structure frameworks.
