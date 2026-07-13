# Computational supplement

This directory accompanies **Square-Lattice Spectral--Curvature Closure in an HSS-Motivated Leading Sector**.

## Files

- `01_ramanujan_eta_jets.sage` — checks the eta-product coefficients through order 50 and derives the normalized jets with Ramanujan's equations; verifies `R1=R5=0` exactly.
- `02_finite_jet_obstruction.sage` — constructs the formal curvature obstruction for unrestricted jets, performs the elimination steps, and verifies the sufficient relations discovered in the calculation.
- `03_generic_x_closure_certificate.sage` — imposes `R1=R5=0` and proves that the generic-`x` residual numerator is the zero polynomial.
- `04_sympy_independent_check.py` — independent exact verification of the square-lattice jet relations, exact constants, and generic-`x` closure from the reduced rational formulas.

## Suggested commands

```bash
sage 01_ramanujan_eta_jets.sage
sage 02_finite_jet_obstruction.sage
sage 03_generic_x_closure_certificate.sage
python 04_sympy_independent_check.py
```

The proof steps use exact rational-function arithmetic. Numerical evaluations printed by the Sage scripts are diagnostics only.
