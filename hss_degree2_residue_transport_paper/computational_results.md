# Computational progression

## 1. Second-order operator envelope

The complete symmetric radial/theta operator envelope through differential order two was tested on 21 radial points. It produced three new quotient directions but did not close:

```text
full rank = 7
augmented rank = 8
```

This ended the search for a missing fixed potential/operator channel.

## 2. Evolving spectral connection

Allowing `alpha_1(X)` and `beta_1(X)` to lie in a Laurent module based on `U=3X+1` closed on 33 points:

```text
rank(A) = rank([A|b]) = 10
```

A matched module based on `W=2X+1` failed:

```text
rank(A) = 10
rank([A|b]) = 11
```

## 3. Fixed geometries

For the primitive product, primitive-plus-double-cover, and `Z2(y)Z2(z)` geometries, a common connection parallel to the leading vector failed, but independent `U`-Laurent functions closed on 50 points. The same minimal support was selected for all three geometries.

## 4. Affine geometry law

The nine connection coefficients were reconstructed as affine functions of `(a,b,c)` and validated on 70 radial points and additional rational geometries. Only `alpha_U` and `alpha_U2` were geometry independent. The geometry-slope map had rank three.

## 5. Exact residue mechanism

The standalone certificate proved:

```text
U = 3X+1 is the volume, metric, and frame-collapse divisor.
ord_U(sigma_bar) = ord_U(D_bar) = -2.
beta_-3 = rho_3 alpha_-3 exactly.
```

The channels are finite Laurent polynomials.

## 6. Full exact formal-CM geometry

The final program rebuilt the level-one and level-two theta jets from Ramanujan's equations, reconstructed the complete order-p^2 geometry, solved the exact seven-equation Laurent system, and checked the four affine rational identities by exact degree-bounded polynomial evaluation.

Final certificate:

```text
intercept Rquad: True (degree <= 13)
slope rA: True (degree <= 10)
slope rB: True (degree <= 10)
slope rC: True (degree <= 10)
```
