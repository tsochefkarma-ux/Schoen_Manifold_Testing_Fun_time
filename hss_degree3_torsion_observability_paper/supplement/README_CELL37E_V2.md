# Cell 37E v2 — parity-block mixing audit

Cell 37E v1 assumed that the classical odd/radial frame diagonalized the leading generalized pencil `(G1,M)`. The local run rejected that assumption: the odd vector is an eigendirection to numerical precision, but the classical radial vector is not.

Cell 37E v2 makes the weaker, correct statement:

- the odd line is orthogonal to and decoupled from the even plane `span{e,r}` at leading order;
- the classical `e,r` basis may mix inside that even plane;
- the degree-three perturbation defines the covector
  `C_o(w)=o^T G3 w` on the even plane;
- `g_or=C_o(r)` is one coordinate of this covector, not a canonical eigenvector-rotation coefficient;
- curvature primitives plus the full odd/even mixing covector retain rank 15.

Copy:

- `37e_degree3_parity_block_mixing_certificate_v2.sage` to `/home/will/research/code/`
- `run_cell37e_parity_block_v2.sh` to `/home/will/research/`

Run:

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell37e_parity_block_v2.sh
bash run_cell37e_parity_block_v2.sh
```

Expected qualitative result:

```text
odd line / even plane invariant: True
classical radial vector is generalized eigenvector: False
curvature-primitives rank: 14 / 15
curvature primitives + full odd/even covector: 15 / 15
```
