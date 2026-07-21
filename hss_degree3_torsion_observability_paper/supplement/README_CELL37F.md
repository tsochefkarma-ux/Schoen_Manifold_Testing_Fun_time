# Cell 37F — degree-three exchange-parity decomposition

This cell converts the Cell 37E v2 observation into a representation-theoretic audit.

It verifies:

- the classical frame parity under exchange of the two HSS factors,
  `P e=e`, `P r=r`, `P o=-o`;
- transpose symmetry of the `(omega,omega)` CM tensors;
- covariance of the mixed orbit,
  `T_(omega2,omega)(a,b)=T_(omega,omega2)(b,a)`;
- numerical vanishing of the odd/even metric block in the diagonal sector;
- rank 14 for the exchange-even curvature-primitives stack;
- one additional quotient direction from the exchange-odd covector
  `C_o(w)=o^T G3 w`, `w in span{e,r}`;
- stable total rank `14+1=15` at tolerances `1e-30`, `1e-40`, and `1e-50`.

This is an information and covariance result. It does not create a new enumerative constraint. Any matrix-valued closure proposed next must first be calibrated on the complete degree-two ambiguity space, whose scalar bare selectivity is already known to be zero.

## Install

Copy:

- `37f_degree3_exchange_parity_decomposition.sage` to `/home/will/research/code/`
- `run_cell37f_exchange_parity.sh` to `/home/will/research/`

Run:

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell37f_exchange_parity.sh
bash run_cell37f_exchange_parity.sh
```

Inputs:

- `results/degree3_cm_jets_cell37a.sobj`
- `results/degree3_information_loss_cell37c.sobj`

Outputs:

- `results/degree3_exchange_parity_cell37f.sobj`
- `results/degree3_exchange_parity_cell37f.txt`
- `results/degree3_matrix_transport_template_cell37f.sobj`
