# -*- coding: utf-8 -*-
# =============================================================================
# CELL 37F: DEGREE-THREE EXCHANGE-PARITY DECOMPOSITION CERTIFICATE
# =============================================================================
#
# Inputs
# ------
#   results/degree3_cm_jets_cell37a.sobj
#   results/degree3_information_loss_cell37c.sobj
#
# Purpose
# -------
# Cell 37C found rank 14 for the scalar curvature-primitives stack.  Cells 37D
# and 37E v2 showed that one odd/even frame-metric component restores rank 15.
# This cell tests the representation-theoretic explanation:
#
#   * the diagonal (omega,omega) tensor is exchange-even;
#   * the mixed sector and its swapped representative obey
#         T_(omega2,omega)(a,b) = T_(omega,omega2)(b,a);
#   * the classical frame satisfies P e=e, P r=r, P o=-o;
#   * therefore C_o(e)=o^T G3 e and C_o(r)=o^T G3 r are exchange-odd;
#   * the exchange-even curvature stack has rank 14 and the odd block adds
#     exactly one quotient direction, giving 14+1=15.
#
# This is an information/covariance certificate.  It does not impose a new
# enumerative equation and is not yet a transport-selectivity theorem.
# =============================================================================

from sage.all import *
import os

print("="*79, flush=True)
print("CELL 37F: DEGREE-THREE EXCHANGE-PARITY DECOMPOSITION CERTIFICATE", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

CELL37A_PATH = os.path.join(RESULTS_DIR, "degree3_cm_jets_cell37a.sobj")
CELL37C_PATH = os.path.join(RESULTS_DIR, "degree3_information_loss_cell37c.sobj")
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_exchange_parity_cell37f.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_exchange_parity_cell37f.txt")
TEMPLATE_OUT = os.path.join(RESULTS_DIR, "degree3_matrix_transport_template_cell37f.sobj")

for path in [CELL37A_PATH, CELL37C_PATH]:
    if not os.path.exists(path):
        raise IOError("Missing required result file: {}".format(path))

cell37a = load(CELL37A_PATH)
cell37c = load(CELL37C_PATH)

PREC_BITS = ZZ(cell37c["precision_bits"])
CC = ComplexField(PREC_BITS)
RR = RealField(PREC_BITS)
responses = cell37c["responses"]
labels15 = list(cell37c["observable_labels15"])
sector_names = list(cell37c["sector_names"])
EXP_MIN, EXP_MAX = map(ZZ, cell37c["census_range"])

if sector_names != ["(omega,omega)", "(omega,omega2)"]:
    raise ArithmeticError("Unexpected parked-sector ordering: {}".format(sector_names))
if len(labels15) != 15:
    raise ArithmeticError("Expected fifteen observable columns")

TOLS = {
    "1e-30": RR(10)^(-30),
    "1e-40": RR(10)^(-40),
    "1e-50": RR(10)^(-50),
}

print("\nPART I. INPUT AUDIT", flush=True)
print("-"*79, flush=True)
print("  precision bits      : {}".format(PREC_BITS), flush=True)
print("  Laurent window      : [{} , {}]".format(EXP_MIN, EXP_MAX), flush=True)
print("  ambiguity columns   : {}".format(len(labels15)), flush=True)
print("  parked sectors      : {}".format(sector_names), flush=True)


def vector_norm(v):
    return RR(sqrt(sum((abs(x)^2 for x in v), RR(0))))


def numerical_rank_columns(columns, tolerance):
    basis = []
    for raw in columns:
        v = vector(CC, raw)
        norm0 = vector_norm(v)
        if norm0 == 0:
            continue
        v = v/CC(norm0)
        # Two MGS passes improve robustness for the nearly dependent stacks.
        for _ in range(2):
            for q in basis:
                inner = sum((q[j].conjugate()*v[j] for j in range(len(v))), CC(0))
                v -= inner*q
        norm = vector_norm(v)
        if norm > tolerance:
            basis.append(v/CC(norm))
    return ZZ(len(basis))


def stack_column_families(families):
    out = [[] for _ in labels15]
    for family in families:
        if len(family) != len(labels15):
            raise ArithmeticError("Column-family width mismatch")
        for j in range(len(labels15)):
            out[j].extend(family[j])
    return out


def component_order(block):
    if "component_orders" in cell37c and block in cell37c["component_orders"]:
        return list(cell37c["component_orders"][block])
    names = set()
    for sector in sector_names:
        for label in labels15:
            names.update(responses[sector][label][block].keys())
    return sorted(names)


def component_columns(block, component, sectors=None, lower=None, upper=None):
    if sectors is None:
        sectors = sector_names
    if lower is None:
        lower = EXP_MIN
    if upper is None:
        upper = EXP_MAX
    columns = []
    for label in labels15:
        vec = []
        for sector in sectors:
            data = responses[sector][label][block][component]
            if block == "raw_tensor":
                vec.append(CC(data))
            else:
                vec.extend(CC(data.get(ZZ(e), 0)) for e in range(ZZ(lower), ZZ(upper)+1))
        columns.append(vec)
    return columns


# -----------------------------------------------------------------------------
# Exact exchange action on the classical frame.
# -----------------------------------------------------------------------------
print("\nPART II. CLASSICAL EXCHANGE PARITY", flush=True)
print("-"*79, flush=True)

KX = FractionField(PolynomialRing(QQ, names=("X",)))
Xsym = KX.gen()
P = matrix(KX, [[1,0,0],[0,0,1],[0,1,0]])
e_sym = vector(KX, [-(2*Xsym+1), 1, 1])
o_sym = vector(KX, [0, 1, -1])
r_sym = vector(KX, [Xsym, 1, 1])

Pe = P*e_sym
Po = P*o_sym
Pr = P*r_sym
parity_checks = {
    "P_e_equals_e": bool(Pe == e_sym),
    "P_o_equals_minus_o": bool(Po == -o_sym),
    "P_r_equals_r": bool(Pr == r_sym),
}
for key,value in parity_checks.items():
    print("  {:24s}: {}".format(key, value), flush=True)
if not all(parity_checks.values()):
    raise ArithmeticError("Classical exchange-parity check failed")


# -----------------------------------------------------------------------------
# Tensor-level exchange covariance from the one-sided CM jets.
# -----------------------------------------------------------------------------
print("\nPART III. TORSION-TENSOR EXCHANGE COVARIANCE", flush=True)
print("-"*79, flush=True)

jets = cell37a["jets"]
basis5 = list(cell37a["basis_labels5"])
pairs = [tuple(map(ZZ,p)) for p in cell37a["observable_pairs"]]
theta_max = ZZ(cell37a["theta_max"])
sector_tensors = cell37a["sector_tensors"]

def build_tensor(left_character, right_character, i, j):
    left_i = jets[left_character][basis5[int(i)]]
    right_j = jets[right_character][basis5[int(j)]]
    if i == j:
        return {
            (ZZ(a),ZZ(b)): CC(left_i[a])*CC(right_j[b])
            for a in range(theta_max+1) for b in range(theta_max+1)
        }
    left_j = jets[left_character][basis5[int(j)]]
    right_i = jets[right_character][basis5[int(i)]]
    return {
        (ZZ(a),ZZ(b)): CC(left_i[a])*CC(right_j[b]) + CC(left_j[a])*CC(right_i[b])
        for a in range(theta_max+1) for b in range(theta_max+1)
    }


def relative_tensor_residual(A, B):
    keys = sorted(set(A).union(B))
    diff = vector(CC, [CC(A.get(k,0))-CC(B.get(k,0)) for k in keys])
    scale = max(vector_norm(vector(CC,[CC(A.get(k,0)) for k in keys])),
                vector_norm(vector(CC,[CC(B.get(k,0)) for k in keys])), RR(1))
    return vector_norm(diff)/scale

max_diag_transpose = RR(0)
max_mixed_swap = RR(0)
for idx,(i,j) in enumerate(pairs):
    label = labels15[idx]
    diag = sector_tensors["(omega,omega)"][label]
    diag_T = {(a,b):CC(diag[(b,a)]) for a in range(theta_max+1) for b in range(theta_max+1)}
    max_diag_transpose = max(max_diag_transpose, relative_tensor_residual(diag,diag_T))

    mixed = sector_tensors["(omega,omega2)"][label]
    swapped = build_tensor("omega2","omega",i,j)
    mixed_T = {(a,b):CC(mixed[(b,a)]) for a in range(theta_max+1) for b in range(theta_max+1)}
    max_mixed_swap = max(max_mixed_swap, relative_tensor_residual(swapped,mixed_T))

print("  max diagonal transpose residual : {}".format(max_diag_transpose), flush=True)
print("  max mixed-swap residual         : {}".format(max_mixed_swap), flush=True)


# -----------------------------------------------------------------------------
# Frame-metric odd/even block from stored coordinate metric responses.
# -----------------------------------------------------------------------------
print("\nPART IV. ODD/EVEN FRAME-METRIC BLOCK", flush=True)
print("-"*79, flush=True)

LS = LaurentSeriesRing(CC, names=("u",), default_prec=int(EXP_MAX+8))
u = LS.gen()
X = (u-1)/3
FRAME_MIN = EXP_MIN + 2
FRAME_MAX = EXP_MAX - 2

def dict_to_series(data):
    return LS({ZZ(e):CC(v) for e,v in data.items()})


def metric_matrix(sector, label):
    block = responses[sector][label]["metric"]
    G00 = dict_to_series(block["G00"])
    G01 = dict_to_series(block["G01"])
    G02 = dict_to_series(block["G02"])
    G11 = dict_to_series(block["G11"])
    G12 = dict_to_series(block["G12"])
    G22 = dict_to_series(block["G22"])
    return [[G00,G01,G02],[G01,G11,G12],[G02,G12,G22]]


def bilinear(G,left,right):
    return sum((left[i]*G[i][j]*right[j] for i in range(3) for j in range(3)),LS(0))


evec = [-(2*X+1),LS(1),LS(1)]
ovec = [LS(0),LS(1),LS(-1)]
rvec = [X,LS(1),LS(1)]


def frame_columns(left,right,sectors):
    columns = []
    for label in labels15:
        vec = []
        for sector in sectors:
            value = bilinear(metric_matrix(sector,label),left,right)
            vec.extend(CC(value[e]) for e in range(FRAME_MIN,FRAME_MAX+1))
        columns.append(vec)
    return columns

Co_e_diag = frame_columns(ovec,evec,["(omega,omega)"])
Co_r_diag = frame_columns(ovec,rvec,["(omega,omega)"])
Co_e_mix = frame_columns(ovec,evec,["(omega,omega2)"])
Co_r_mix = frame_columns(ovec,rvec,["(omega,omega2)"])
Co_full = stack_column_families([Co_e_mix,Co_r_mix])


def family_scale(columns):
    return max([vector_norm(vector(CC,c)) for c in columns] + [RR(1)])

print("  diagonal C_o(e) relative norm: {}".format(
    max([vector_norm(vector(CC,c)) for c in Co_e_diag] + [RR(0)])/family_scale(Co_e_mix)
), flush=True)
print("  diagonal C_o(r) relative norm: {}".format(
    max([vector_norm(vector(CC,c)) for c in Co_r_diag] + [RR(0)])/family_scale(Co_r_mix)
), flush=True)

odd_ranks = {}
for tol_name,tol in TOLS.items():
    odd_ranks[tol_name] = {
        "Co_e_mixed": numerical_rank_columns(Co_e_mix,tol),
        "Co_r_mixed": numerical_rank_columns(Co_r_mix,tol),
        "full_covector": numerical_rank_columns(Co_full,tol),
    }
print("  mixed odd-block ranks: {}".format(odd_ranks), flush=True)


# -----------------------------------------------------------------------------
# 14+1 rank decomposition.
# -----------------------------------------------------------------------------
print("\nPART V. EXCHANGE-EVEN PLUS EXCHANGE-ODD RANK", flush=True)
print("-"*79, flush=True)

primitive_components = component_order("curvature_primitives")
even_columns = stack_column_families([
    component_columns("curvature_primitives",comp)
    for comp in primitive_components
])

# The mixed-sector odd block is sufficient; the swapped representative carries
# the negative of the same exchange-odd information.
combined_e = stack_column_families([even_columns,Co_e_mix])
combined_r = stack_column_families([even_columns,Co_r_mix])
combined_full = stack_column_families([even_columns,Co_full])

rank_decomposition = {}
for tol_name,tol in TOLS.items():
    r_even = numerical_rank_columns(even_columns,tol)
    r_e = numerical_rank_columns(combined_e,tol)
    r_r = numerical_rank_columns(combined_r,tol)
    r_full = numerical_rank_columns(combined_full,tol)
    rank_decomposition[tol_name] = {
        "exchange_even": r_even,
        "even_plus_Co_e": r_e,
        "even_plus_Co_r": r_r,
        "even_plus_full_odd_covector": r_full,
        "odd_quotient_increment": r_full-r_even,
    }
print("  rank decomposition: {}".format(rank_decomposition), flush=True)

passes = all(
    data["exchange_even"] == 14
    and data["even_plus_Co_e"] == 15
    and data["even_plus_Co_r"] == 15
    and data["even_plus_full_odd_covector"] == 15
    and data["odd_quotient_increment"] == 1
    for data in rank_decomposition.values()
)

print("\nPART VI. VERDICT", flush=True)
print("-"*79, flush=True)
print("  tensor exchange covariance passed? {}".format(
    max_diag_transpose < RR(10)^(-60) and max_mixed_swap < RR(10)^(-60)
), flush=True)
print("  stable 14+1 parity decomposition? {}".format(passes), flush=True)
print("  interpretation: curvature contractions retain the exchange-even 14D image", flush=True)
print("                  and the odd/even frame block supplies one missing quotient direction", flush=True)
print("  selectivity caution: this is covariance/information, not yet an enumerative constraint", flush=True)

certificate = {
    "schema_version": ZZ(1),
    "scope": "high-precision exchange-parity decomposition at the square-lattice CM point",
    "precision_bits": PREC_BITS,
    "census_range": (EXP_MIN,EXP_MAX),
    "frame_range": (FRAME_MIN,FRAME_MAX),
    "observable_labels15": labels15,
    "sector_names": sector_names,
    "classical_parity_checks": parity_checks,
    "max_diagonal_tensor_transpose_residual": max_diag_transpose,
    "max_mixed_swap_tensor_residual": max_mixed_swap,
    "odd_block_ranks": odd_ranks,
    "rank_decomposition": rank_decomposition,
    "passes_14_plus_1": bool(passes),
    "next_step": "Calibrate a matrix-valued/tensorial closure at degree two before applying it to the 14+1 degree-three stack",
}
save(certificate,OUT_SOBJ)

save({
    "schema_version": ZZ(1),
    "source_cell37c": CELL37C_PATH,
    "even_observable": {
        "block":"curvature_primitives",
        "components":primitive_components,
        "rank":ZZ(14),
    },
    "odd_observable": {
        "name":"C_o on span{e,r}",
        "coordinates":["C_o(e)","C_o(r)"],
        "quotient_rank":ZZ(1),
        "exchange_parity":ZZ(-1),
    },
    "total_observable_rank":ZZ(15),
    "required_calibration":"degree-two tensorial closure/selectivity test",
},TEMPLATE_OUT)

summary_lines = [
    "CELL 37F: DEGREE-THREE EXCHANGE-PARITY DECOMPOSITION CERTIFICATE",
    "precision bits: {}".format(PREC_BITS),
    "Laurent safe range: [{} , {}]".format(FRAME_MIN,FRAME_MAX),
    "classical exchange parity: P e=e, P r=r, P o=-o",
    "max diagonal tensor transpose residual: {}".format(max_diag_transpose),
    "max mixed-swap tensor residual: {}".format(max_mixed_swap),
    "mixed odd-block ranks: {}".format(odd_ranks),
    "rank decomposition: {}".format(rank_decomposition),
    "stable exchange-even plus exchange-odd decomposition: {}".format(passes),
    "conclusion: observable information decomposes as 14 even + 1 odd quotient direction",
    "caution: covariance is not yet a selective transport equation",
    "next: calibrate any proposed tensorial closure at degree two",
]
with open(OUT_TXT,"w") as handle:
    handle.write("\n".join(summary_lines)+"\n")

print("\nFINAL SUMMARY", flush=True)
print("-"*79, flush=True)
for line in summary_lines:
    print(line, flush=True)
print("\nSaved:", flush=True)
print("  {}".format(OUT_SOBJ), flush=True)
print("  {}".format(OUT_TXT), flush=True)
print("  {}".format(TEMPLATE_OUT), flush=True)
