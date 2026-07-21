# -*- coding: utf-8 -*-
# =============================================================================
# CELL 37B: DEGREE-THREE LINEAR CURVATURE RESPONSE AND SUPPORT CENSUS
# =============================================================================
#
# Input
# -----
#   results/degree3_cm_jets_cell37a.sobj
#   results/degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj
#
# Purpose
# -------
# 1. Reconstruct the common classical/degree-one HSS geometry as a Laurent
#    jet at the frame divisor U=3X+1=0.
# 2. Insert each of the fifteen observable degree-three symmetric CM-jet
#    tensors, separately in the two parked torsion sectors.
# 3. Compute the linear p^3 curvature-defect response
#
#       G_j(U) = Delta S_curv,3 - alpha0 sigma_3 - beta0 D_3.
#
# 4. Record the Laurent support and the numerical rank of the fifteen shared
#    ambiguity columns.
# 5. Build a *maximal support-matched* convolution module from the two known
#    leading channels.  This is only an absorption diagnostic, not yet the
#    physically admissible degree-three connection.  Cell 37C must derive the
#    actual pole/support restrictions before inverse claims are made.
#
# The p^3 linear response depends only on K0, K1 and the varied K3 column.
# Fixed K2 terms cancel from the derivative with respect to the degree-three
# ambiguity.  In particular,
#
#   Delta G^{-1}_3 = -G0^{-1} G3 G0^{-1},
#
# and the only denominator variation is
#
#   Delta Den_3 = ge3*go0 + ge0*go3.
# =============================================================================

from sage.all import *
from math import factorial, comb
import os
import sys
import gc

print("="*79, flush=True)
print("CELL 37B: DEGREE-THREE LINEAR CURVATURE RESPONSE", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

PREC_BITS = ZZ(sys.argv[1]) if len(sys.argv) > 1 else ZZ(256)
LAURENT_PREC = ZZ(sys.argv[2]) if len(sys.argv) > 2 else ZZ(32)
DISPLAY_TOL_DIGITS = ZZ(sys.argv[3]) if len(sys.argv) > 3 else ZZ(50)

if PREC_BITS < 128:
    raise ValueError("Use at least 128 bits of precision")
if LAURENT_PREC < 20:
    raise ValueError("Use Laurent precision at least 20")

CELL37A_PATH = os.path.join(RESULTS_DIR, "degree3_cm_jets_cell37a.sobj")
CELL36_PATH = os.path.join(
    RESULTS_DIR,
    "degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj",
)
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_linear_curvature_response_cell37b.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_linear_curvature_response_cell37b.txt")
TEMPLATE_OUT = os.path.join(RESULTS_DIR, "degree3_transport_selectivity_input_cell37b.sobj")
PARTIAL_OUT = os.path.join(RESULTS_DIR, "degree3_linear_curvature_response_cell37b_partial.sobj")

for path in [CELL37A_PATH, CELL36_PATH]:
    if not os.path.exists(path):
        raise IOError("Missing required result file: {}".format(path))

cell37a = load(CELL37A_PATH)
cell36 = load(CELL36_PATH)

if ZZ(cell37a["theta_max"]) < 6:
    raise ArithmeticError("Cell 37A does not contain enough theta jets")
if len(cell37a["observable_labels15"]) != 15:
    raise ArithmeticError("Cell 37A does not contain fifteen observable tensors")
if ZZ(cell36["degree3"]["observable_symmetric_dimension"]) != 15:
    raise ArithmeticError("Cell 36 observable dimension is not fifteen")

CC = ComplexField(PREC_BITS)
RR = RealField(PREC_BITS)
LS = LaurentSeriesRing(CC, names=("u",), default_prec=LAURENT_PREC)
u = LS.gen()

L = CC(cell37a["L"])
X = LS((u-1)/3)
ORDER = 4

leading = cell37a["leading_channel_data"]
B_jets = [CC(value) for value in leading["B_jets"]]
alpha0 = CC(leading["alpha0"])
beta0 = CC(leading["beta0"])

labels15 = list(cell37a["observable_labels15"])
sector_names = ["(omega,omega)", "(omega,omega2)"]
for sector in sector_names:
    if sector not in cell37a["sector_tensors"]:
        raise ArithmeticError("Missing sector tensor block {}".format(sector))

print("\nPART I. INPUT AUDIT", flush=True)
print("-"*79, flush=True)
print("  precision bits          : {}".format(PREC_BITS), flush=True)
print("  Laurent precision       : O(U^{})".format(LAURENT_PREC), flush=True)
print("  observable columns      : {}".format(len(labels15)), flush=True)
print("  parked sectors          : {}".format(sector_names), flush=True)
print("  alpha0                  : {}".format(alpha0), flush=True)
print("  beta0                   : {}".format(beta0), flush=True)


def total_degree(index):
    return int(index[0]) + int(index[1]) + int(index[2])


class Jet:
    def __init__(self, coeffs=None):
        self.c = {}
        if coeffs:
            for key, value in coeffs.items():
                kk = tuple(int(v) for v in key)
                vv = LS(value)
                if total_degree(kk) <= ORDER and vv != 0:
                    self.c[kk] = vv

    @staticmethod
    def const(value):
        return Jet({(0,0,0): LS(value)})

    @staticmethod
    def var(index, base):
        key = [0,0,0]
        key[index] = 1
        return Jet({(0,0,0): LS(base), tuple(key): LS(1)})

    def scale(self, value):
        value = LS(value)
        return Jet({key:value*entry for key,entry in self.c.items()})

    def __add__(self, other):
        other = tojet(other)
        out = dict(self.c)
        for key,value in other.c.items():
            out[key] = out.get(key, LS(0)) + value
            if out[key] == 0:
                del out[key]
        return Jet(out)

    __radd__ = __add__

    def __neg__(self):
        return Jet({key:-value for key,value in self.c.items()})

    def __sub__(self, other):
        return self + (-tojet(other))

    def __rsub__(self, other):
        return tojet(other) - self

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return self.scale(other)
        out = {}
        for a,va in self.c.items():
            for b,vb in other.c.items():
                key = (a[0]+b[0], a[1]+b[1], a[2]+b[2])
                if total_degree(key) <= ORDER:
                    out[key] = out.get(key, LS(0)) + va*vb
        return Jet(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Jet):
            return self*other.inv()
        return self.scale(LS(1)/LS(other))

    def __rtruediv__(self, other):
        return tojet(other)*self.inv()

    def __pow__(self, exponent):
        exponent = int(exponent)
        if exponent < 0:
            return self.inv()^(-exponent)
        if exponent == 0:
            return Jet.const(1)
        result = Jet.const(1)
        base = self
        n = exponent
        while n:
            if n & 1:
                result = result*base
            base = base*base
            n >>= 1
        return result

    def inv(self):
        a0 = self.c.get((0,0,0), LS(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant jet coefficient")
        relative = (self-Jet.const(a0)).scale(LS(1)/a0)
        result = Jet.const(0)
        term = Jet.const(1)
        sign = LS(1)
        for _ in range(ORDER+1):
            result += term.scale(sign)
            term = term*relative
            sign = -sign
        return result.scale(LS(1)/a0)

    def log(self):
        a0 = self.c.get((0,0,0), LS(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant jet coefficient")
        relative = (self-Jet.const(a0)).scale(LS(1)/a0)
        result = Jet.const(0)
        term = Jet.const(1)
        for power in range(1,ORDER+1):
            term = term*relative
            result += term.scale(CC((-1)^(power+1))/CC(power))
        return result

    def deriv_value(self, counts):
        counts = tuple(int(v) for v in counts)
        coefficient = self.c.get(counts, LS(0))
        multiplier = ZZ(1)
        for count in counts:
            multiplier *= factorial(count)
        return LS(multiplier)*coefficient


def tojet(value):
    return value if isinstance(value,Jet) else Jet.const(value)


def mat_inv_3(M):
    a,b,c = M[0]
    d,e,f = M[1]
    g,h,i = M[2]
    determinant = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    return [
        [(e*i-f*h)/determinant, (c*h-b*i)/determinant, (b*f-c*e)/determinant],
        [(f*g-d*i)/determinant, (a*i-c*g)/determinant, (c*d-a*f)/determinant],
        [(d*h-e*g)/determinant, (b*g-a*h)/determinant, (a*e-b*d)/determinant],
    ]


def matmul(A,B):
    return [[sum((A[r][k]*B[k][c] for k in range(3)),LS(0))
             for c in range(3)] for r in range(3)]


def matscale(A,scalar):
    return [[scalar*A[r][c] for c in range(3)] for r in range(3)]


def dot(A,left,right=None):
    if right is None:
        right = left
    return sum((left[r]*A[r][c]*right[c]
                for r in range(3) for c in range(3)),LS(0))


def normal_deriv_value(jet, indices):
    return jet.deriv_value((indices.count(0),indices.count(1),indices.count(2)))


def twisted_deriv_value(jet, indices, p_degree):
    aa = indices.count(0)
    bb = indices.count(1)
    cc = indices.count(2)
    return sum((
        LS(comb(aa,power))*(-LS(p_degree)*LS(L))^(aa-power)
        *jet.deriv_value((power,bb,cc))
        for power in range(aa+1)
    ),LS(0))


def one_variable_jet(jets, coordinate_index, theta_shift=0):
    out = Jet.const(0)
    key = [0,0,0]
    key[coordinate_index] = 1
    displacement = Jet({tuple(key):LS(1)})
    for m in range(ORDER+1):
        out += ((-LS(L))^m*LS(jets[m+theta_shift])/LS(factorial(m)))*(displacement^m)
    return out


def tensor_jet(tensor, shift_y=0, shift_z=0):
    """Taylor jet from theta_y/theta_z tensor entries."""
    dy = Jet({(0,1,0):LS(1)})
    dz = Jet({(0,0,1):LS(1)})
    out = Jet.const(0)
    for a in range(ORDER+1):
        for b in range(ORDER+1-a):
            key = (ZZ(a+shift_y), ZZ(b+shift_z))
            if key not in tensor:
                raise ArithmeticError("Missing tensor jet {}".format(key))
            coefficient = (
                (-LS(L))^(a+b)*LS(CC(tensor[key]))
                /LS(factorial(a)*factorial(b))
            )
            out += coefficient*(dy^a)*(dz^b)
    return out


print("\nPART II. COMMON K0/K1 LAURENT GEOMETRY", flush=True)
print("-"*79, flush=True)

xj = Jet.var(0, X)
yj = Jet.var(1, LS(1))
zj = Jet.var(2, LS(1))

Vjet = 9*xj*yj*zj + LS(3)/2*yj*yj*zj + LS(3)/2*yj*zj*zj
K0jet = -Vjet.log()

By = one_variable_jet(B_jets, 1, 0)
Bz = one_variable_jet(B_jets, 2, 0)
Ty = one_variable_jet(B_jets, 1, 1)
Tz = one_variable_jet(B_jets, 2, 1)

K1bar = (
    By*Bz*(1+LS(L)*xj)
    + LS(L)*yj*Ty*Bz
    + LS(L)*zj*By*Tz
)/(2*LS(L)^3*Vjet)

V0 = Vjet.deriv_value((0,0,0))
G0 = [[normal_deriv_value(K0jet,[r,c]) for c in range(3)] for r in range(3)]
G1 = [[twisted_deriv_value(K1bar,[r,c],1) for c in range(3)] for r in range(3)]
M = [[normal_deriv_value(Vjet,[r,c]) for c in range(3)] for r in range(3)]

edir = [-(2*X+1),LS(1),LS(1)]
odir = [LS(0),LS(1),LS(-1)]
rdir = [X,LS(1),LS(1)]

mu_e_1 = dot(G1,edir)/dot(M,edir)
mu_o_1 = dot(G1,odir)/dot(M,odir)
mu_r_1 = dot(G1,rdir)/dot(M,rdir)
sigma1 = mu_o_1-mu_e_1
D1 = (mu_e_1+mu_o_1)/2 + 2*mu_r_1

A0 = []
for q in range(3):
    value = LS(0)
    for index in range(3):
        value += edir[index]*(
            normal_deriv_value(K0jet,[index,1,q])
            -normal_deriv_value(K0jet,[index,2,q])
        )
    A0.append(value)

K4_0 = LS(0)
for first in range(3):
    for second in range(3):
        pattern = [first,second]
        K4_0 += edir[first]*edir[second]*(
            normal_deriv_value(K0jet,pattern+[1,1])
            -2*normal_deriv_value(K0jet,pattern+[1,2])
            +normal_deriv_value(K0jet,pattern+[2,2])
        )

Ginv0 = mat_inv_3(G0)
Curv0 = -K4_0 + dot(Ginv0,A0)
ge0 = dot(G0,edir)
go0 = dot(G0,odir)
Den0 = ge0*go0

print("  common geometry built at U=0", flush=True)
print("  ord_U(V)      = {}".format(V0.valuation()), flush=True)
print("  ord_U(sigma1) = {}".format(sigma1.valuation()), flush=True)
print("  ord_U(D1)     = {}".format(D1.valuation()), flush=True)


def make_K3(tensor):
    F3 = tensor_jet(tensor, 0, 0)
    thetaY = tensor_jet(tensor, 1, 0)
    thetaZ = tensor_jet(tensor, 0, 1)
    return (
        F3*(1+3*LS(L)*xj)
        + LS(L)*yj*thetaY
        + LS(L)*zj*thetaZ
    )/(2*LS(L)^3*Vjet)


def order3_linear_data(K3bar):
    G3 = [[twisted_deriv_value(K3bar,[r,c],3) for c in range(3)] for r in range(3)]

    mu_e_3 = dot(G3,edir)/dot(M,edir)
    mu_o_3 = dot(G3,odir)/dot(M,odir)
    mu_r_3 = dot(G3,rdir)/dot(M,rdir)
    sigma3 = mu_o_3-mu_e_3
    D3 = (mu_e_3+mu_o_3)/2 + 2*mu_r_3

    A3 = []
    for q in range(3):
        value = LS(0)
        for index in range(3):
            value += edir[index]*(
                twisted_deriv_value(K3bar,[index,1,q],3)
                -twisted_deriv_value(K3bar,[index,2,q],3)
            )
        A3.append(value)

    K4_3 = LS(0)
    for first in range(3):
        for second in range(3):
            pattern = [first,second]
            K4_3 += edir[first]*edir[second]*(
                twisted_deriv_value(K3bar,pattern+[1,1],3)
                -2*twisted_deriv_value(K3bar,pattern+[1,2],3)
                +twisted_deriv_value(K3bar,pattern+[2,2],3)
            )

    delta_Ginv3 = matscale(matmul(matmul(Ginv0,G3),Ginv0), LS(-1))
    delta_Curv3 = (
        -K4_3
        + dot(delta_Ginv3,A0)
        + dot(Ginv0,A3,A0)
        + dot(Ginv0,A0,A3)
    )

    ge3 = dot(G3,edir)
    go3 = dot(G3,odir)
    delta_Den3 = ge3*go0 + ge0*go3

    delta_Bcurv3 = delta_Curv3/Den0 - Curv0*delta_Den3/Den0^2
    delta_Scurv3 = -delta_Bcurv3/V0
    defect = delta_Scurv3 - LS(alpha0)*sigma3 - LS(beta0)*D3

    return {
        "Scurv3": delta_Scurv3,
        "sigma3": sigma3,
        "D3": D3,
        "defect": defect,
    }


def coefficient_dict(series, exp_min, exp_max):
    return {ZZ(e): CC(series[e]) for e in range(ZZ(exp_min), ZZ(exp_max)+1)}


def numerical_support(series, exp_min, exp_max, tol):
    nonzero = []
    for exponent in range(ZZ(exp_min), ZZ(exp_max)+1):
        if abs(CC(series[exponent])) > tol:
            nonzero.append(ZZ(exponent))
    if not nonzero:
        return None
    return (min(nonzero), max(nonzero))


def numerical_rank_columns(columns, tolerance):
    """Modified Gram-Schmidt rank after independent column normalization."""
    if not columns:
        return ZZ(0), []
    vectors = []
    accepted_norms = []
    for raw in columns:
        v = vector(CC, raw)
        norm0 = sqrt(sum((abs(entry)^2 for entry in v), RR(0)))
        if norm0 == 0:
            continue
        v = v/CC(norm0)
        for q in vectors:
            inner = sum((q[j].conjugate()*v[j] for j in range(len(v))), CC(0))
            v -= inner*q
        norm = sqrt(sum((abs(entry)^2 for entry in v), RR(0)))
        if norm > tolerance:
            vectors.append(v/CC(norm))
            accepted_norms.append(RR(norm))
    return ZZ(len(vectors)), accepted_norms


print("\nPART III. FIFTEEN LINEAR RESPONSE COLUMNS", flush=True)
print("-"*79, flush=True)

# Wide census range.  The upper edge is kept below the Laurent precision to
# avoid reading coefficients whose precision may have been lost.
EXP_MIN = ZZ(-10)
EXP_MAX = ZZ(min(12, LAURENT_PREC-8))
zero_tol = RR(10)^(-min(ZZ(DISPLAY_TOL_DIGITS), ZZ(PREC_BITS//5)))
rank_tol = RR(10)^(-min(ZZ(40), ZZ(PREC_BITS//6)))

responses = {sector:{} for sector in sector_names}
if os.path.exists(PARTIAL_OUT):
    partial = load(PARTIAL_OUT)
    compatible = (
        ZZ(partial.get("precision_bits", -1)) == PREC_BITS
        and ZZ(partial.get("laurent_precision", -1)) == LAURENT_PREC
        and tuple(partial.get("census_range", (None,None))) == (EXP_MIN,EXP_MAX)
        and list(partial.get("observable_labels15", [])) == labels15
    )
    if compatible:
        stored = partial.get("responses", {})
        for sector in sector_names:
            responses[sector].update(stored.get(sector, {}))
        completed = sum(len(responses[sector]) for sector in sector_names)
        print("  resumed {} / 30 completed response columns".format(completed), flush=True)
    else:
        print("  ignoring incompatible partial checkpoint", flush=True)

all_supports = []
for sector in sector_names:
    print("  sector {}".format(sector), flush=True)
    for idx,label in enumerate(labels15):
        if label in responses[sector]:
            support = responses[sector][label]["support"]
            if support is not None:
                all_supports.append(tuple(support))
            print("    {:2d}/15 {:45s} support {} [cached]".format(idx+1,label,support), flush=True)
            continue
        tensor = cell37a["sector_tensors"][sector][label]
        K3 = make_K3(tensor)
        data = order3_linear_data(K3)
        support = numerical_support(data["defect"], EXP_MIN, EXP_MAX, zero_tol)
        if support is not None:
            all_supports.append(support)
        responses[sector][label] = {
            "defect": coefficient_dict(data["defect"], EXP_MIN, EXP_MAX),
            "sigma3": coefficient_dict(data["sigma3"], EXP_MIN, EXP_MAX),
            "D3": coefficient_dict(data["D3"], EXP_MIN, EXP_MAX),
            "support": support,
        }
        save({
            "schema_version": ZZ(1),
            "precision_bits": PREC_BITS,
            "laurent_precision": LAURENT_PREC,
            "census_range": (EXP_MIN,EXP_MAX),
            "observable_labels15": labels15,
            "responses": responses,
        }, PARTIAL_OUT)
        print("    {:2d}/15 {:45s} support {} [saved]".format(idx+1,label,support), flush=True)
        del K3, data
        gc.collect()

if not all_supports:
    raise ArithmeticError("All degree-three response columns vanished")
global_min = min(item[0] for item in all_supports)
global_max = max(item[1] for item in all_supports)
row_exponents = list(range(global_min, global_max+1))
row_labels = [(sector, ZZ(e)) for sector in sector_names for e in row_exponents]

G_columns = []
for label in labels15:
    column = []
    for sector in sector_names:
        column.extend(responses[sector][label]["defect"][ZZ(e)] for e in row_exponents)
    G_columns.append(column)

rank_G, rank_G_residuals = numerical_rank_columns(G_columns, rank_tol)
print("\n  common numerical support window: [{} , {}]".format(global_min,global_max), flush=True)
print("  stacked response rows          : {}".format(len(row_labels)), flush=True)
print("  numerical rank of G            : {} / 15".format(rank_G), flush=True)

print("\nPART IV. MAXIMAL SUPPORT-MATCHED ABSORPTION DIAGNOSTIC", flush=True)
print("-"*79, flush=True)

# This deliberately large transfer module is a diagnostic only.  It asks
# whether arbitrary Laurent rows in the observed support can be generated by
# shifts of the two leading channels.  It is not claimed to be the physical
# p^3 connection module.
channel_sigma = coefficient_dict(sigma1, -2, 1)
channel_D = coefficient_dict(D1, -2, -1)
alpha_powers = list(range(global_min+2, global_max))   # sigma support [-2,1]
beta_powers = list(range(global_min+2, global_max+2)) # D support [-2,-1]

T_columns = []
T_labels = []
for sector in sector_names:
    sector_row_offset = sector_names.index(sector)*len(row_exponents)
    for power in alpha_powers:
        col = [CC(0)]*len(row_labels)
        for row_index,e in enumerate(row_exponents):
            coefficient = channel_sigma.get(ZZ(e-power), CC(0))
            col[sector_row_offset+row_index] = coefficient
        T_columns.append(col)
        T_labels.append((sector,"alpha",ZZ(power)))
    for power in beta_powers:
        col = [CC(0)]*len(row_labels)
        for row_index,e in enumerate(row_exponents):
            coefficient = channel_D.get(ZZ(e-power), CC(0))
            col[sector_row_offset+row_index] = coefficient
        T_columns.append(col)
        T_labels.append((sector,"beta",ZZ(power)))

rank_T, _ = numerical_rank_columns(T_columns, rank_tol)
rank_TG, _ = numerical_rank_columns(T_columns + G_columns, rank_tol)
maximal_selectivity = ZZ(rank_TG-rank_T)

print("  alpha shift powers      : {}".format(alpha_powers), flush=True)
print("  beta shift powers       : {}".format(beta_powers), flush=True)
print("  maximal-module columns  : {}".format(len(T_columns)), flush=True)
print("  rank(T_max)             : {}".format(rank_T), flush=True)
print("  rank([T_max|G])         : {}".format(rank_TG), flush=True)
print("  support-matched selectivity: {}".format(maximal_selectivity), flush=True)
print("  WARNING: T_max is an absorption diagnostic, not the physical p^3 module.", flush=True)

certificate = {
    "schema_version": ZZ(1),
    "scope": "high-precision degree-three linear curvature-response and Laurent-support census",
    "precision_bits": PREC_BITS,
    "laurent_precision": LAURENT_PREC,
    "response_formula": "Delta Scurv3 - alpha0 sigma3 - beta0 D3",
    "linearization_note": "fixed K2 cancels from variation with respect to K3 ambiguity",
    "observable_labels15": labels15,
    "sector_names": sector_names,
    "response_coefficients": responses,
    "support_window": (global_min,global_max),
    "row_exponents": row_exponents,
    "row_labels": row_labels,
    "G_columns": G_columns,
    "G_rank_numerical": rank_G,
    "rank_tolerance": rank_tol,
    "zero_tolerance": zero_tol,
    "leading_channels": {
        "sigma_support": (-2,1),
        "D_support": (-2,-1),
        "sigma_coefficients": channel_sigma,
        "D_coefficients": channel_D,
    },
    "maximal_support_matched_module": {
        "warning": "diagnostic only; not the physically admissible degree-three connection",
        "alpha_powers": alpha_powers,
        "beta_powers": beta_powers,
        "column_labels": T_labels,
        "columns": T_columns,
        "rank_T": rank_T,
        "rank_TG": rank_TG,
        "selectivity": maximal_selectivity,
    },
    "next_step": "Cell 37C derive admissible p^3 connection support and compute physical selectivity",
}
save(certificate, OUT_SOBJ)

# Populate a standalone selectivity interface without mutating Cell 36.
template = {
    "schema_version": ZZ(1),
    "observable_dimension": ZZ(15),
    "observable_labels": labels15,
    "row_labels": row_labels,
    "G_ambiguity": G_columns,
    "G_rank_numerical": rank_G,
    "T_admissible": None,
    "T_structured": None,
    "selectivity_formula": "rank([T|G])-rank(T)",
    "support_window": (global_min,global_max),
    "leading_channels": certificate["leading_channels"],
    "warning": "T_admissible must be derived geometrically in Cell 37C",
}
save(template, TEMPLATE_OUT)

summary_lines = [
    "CELL 37B: DEGREE-THREE LINEAR CURVATURE RESPONSE",
    "precision bits: {}".format(PREC_BITS),
    "Laurent support census window: [{} , {}]".format(EXP_MIN,EXP_MAX),
    "observed common defect support: [{} , {}]".format(global_min,global_max),
    "stacked parked-sector row count: {}".format(len(row_labels)),
    "numerical rank of fifteen ambiguity responses: {} / 15".format(rank_G),
    "maximal support-matched transfer rank: {}".format(rank_T),
    "rank after adjoining G: {}".format(rank_TG),
    "maximal support-matched selectivity: {}".format(maximal_selectivity),
    "scope warning: maximal module is diagnostic, not the physical p^3 connection",
    "next: derive the admissible p^3 connection support in Cell 37C",
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
print("  {}".format(PARTIAL_OUT), flush=True)
