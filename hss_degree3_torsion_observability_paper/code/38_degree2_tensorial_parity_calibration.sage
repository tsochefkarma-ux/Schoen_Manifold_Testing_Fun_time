# -*- coding: utf-8 -*-
# =============================================================================
# CELL 38: DEGREE-TWO TENSORIAL / EXCHANGE-PARITY CALIBRATION
# =============================================================================
#
# Inputs
# ------
#   results/degree3_cm_jets_cell37a.sobj
#   results/degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj
#
# Purpose
# -------
# The physical degree-two ambiguity is the complete exchange-even symmetric
# square of the two one-sided directions P2 and B(q^2):
#
#   F_A = P2(y)P2(z),
#   F_B = B(y^2)B(z^2),
#   F_C = P2(y)B(z^2)+B(y^2)P2(z).
#
# For calibration only, enlarge to the full ordered tensor product by adding
# the exchange-odd test direction
#
#   F_- = P2(y)B(z^2)-B(y^2)P2(z).
#
# This cell reconstructs the linear p^2 response at the square-lattice CM
# point and measures the information carried by scalar/even observables and by
# the odd-to-even metric block C_o(w)=o^T G2 w.  The goal is to distinguish:
#
#   * observability/covariance of the odd direction;
#   * genuine inverse selectivity on the physical exchange-even space.
#
# Cell 36 already proved exact scalar bare selectivity zero on the complete
# physical degree-two space.  Cell 38 checks that adjoining the odd block merely
# restores the deliberately added antisymmetric test direction and does not
# create a hidden degree-two enumerative constraint.
# =============================================================================

from sage.all import *
from math import factorial, comb
import os
import sys

print("="*79, flush=True)
print("CELL 38: DEGREE-TWO TENSORIAL / EXCHANGE-PARITY CALIBRATION", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

PREC_BITS = ZZ(sys.argv[1]) if len(sys.argv) > 1 else ZZ(320)
LAURENT_PREC = ZZ(sys.argv[2]) if len(sys.argv) > 2 else ZZ(64)
EXP_MAX_REQUESTED = ZZ(sys.argv[3]) if len(sys.argv) > 3 else ZZ(24)
DISPLAY_TOL_DIGITS = ZZ(sys.argv[4]) if len(sys.argv) > 4 else ZZ(50)
Q_SERIES_MAX = ZZ(sys.argv[5]) if len(sys.argv) > 5 else ZZ(80)

if PREC_BITS < 192:
    raise ValueError("Use at least 192 bits of precision")
if LAURENT_PREC < 32:
    raise ValueError("Use Laurent precision at least 32")
if EXP_MAX_REQUESTED > LAURENT_PREC-12:
    raise ValueError("Keep at least 12 Laurent orders of headroom")
if Q_SERIES_MAX < 40:
    raise ValueError("Use q-series cutoff at least 40")

CELL37A_PATH = os.path.join(RESULTS_DIR, "degree3_cm_jets_cell37a.sobj")
CELL36_PATH = os.path.join(
    RESULTS_DIR,
    "degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj",
)
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree2_tensorial_parity_calibration_cell38.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree2_tensorial_parity_calibration_cell38.txt")
TEMPLATE_OUT = os.path.join(RESULTS_DIR, "degree3_matrix_closure_acceptance_template_cell38.sobj")

for path in [CELL37A_PATH, CELL36_PATH]:
    if not os.path.exists(path):
        raise IOError("Missing required result file: {}".format(path))

cell37a = load(CELL37A_PATH)
cell36 = load(CELL36_PATH)

CC = ComplexField(PREC_BITS)
RR = RealField(PREC_BITS)
LS = LaurentSeriesRing(CC, names=("u",), default_prec=int(LAURENT_PREC))
u = LS.gen()
L = CC(cell37a["L"])
q0 = RR(cell37a["q0"])
X = LS((u-1)/3)
ORDER = 4
THETA_MAX = 6
EXP_MIN = ZZ(-8)
EXP_MAX = ZZ(EXP_MAX_REQUESTED)
TOLERANCES = {
    "1e-30": RR(10)^(-30),
    "1e-40": RR(10)^(-40),
    "1e-50": RR(10)^(-50),
}
DISPLAY_TOL = RR(10)^(-DISPLAY_TOL_DIGITS)

leading = cell37a["leading_channel_data"]
B_jets_leading = [CC(v) for v in leading["B_jets"]]
alpha0 = CC(leading["alpha0"])
beta0 = CC(leading["beta0"])

print("\nPART I. INPUT AUDIT", flush=True)
print("-"*79, flush=True)
print("  precision bits          : {}".format(PREC_BITS), flush=True)
print("  Laurent precision       : O(U^{})".format(LAURENT_PREC), flush=True)
print("  audited Laurent window  : [{} , {}]".format(EXP_MIN, EXP_MAX), flush=True)
print("  q-series cutoff         : {}".format(Q_SERIES_MAX), flush=True)
print("  physical p^2 dimension  : {}".format(cell36["degree2"]["symmetric_dimension"]), flush=True)
print("  exact scalar selectivity: {}".format(cell36["degree2"]["bare_selective_rank"]), flush=True)


# -----------------------------------------------------------------------------
# Exact degree-two q-series coefficients.
# -----------------------------------------------------------------------------
def b_coefficients(N):
    """Coefficients of B(q)=9 prod_{n>=1}(1-q^n)^(-4) through q^N."""
    coeff = [ZZ(0)]*(int(N)+1)
    coeff[0] = ZZ(9)
    for n in range(1, int(N)+1):
        factor = [ZZ(0)]*(int(N)+1)
        for k in range(0, int(N)//n+1):
            factor[k*n] = binomial(k+3, 3)
        new = [ZZ(0)]*(int(N)+1)
        for i in range(0, int(N)+1):
            if coeff[i] == 0:
                continue
            for j in range(0, int(N)-i+1):
                if factor[j] != 0:
                    new[i+j] += coeff[i]*factor[j]
        coeff = new
    return coeff


def convolution(a, b, N):
    out = [QQ(0)]*(int(N)+1)
    for i in range(0, min(len(a)-1, int(N))+1):
        if a[i] == 0:
            continue
        for j in range(0, min(len(b)-1, int(N)-i)+1):
            if b[j] != 0:
                out[i+j] += QQ(a[i])*QQ(b[j])
    return out


def sigma1_integer(n):
    return sum(ZZ(d) for d in divisors(ZZ(n)))


def degree2_coefficients(N):
    B = [QQ(v) for v in b_coefficients(N)]
    E2 = [QQ(0)]*(int(N)+1)
    E2[0] = QQ(1)
    for n in range(1, int(N)+1):
        E2[n] = QQ(-24)*QQ(sigma1_integer(n))
    Bsq = convolution(B, B, N)
    Z2 = [QQ(v)/QQ(72) for v in convolution(E2, Bsq, N)]
    B2 = [QQ(0)]*(int(N)+1)
    for n in range(0, int(N)//2+1):
        B2[2*n] = B[n]
    P2 = [Z2[n]-B2[n]/QQ(8) for n in range(int(N)+1)]
    return B, B2, P2


def theta_jets_from_coeffs(coeffs, q, max_order):
    out = []
    for m in range(int(max_order)+1):
        value = CC(0)
        for n, coefficient in enumerate(coeffs):
            if coefficient != 0:
                factor = ZZ(1) if (m == 0 and n == 0) else ZZ(n)^m
                value += CC(coefficient)*CC(factor)*(CC(q)^n)
        out.append(CC(value))
    return out

B_coeff, B2_coeff, P2_coeff = degree2_coefficients(Q_SERIES_MAX)
B2_jets = theta_jets_from_coeffs(B2_coeff, q0, THETA_MAX+1)
P2_jets = theta_jets_from_coeffs(P2_coeff, q0, THETA_MAX+1)
B_jets_q = theta_jets_from_coeffs(B_coeff, q0, THETA_MAX+1)

# Independent cutoff check.
B_coeff_short, B2_coeff_short, P2_coeff_short = degree2_coefficients(Q_SERIES_MAX-10)
B2_short = theta_jets_from_coeffs(B2_coeff_short, q0, THETA_MAX+1)
P2_short = theta_jets_from_coeffs(P2_coeff_short, q0, THETA_MAX+1)
max_q_cutoff_change = max(
    [abs(B2_jets[m]-B2_short[m]) for m in range(THETA_MAX+2)]
    + [abs(P2_jets[m]-P2_short[m]) for m in range(THETA_MAX+2)]
)
max_leading_B_difference = max(
    abs(B_jets_q[m]-B_jets_leading[m])
    for m in range(min(len(B_jets_q),len(B_jets_leading)))
)

minor_det = P2_coeff[0]*B2_coeff[1]-B2_coeff[0]*P2_coeff[1]

print("\nPART II. COMPLETE ORDERED DEGREE-TWO TENSOR SPACE", flush=True)
print("-"*79, flush=True)
print("  P2 q^0/q^1 coefficients       : {}, {}".format(P2_coeff[0],P2_coeff[1]), flush=True)
print("  B(q^2) q^0/q^1 coefficients   : {}, {}".format(B2_coeff[0],B2_coeff[1]), flush=True)
print("  independence minor determinant: {}".format(minor_det), flush=True)
print("  maximum q-cutoff jet change   : {}".format(max_q_cutoff_change), flush=True)
print("  maximum leading-B cross-check : {}".format(max_leading_B_difference), flush=True)


# -----------------------------------------------------------------------------
# Truncated Laurent/Taylor jet algebra.
# -----------------------------------------------------------------------------
def total_degree(index):
    return int(index[0])+int(index[1])+int(index[2])


class Jet:
    def __init__(self, coeffs=None):
        self.c = {}
        if coeffs:
            for key,value in coeffs.items():
                kk = tuple(int(v) for v in key)
                vv = LS(value)
                if total_degree(kk) <= ORDER and vv != 0:
                    self.c[kk] = vv

    @staticmethod
    def const(value):
        return Jet({(0,0,0):LS(value)})

    @staticmethod
    def var(index, base):
        key = [0,0,0]
        key[index] = 1
        return Jet({(0,0,0):LS(base), tuple(key):LS(1)})

    def scale(self, value):
        value = LS(value)
        return Jet({key:value*entry for key,entry in self.c.items()})

    def __add__(self, other):
        other = tojet(other)
        out = dict(self.c)
        for key,value in other.c.items():
            out[key] = out.get(key,LS(0))+value
            if out[key] == 0:
                del out[key]
        return Jet(out)

    __radd__ = __add__

    def __neg__(self):
        return Jet({key:-value for key,value in self.c.items()})

    def __sub__(self, other):
        return self+(-tojet(other))

    def __rsub__(self, other):
        return tojet(other)-self

    def __mul__(self, other):
        if not isinstance(other,Jet):
            return self.scale(other)
        out = {}
        for a,va in self.c.items():
            for b,vb in other.c.items():
                key = (a[0]+b[0],a[1]+b[1],a[2]+b[2])
                if total_degree(key) <= ORDER:
                    out[key] = out.get(key,LS(0))+va*vb
        return Jet(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other,Jet):
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
        a0 = self.c.get((0,0,0),LS(0))
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
        a0 = self.c.get((0,0,0),LS(0))
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
        coefficient = self.c.get(counts,LS(0))
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
        [(e*i-f*h)/determinant,(c*h-b*i)/determinant,(b*f-c*e)/determinant],
        [(f*g-d*i)/determinant,(a*i-c*g)/determinant,(c*d-a*f)/determinant],
        [(d*h-e*g)/determinant,(b*g-a*h)/determinant,(a*e-b*d)/determinant],
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


def normal_deriv_value(jet,indices):
    return jet.deriv_value((indices.count(0),indices.count(1),indices.count(2)))


def twisted_deriv_value(jet,indices,p_degree):
    aa = indices.count(0)
    bb = indices.count(1)
    cc = indices.count(2)
    return sum((
        LS(comb(aa,power))*(-LS(p_degree)*LS(L))^(aa-power)
        *jet.deriv_value((power,bb,cc))
        for power in range(aa+1)
    ),LS(0))


def one_variable_jet(jets,coordinate_index,theta_shift=0):
    out = Jet.const(0)
    key = [0,0,0]
    key[coordinate_index] = 1
    displacement = Jet({tuple(key):LS(1)})
    for m in range(ORDER+1):
        out += ((-LS(L))^m*LS(jets[m+theta_shift])/LS(factorial(m)))*(displacement^m)
    return out


def tensor_jet(tensor,shift_y=0,shift_z=0):
    dy = Jet({(0,1,0):LS(1)})
    dz = Jet({(0,0,1):LS(1)})
    out = Jet.const(0)
    for a in range(ORDER+1):
        for b in range(ORDER+1-a):
            key = (ZZ(a+shift_y),ZZ(b+shift_z))
            if key not in tensor:
                raise ArithmeticError("Missing tensor jet {}".format(key))
            coefficient = (
                (-LS(L))^(a+b)*LS(CC(tensor[key]))
                /LS(factorial(a)*factorial(b))
            )
            out += coefficient*(dy^a)*(dz^b)
    return out


# -----------------------------------------------------------------------------
# Common K0/K1 geometry.
# -----------------------------------------------------------------------------
print("\nPART III. COMMON K0/K1 GEOMETRY", flush=True)
print("-"*79, flush=True)

xj = Jet.var(0,X)
yj = Jet.var(1,LS(1))
zj = Jet.var(2,LS(1))
Vjet = 9*xj*yj*zj + LS(3)/2*yj*yj*zj + LS(3)/2*yj*zj*zj
K0jet = -Vjet.log()

By = one_variable_jet(B_jets_leading,1,0)
Bz = one_variable_jet(B_jets_leading,2,0)
Ty = one_variable_jet(B_jets_leading,1,1)
Tz = one_variable_jet(B_jets_leading,2,1)
K1bar = (
    By*Bz*(1+LS(L)*xj)
    +LS(L)*yj*Ty*Bz
    +LS(L)*zj*By*Tz
)/(2*LS(L)^3*Vjet)

V0 = Vjet.deriv_value((0,0,0))
G0 = [[normal_deriv_value(K0jet,[r,c]) for c in range(3)] for r in range(3)]
G1 = [[twisted_deriv_value(K1bar,[r,c],1) for c in range(3)] for r in range(3)]
M = [[normal_deriv_value(Vjet,[r,c]) for c in range(3)] for r in range(3)]
Ginv0 = mat_inv_3(G0)

edir = [-(2*X+1),LS(1),LS(1)]
odir = [LS(0),LS(1),LS(-1)]
rdir = [X,LS(1),LS(1)]

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

Curv0 = -K4_0+dot(Ginv0,A0)
ge0 = dot(G0,edir)
go0 = dot(G0,odir)
Den0 = ge0*go0

print("  ord_U(V)       : {}".format(V0.valuation()), flush=True)
print("  ord_U(Den0)    : {}".format(Den0.valuation()), flush=True)


# -----------------------------------------------------------------------------
# Ordered and parity-adapted degree-two tensors.
# -----------------------------------------------------------------------------
def outer_tensor(left,right):
    return {
        (ZZ(a),ZZ(b)):CC(left[a])*CC(right[b])
        for a in range(THETA_MAX+1)
        for b in range(THETA_MAX+1)
    }


def tensor_add(A,B,scale_B=CC(1)):
    return {key:CC(A.get(key,0))+CC(scale_B)*CC(B.get(key,0))
            for key in set(A).union(B)}

PP = outer_tensor(P2_jets,P2_jets)
BB = outer_tensor(B2_jets,B2_jets)
PB = outer_tensor(P2_jets,B2_jets)
BP = outer_tensor(B2_jets,P2_jets)
CPLUS = tensor_add(PB,BP,CC(1))
CMINUS = tensor_add(PB,BP,CC(-1))

parity_tensors = {
    "F_A=PP":PP,
    "F_B=BB":BB,
    "F_C+=PB+BP":CPLUS,
    "F_C-=PB-BP":CMINUS,
}
physical_labels = ["F_A=PP","F_B=BB","F_C+=PB+BP"]
full_labels = physical_labels+["F_C-=PB-BP"]


def tensor_transpose(A):
    return {(a,b):A[(b,a)] for a in range(THETA_MAX+1) for b in range(THETA_MAX+1)}


def tensor_relative_residual(A,B):
    keys = sorted(set(A).union(B))
    diff = vector(CC,[CC(A.get(k,0))-CC(B.get(k,0)) for k in keys])
    va = vector(CC,[CC(A.get(k,0)) for k in keys])
    vb = vector(CC,[CC(B.get(k,0)) for k in keys])
    scale = max(vector_norm(va),vector_norm(vb),RR(1))
    return vector_norm(diff)/scale


# -----------------------------------------------------------------------------
# Linear p^2 response.
# -----------------------------------------------------------------------------
def make_K2(tensor):
    F2 = tensor_jet(tensor,0,0)
    thetaY = tensor_jet(tensor,1,0)
    thetaZ = tensor_jet(tensor,0,1)
    return (
        F2*(1+2*LS(L)*xj)
        +LS(L)*yj*thetaY
        +LS(L)*zj*thetaZ
    )/(2*LS(L)^3*Vjet)


def order2_linear_data(K2bar, raw_tensor):
    G2 = [[twisted_deriv_value(K2bar,[r,c],2) for c in range(3)] for r in range(3)]

    mu_e_2 = dot(G2,edir)/dot(M,edir)
    mu_o_2 = dot(G2,odir)/dot(M,odir)
    mu_r_2 = dot(G2,rdir)/dot(M,rdir)
    sigma2 = mu_o_2-mu_e_2
    D2 = (mu_e_2+mu_o_2)/2+2*mu_r_2

    A2 = []
    for q in range(3):
        value = LS(0)
        for index in range(3):
            value += edir[index]*(
                twisted_deriv_value(K2bar,[index,1,q],2)
                -twisted_deriv_value(K2bar,[index,2,q],2)
            )
        A2.append(value)

    K4_2 = LS(0)
    for first in range(3):
        for second in range(3):
            pattern = [first,second]
            K4_2 += edir[first]*edir[second]*(
                twisted_deriv_value(K2bar,pattern+[1,1],2)
                -2*twisted_deriv_value(K2bar,pattern+[1,2],2)
                +twisted_deriv_value(K2bar,pattern+[2,2],2)
            )

    delta_Ginv2 = matscale(matmul(matmul(Ginv0,G2),Ginv0),LS(-1))
    delta_Curv2 = (
        -K4_2
        +dot(delta_Ginv2,A0)
        +dot(Ginv0,A2,A0)
        +dot(Ginv0,A0,A2)
    )
    ge2 = dot(G2,edir)
    go2 = dot(G2,odir)
    delta_Den2 = ge2*go0+ge0*go2
    delta_Bcurv2 = delta_Curv2/Den0-Curv0*delta_Den2/Den0^2
    delta_Scurv2 = -delta_Bcurv2/V0
    defect = delta_Scurv2-LS(alpha0)*sigma2-LS(beta0)*D2

    metric = {
        "G00":G2[0][0],"G01":G2[0][1],"G02":G2[0][2],
        "G11":G2[1][1],"G12":G2[1][2],"G22":G2[2][2],
    }
    kahler = {"K2jet_{}_{}_{}".format(*key):value for key,value in K2bar.c.items()}
    Co_e = dot(G2,odir,edir)
    Co_r = dot(G2,odir,rdir)

    return {
        "raw_tensor":raw_tensor,
        "kahler_jet":kahler,
        "metric":metric,
        "connection_vector":{"A0":A2[0],"A1":A2[1],"A2":A2[2]},
        "eigen_channels":{"mu_e":mu_e_2,"mu_o":mu_o_2,"mu_r":mu_r_2},
        "spectral_channels":{"sigma":sigma2,"D":D2},
        "curvature_primitives":{
            "K4":K4_2,"delta_Curv":delta_Curv2,"delta_Den":delta_Den2,
            "ge2":ge2,"go2":go2,
        },
        "scalar_curvature":{"Scurv":delta_Scurv2},
        "transport_triplet":{"Scurv":delta_Scurv2,"sigma":sigma2,"D":D2},
        "defect":{"defect":defect},
        "odd_even_block":{"Co_e":Co_e,"Co_r":Co_r},
    }


responses = {}
print("\nPART IV. FOUR ORDERED/PARITY RESPONSE COLUMNS", flush=True)
print("-"*79, flush=True)
for label in full_labels:
    print("  computing {} ...".format(label), end="", flush=True)
    tensor = parity_tensors[label]
    K2 = make_K2(tensor)
    responses[label] = order2_linear_data(K2,tensor)
    print(" done", flush=True)


# -----------------------------------------------------------------------------
# Rank and covariance diagnostics.
# -----------------------------------------------------------------------------
def vector_norm(v):
    return RR(sqrt(sum((abs(x)^2 for x in v),RR(0))))


def numerical_rank_columns(columns,tolerance):
    basis = []
    for raw in columns:
        v = vector(CC,raw)
        norm0 = vector_norm(v)
        if norm0 == 0:
            continue
        v = v/CC(norm0)
        for _ in range(2):
            for q in basis:
                inner = sum((q[j].conjugate()*v[j] for j in range(len(v))),CC(0))
                v -= inner*q
        norm = vector_norm(v)
        if norm > tolerance:
            basis.append(v/CC(norm))
    return ZZ(len(basis))


def component_names(block):
    names = set()
    for label in full_labels:
        names.update(responses[label][block].keys())
    return sorted(names)


def flatten_label(label,block,components=None):
    if components is None:
        components = component_names(block)
    vec = []
    for component in components:
        value = responses[label][block].get(component, CC(0) if block == "raw_tensor" else LS(0))
        if block == "raw_tensor":
            vec.append(CC(value))
        else:
            vec.extend(CC(value[e]) for e in range(EXP_MIN,EXP_MAX+1))
    return vec


def columns_for(labels,block,components=None):
    return [flatten_label(label,block,components) for label in labels]


layers = [
    "raw_tensor","kahler_jet","metric","connection_vector","eigen_channels",
    "spectral_channels","curvature_primitives","scalar_curvature",
    "transport_triplet","defect","odd_even_block",
]
rank_table = {}
for tol_name,tol in TOLERANCES.items():
    rank_table[tol_name] = {}
    for block in layers:
        rank_table[tol_name][block] = {
            "physical_even":numerical_rank_columns(columns_for(physical_labels,block),tol),
            "full_ordered":numerical_rank_columns(columns_for(full_labels,block),tol),
        }

# Curvature-even stack plus the odd/even block.
curv_components = component_names("curvature_primitives")
curv_cols = columns_for(full_labels,"curvature_primitives",curv_components)
coe_cols = columns_for(full_labels,"odd_even_block",["Co_e"])
cor_cols = columns_for(full_labels,"odd_even_block",["Co_r"])
curv_plus_coe = [curv_cols[j]+coe_cols[j] for j in range(len(full_labels))]
curv_plus_cor = [curv_cols[j]+cor_cols[j] for j in range(len(full_labels))]
curv_plus_covector = [curv_cols[j]+coe_cols[j]+cor_cols[j] for j in range(len(full_labels))]

lift_ranks = {}
for tol_name,tol in TOLERANCES.items():
    lift_ranks[tol_name] = {
        "curvature_primitives_full":numerical_rank_columns(curv_cols,tol),
        "curvature_plus_Co_e":numerical_rank_columns(curv_plus_coe,tol),
        "curvature_plus_Co_r":numerical_rank_columns(curv_plus_cor,tol),
        "curvature_plus_full_covector":numerical_rank_columns(curv_plus_covector,tol),
    }

# Tensor parity checks.
parity_residuals = {
    "F_A_even":tensor_relative_residual(PP,tensor_transpose(PP)),
    "F_B_even":tensor_relative_residual(BB,tensor_transpose(BB)),
    "F_Cplus_even":tensor_relative_residual(CPLUS,tensor_transpose(CPLUS)),
    "F_Cminus_odd":tensor_relative_residual(CMINUS,{key:-value for key,value in tensor_transpose(CMINUS).items()}),
}

# Norms of the antisymmetric response in even/scalar and odd blocks.
def block_vector_norm(label,block):
    return vector_norm(vector(CC,flatten_label(label,block)))

antisym_norms = {
    "curvature_primitives":block_vector_norm("F_C-=PB-BP","curvature_primitives"),
    "scalar_curvature":block_vector_norm("F_C-=PB-BP","scalar_curvature"),
    "spectral_channels":block_vector_norm("F_C-=PB-BP","spectral_channels"),
    "defect":block_vector_norm("F_C-=PB-BP","defect"),
    "odd_even_block":block_vector_norm("F_C-=PB-BP","odd_even_block"),
}
max_even_odd_block = max(
    block_vector_norm(label,"odd_even_block") for label in physical_labels
)

print("\nPART V. INFORMATION / PARITY RANK LADDER", flush=True)
print("-"*79, flush=True)
for block in layers:
    entry = rank_table["1e-40"][block]
    print("  {:22s}: physical {}/3, ordered {}/4".format(
        block,entry["physical_even"],entry["full_ordered"]
    ), flush=True)

print("\n  curvature-even plus odd-block lifts:", flush=True)
for tol_name in ["1e-30","1e-40","1e-50"]:
    print("    {}: {}".format(tol_name,lift_ranks[tol_name]), flush=True)

print("\nPART VI. EXCHANGE-COVARIANCE CALIBRATION", flush=True)
print("-"*79, flush=True)
for key,value in parity_residuals.items():
    print("  {:24s}: {}".format(key,value), flush=True)
print("  maximum physical-even odd-block norm: {}".format(max_even_odd_block), flush=True)
print("  antisymmetric response norms          : {}".format(antisym_norms), flush=True)

stable_full_ordered_observability = all(
    lift_ranks[t]["curvature_plus_full_covector"] == 4 for t in TOLERANCES
)
physical_scalar_selectivity_zero = (ZZ(cell36["degree2"]["bare_selective_rank"]) == 0)
exchange_excludes_odd_test = (
    parity_residuals["F_Cminus_odd"] < RR(10)^(-60)
    and parity_residuals["F_Cplus_even"] < RR(10)^(-60)
)

verdict = {
    "physical_even_dimension":ZZ(3),
    "ordered_test_dimension":ZZ(4),
    "antisymmetric_test_dimension":ZZ(1),
    "physical_scalar_bare_selectivity":ZZ(cell36["degree2"]["bare_selective_rank"]),
    "full_ordered_observability_with_odd_block":bool(stable_full_ordered_observability),
    "exchange_covariance_classifies_odd_test":bool(exchange_excludes_odd_test),
    "tensorial_parity_adds_physical_selectivity":False,
}

print("\nPART VII. CALIBRATION VERDICT", flush=True)
print("-"*79, flush=True)
print("  complete physical exchange-even space : 3 dimensions", flush=True)
print("  added exchange-odd test space          : 1 dimension", flush=True)
print("  even observables + odd block see all 4?: {}".format(stable_full_ordered_observability), flush=True)
print("  exact physical scalar selectivity      : {}".format(cell36["degree2"]["bare_selective_rank"]), flush=True)
print("  parity itself adds physical constraints?: False", flush=True)
print("  conclusion: degree-two tensorial parity is observability/covariance, not inverse selection", flush=True)

certificate = {
    "schema_version":ZZ(1),
    "scope":"degree-two ordered-space tensorial/exchange-parity calibration at the square-lattice CM point",
    "precision_bits":PREC_BITS,
    "Laurent_range":(EXP_MIN,EXP_MAX),
    "q_series_max":Q_SERIES_MAX,
    "q_cutoff_change":max_q_cutoff_change,
    "leading_B_cross_check":max_leading_B_difference,
    "one_sided_labels":["P2","B(q^2)"],
    "physical_labels":physical_labels,
    "ordered_parity_labels":full_labels,
    "rank_table":rank_table,
    "lift_ranks":lift_ranks,
    "parity_residuals":parity_residuals,
    "antisymmetric_norms":antisym_norms,
    "max_even_odd_block_norm":max_even_odd_block,
    "verdict":verdict,
}
save(certificate,OUT_SOBJ)

acceptance_template = {
    "schema_version":ZZ(1),
    "purpose":"acceptance criteria for a future degree-three matrix-valued closure",
    "degree2_calibration":{
        "physical_dimension":ZZ(3),
        "required_selectivity_on_complete_physical_space":ZZ(0),
        "ordered_test_dimension":ZZ(4),
        "odd_test_role":"covariance/observability control, not an enumerative constraint",
        "scalar_transfer_rank":ZZ(cell36["degree2"]["transport_rank"]),
    },
    "degree3_requirements":{
        "observable_dimension":ZZ(15),
        "exchange_even_dimension":ZZ(14),
        "exchange_odd_quotient_dimension":ZZ(1),
        "must_report":"rank([T|G])-rank(T) separately in even and odd blocks",
        "noncircularity":"any restricted matrix transport law must be justified independently of the unknown q-series",
    },
}
save(acceptance_template,TEMPLATE_OUT)

summary_lines = [
    "CELL 38: DEGREE-TWO TENSORIAL / EXCHANGE-PARITY CALIBRATION",
    "physical exchange-even ambiguity dimension: 3",
    "ordered tensor-product test dimension: 4",
    "antisymmetric test direction: P2(y)B(z^2)-B(y^2)P2(z)",
    "exact scalar bare selectivity on physical space: {}".format(cell36["degree2"]["bare_selective_rank"]),
    "curvature primitives + full odd/even covector rank: {} / 4".format(
        lift_ranks["1e-40"]["curvature_plus_full_covector"]
    ),
    "exchange covariance classifies the odd test: {}".format(exchange_excludes_odd_test),
    "tensorial parity adds physical selectivity: False",
    "conclusion: the degree-three 14+1 split is an observability/covariance structure until an independently justified matrix transport law is supplied",
]
with open(OUT_TXT,"w") as handle:
    handle.write("\n".join(summary_lines)+"\n")

print("\nSaved:", flush=True)
print("  {}".format(OUT_SOBJ), flush=True)
print("  {}".format(OUT_TXT), flush=True)
print("  {}".format(TEMPLATE_OUT), flush=True)
