# -*- coding: utf-8 -*-
# =============================================================================
# CELL 37E: DEGREE-THREE SPECTRAL-FRAME ROTATION CERTIFICATE
# =============================================================================
#
# Inputs
# ------
#   results/degree3_cm_jets_cell37a.sobj
#   results/degree3_information_loss_cell37c.sobj
#   results/degree3_minimal_observable_lift_cell37d.sobj
#
# Purpose
# -------
# Cell 37D found that one frame-metric component lifts the rank-fourteen
# curvature-primitives stack to full rank fifteen.  This cell gives the
# preferred odd/radial component a canonical spectral meaning.
#
# The leading generalized spectral pencil is
#
#       G1 v_i = mu_i M v_i,
#
# in the classical frame (e,o,r).  For a p^3 perturbation G3, the first-order
# rotation of the r-vector toward o is
#
#       Omega_or = (o^T G3 r) / ((mu_r-mu_o) o^T M o).
#
# It is characterized by the projected linearized eigenvector equation
#
#       o^T[ G3 r + Omega_or (G1-mu_r M)o ] = 0.
#
# This is a tensorial covariance identity, not an additional enumerative
# constraint.  It identifies the missing local observable that a future
# matrix-valued transport law must carry.
# =============================================================================

from sage.all import *
from math import factorial, comb
import os

print("="*79, flush=True)
print("CELL 37E: DEGREE-THREE SPECTRAL-FRAME ROTATION CERTIFICATE", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

CELL37A_PATH = os.path.join(RESULTS_DIR, "degree3_cm_jets_cell37a.sobj")
CELL37C_PATH = os.path.join(RESULTS_DIR, "degree3_information_loss_cell37c.sobj")
CELL37D_PATH = os.path.join(RESULTS_DIR, "degree3_minimal_observable_lift_cell37d.sobj")
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_spectral_frame_rotation_cell37e.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_spectral_frame_rotation_cell37e.txt")
TEMPLATE_OUT = os.path.join(RESULTS_DIR, "degree3_matrix_transport_input_cell37e.sobj")

for path in [CELL37A_PATH, CELL37C_PATH, CELL37D_PATH]:
    if not os.path.exists(path):
        raise IOError("Missing required result file: {}".format(path))

cell37a = load(CELL37A_PATH)
cell37c = load(CELL37C_PATH)
cell37d = load(CELL37D_PATH)

PREC_BITS = ZZ(cell37c["precision_bits"])
CC = ComplexField(PREC_BITS)
RR = RealField(PREC_BITS)
EXP_MIN, EXP_MAX = map(ZZ, cell37c["census_range"])
LAURENT_PREC = max(ZZ(cell37c.get("laurent_precision", EXP_MAX+16)), EXP_MAX+16)
LS = LaurentSeriesRing(CC, names=("u",), default_prec=int(LAURENT_PREC))
u = LS.gen()
X = (u-1)/3
L = CC(cell37a["L"])
responses = cell37c["responses"]
labels15 = list(cell37c["observable_labels15"])
sector_names = list(cell37c["sector_names"])
B_jets = [CC(x) for x in cell37a["leading_channel_data"]["B_jets"]]

if len(labels15) != 15 or len(sector_names) != 2:
    raise ArithmeticError("Unexpected Cell 37C dimensions")
if len(B_jets) < 4:
    raise ArithmeticError("Need at least four leading B-channel jets")

TOLS = {
    "1e-30": RR(10)^(-30),
    "1e-40": RR(10)^(-40),
    "1e-50": RR(10)^(-50),
}

print("\nPART I. INPUT AUDIT", flush=True)
print("-"*79, flush=True)
print("  precision bits        : {}".format(PREC_BITS), flush=True)
print("  Laurent range         : [{} , {}]".format(EXP_MIN, EXP_MAX), flush=True)
print("  ambiguity columns     : {}".format(len(labels15)), flush=True)
print("  parked sectors        : {}".format(sector_names), flush=True)
print("  Cell 37D preferred    : {}".format(cell37d.get("preferred_lift")), flush=True)


def total_degree(index):
    return int(index[0])+int(index[1])+int(index[2])


ORDER = 2


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
        return Jet({k:value*v for k,v in self.c.items()})

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
        return Jet({k:-v for k,v in self.c.items()})

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
        rel = (self-Jet.const(a0)).scale(LS(1)/a0)
        result = Jet.const(0)
        term = Jet.const(1)
        sign = LS(1)
        for _ in range(ORDER+1):
            result += term.scale(sign)
            term = term*rel
            sign = -sign
        return result.scale(LS(1)/a0)

    def log(self):
        a0 = self.c.get((0,0,0),LS(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant jet coefficient")
        rel = (self-Jet.const(a0)).scale(LS(1)/a0)
        result = Jet.const(0)
        term = Jet.const(1)
        for power in range(1,ORDER+1):
            term = term*rel
            result += term.scale(CC((-1)^(power+1))/CC(power))
        return result

    def deriv_value(self, counts):
        counts = tuple(int(v) for v in counts)
        multiplier = ZZ(1)
        for count in counts:
            multiplier *= factorial(count)
        return LS(multiplier)*self.c.get(counts,LS(0))


def tojet(value):
    return value if isinstance(value,Jet) else Jet.const(value)


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


def dot(A,left,right=None):
    if right is None:
        right = left
    return sum((left[r]*A[r][c]*right[c]
                for r in range(3) for c in range(3)),LS(0))


def matvec(A,v):
    return [sum((A[r][c]*v[c] for c in range(3)),LS(0)) for r in range(3)]


def vec_sub(a,b):
    return [a[i]-b[i] for i in range(len(a))]


def vec_scale(a,s):
    return [s*x for x in a]


def max_coeff_abs(series, lower, upper):
    return max([abs(CC(series[e])) for e in range(ZZ(lower),ZZ(upper)+1)]+[RR(0)])


# -----------------------------------------------------------------------------
# Leading generalized spectral pencil.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Leading parity decomposition of the generalized pencil.
# -----------------------------------------------------------------------------
print("\nPART II. LEADING PARITY-BLOCK AUDIT", flush=True)
print("-"*79, flush=True)

xj = Jet.var(0,X)
yj = Jet.var(1,LS(1))
zj = Jet.var(2,LS(1))
Vjet = 9*xj*yj*zj + LS(3)/2*yj*yj*zj + LS(3)/2*yj*zj*zj
K0jet = -Vjet.log()
By = one_variable_jet(B_jets,1,0)
Bz = one_variable_jet(B_jets,2,0)
Ty = one_variable_jet(B_jets,1,1)
Tz = one_variable_jet(B_jets,2,1)
K1bar = (
    By*Bz*(1+LS(L)*xj)
    +LS(L)*yj*Ty*Bz
    +LS(L)*zj*By*Tz
)/(2*LS(L)^3*Vjet)

G1 = [[twisted_deriv_value(K1bar,[r,c],1) for c in range(3)] for r in range(3)]
M = [[normal_deriv_value(Vjet,[r,c]) for c in range(3)] for r in range(3)]

edir = [-(2*X+1),LS(1),LS(1)]
odir = [LS(0),LS(1),LS(-1)]
rdir = [X,LS(1),LS(1)]

SAFE_LOWER = EXP_MIN+2
SAFE_UPPER = EXP_MAX-2
FRAME_TOL = RR(10)^(-60)


def series_size(series):
    return max_coeff_abs(series,SAFE_LOWER,SAFE_UPPER)


def vector_size(vec):
    return max(series_size(x) for x in vec)

Me = dot(M,edir)
Mo = dot(M,odir)
Mr = dot(M,rdir)
Moe = dot(M,odir,edir)
Mor = dot(M,odir,rdir)
Mer = dot(M,edir,rdir)
G1oe = dot(G1,odir,edir)
G1or = dot(G1,odir,rdir)
G1er = dot(G1,edir,rdir)

mu_o = dot(G1,odir)/Mo
res_o = vec_sub(matvec(G1,odir),vec_scale(matvec(M,odir),mu_o))
res_o_size = vector_size(res_o)

# The classical radial vector is deliberately audited rather than assumed to
# be a generalized eigenvector.
mu_r_rayleigh = dot(G1,rdir)/Mr
res_r = vec_sub(matvec(G1,rdir),vec_scale(matvec(M,rdir),mu_r_rayleigh))
res_r_size = vector_size(res_r)

odd_M_decoupling = max(series_size(Moe),series_size(Mor))
odd_G1_decoupling = max(series_size(G1oe),series_size(G1or))
parity_block_pass = bool(
    odd_M_decoupling < FRAME_TOL
    and odd_G1_decoupling < FRAME_TOL
    and res_o_size < FRAME_TOL
)
radial_is_eigenvector = bool(res_r_size < FRAME_TOL)

print("  ord_U(o^T M o)                    : {}".format(Mo.valuation()), flush=True)
print("  max odd/even M coupling           : {}".format(odd_M_decoupling), flush=True)
print("  max odd/even G1 coupling          : {}".format(odd_G1_decoupling), flush=True)
print("  odd generalized-eigen residual    : {}".format(res_o_size), flush=True)
print("  radial Rayleigh residual          : {}".format(res_r_size), flush=True)
print("  odd line / even plane invariant?  : {}".format(parity_block_pass), flush=True)
print("  classical r is an eigenvector?    : {}".format(radial_is_eigenvector), flush=True)
print("  even-block coupling o-free G1(e,r): ord {}".format(G1er.valuation()), flush=True)

if not parity_block_pass:
    raise ArithmeticError("Odd/even parity block failed at leading order")

# The failed v1 interpretation is an expected audit result: the even plane is
# invariant, but the classical e/r basis need not diagonalize its 2x2 block.


# -----------------------------------------------------------------------------
# Stored p^3 metric columns and odd-to-even mixing covector.
# -----------------------------------------------------------------------------
def dict_to_series(data):
    return LS({ZZ(e):CC(v) for e,v in data.items()})


def metric_matrix(sector,label):
    block = responses[sector][label]["metric"]
    G00 = dict_to_series(block["G00"])
    G01 = dict_to_series(block["G01"])
    G02 = dict_to_series(block["G02"])
    G11 = dict_to_series(block["G11"])
    G12 = dict_to_series(block["G12"])
    G22 = dict_to_series(block["G22"])
    return [[G00,G01,G02],[G01,G11,G12],[G02,G12,G22]]


def coefficient_vector(series,lower,upper):
    return [CC(series[e]) for e in range(ZZ(lower),ZZ(upper)+1)]


def family_columns(series_map,lower,upper):
    columns = []
    for label in labels15:
        vec = []
        for sector in sector_names:
            vec.extend(coefficient_vector(series_map[sector][label],lower,upper))
        columns.append(vec)
    return columns


def stack_families(families):
    out = [[] for _ in labels15]
    for family in families:
        for j in range(15):
            out[j].extend(family[j])
    return out


def numerical_rank(columns,tol):
    basis=[]
    for raw in columns:
        v=vector(CC,raw)
        n0=RR(sqrt(sum((abs(x)^2 for x in v),RR(0))))
        if n0 == 0:
            continue
        v=v/CC(n0)
        # Two-pass modified Gram-Schmidt for better conditioning.
        for _ in range(2):
            for q in basis:
                ip=sum((q[k].conjugate()*v[k] for k in range(len(v))),CC(0))
                v-=ip*q
        n=RR(sqrt(sum((abs(x)^2 for x in v),RR(0))))
        if n>tol:
            basis.append(v/CC(n))
    return ZZ(len(basis))


print("\nPART III. ODD-TO-EVEN MIXING COVECTOR", flush=True)
print("-"*79, flush=True)

# C_o(w)=o^T G3 w is a covector on the invariant even plane span{e,r}.
# Dividing by o^T M o fixes the normalization of the odd vector but does not
# choose a preferred basis in the even plane.
c_oe = {sector:{} for sector in sector_names}
c_or = {sector:{} for sector in sector_names}
xi_oe = {sector:{} for sector in sector_names}
xi_or = {sector:{} for sector in sector_names}
covariance_residuals=[]

for sector in sector_names:
    for label in labels15:
        G3 = metric_matrix(sector,label)
        ce = dot(G3,odir,edir)
        cr = dot(G3,odir,rdir)
        c_oe[sector][label]=ce
        c_or[sector][label]=cr
        xi_oe[sector][label]=ce/Mo
        xi_or[sector][label]=cr/Mo

        # Explicit covariance check under e' = e+r, r' = r.
        ep = [edir[i]+rdir[i] for i in range(3)]
        direct = dot(G3,odir,ep)
        covariance_residuals.append(series_size(direct-(ce+cr)))

print("  maximum even-basis covariance residual: {}".format(max(covariance_residuals)), flush=True)
print("  interpretation: (C_o(e),C_o(r)) is an even-plane covector", flush=True)
print("  caution: C_o(r) alone is a coordinate, not a canonical eigen-rotation", flush=True)


# -----------------------------------------------------------------------------
# Rank comparison.
# -----------------------------------------------------------------------------
print("\nPART IV. INFORMATION RANK", flush=True)
print("-"*79, flush=True)

primitive_components = list(cell37c["component_orders"]["curvature_primitives"])
primitive_families=[]
for component in primitive_components:
    cols=[]
    for label in labels15:
        vec=[]
        for sector in sector_names:
            data=responses[sector][label]["curvature_primitives"][component]
            vec.extend(CC(data.get(ZZ(e),0)) for e in range(EXP_MIN,EXP_MAX+1))
        cols.append(vec)
    primitive_families.append(cols)
primitive_columns=stack_families(primitive_families)

mix_lower=SAFE_LOWER
mix_upper=SAFE_UPPER
coe_columns=family_columns(c_oe,mix_lower,mix_upper)
cor_columns=family_columns(c_or,mix_lower,mix_upper)
xie_columns=family_columns(xi_oe,mix_lower,mix_upper)
xir_columns=family_columns(xi_or,mix_lower,mix_upper)
full_covector_columns=stack_families([xie_columns,xir_columns])

rank_families={
    "curvature_primitives":primitive_columns,
    "C_o(e)":coe_columns,
    "C_o(r)":cor_columns,
    "Xi_o(e)":xie_columns,
    "Xi_o(r)":xir_columns,
    "odd_even_covector":full_covector_columns,
    "primitives_plus_C_o(e)":stack_families([primitive_columns,coe_columns]),
    "primitives_plus_C_o(r)":stack_families([primitive_columns,cor_columns]),
    "primitives_plus_odd_even_covector":stack_families([primitive_columns,full_covector_columns]),
}
ranks={}
for name,cols in rank_families.items():
    ranks[name]={key:numerical_rank(cols,tol) for key,tol in TOLS.items()}
    print("  {:36s}: {}".format(name,ranks[name]), flush=True)

for key in TOLS:
    if ranks["curvature_primitives"][key] != 14:
        raise ArithmeticError("Curvature primitive rank is not stably fourteen")
    if ranks["primitives_plus_odd_even_covector"][key] != 15:
        raise ArithmeticError("Odd-even mixing covector does not lift to rank fifteen")


# -----------------------------------------------------------------------------
# Verdict and next interface.
# -----------------------------------------------------------------------------
print("\nPART V. VERDICT", flush=True)
print("-"*79, flush=True)
print("  v1 spectral-rotation assumption rejected: {}".format(not radial_is_eigenvector), flush=True)
print("  leading pencil has an odd line plus an invariant even plane: PASS", flush=True)
print("  curvature primitives + odd/even mixing covector: 15 / 15", flush=True)
print("  physical next object: a parity-block matrix transport law", flush=True)
print("  no unique eigen-rotation coordinate has yet been selected", flush=True)

certificate={
    "schema_version":ZZ(2),
    "scope":"high-precision parity-block and odd-to-even mixing-covector certificate at the square-lattice CM point",
    "source_cell37a":CELL37A_PATH,
    "source_cell37c":CELL37C_PATH,
    "source_cell37d":CELL37D_PATH,
    "precision_bits":PREC_BITS,
    "census_range":(EXP_MIN,EXP_MAX),
    "safe_range":(SAFE_LOWER,SAFE_UPPER),
    "observable_labels15":labels15,
    "sector_names":sector_names,
    "leading_parity_audit":{
        "Me":Me,"Mo":Mo,"Mr":Mr,
        "Moe":Moe,"Mor":Mor,"Mer":Mer,
        "G1oe":G1oe,"G1or":G1or,"G1er":G1er,
        "mu_o":mu_o,
        "odd_eigen_residual":res_o_size,
        "radial_rayleigh_residual":res_r_size,
        "parity_block_pass":parity_block_pass,
        "classical_r_is_eigenvector":radial_is_eigenvector,
    },
    "odd_even_mixing":{
        "C_o_e":c_oe,
        "C_o_r":c_or,
        "Xi_o_e":xi_oe,
        "Xi_o_r":xi_or,
        "definition":"Xi_o(w)=(o^T G3 w)/(o^T M o), w in span{e,r}",
        "basis_covariance_residual_max":max(covariance_residuals),
    },
    "ranks":ranks,
}
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_parity_block_mixing_cell37e_v2.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_parity_block_mixing_cell37e_v2.txt")
TEMPLATE_OUT = os.path.join(RESULTS_DIR, "degree3_matrix_transport_input_cell37e_v2.sobj")
save(certificate,OUT_SOBJ)

save({
    "schema_version":ZZ(2),
    "source_cell37e_v2":OUT_SOBJ,
    "observable_labels15":labels15,
    "sector_names":sector_names,
    "diagonal_curvature_stack":{
        "block":"curvature_primitives",
        "components":primitive_components,
    },
    "off_diagonal_parity_block":{
        "name":"Xi_o",
        "definition":"even-plane covector w -> (o^T G3 w)/(o^T M o)",
        "components_in_classical_basis":{"e":xi_oe,"r":xi_or},
    },
    "full_information_rank":ZZ(15),
    "next_required_step":"construct and calibrate a parity-block matrix transport/covariance law; do not assume the classical radial vector is a generalized eigenvector",
},TEMPLATE_OUT)

summary_lines=[
    "CELL 37E v2: DEGREE-THREE PARITY-BLOCK MIXING CERTIFICATE",
    "precision bits: {}".format(PREC_BITS),
    "Laurent safe range: [{} , {}]".format(SAFE_LOWER,SAFE_UPPER),
    "odd line / even plane invariant: {}".format(parity_block_pass),
    "classical radial vector is generalized eigenvector: {}".format(radial_is_eigenvector),
    "odd generalized-eigen residual: {}".format(res_o_size),
    "radial Rayleigh residual: {}".format(res_r_size),
    "curvature-primitives rank: {} / 15".format(ranks["curvature_primitives"]["1e-40"]),
    "curvature primitives + C_o(e): {} / 15".format(ranks["primitives_plus_C_o(e)"]["1e-40"]),
    "curvature primitives + C_o(r): {} / 15".format(ranks["primitives_plus_C_o(r)"]["1e-40"]),
    "curvature primitives + full odd/even covector: {} / 15".format(ranks["primitives_plus_odd_even_covector"]["1e-40"]),
    "conclusion: the missing information lies in the odd-to-even perturbation block",
    "caution: g_or is one classical-frame coordinate, not yet a canonical spectral rotation",
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
