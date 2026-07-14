from sage.all import *
from math import factorial, comb
import gc

# ============================================================================
# EXACT DEGREE-TWO VOLUME-DIVISOR DEFECT CERTIFICATE
# ============================================================================
#
# Standalone Sage/SageCell program.
#
# Goal
# ----
# Derive, rather than assume,
#
#   * the Laurent expansions of sigma_1 and D_1 at U=3X+1=0;
#   * the exact degree-two defect components Rquad,rA,rB,rC;
#   * the seven polar/constant connection coefficients after fixing the
#     universal drift
#
#       alpha_U   = alpha0*C0^2/(81 L),
#       alpha_U^2 = -alpha0*C0^2/243;
#
#   * and finally verify the four rational-function identities
#
#       Rquad = response(p0) + universal drift,
#       rA    = response(pA),
#       rB    = response(pB),
#       rC    = response(pC).
#
# If all four cleared numerators vanish, then for arbitrary factorized geometry
# (a,b,c) the exact affine connection
#
#       p(a,b,c) = p0 + a*pA + b*pB + c*pC
#
# closes identically as a rational function of X.  The level-two CM constants
# E2(2i),E4(2i),E6(2i),B(q0^2) are kept formal, while all their theta jets are
# generated exactly by Ramanujan's differential equations.  Thus a successful
# result is stronger than a numerical CM specialization.
#
# Runtime notes
# -------------
# This is an exact symbolic calculation and can be demanding.  SageCell users
# should leave VERIFY_FULL_IDENTITIES=True for the decisive run, but may set it
# False for a quicker first pass through the Laurent solve.
# ============================================================================

print("="*78)
print("EXACT DEGREE-TWO VOLUME-DIVISOR DEFECT CERTIFICATE")
print("Ramanujan CM jets, exact p^2 geometry, Laurent solve, rational identity")
print("="*78)
print()

VERIFY_FULL_IDENTITIES = True
PRINT_FULL_CONNECTION_MAP = False
WRITE_CHECKPOINT = True
CHECKPOINT_NAME = "exact_degree2_connection_certificate.sobj"

ORDER = 4
LAURENT_PREC = 8          # coefficients through U^7 are available
MATCH_EXPONENTS = list(range(-5, 2))  # -5,-4,...,1: seven equations

# ----------------------------------------------------------------------------
# PART I. Coefficient field and Ramanujan-generated CM jets.
# ----------------------------------------------------------------------------

print("PART I. FORMAL CM FIELD AND EXACT THETA JETS")
print("-"*78)

parameter_names = [
    "L", "E4", "C0",
    "C2", "E22", "E42", "E62",
]
Rpar = PolynomialRing(QQ, names=parameter_names)
K = Rpar.fraction_field()
pg = dict(zip(parameter_names, K.gens()))

L = pg["L"]
E4 = pg["E4"]
C0 = pg["C0"]
C2 = pg["C2"]
E22 = pg["E22"]
E42 = pg["E42"]
E62 = pg["E62"]

# Radial rational-function field K(X).
RX = PolynomialRing(K, names=("X",))
X = RX.gen()
FX = RX.fraction_field()

U = FX(3*X + 1)
Ycm = K(L^2*E4)
alpha0 = K((Ycm + 12)/72)
beta0 = K(-(7*Ycm^2 - 552*Ycm - 13392)/(288*(Ycm + 36)))

# Ramanujan differential ring for theta=q d/dq.
Rram = PolynomialRing(QQ, names=("e2", "e4", "e6", "bb"))
e2r, e4r, e6r, bbr = Rram.gens()
ram_derivatives = {
    e2r: (e2r^2 - e4r)/12,
    e4r: (e2r*e4r - e6r)/3,
    e6r: (e2r*e6r - e4r^2)/2,
    bbr: ((1-e2r)/6)*bbr,
}


def theta_ram(poly):
    poly = Rram(poly)
    return Rram(sum(
        poly.derivative(generator)*ram_derivatives[generator]
        for generator in Rram.gens()
    ))


def ram_jets(seed, max_order):
    out = []
    current = Rram(seed)
    for _ in range(max_order + 1):
        out.append(current)
        current = theta_ram(current)
    return out


B_ram = ram_jets(bbr, ORDER + 1)
Z2_ram = ram_jets(e2r*bbr^2/72, ORDER + 1)

square_subs = {
    e2r: K(6/L),
    e4r: E4,
    e6r: K(0),
    bbr: C0,
}
level2_subs = {
    e2r: E22,
    e4r: E42,
    e6r: E62,
    bbr: C2,
}


def coerce_substitution(poly, substitutions):
    """Evaluate a Ramanujan polynomial in K without cross-parent subs()."""
    values = [K(substitutions[generator]) for generator in Rram.gens()]
    result = K(0)
    for exponent_tuple, coefficient in Rram(poly).dict().items():
        term = K(coefficient)
        for index, exponent in enumerate(exponent_tuple):
            if exponent:
                term *= values[index]^exponent
        result += term
    return K(result)


Bjets = [coerce_substitution(poly, square_subs) for poly in B_ram]
Z2jets = [coerce_substitution(poly, square_subs) for poly in Z2_ram]

# For f(q^2), theta_q^m f(q^2)=2^m (theta^m f)(q^2).
B2jets = [K(2^m)*coerce_substitution(B_ram[m], level2_subs)
          for m in range(ORDER + 2)]
P2jets = [K(Z2jets[m] - B2jets[m]/8) for m in range(ORDER + 2)]

print("  B and Z2 jets generated through theta^{}".format(ORDER+1))
print("  level-two B(q^2) jets generated formally from E2(2i),E4(2i),E6(2i)")
print("  alpha0 = {}".format(alpha0))
print("  beta0  = {}".format(beta0))
print()

# ----------------------------------------------------------------------------
# PART II. Exact truncated Taylor-jet algebra over K(X).
# ----------------------------------------------------------------------------

print("PART II. EXACT TRUNCATED GEOMETRY")
print("-"*78)


def total_degree(index):
    return int(index[0]) + int(index[1]) + int(index[2])


class Jet:
    def __init__(self, coeffs=None):
        self.c = {}
        if coeffs:
            for key, value in coeffs.items():
                kk = tuple(int(v) for v in key)
                vv = FX(value)
                if total_degree(kk) <= ORDER and vv != 0:
                    self.c[kk] = vv

    @staticmethod
    def const(value):
        return Jet({(0,0,0): FX(value)})

    @staticmethod
    def var(index, base):
        key = [0,0,0]
        key[index] = 1
        return Jet({(0,0,0): FX(base), tuple(key): FX(1)})

    def scale(self, value):
        value = FX(value)
        return Jet({key:value*entry for key,entry in self.c.items()})

    def __add__(self, other):
        other = tojet(other)
        out = dict(self.c)
        for key,value in other.c.items():
            out[key] = out.get(key, FX(0)) + value
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
                    out[key] = out.get(key, FX(0)) + va*vb
        return Jet(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Jet):
            return self*other.inv()
        return self.scale(FX(1)/FX(other))

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
        a0 = self.c.get((0,0,0), FX(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.inv")
        relative = (self-Jet.const(a0)).scale(FX(1)/a0)
        result = Jet.const(0)
        term = Jet.const(1)
        sign = FX(1)
        for _ in range(ORDER+1):
            result += term.scale(sign)
            term = term*relative
            sign = -sign
        return result.scale(FX(1)/a0)

    def log(self):
        a0 = self.c.get((0,0,0), FX(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.log")
        relative = (self-Jet.const(a0)).scale(FX(1)/a0)
        result = Jet.const(0)
        term = Jet.const(1)
        for power in range(1,ORDER+1):
            term = term*relative
            result += term.scale(FX(QQ((-1)^(power+1))/power))
        return result

    def deriv_value(self, counts):
        counts = tuple(int(v) for v in counts)
        coefficient = self.c.get(counts, FX(0))
        multiplier = ZZ(1)
        for count in counts:
            multiplier *= factorial(count)
        return FX(multiplier)*coefficient


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
    return [[sum((A[r][k]*B[k][c] for k in range(3)),FX(0))
             for c in range(3)] for r in range(3)]


def matadd(A,B):
    return [[A[r][c]+B[r][c] for c in range(3)] for r in range(3)]


def matscale(A,scalar):
    return [[scalar*A[r][c] for c in range(3)] for r in range(3)]


def dot(A,left,right=None):
    if right is None:
        right = left
    return sum((left[r]*A[r][c]*right[c]
                for r in range(3) for c in range(3)),FX(0))


def normal_deriv_value(jet, indices):
    return jet.deriv_value((indices.count(0),indices.count(1),indices.count(2)))


def twisted_deriv_value(jet, indices, p_degree):
    aa = indices.count(0)
    bb = indices.count(1)
    cc = indices.count(2)
    return sum((
        FX(comb(aa,power))*(-FX(p_degree)*FX(L))^(aa-power)
        *jet.deriv_value((power,bb,cc))
        for power in range(aa+1)
    ),FX(0))


def one_variable_jet(jets, coordinate_index):
    out = Jet.const(0)
    key = [0,0,0]
    key[coordinate_index] = 1
    displacement = Jet({tuple(key):FX(1)})
    for m in range(ORDER+1):
        out += ((-FX(L))^m*FX(jets[m])/FX(factorial(m)))*(displacement^m)
    return out


def one_variable_theta_jet(jets, coordinate_index):
    out = Jet.const(0)
    key = [0,0,0]
    key[coordinate_index] = 1
    displacement = Jet({tuple(key):FX(1)})
    for m in range(ORDER+1):
        out += ((-FX(L))^m*FX(jets[m+1])/FX(factorial(m)))*(displacement^m)
    return out


print("  building common K0/K1 geometry...",flush=True)

xj = Jet.var(0,FX(X))
yj = Jet.var(1,FX(1))
zj = Jet.var(2,FX(1))

Vjet = 9*xj*yj*zj + FX(3)/2*yj*yj*zj + FX(3)/2*yj*zj*zj
K0jet = -Vjet.log()

By = one_variable_jet(Bjets,1)
Bz = one_variable_jet(Bjets,2)
Ty = one_variable_theta_jet(Bjets,1)
Tz = one_variable_theta_jet(Bjets,2)

K1bar = (
    By*Bz*(1+FX(L)*xj)
    + FX(L)*yj*Ty*Bz
    + FX(L)*zj*By*Tz
)/(2*FX(L)^3*Vjet)

V0 = Vjet.deriv_value((0,0,0))
G0 = [[normal_deriv_value(K0jet,[r,c]) for c in range(3)] for r in range(3)]
G1 = [[twisted_deriv_value(K1bar,[r,c],1) for c in range(3)] for r in range(3)]
M = [[normal_deriv_value(Vjet,[r,c]) for c in range(3)] for r in range(3)]

edir = [-(2*FX(X)+1),FX(1),FX(1)]
odir = [FX(0),FX(1),FX(-1)]
rdir = [FX(X),FX(1),FX(1)]

mu_e_1 = dot(G1,edir)/dot(M,edir)
mu_o_1 = dot(G1,odir)/dot(M,odir)
mu_r_1 = dot(G1,rdir)/dot(M,rdir)
sigma1 = FX(mu_o_1-mu_e_1)
D1 = FX((mu_e_1+mu_o_1)/2 + 2*mu_r_1)

A0 = []
A1 = []
for q in range(3):
    a0q = FX(0)
    a1q = FX(0)
    for index in range(3):
        a0q += edir[index]*(
            normal_deriv_value(K0jet,[index,1,q])
            -normal_deriv_value(K0jet,[index,2,q])
        )
        a1q += edir[index]*(
            twisted_deriv_value(K1bar,[index,1,q],1)
            -twisted_deriv_value(K1bar,[index,2,q],1)
        )
    A0.append(FX(a0q))
    A1.append(FX(a1q))

K4_0 = FX(0)
K4_1 = FX(0)
for first in range(3):
    for second in range(3):
        pattern = [first,second]
        K4_0 += edir[first]*edir[second]*(
            normal_deriv_value(K0jet,pattern+[1,1])
            -2*normal_deriv_value(K0jet,pattern+[1,2])
            +normal_deriv_value(K0jet,pattern+[2,2])
        )
        K4_1 += edir[first]*edir[second]*(
            twisted_deriv_value(K1bar,pattern+[1,1],1)
            -2*twisted_deriv_value(K1bar,pattern+[1,2],1)
            +twisted_deriv_value(K1bar,pattern+[2,2],1)
        )

Ginv0 = mat_inv_3(G0)
Ginv1 = matscale(matmul(matmul(Ginv0,G1),Ginv0),FX(-1))
Curv0 = FX(-K4_0 + dot(Ginv0,A0))
Curv1 = FX(
    -K4_1 + dot(Ginv1,A0)
    + dot(Ginv0,A1,A0) + dot(Ginv0,A0,A1)
)
ge0,ge1 = dot(G0,edir),dot(G1,edir)
go0,go1 = dot(G0,odir),dot(G1,odir)
Den0 = FX(ge0*go0)
Den1 = FX(ge1*go0+ge0*go1)

print("  common geometry complete",flush=True)


def order2_data(K2bar):
    G2 = [[twisted_deriv_value(K2bar,[r,c],2) for c in range(3)] for r in range(3)]

    mu_e_2 = dot(G2,edir)/dot(M,edir)
    mu_o_2 = dot(G2,odir)/dot(M,odir)
    mu_r_2 = dot(G2,rdir)/dot(M,rdir)
    sigma2 = FX(mu_o_2-mu_e_2)
    D2 = FX((mu_e_2+mu_o_2)/2+2*mu_r_2)

    A2 = []
    for q in range(3):
        value = FX(0)
        for index in range(3):
            value += edir[index]*(
                twisted_deriv_value(K2bar,[index,1,q],2)
                -twisted_deriv_value(K2bar,[index,2,q],2)
            )
        A2.append(FX(value))

    K4_2 = FX(0)
    for first in range(3):
        for second in range(3):
            pattern = [first,second]
            K4_2 += edir[first]*edir[second]*(
                twisted_deriv_value(K2bar,pattern+[1,1],2)
                -2*twisted_deriv_value(K2bar,pattern+[1,2],2)
                +twisted_deriv_value(K2bar,pattern+[2,2],2)
            )

    Ginv2 = matadd(
        matmul(matmul(matmul(matmul(Ginv0,G1),Ginv0),G1),Ginv0),
        matscale(matmul(matmul(Ginv0,G2),Ginv0),FX(-1)),
    )
    Curv2 = FX(
        -K4_2 + dot(Ginv2,A0)
        + dot(Ginv1,A1,A0) + dot(Ginv1,A0,A1)
        + dot(Ginv0,A2,A0) + dot(Ginv0,A0,A2)
        + dot(Ginv0,A1,A1)
    )

    ge2 = dot(G2,edir)
    go2 = dot(G2,odir)
    Den2 = FX(ge2*go0+ge1*go1+ge0*go2)
    Bcurv2 = FX(
        Curv2/Den0 - Curv1*Den1/Den0^2
        + Curv0*(Den1^2/Den0^3-Den2/Den0^2)
    )
    Scurv2 = FX(-Bcurv2/V0)
    return {"Scurv2":Scurv2,"sigma2":sigma2,"D2":D2}


# Build degree-two one-variable jets once.
Py = one_variable_jet(P2jets,1)
Pz = one_variable_jet(P2jets,2)
TPy = one_variable_theta_jet(P2jets,1)
TPz = one_variable_theta_jet(P2jets,2)
B2y = one_variable_jet(B2jets,1)
B2z = one_variable_jet(B2jets,2)
TB2y = one_variable_theta_jet(B2jets,1)
TB2z = one_variable_theta_jet(B2jets,2)


def make_K2(F2,thetaY,thetaZ):
    return (
        F2*(1+2*FX(L)*xj)
        + FX(L)*yj*thetaY + FX(L)*zj*thetaZ
    )/(2*FX(L)^3*Vjet)


K2_zero = Jet.const(0)
K2_A = make_K2(Py*Pz,TPy*Pz,Py*TPz)
K2_B = make_K2(B2y*B2z,TB2y*B2z,B2y*TB2z)
K2_C = make_K2(Py*B2z+B2y*Pz,
               TPy*B2z+TB2y*Pz,
               Py*TB2z+B2y*TPz)

print("  computing base quadratic K1 x K1 response...",flush=True)
base2 = order2_data(K2_zero)
print("  computing F_A response...",flush=True)
dataA = order2_data(K2_A)
gc.collect()
print("  computing F_B response...",flush=True)
dataB = order2_data(K2_B)
gc.collect()
print("  computing F_C response...",flush=True)
dataC = order2_data(K2_C)
gc.collect()


def corrected_column(total):
    return FX(
        total["Scurv2"]-base2["Scurv2"]
        -FX(alpha0)*total["sigma2"]-FX(beta0)*total["D2"]
    )


Rquad = FX(base2["Scurv2"])
rA = corrected_column(dataA)
rB = corrected_column(dataB)
rC = corrected_column(dataC)

print("  exact defect components constructed",flush=True)
print()

# ----------------------------------------------------------------------------
# PART III. Laurent extraction and direct channel verification.
# ----------------------------------------------------------------------------

print("PART III. DIRECT LAURENT EXTRACTION AT U=0")
print("-"*78)

LS = LaurentSeriesRing(K,names=("u",),default_prec=LAURENT_PREC)
u = LS.gen()
X_at_u = (u-1)/3


def eval_polynomial_series(poly,arg):
    coefficients = poly.list()
    result = LS(0)
    for coefficient in reversed(coefficients):
        result = result*arg + LS(coefficient)
    return result


def to_laurent(expr,precision=LAURENT_PREC):
    expr = FX(expr)
    numerator = RX(expr.numerator())
    denominator = RX(expr.denominator())
    num_series = eval_polynomial_series(numerator,X_at_u)
    den_series = eval_polynomial_series(denominator,X_at_u)
    return (num_series/den_series).add_bigoh(precision)


def series_coefficient(series,exponent):
    try:
        return K(series[exponent])
    except (IndexError,ValueError):
        return K(0)


def is_zero_K(expr):
    value = K(expr)
    return value.numerator() == 0


def is_zero_FX(expr):
    value = FX(expr)
    return value.numerator() == 0


sigma_series = to_laurent(sigma1)
D_series = to_laurent(D1)

s = {k:series_coefficient(sigma_series,k) for k in range(-2,5)}
d = {k:series_coefficient(D_series,k) for k in range(-2,5)}

expected_s = {
    -2:K(-C0^2*(L^2*E4+108)/(972*L^3)),
    -1:K(-C0^2*(L^2*E4+84)/(1944*L^2)),
     0:K(-C0^2/(243*L)),
}
expected_d = {
    -2:K(C0^2*(L^2*E4+36)/(648*L^3)),
    -1:K(C0^2*(L^2*E4+36)/(1944*L^2)),
     0:K(0),
}

for exponent,value in expected_s.items():
    if not is_zero_K(s[exponent]-value):
        raise ArithmeticError("sigma Laurent derivation failed at U^{}".format(exponent))
for exponent,value in expected_d.items():
    if not is_zero_K(d[exponent]-value):
        raise ArithmeticError("D Laurent derivation failed at U^{}".format(exponent))

rho3 = K(-s[-2]/d[-2])
rho3_expected = K(2*(L^2*E4+108)/(3*(L^2*E4+36)))
if not is_zero_K(rho3-rho3_expected):
    raise ArithmeticError("rho3 recognition failed")

print("  sigma and D Laurent coefficients derived from the original geometry")
print("  ord_U(sigma)=ord_U(D)=-2")
print("  rho3 = {}".format(rho3))
print("  s_1,s_2,s_3,s_4 and d_1,d_2,d_3,d_4 are now exact derived values")
print()

# Defect components and their Laurent coefficients.
defect_functions = [Rquad,rA,rB,rC]
defect_names = ["intercept Rquad","slope rA","slope rB","slope rC"]
defect_series = [to_laurent(item) for item in defect_functions]
F = [{k:series_coefficient(series,k) for k in MATCH_EXPONENTS}
     for series in defect_series]

for name,coefficients in zip(defect_names,F):
    valuation = None
    for exponent in MATCH_EXPONENTS:
        if not is_zero_K(coefficients[exponent]):
            valuation = exponent
            break
    print("  {} Laurent valuation in tested range: {}".format(name,valuation))
print()

# ----------------------------------------------------------------------------
# PART IV. Exact seven-equation Laurent solve for the affine connection map.
# ----------------------------------------------------------------------------

print("PART IV. EXACT LAURENT SOLVE")
print("-"*78)

# Unknown polar/constant vector:
#   [alpha_-3,beta_-3,alpha_-2,beta_-2,alpha_-1,alpha_0,beta_0]
# The symbols alpha0,beta0 above denote the leading closure coefficients; the
# last two entries here are degree-two connection constants.
unknown_names = [
    "alpha_U_-3","beta_U_-3","alpha_U_-2","beta_U_-2",
    "alpha_U_-1","alpha_const","beta_const"
]


def coeff(mapping,index):
    return mapping.get(index,K(0))


connection_matrix_rows = []
for exponent in MATCH_EXPONENTS:
    connection_matrix_rows.append([
        coeff(s,exponent+3),
        coeff(d,exponent+3),
        coeff(s,exponent+2),
        coeff(d,exponent+2),
        coeff(s,exponent+1),
        coeff(s,exponent),
        coeff(d,exponent),
    ])
connection_matrix = matrix(K,connection_matrix_rows)

universal_a1 = K(alpha0*C0^2/(81*L))
universal_a2 = K(-alpha0*C0^2/243)


def universal_coefficient(exponent):
    return K(universal_a1*coeff(s,exponent-1)
             + universal_a2*coeff(s,exponent-2))


rhs_columns = []
for component_index in range(4):
    column = []
    for exponent in MATCH_EXPONENTS:
        value = F[component_index][exponent]
        if component_index == 0:
            value -= universal_coefficient(exponent)
        column.append(K(value))
    rhs_columns.append(column)

rhs_matrix = matrix(K,7,4,[rhs_columns[col][row]
                           for row in range(7) for col in range(4)])

if connection_matrix.rank() != 7:
    raise ArithmeticError("The seven-equation Laurent connection matrix is singular")

print("  seven-channel Laurent matrix rank = 7")
print("  solving intercept and three geometry slopes simultaneously...",flush=True)
solution_matrix = connection_matrix.solve_right(rhs_matrix)

# solution_matrix[row,column]: row=connection coefficient, column=0,A,B,C.
AFFINE_CONNECTION_MAP = {
    unknown_names[row]: [K(solution_matrix[row,col]) for col in range(4)]
    for row in range(7)
}
AFFINE_CONNECTION_MAP["alpha_U_1"] = [universal_a1,K(0),K(0),K(0)]
AFFINE_CONNECTION_MAP["alpha_U_2"] = [universal_a2,K(0),K(0),K(0)]

# Exact deepest-pole alignment, including all four affine columns.
for col in range(4):
    am3_value = solution_matrix[0,col]
    bm3_value = solution_matrix[1,col]
    if not is_zero_K(bm3_value-rho3*am3_value):
        raise ArithmeticError("beta_-3=rho3 alpha_-3 failed in affine column {}".format(col))

print("  exact affine connection map solved")
print("  beta_-3=rho3*alpha_-3 certified for intercept and all three slopes")
print("  universal positive-U drift has zero geometry slopes by construction")
print()

# ----------------------------------------------------------------------------
# PART V. CoCalc-optimized exact rational-function certificate.
# ----------------------------------------------------------------------------

print("PART V. COCALC-OPTIMIZED EXACT RATIONAL-FUNCTION CERTIFICATE")
print("-"*78)
print("  saving the complete Part-IV state before verification...", flush=True)

PART4_CHECKPOINT = "exact_degree2_part4_checkpoint.sobj"
FINAL_CHECKPOINT = "exact_degree2_fast_certificate.sobj"
TEXT_SUMMARY = "exact_degree2_fast_certificate_summary.txt"


def num_den_RX(value):
    value = FX(value)
    return RX(value.numerator()), RX(value.denominator())


def degree_or_minus_one(poly):
    poly = RX(poly)
    return -1 if poly == 0 else int(poly.degree())


def horner_eval(poly, point):
    poly = RX(poly)
    point = K(point)
    result = K(0)
    for coefficient in reversed(poly.list()):
        result = result*point + K(coefficient)
    return result


sigma_num, sigma_den = num_den_RX(sigma1)
D_num, D_den = num_den_RX(D1)
defect_num_den = [num_den_RX(item) for item in defect_functions]
solution_entries = [
    [K(solution_matrix[row,col]) for col in range(solution_matrix.ncols())]
    for row in range(solution_matrix.nrows())
]

PART4_DATA = {
    "sigma_num": sigma_num,
    "sigma_den": sigma_den,
    "D_num": D_num,
    "D_den": D_den,
    "defect_num_den": defect_num_den,
    "defect_names": list(defect_names),
    "solution_entries": solution_entries,
    "universal_a1": K(universal_a1),
    "universal_a2": K(universal_a2),
    "rho3": K(rho3),
    "connection_map": AFFINE_CONNECTION_MAP,
    "channel_laurent": {"sigma":s,"D":d},
    "defect_laurent": F,
}
save(PART4_DATA, PART4_CHECKPOINT)
print("  Part-IV checkpoint written: {}".format(PART4_CHECKPOINT), flush=True)

# Compact rho3 recognition, avoiding the enormous expanded display.
rho3_compact = K(2*(L^2*E4+108)/(3*(L^2*E4+36)))
if not is_zero_K(K(rho3)-rho3_compact):
    raise ArithmeticError("compact rho3 recognition failed")
print("  rho3 = 2*(L^2*E4+108)/(3*(L^2*E4+36))", flush=True)

# Clear U^3 structurally.  We prove the resulting polynomial N(X) is zero by
# checking D+1 distinct exact rational X-values, where D is a rigorous degree
# upper bound.  This avoids global multivariate fraction normalization.
U_poly = RX(3*X+1)


def connection_polynomials(column_index):
    am3,bm3,am2,bm2,am1,aconst,bconst = [
        K(solution_entries[row][column_index]) for row in range(7)
    ]
    Atilde = RX(am3 + am2*U_poly + am1*U_poly^2 + aconst*U_poly^3)
    Btilde = RX(bm3 + bm2*U_poly + bconst*U_poly^3)
    if column_index == 0:
        Atilde += RX(K(universal_a1)*U_poly^4 + K(universal_a2)*U_poly^5)
    return RX(Atilde), RX(Btilde)


def identity_degree_bound(index,Atilde,Btilde):
    nF,dF = defect_num_den[index]
    bounds = [
        degree_or_minus_one(nF)+degree_or_minus_one(sigma_den)
        +degree_or_minus_one(D_den)+3,
        degree_or_minus_one(dF)+degree_or_minus_one(Atilde)
        +degree_or_minus_one(sigma_num)+degree_or_minus_one(D_den),
        degree_or_minus_one(dF)+degree_or_minus_one(Btilde)
        +degree_or_minus_one(D_num)+degree_or_minus_one(sigma_den),
    ]
    return max(bounds), bounds


def cleared_value(index,Atilde,Btilde,point):
    nF,dF = defect_num_den[index]
    x0 = K(point)
    u0 = K(3*x0+1)
    nF0,dF0 = horner_eval(nF,x0),horner_eval(dF,x0)
    nS0,dS0 = horner_eval(sigma_num,x0),horner_eval(sigma_den,x0)
    nD0,dD0 = horner_eval(D_num,x0),horner_eval(D_den,x0)
    A0,B0 = horner_eval(Atilde,x0),horner_eval(Btilde,x0)
    return K(nF0*dS0*dD0*u0^3
             - dF0*(A0*nS0*dD0+B0*nD0*dS0))


def rational_test_points(count):
    points=[]
    shell=0
    while len(points)<count:
        candidates=[QQ(0)] if shell==0 else [
            QQ(shell),QQ(-shell),(QQ(2*shell-1)/QQ(2)),(QQ(-(2*shell-1))/QQ(2))
        ]
        for value in candidates:
            if value not in points:
                points.append(value)
                if len(points)>=count:
                    break
        shell+=1
    return points


FULL_IDENTITY_RESULTS=[]
IDENTITY_WITNESSES=[]
IDENTITY_DEGREE_BOUNDS=[]

for index,name in enumerate(defect_names):
    print("-"*78,flush=True)
    print("  verifying {}".format(name),flush=True)
    Atilde,Btilde=connection_polynomials(index)
    bound,term_bounds=identity_degree_bound(index,Atilde,Btilde)
    points=rational_test_points(bound+1)
    print("    term degree bounds: {}".format(term_bounds),flush=True)
    print("    cleared numerator degree <= {}".format(bound),flush=True)
    print("    exact evaluation count: {}".format(bound+1),flush=True)

    passed=True
    witness=None
    start_time=cputime()
    for j,point in enumerate(points):
        value=cleared_value(index,Atilde,Btilde,point)
        if not is_zero_K(value):
            passed=False
            witness={"point":point,"value":value}
            print("    NONZERO exact witness at X={}".format(point),flush=True)
            break
        if (j+1)%5==0 or j+1==len(points):
            print("    passed {}/{} evaluations".format(j+1,len(points)),flush=True)
            gc.collect()
    print("    exact identity certified? {}".format(passed),flush=True)
    print("    CPU seconds: {}".format(cputime(start_time)),flush=True)
    FULL_IDENTITY_RESULTS.append(passed)
    IDENTITY_WITNESSES.append(witness)
    IDENTITY_DEGREE_BOUNDS.append(bound)

CERTIFICATE = dict(PART4_DATA)
CERTIFICATE.update({
    "verification_method":"exact degree bound plus D+1 rational evaluations",
    "full_identity_results":FULL_IDENTITY_RESULTS,
    "identity_witnesses":IDENTITY_WITNESSES,
    "identity_degree_bounds":IDENTITY_DEGREE_BOUNDS,
})
save(CERTIFICATE,FINAL_CHECKPOINT)

with open(TEXT_SUMMARY,"w") as handle:
    handle.write("EXACT DEGREE-TWO COCALC CERTIFICATE\n")
    handle.write("Method: degree bound plus exact rational evaluations\n")
    for name,result,bound,witness in zip(
        defect_names,FULL_IDENTITY_RESULTS,IDENTITY_DEGREE_BOUNDS,IDENTITY_WITNESSES
    ):
        handle.write("{}: {} (degree <= {})\n".format(name,result,bound))
        if witness is not None:
            handle.write("  witness X={}\n".format(witness["point"]))
    handle.write("rho3 = 2*(L^2*E4+108)/(3*(L^2*E4+36))\n")

print()
print("="*78)
print("EXACT DEGREE-TWO COCALC CERTIFICATE SUMMARY")
print("  Part-IV checkpoint: {}".format(PART4_CHECKPOINT))
for name,result,bound in zip(defect_names,FULL_IDENTITY_RESULTS,IDENTITY_DEGREE_BOUNDS):
    print("  {:18s}: {}  (degree <= {})".format(name,result,bound))
print("  all four identities exact? {}".format(all(FULL_IDENTITY_RESULTS)))
print("  final checkpoint: {}".format(FINAL_CHECKPOINT))
print("  text summary: {}".format(TEXT_SUMMARY))
print("="*78)
if all(FULL_IDENTITY_RESULTS):
    print("SUCCESS")
    print("The affine degree-two spectral connection is an exact rational identity")
    print("over the formal square-lattice/level-two CM field.")
else:
    print("PARTIAL SUCCESS")
    print("The first nonzero exact evaluation is saved as a diagnostic witness.")
print("="*78)
