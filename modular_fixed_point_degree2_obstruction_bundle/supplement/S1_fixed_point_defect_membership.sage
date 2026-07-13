from sage.all import *
from math import factorial, comb

# ============================================================================
# FIXED-POINT DEFECT MEMBERSHIP TEST
# ============================================================================
#
# Purpose
# -------
# Keep the modular data formal and test whether the cleared two-channel
# curvature-closure defect belongs to the square-lattice fixed-point ideal
#
#       I_i = < delta, E6 >,        delta := L*E2 - 6.
#
# The script uses the natural off-fixed-point extensions
#
#   alpha(E4) = (L^2 E4 + 12)/72,
#
#   beta(E4)  = -(7 L^4 E4^2 - 552 L^2 E4 - 13392)
#               / (288 (L^2 E4 + 36)),
#
# which reduce to the exact square-lattice coefficients at tau=i.
#
# It constructs the normalized eta-product jets from Ramanujan's equations,
# builds the same metric/spectral/curvature channels as the previous scripts,
# clears the defect denominator, and proves an exact polynomial certificate
#
#       N = delta*A + E6*B.
#
# IMPORTANT
# ---------
# This is a formal algebraic fixed-point-ideal test.  It does not yet prove
# that the defect is modular-covariant away from the imaginary HSS slice.
# ============================================================================

ORDER = 4
PRINT_FULL_CERTIFICATE = False   # Set True only if you really want huge output.
PRINT_LINEAR_FACTORS = False

print("="*78)
print("FIXED-POINT DEFECT MEMBERSHIP TEST")
print("Testing I_i = <delta, E6>, where delta = L*E2 - 6")
print("="*78)

# ----------------------------------------------------------------------------
# Exact coefficient field
# ----------------------------------------------------------------------------
# X     = v0 on the symmetric slice (X,1,1)
# L     = 2*pi, kept formal
# E4    = formal weight-four Eisenstein value
# delta = L*E2 - 6, the first fixed-point generator
# E6    = formal weight-six Eisenstein value
# C0    = B itself; all higher C_m are C0*u_m
#
# We use a fraction field because the channels are rational functions.

names = ["X", "L", "E4", "delta", "E6", "C0"]
R = PolynomialRing(QQ, names)
K = R.fraction_field()

rX, rL, rE4, rdelta, rE6, rC0 = R.gens()
X, L, E4, delta, E6, C0 = [K(g) for g in R.gens()]

E2 = (delta + 6)/L

print("Ring variables:", names)
print("E2 represented as (delta+6)/L")
print()

# ----------------------------------------------------------------------------
# Rational differentiation and Ramanujan derivation
# ----------------------------------------------------------------------------

def dK(expr, rvar):
    """Differentiate a fraction-field expression with respect to rvar."""
    expr = K(expr)
    num = R(expr.numerator())
    den = R(expr.denominator())
    return K(num.derivative(rvar)*den - num*den.derivative(rvar)) / K(den)**2


# Since delta=L*E2-6 and theta(L)=0,
#
#   theta(delta) = L*(E2^2-E4)/12.
#
# Ramanujan:
#   theta(E4) = (E2*E4-E6)/3,
#   theta(E6) = (E2*E6-E4^2)/2.

theta_delta = L*(E2**2 - E4)/12
theta_E4 = (E2*E4 - E6)/3
theta_E6 = (E2*E6 - E4**2)/2


def theta(expr):
    """Ramanujan theta=q*d/dq in variables delta,E4,E6."""
    return (
        dK(expr, rdelta) * theta_delta
        + dK(expr, rE4) * theta_E4
        + dK(expr, rE6) * theta_E6
    )


# For B=9*q^(1/6)*eta^(-4), theta(log B)=(1-E2)/6.
Alog = (1 - E2)/6

# u_m = theta^m B / B, with u_{m+1}=theta(u_m)+Alog*u_m.
u = [K(1)]
for m in range(5):
    u.append(K(theta(u[-1]) + Alog*u[-1]))

C = [K(C0*u[m]) for m in range(6)]

print("Built normalized Ramanujan jets u_0,...,u_5.")
print("Quick jet-ideal checks:")

R1_norm = K(-L + 6*L*u[1] + 6)
R3_norm = K((L - 9) - 54*L*u[2] + 108*L*u[3])
R5_norm = K(
    -54*L**2*u[2]**2
    + 12*L**2*u[3]
    + 108*L**2*u[2]*u[3]
    - 45*L**2*u[4]
    + 54*L**2*u[5]
    + L
    - 72*L*u[2]
    + 216*L*u[3]
    - 270*L*u[4]
    - 8
)


def specialize_fixed(expr):
    """Set delta=0 and E6=0 in a fraction-field expression."""
    expr = K(expr)
    num = R(expr.numerator()).subs({rdelta: 0, rE6: 0})
    den = R(expr.denominator()).subs({rdelta: 0, rE6: 0})
    return K(num)/K(den)


# Explicit small fixed-point-ideal certificates for the jet relations.
R5_A = (
    15*L**2*E4**2 + 3*L**2*E4 + 5*L**2 + 120*L + 108
    - 12*L*delta + 36*delta + 3*delta**2
)/(144*L)
R5_B = -(
    6*L**2*E4 - 2*L*delta + 3*delta**2 + 6*delta - 72
)/72

print("  R1 = -delta ?", bool(R1_norm + delta == 0))
print("  R3 = (3*delta-L*E6)/2 ?",
      bool(R3_norm - (3*delta-L*E6)/2 == 0))
print("  R5 = delta*R5_A + E6*R5_B ?",
      bool(R5_norm - delta*R5_A - E6*R5_B == 0))
print("  R1|fixed = 0 ?", bool(specialize_fixed(R1_norm) == 0))
print("  R3|fixed = 0 ?", bool(specialize_fixed(R3_norm) == 0))
print("  R5|fixed = 0 ?", bool(specialize_fixed(R5_norm) == 0))
print()

# ============================================================================
# Truncated multivariate Taylor-jet algebra over K
# ============================================================================

def total_degree(k):
    return int(k[0]) + int(k[1]) + int(k[2])


class Jet:
    def __init__(self, coeffs=None):
        self.c = {}
        if coeffs:
            for k, value in coeffs.items():
                kk = tuple(int(a) for a in k)
                vv = K(value)
                if total_degree(kk) <= ORDER and vv != 0:
                    self.c[kk] = vv

    @staticmethod
    def const(a):
        return Jet({(0, 0, 0): K(a)})

    @staticmethod
    def var(idx, base):
        k = [0, 0, 0]
        k[idx] = 1
        return Jet({(0, 0, 0): K(base), tuple(k): K(1)})

    def scale(self, a):
        a = K(a)
        return Jet({k: a*v for k, v in self.c.items()})

    def __add__(self, other):
        other = tojet(other)
        out = dict(self.c)
        for k, value in other.c.items():
            out[k] = out.get(k, K(0)) + value
        return Jet(out)

    __radd__ = __add__

    def __neg__(self):
        return Jet({k: -v for k, v in self.c.items()})

    def __sub__(self, other):
        return self + (-tojet(other))

    def __rsub__(self, other):
        return tojet(other) + (-self)

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return self.scale(other)
        out = {}
        for a, va in self.c.items():
            for b, vb in other.c.items():
                k = (a[0]+b[0], a[1]+b[1], a[2]+b[2])
                if total_degree(k) <= ORDER:
                    out[k] = out.get(k, K(0)) + va*vb
        return Jet(out)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, Jet):
            return self.scale(K(1)/K(other))
        return self * other.inv()

    def __rtruediv__(self, other):
        return tojet(other) * self.inv()

    def __pow__(self, n):
        n = int(n)
        if n < 0:
            return self.inv().__pow__(-n)
        if n == 0:
            return Jet.const(1)
        out = Jet.const(1)
        base = self
        m = n
        while m > 0:
            if m % 2 == 1:
                out = out * base
            base = base * base
            m //= 2
        return out

    def inv(self):
        a0 = self.c.get((0, 0, 0), K(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.inv()")
        r = (self - Jet.const(a0)).scale(K(1)/a0)
        out = Jet.const(0)
        term = Jet.const(1)
        sign = K(1)
        for _ in range(ORDER + 1):
            out = out + term.scale(sign)
            term = term * r
            sign = -sign
        return out.scale(K(1)/a0)

    def log(self):
        """Formal log; the irrelevant constant log(a0) is omitted."""
        a0 = self.c.get((0, 0, 0), K(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.log()")
        r = (self - Jet.const(a0)).scale(K(1)/a0)
        out = Jet.const(0)
        term = Jet.const(1)
        for k in range(1, ORDER + 1):
            term = term * r
            out = out + term.scale(K((-1)**(k+1))/K(k))
        return out

    def deriv_value(self, counts):
        counts = tuple(int(a) for a in counts)
        coeff = self.c.get(counts, K(0))
        mult = K(1)
        for a in counts:
            mult *= factorial(a)
        return coeff*mult


def tojet(value):
    if isinstance(value, Jet):
        return value
    return Jet.const(value)

# ----------------------------------------------------------------------------
# Small matrix helpers
# ----------------------------------------------------------------------------

def mat_inv_3(M):
    a,b,c = M[0]
    d,e,f = M[1]
    g,h,i = M[2]
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    return [
        [(e*i-f*h)/det, (c*h-b*i)/det, (b*f-c*e)/det],
        [(f*g-d*i)/det, (a*i-c*g)/det, (c*d-a*f)/det],
        [(d*h-e*g)/det, (b*g-a*h)/det, (a*e-b*d)/det],
    ]


def matmul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]


def mat_neg(A):
    return [[-A[i][j] for j in range(3)] for i in range(3)]


def dot(A, v, w=None):
    if w is None:
        w = v
    return sum(v[i]*A[i][j]*w[j] for i in range(3) for j in range(3))

# ----------------------------------------------------------------------------
# Derivative helpers and one-variable B jets
# ----------------------------------------------------------------------------

def normal_deriv_value(jet, inds):
    return jet.deriv_value((inds.count(0), inds.count(1), inds.count(2)))


def pfree_deriv_value(jet, inds):
    """D_x=partial_x-L; D_y=partial_y; D_z=partial_z."""
    a = inds.count(0)
    b = inds.count(1)
    c = inds.count(2)
    total = K(0)
    for j in range(a+1):
        total += K(comb(a, j))*(-L)**(a-j)*jet.deriv_value((j,b,c))
    return total


def Bjet_y():
    out = Jet.const(0)
    dy = Jet({(0,1,0): K(1)})
    for m in range(ORDER+1):
        out += ((-L)**m*C[m]/K(factorial(m)))*(dy**m)
    return out


def Bjet_z():
    out = Jet.const(0)
    dz = Jet({(0,0,1): K(1)})
    for m in range(ORDER+1):
        out += ((-L)**m*C[m]/K(factorial(m)))*(dz**m)
    return out


def Thetajet_y():
    out = Jet.const(0)
    dy = Jet({(0,1,0): K(1)})
    for m in range(ORDER+1):
        out += ((-L)**m*C[m+1]/K(factorial(m)))*(dy**m)
    return out


def Thetajet_z():
    out = Jet.const(0)
    dz = Jet({(0,0,1): K(1)})
    for m in range(ORDER+1):
        out += ((-L)**m*C[m+1]/K(factorial(m)))*(dz**m)
    return out


def primitive_numerator(expr):
    expr = K(expr)
    num = R(expr.numerator())
    try:
        return num.primitive_part()
    except Exception:
        content = num.content()
        return R(num/content) if content != 0 else num

# ============================================================================
# Build the three p-free channels at generic X
# ============================================================================

def channels_at_generic_X():
    x = Jet.var(0, X)
    y = Jet.var(1, K(1))
    z = Jet.var(2, K(1))

    # Restricted HSS cubic volume.
    V = 9*x*y*z + K(QQ(3)/2)*y*y*z + K(QQ(3)/2)*y*z*z
    K0 = -V.log()

    By = Bjet_y()
    Bz = Bjet_z()
    Ty = Thetajet_y()
    Tz = Thetajet_z()

    # p-free leading correction.
    K1bar = (
        By*Bz*(1 + L*x)
        + L*y*Ty*Bz
        + L*z*By*Tz
    )/(2*L**3*V)

    V0 = V.deriv_value((0,0,0))

    G0 = [[normal_deriv_value(K0, [i,j]) for j in range(3)]
          for i in range(3)]
    G1 = [[pfree_deriv_value(K1bar, [i,j]) for j in range(3)]
          for i in range(3)]
    M = [[normal_deriv_value(V, [i,j]) for j in range(3)]
         for i in range(3)]

    # Exchange-even, exchange-odd, and radial directions.
    e = [-(2*X+1), K(1), K(1)]
    o = [K(0), K(1), K(-1)]
    r = [X, K(1), K(1)]

    mu_e = dot(G1,e)/dot(M,e)
    mu_o = dot(G1,o)/dot(M,o)
    mu_r = dot(G1,r)/dot(M,r)

    sigma = mu_o - mu_e
    Dchan = (mu_e+mu_o)/2 + 2*mu_r

    A0 = []
    A1 = []
    for q in range(3):
        a0q = K(0)
        a1q = K(0)
        for i in range(3):
            a0q += e[i]*(
                normal_deriv_value(K0,[i,1,q])
                - normal_deriv_value(K0,[i,2,q])
            )
            a1q += e[i]*(
                pfree_deriv_value(K1bar,[i,1,q])
                - pfree_deriv_value(K1bar,[i,2,q])
            )
        A0.append(a0q)
        A1.append(a1q)

    K4_0 = K(0)
    K4_1 = K(0)
    for i in range(3):
        for j in range(3):
            K4_0 += e[i]*e[j]*(
                normal_deriv_value(K0,[i,j,1,1])
                - 2*normal_deriv_value(K0,[i,j,1,2])
                + normal_deriv_value(K0,[i,j,2,2])
            )
            K4_1 += e[i]*e[j]*(
                pfree_deriv_value(K1bar,[i,j,1,1])
                - 2*pfree_deriv_value(K1bar,[i,j,1,2])
                + pfree_deriv_value(K1bar,[i,j,2,2])
            )

    Ginv0 = mat_inv_3(G0)
    Ginv1 = mat_neg(matmul(matmul(Ginv0,G1),Ginv0))

    Curv0 = -K4_0 + dot(Ginv0,A0)
    Curv1 = (
        -K4_1
        + dot(Ginv1,A0)
        + dot(Ginv0,A1,A0)
        + dot(Ginv0,A0,A1)
    )

    Den0 = dot(G0,e)*dot(G0,o)
    Den1 = dot(G1,e)*dot(G0,o) + dot(G0,e)*dot(G1,o)

    B0 = Curv0/Den0
    B1_num = Curv1/Den0
    B1_den = -Curv0*Den1/Den0**2
    B1 = B1_num + B1_den

    Scurv = -B1/V0
    denominator_identity = (-B1_den/V0) - (mu_e+mu_o)/3

    return {
        "Scurv": K(Scurv),
        "sigma": K(sigma),
        "D": K(Dchan),
        "B0": K(B0),
        "denominator_identity": K(denominator_identity),
    }


print("Building generic-X metric, spectral, and curvature channels...")
ch = channels_at_generic_X()
print("Channels built.")
print()

# Sanity checks inherited from the earlier formal calculation.
B0_test = primitive_numerator(ch["B0"] + K(QQ(1)/3))
den_test = primitive_numerator(ch["denominator_identity"])

print("Sanity checks:")
print("  B0 = -1/3 ?", bool(B0_test == 0))
print("  denominator-response identity ?", bool(den_test == 0))
print()

# Exact square-lattice coefficient functions, naturally extended in E4.
alpha = (L**2*E4 + 12)/72
beta = -(7*L**4*E4**2 - 552*L**2*E4 - 13392) / (
    288*(L**2*E4 + 36)
)

print("Constructing the generic modular closure defect...")
defect = K(ch["Scurv"] - alpha*ch["sigma"] - beta*ch["D"])

N = primitive_numerator(defect)
Dden = R(defect.denominator())

print("Cleared numerator built.")
print("  numerator term count:", len(N.dict()))
print("  denominator term count:", len(Dden.dict()))
print()

# Ensure the rational expression is defined generically on the fixed locus.
D_fixed = R(Dden.subs({rdelta:0, rE6:0}))
print("Denominator remains nonzero as a polynomial on the fixed locus?",
      bool(D_fixed != 0))

# ----------------------------------------------------------------------------
# Exact ideal-membership certificate N = delta*A + E6*B
# ----------------------------------------------------------------------------

N_delta0 = R(N.subs({rdelta:0}))
N_fixed = R(N_delta0.subs({rE6:0}))

print("Cleared defect numerator vanishes at delta=E6=0?", bool(N_fixed == 0))

A_certificate, A_remainder = (N - N_delta0).quo_rem(rdelta)
B_certificate, B_remainder = N_delta0.quo_rem(rE6)

certificate_ok = bool(
    A_remainder == 0
    and B_remainder == 0
    and N == rdelta*A_certificate + rE6*B_certificate
)

print("Exact certificate N = delta*A + E6*B verified?", certificate_ok)
print()

if not certificate_ok:
    print("FAILURE DETAILS")
    print("  remainder after division by delta:", A_remainder)
    print("  remainder after division by E6:", B_remainder)
    raise ArithmeticError("Fixed-point ideal certificate failed.")

print("Certificate sizes:")
print("  A term count:", len(A_certificate.dict()))
print("  B term count:", len(B_certificate.dict()))
print()

# First-order transverse defect coefficients.
A_linear = R(A_certificate.subs({rdelta:0, rE6:0}))
B_linear = R(B_certificate.subs({rdelta:0, rE6:0}))

print("First-order transverse data:")
print("  coefficient of delta is identically zero?", bool(A_linear == 0))
print("  coefficient of E6 is identically zero?", bool(B_linear == 0))
print("  delta-linear term count:", len(A_linear.dict()))
print("  E6-linear term count:", len(B_linear.dict()))
print()

if PRINT_LINEAR_FACTORS:
    print("Factored first-order delta coefficient:")
    try:
        print(A_linear.factor())
    except Exception as err:
        print("  factorization failed:", err)
        print(A_linear)
    print()

    print("Factored first-order E6 coefficient:")
    try:
        print(B_linear.factor())
    except Exception as err:
        print("  factorization failed:", err)
        print(B_linear)
    print()

if PRINT_FULL_CERTIFICATE:
    print("A certificate =")
    print(A_certificate)
    print()
    print("B certificate =")
    print(B_certificate)
    print()

print("="*78)
print("SUCCESS")
print("The cleared closure defect lies in <L*E2-6, E6>.")
print("Please copy the complete printed output back into the chat.")
print("="*78)
