from sage.all import *
from math import factorial, comb

# ============================================================
# GENERIC-v0 CONDITIONAL FORMAL PROOF
#
# Goal:
#   Prove, as a rational identity in X=v0, that
#
#       Scurv_bar(X) = alpha*sigma_bar(X) + beta*Dchan_bar(X)
#
# after imposing the two HSS/CM finite-jet identities:
#
#   R1 = -L*C0 + 6*L*C1 + 6*C0 = 0
#
# and
#
#   R5 =
#   -54L^2C2^2 + 12L^2C0C3 + 108L^2C2C3
#   -45L^2C0C4 + 54L^2C0C5
#   +LC0^2 -72LC0C2 +216LC0C3 -270LC0C4 -8C0^2 = 0.
#
# We solve R1 for C1 and R5 for C5 before constructing the channels.
#
# If the final residual numerator is zero, then the closure is proven
# formally for arbitrary v0 on the HSS CM line.
# ============================================================

ORDER = 4

X1 = QQ(6)/5
X2 = QQ(2)

DO_FACTOR_IF_NONZERO = False

print("="*70)
print("GENERIC-v0 CONDITIONAL FORMAL PROOF")
print("="*70)

# ------------------------------------------------------------
# Ring setup
# ------------------------------------------------------------

# Variables:
#   X  = symbolic v0
#   L  = 2*pi
#   C0,C2,C3,C4 are independent after imposing R1,R5.
names = ["X", "L", "C0", "C2", "C3", "C4"]
R = PolynomialRing(QQ, names)
K = R.fraction_field()

X, L, C0, C2, C3, C4 = [K(g) for g in R.gens()]

# R1 gives:
#   C1 = (L-6)C0/(6L)
C1 = ((L - 6) * C0) / (6 * L)

# R5 gives:
#   C5 = ...
C5 = (
    54*L**2*C2**2
    - 12*L**2*C0*C3
    - 108*L**2*C2*C3
    + 45*L**2*C0*C4
    - L*C0**2
    + 72*L*C0*C2
    - 216*L*C0*C3
    + 270*L*C0*C4
    + 8*C0**2
) / (54*L**2*C0)

C = [C0, C1, C2, C3, C4, C5]

print("Reduced variables:", names)
print("C1 imposed from R1.")
print("C5 imposed from R5.")
print()

# ============================================================
# Tiny formal Taylor-jet algebra over K
# ============================================================

def total_degree(k):
    return int(k[0]) + int(k[1]) + int(k[2])


class Jet:
    def __init__(self, coeffs=None):
        self.c = {}
        if coeffs:
            for k, v in coeffs.items():
                kk = tuple(int(a) for a in k)
                vv = K(v)
                if total_degree(kk) <= ORDER and vv != 0:
                    self.c[kk] = vv

    @staticmethod
    def const(a):
        return Jet({(0, 0, 0): K(a)})

    @staticmethod
    def var(idx, base):
        k = [0, 0, 0]
        k[idx] = 1
        return Jet({
            (0, 0, 0): K(base),
            tuple(k): K(1),
        })

    def scale(self, a):
        a = K(a)
        return Jet({k: a*v for k, v in self.c.items()})

    def __add__(self, other):
        other = tojet(other)
        out = dict(self.c)
        for k, v in other.c.items():
            out[k] = out.get(k, K(0)) + v
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
        """
        Formal log. The constant log(a0) is irrelevant for all derivatives
        of positive order, so we omit it.
        """
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

        return coeff * mult


def tojet(x):
    if isinstance(x, Jet):
        return x
    return Jet.const(x)

# ============================================================
# Matrix helpers
# ============================================================

def mat_inv_3(A):
    a,b,c = A[0]
    d,e,f = A[1]
    g,h,i = A[2]

    det = (
        a*(e*i - f*h)
        - b*(d*i - f*g)
        + c*(d*h - e*g)
    )

    return [
        [(e*i - f*h)/det, (c*h - b*i)/det, (b*f - c*e)/det],
        [(f*g - d*i)/det, (a*i - c*g)/det, (c*d - a*f)/det],
        [(d*h - e*g)/det, (b*g - a*h)/det, (a*e - b*d)/det],
    ]


def matmul(A, B):
    return [
        [sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def mat_neg(A):
    return [[-A[i][j] for j in range(3)] for i in range(3)]


def dot(A, u, v=None):
    if v is None:
        v = u
    return sum(u[i]*A[i][j]*v[j] for i in range(3) for j in range(3))

# ============================================================
# Derivative helpers
# ============================================================

def normal_deriv_value(jet, inds):
    return jet.deriv_value((inds.count(0), inds.count(1), inds.count(2)))


def pfree_deriv_value(jet, inds):
    """
    D_x = partial_x - L,
    D_y = partial_y,
    D_z = partial_z.
    """
    a = inds.count(0)
    b = inds.count(1)
    c = inds.count(2)

    total = K(0)

    for j in range(a + 1):
        total += K(comb(a, j)) * ((-L)**(a-j)) * jet.deriv_value((j, b, c))

    return total


def Bjet_y():
    out = Jet.const(0)
    dy = Jet({(0, 1, 0): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m] / K(factorial(m))) * (dy**m)

    return out


def Bjet_z():
    out = Jet.const(0)
    dz = Jet({(0, 0, 1): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m] / K(factorial(m))) * (dz**m)

    return out


def Thetajet_y():
    out = Jet.const(0)
    dy = Jet({(0, 1, 0): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m+1] / K(factorial(m))) * (dy**m)

    return out


def Thetajet_z():
    out = Jet.const(0)
    dz = Jet({(0, 0, 1): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m+1] / K(factorial(m))) * (dz**m)

    return out


def primitive_numerator(expr):
    expr = K(expr)
    num = R(expr.numerator())

    try:
        return num.primitive_part()
    except Exception:
        cont = num.content()
        if cont != 0:
            return R(num/cont)
        return num

# ============================================================
# Channel construction at arbitrary x0
# ============================================================

def channels_at_x(x0):
    """
    x0 may be a rational number or the symbolic variable X.
    Returns:
      F = Scurv_bar
      S = sigma_bar
      D = Dchan_bar
      B0
      den_check
    """
    x0 = K(x0)

    x = Jet.var(0, x0)
    y = Jet.var(1, K(1))
    z = Jet.var(2, K(1))

    # HSS restricted cubic.
    V = 9*x*y*z + K(QQ(3)/2)*y*y*z + K(QQ(3)/2)*y*z*z

    K0 = -V.log()

    By = Bjet_y()
    Bz = Bjet_z()
    Ty = Thetajet_y()
    Tz = Thetajet_z()

    # p-free leading correction K1bar.
    K1bar = (
        By*Bz*(1 + L*x)
        + L*y*Ty*Bz
        + L*z*By*Tz
    ) / (2 * L**3 * V)

    V0 = V.deriv_value((0, 0, 0))

    G0 = [
        [normal_deriv_value(K0, [i, j]) for j in range(3)]
        for i in range(3)
    ]

    G1bar = [
        [pfree_deriv_value(K1bar, [i, j]) for j in range(3)]
        for i in range(3)
    ]

    M = [
        [normal_deriv_value(V, [i, j]) for j in range(3)]
        for i in range(3)
    ]

    # Even, odd, radial directions at v=(x0,1,1).
    e = [-(2*x0 + 1), K(1), K(1)]
    o = [K(0), K(1), K(-1)]
    r = [x0, K(1), K(1)]

    mu_e = dot(G1bar, e) / dot(M, e)
    mu_o = dot(G1bar, o) / dot(M, o)
    mu_r = dot(G1bar, r) / dot(M, r)

    sigma = mu_o - mu_e
    Dchan = (mu_e + mu_o)/2 + 2*mu_r

    # Curvature pieces.
    A0 = []
    A1 = []

    for q in range(3):
        A0q = K(0)
        A1q = K(0)

        for i in range(3):
            A0q += e[i] * (
                normal_deriv_value(K0, [i, 1, q])
                - normal_deriv_value(K0, [i, 2, q])
            )

            A1q += e[i] * (
                pfree_deriv_value(K1bar, [i, 1, q])
                - pfree_deriv_value(K1bar, [i, 2, q])
            )

        A0.append(A0q)
        A1.append(A1q)

    K4_0 = K(0)
    K4_1 = K(0)

    for i in range(3):
        for j in range(3):
            K4_0 += e[i]*e[j] * (
                normal_deriv_value(K0, [i, j, 1, 1])
                - 2*normal_deriv_value(K0, [i, j, 1, 2])
                + normal_deriv_value(K0, [i, j, 2, 2])
            )

            K4_1 += e[i]*e[j] * (
                pfree_deriv_value(K1bar, [i, j, 1, 1])
                - 2*pfree_deriv_value(K1bar, [i, j, 1, 2])
                + pfree_deriv_value(K1bar, [i, j, 2, 2])
            )

    Ginv0 = mat_inv_3(G0)
    Ginv1 = mat_neg(matmul(matmul(Ginv0, G1bar), Ginv0))

    R0 = -K4_0 + dot(Ginv0, A0)

    R1 = (
        -K4_1
        + dot(Ginv1, A0)
        + dot(Ginv0, A1, A0)
        + dot(Ginv0, A0, A1)
    )

    Den0 = dot(G0, e) * dot(G0, o)
    Den1 = dot(G1bar, e) * dot(G0, o) + dot(G0, e) * dot(G1bar, o)

    B0 = R0 / Den0

    B1_num = R1 / Den0
    B1_den = -R0 * Den1 / Den0**2
    B1 = B1_num + B1_den

    Scurv = -B1 / V0

    den_check = (-B1_den / V0) - (mu_e + mu_o)/3

    return {
        "F": Scurv,
        "S": sigma,
        "D": Dchan,
        "B0": B0,
        "den_check": den_check,
    }

# ============================================================
# Main proof calculation
# ============================================================

print("Building channels at x1 =", X1)
ch1 = channels_at_x(K(X1))

print("Building channels at x2 =", X2)
ch2 = channels_at_x(K(X2))

print("Building channels at generic X")
chX = channels_at_x(X)

print()
print("Sanity checks:")
B0_num = primitive_numerator(chX["B0"] + K(QQ(1)/3))
den_num = primitive_numerator(chX["den_check"])

print("  generic B0 + 1/3 numerator zero?", bool(B0_num == 0))
print("  generic denominator identity numerator zero?", bool(den_num == 0))

# Fit alpha,beta from x1,x2.
F1 = ch1["F"]
F2 = ch2["F"]

S1 = ch1["S"]
S2 = ch2["S"]

D1 = ch1["D"]
D2 = ch2["D"]

det = S1*D2 - S2*D1

alpha = (F1*D2 - F2*D1) / det
beta  = (S1*F2 - S2*F1) / det

print()
print("alpha and beta built as rational functions in reduced jet variables.")
print("alpha numerator terms:", len(R(alpha.numerator()).dict()))
print("alpha denominator terms:", len(R(alpha.denominator()).dict()))
print("beta numerator terms:", len(R(beta.numerator()).dict()))
print("beta denominator terms:", len(R(beta.denominator()).dict()))

# Generic residual.
print()
print("Building generic residual...")
resid = chX["F"] - alpha*chX["S"] - beta*chX["D"]

num = primitive_numerator(resid)

print()
print("="*70)
print("GENERIC RESIDUAL RESULT")
print("="*70)
print("residual numerator zero?", bool(num == 0))

if num == 0:
    print()
    print("SUCCESS:")
    print("After imposing R1 and R5, the closure holds as a rational")
    print("identity in the generic variable X=v0.")
else:
    print()
    print("NONZERO:")
    print("The residual did not vanish identically.")
    print("degree:", num.degree())
    print("term count:", len(num.dict()))
    for nm, gen in zip(names, R.gens()):
        print("degree in {}: {}".format(nm, num.degree(gen)))

    if DO_FACTOR_IF_NONZERO:
        print()
        print("Trying to factor the residual numerator...")
        try:
            print(num.factor())
        except Exception as err:
            print("factor failed:", err)

print()
print("Done.")
