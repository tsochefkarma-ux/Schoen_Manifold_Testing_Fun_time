from sage.all import *
from math import factorial, comb

# ============================================================================
# DEGREE-TWO SCHOEN FACTORIZED-GLUING TEST IN THE FULL p^2 GEOMETRY
# ============================================================================
#
# PURPOSE
# -------
# The exact Schoen degeneration formula at base degree k=2 contains the
# contact partitions (2) and (1,1).  The primary literature proves the
# E8 x E8 quasi-bi-Jacobi structure, but it does not print a closed HSS-
# specialized formula for those two contact contributions.
#
# This cell therefore tests the strongest natural factorized degree-two span
# that can be built from the exact rational-elliptic-surface data already
# reconstructed:
#
#   B(q)       = 9 prod_{n>=1}(1-q^n)^(-4),
#   Z2(q)      = (E2(q)/72) B(q)^2,
#   P2(q)      = Z2(q) - (1/8) B(q^2),
#   B2(q)      = B(q^2).
#
# On the symmetric two-surface slice we use the three-dimensional span
#
#   F_A(y,z) = P2(y) P2(z),
#   F_B(y,z) = B2(y) B2(z),
#   F_C(y,z) = P2(y) B2(z) + B2(y) P2(z).
#
# Every product made from Z2 and B2 lies in this span.  In particular:
#
#   Z2(y)Z2(z) = F_A + (1/8)F_C + (1/64)F_B,
#
# and the simplest primitive-plus-multiple-cover candidate is
#
#   F_A + (1/8)F_B.
#
# IMPORTANT IMPROVEMENT OVER THE LINEARIZED JET TEST
# --------------------------------------------------
# This cell computes the actual p^2 coefficient of the metric and curvature.
# It includes the quadratic p^1 x p^1 terms from:
#
#   - the inverse metric,
#   - the curvature numerator,
#   - and the bisectional-curvature denominator.
#
# The p^2 potential is modelled by the same HSS differential operator as the
# p term, with the correct p^2 twist D_x = partial_x - 2L and the factor
# (1 + 2Lx).  This is the natural coefficientwise extension of the leading
# HSS Hessian model.
#
# TESTS
# -----
# 1. Determine whether ANY combination a*F_A+b*F_B+c*F_C can restore the
#    two-channel closure at order p^2 after allowing alpha and beta to acquire
#    first p-corrections alpha_1,beta_1.
#
# 2. Test the specific candidates:
#       - primitive product F_A,
#       - primitive + GV double cover F_A + F_B/8,
#       - uninserted product Z2(y)Z2(z).
#
# The calculation is high-precision numerical, with independent validation
# points.  It is a feasibility/obstruction test for the natural factorized
# gluing span, not a substitute for the missing exact relative-contact formula.
# ============================================================================

BITS = 320
SERIES_PREC = 70
ORDER = 4

TRAIN_X = [
    QQ(1),
    QQ(5)/4,
    QQ(3)/2,
    QQ(2),
    QQ(5)/2,
]

VALIDATE_X = [
    QQ(3)/4,
    QQ(9)/10,
    QQ(11)/10,
    QQ(4)/3,
    QQ(7)/4,
    QQ(9)/4,
    QQ(3),
    QQ(4),
]

CANDIDATE_INTERPOLATION_X = [QQ(5)/4, QQ(3)/2]

# Absolute numerical threshold.  With 320-bit arithmetic, successful exact-
# structure matches should normally be vastly below this.
ZERO_TOL = RealField(BITS)(10)**(-70)

RF = RealField(BITS)

print("="*78, flush=True)
print("FULL p^2 FACTORIZED SCHOEN-GLUING TEST", flush=True)
print("Precision: {} bits; q-series through q^{}".format(BITS, SERIES_PREC), flush=True)
print("="*78, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Physical square-lattice constants
# ----------------------------------------------------------------------------

pi = RF.pi()
L = 2*pi
q0 = exp(-L)

# Sage versions differ on whether RealField exposes a .gamma() method.
# Evaluate Gamma(1/4) through Sage's global symbolic gamma function, then
# coerce the result into the chosen high-precision real field.
gamma_quarter = RF(gamma(QQ(1)/4).n(prec=BITS))
E4i = 3*gamma_quarter**8/(2*pi)**6
Y = L**2*E4i

alpha0 = (Y + 12)/72
beta0 = -(7*Y**2 - 552*Y - 13392)/(288*(Y + 36))

print("Square-lattice constants:", flush=True)
print("  q0       =", q0, flush=True)
print("  L        =", L, flush=True)
print("  E4(i)    =", E4i, flush=True)
print("  alpha_0  =", alpha0, flush=True)
print("  beta_0   =", beta0, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Exact q-series for B, E2, Z2, B(q^2), and P2
# ----------------------------------------------------------------------------

PSQ = PowerSeriesRing(QQ, "q13", default_prec=SERIES_PREC + 1)
q13 = PSQ.gen()

B_series = PSQ(9)
for n in range(1, SERIES_PREC + 1):
    B_series *= (1 - q13**n)**(-4)
B_series = B_series.add_bigoh(SERIES_PREC + 1)

E2_series = PSQ(1)
for n in range(1, SERIES_PREC + 1):
    E2_series += -24*sigma(n, 1)*q13**n
E2_series = E2_series.add_bigoh(SERIES_PREC + 1)

Z2_series = (E2_series*B_series**2/72).add_bigoh(SERIES_PREC + 1)
B2_series = B_series(q13**2).add_bigoh(SERIES_PREC + 1)
P2_series = (Z2_series - B2_series/8).add_bigoh(SERIES_PREC + 1)


def theta_jets_from_series(series, max_order):
    """Return [theta^m series(q0)] for m=0,...,max_order."""
    coeffs = [QQ(series[n]) for n in range(SERIES_PREC + 1)]
    output = []

    for m in range(max_order + 1):
        total = RF(0)
        for n, coefficient in enumerate(coeffs):
            if coefficient == 0:
                continue
            if m == 0:
                weight = ZZ(1)
            elif n == 0:
                weight = ZZ(0)
            else:
                weight = ZZ(n)**m
            total += RF(coefficient)*RF(weight)*q0**n
        output.append(total)

    return output


Bjets = theta_jets_from_series(B_series, ORDER + 1)
Z2jets = theta_jets_from_series(Z2_series, ORDER + 1)
B2jets = theta_jets_from_series(B2_series, ORDER + 1)
P2jets = theta_jets_from_series(P2_series, ORDER + 1)

print("Series sanity checks:", flush=True)
print("  B(q0)        =", Bjets[0], flush=True)
print("  Z2(q0)       =", Z2jets[0], flush=True)
print("  P2(q0)       =", P2jets[0], flush=True)
print("  B(q0^2)      =", B2jets[0], flush=True)
print("  Z2=P2+B2/8 ? =", abs(Z2jets[0] - P2jets[0] - B2jets[0]/8) < ZERO_TOL, flush=True)
assert abs(Z2jets[0] - P2jets[0] - B2jets[0]/8) < ZERO_TOL
print(flush=True)

# ----------------------------------------------------------------------------
# Truncated Taylor-jet algebra in (x,y,z)
# ----------------------------------------------------------------------------


def total_degree(index):
    return int(index[0]) + int(index[1]) + int(index[2])


class Jet:
    def __init__(self, coeffs=None):
        self.c = {}
        if coeffs:
            for key, value in coeffs.items():
                kk = tuple(int(a) for a in key)
                vv = RF(value)
                if total_degree(kk) <= ORDER and vv != 0:
                    self.c[kk] = vv

    @staticmethod
    def const(value):
        return Jet({(0, 0, 0): RF(value)})

    @staticmethod
    def var(index, base):
        key = [0, 0, 0]
        key[index] = 1
        return Jet({(0, 0, 0): RF(base), tuple(key): RF(1)})

    def scale(self, value):
        value = RF(value)
        return Jet({key: value*entry for key, entry in self.c.items()})

    def __add__(self, other):
        other = tojet(other)
        out = dict(self.c)
        for key, value in other.c.items():
            out[key] = out.get(key, RF(0)) + value
        return Jet(out)

    __radd__ = __add__

    def __neg__(self):
        return Jet({key: -value for key, value in self.c.items()})

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
                key = (a[0]+b[0], a[1]+b[1], a[2]+b[2])
                if total_degree(key) <= ORDER:
                    out[key] = out.get(key, RF(0)) + va*vb
        return Jet(out)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, Jet):
            return self.scale(RF(1)/RF(other))
        return self*other.inv()

    def __rtruediv__(self, other):
        return tojet(other)*self.inv()

    def __pow__(self, exponent):
        exponent = int(exponent)
        if exponent < 0:
            return self.inv().__pow__(-exponent)
        if exponent == 0:
            return Jet.const(1)

        out = Jet.const(1)
        base = self
        n = exponent
        while n > 0:
            if n % 2 == 1:
                out = out*base
            base = base*base
            n //= 2
        return out

    def inv(self):
        a0 = self.c.get((0, 0, 0), RF(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.inv()")

        relative = (self - Jet.const(a0)).scale(RF(1)/a0)
        out = Jet.const(0)
        term = Jet.const(1)
        sign = RF(1)

        for _ in range(ORDER + 1):
            out = out + term.scale(sign)
            term = term*relative
            sign = -sign

        return out.scale(RF(1)/a0)

    def log(self):
        a0 = self.c.get((0, 0, 0), RF(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.log()")

        relative = (self - Jet.const(a0)).scale(RF(1)/a0)
        out = Jet.const(0)
        term = Jet.const(1)

        for power in range(1, ORDER + 1):
            term = term*relative
            out = out + term.scale(RF(QQ((-1)**(power+1))/power))

        return out

    def deriv_value(self, counts):
        counts = tuple(int(a) for a in counts)
        coefficient = self.c.get(counts, RF(0))
        multiplier = ZZ(1)
        for count in counts:
            multiplier *= factorial(count)
        return coefficient*RF(multiplier)


def tojet(value):
    if isinstance(value, Jet):
        return value
    return Jet.const(value)


# ----------------------------------------------------------------------------
# Matrix helpers
# ----------------------------------------------------------------------------


def mat_inv_3(M):
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]

    determinant = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)

    return [
        [(e*i-f*h)/determinant, (c*h-b*i)/determinant,
         (b*f-c*e)/determinant],
        [(f*g-d*i)/determinant, (a*i-c*g)/determinant,
         (c*d-a*f)/determinant],
        [(d*h-e*g)/determinant, (b*g-a*h)/determinant,
         (a*e-b*d)/determinant],
    ]


def matmul(A, B):
    return [
        [
            sum(A[row][k]*B[k][column] for k in range(3))
            for column in range(3)
        ]
        for row in range(3)
    ]


def matadd(A, B):
    return [
        [A[row][column] + B[row][column] for column in range(3)]
        for row in range(3)
    ]


def matscale(A, scalar):
    return [
        [scalar*A[row][column] for column in range(3)]
        for row in range(3)
    ]


def dot(A, vector_left, vector_right=None):
    if vector_right is None:
        vector_right = vector_left

    return sum(
        vector_left[row]*A[row][column]*vector_right[column]
        for row in range(3)
        for column in range(3)
    )


# ----------------------------------------------------------------------------
# Differential helpers
# ----------------------------------------------------------------------------


def normal_deriv_value(jet, indices):
    return jet.deriv_value((
        indices.count(0),
        indices.count(1),
        indices.count(2),
    ))


def twisted_deriv_value(jet, indices, p_degree):
    """Apply (partial_x-p_degree*L)^a partial_y^b partial_z^c."""
    a = indices.count(0)
    b = indices.count(1)
    c = indices.count(2)

    total = RF(0)
    for power in range(a + 1):
        total += (
            RF(comb(a, power))
            *(-RF(p_degree)*L)**(a-power)
            *jet.deriv_value((power, b, c))
        )
    return total


def one_variable_jet(jets, coordinate_index):
    out = Jet.const(0)
    key = [0, 0, 0]
    key[coordinate_index] = 1
    displacement = Jet({tuple(key): RF(1)})

    for m in range(ORDER + 1):
        out += ((-L)**m*jets[m]/RF(factorial(m)))*(displacement**m)

    return out


def one_variable_theta_jet(jets, coordinate_index):
    out = Jet.const(0)
    key = [0, 0, 0]
    key[coordinate_index] = 1
    displacement = Jet({tuple(key): RF(1)})

    for m in range(ORDER + 1):
        out += ((-L)**m*jets[m+1]/RF(factorial(m)))*(displacement**m)

    return out


# ----------------------------------------------------------------------------
# Build all order-p and order-p^2 data at one radial point X
# ----------------------------------------------------------------------------


def component_data_at_X(X_value, basis_name=None):
    """
    Return order-p^2 closure data.

    basis_name=None gives K2=0 and hence the fixed quadratic K1 x K1 term.
    basis_name in {'A','B','C'} inserts the corresponding F2 basis element.
    """
    Xv = RF(X_value)

    x = Jet.var(0, Xv)
    y = Jet.var(1, RF(1))
    z = Jet.var(2, RF(1))

    V = 9*x*y*z + RF(3)/2*y*y*z + RF(3)/2*y*z*z
    K0 = -V.log()

    By = one_variable_jet(Bjets, 1)
    Bz = one_variable_jet(Bjets, 2)
    Ty = one_variable_theta_jet(Bjets, 1)
    Tz = one_variable_theta_jet(Bjets, 2)

    K1bar = (
        By*Bz*(1 + L*x)
        + L*y*Ty*Bz
        + L*z*By*Tz
    )/(2*L**3*V)

    if basis_name is None:
        K2bar = Jet.const(0)
    else:
        Py = one_variable_jet(P2jets, 1)
        Pz = one_variable_jet(P2jets, 2)
        TPy = one_variable_theta_jet(P2jets, 1)
        TPz = one_variable_theta_jet(P2jets, 2)

        B2y = one_variable_jet(B2jets, 1)
        B2z = one_variable_jet(B2jets, 2)
        TB2y = one_variable_theta_jet(B2jets, 1)
        TB2z = one_variable_theta_jet(B2jets, 2)

        if basis_name == "A":
            F2 = Py*Pz
            ThetaYF2 = TPy*Pz
            ThetaZF2 = Py*TPz
        elif basis_name == "B":
            F2 = B2y*B2z
            ThetaYF2 = TB2y*B2z
            ThetaZF2 = B2y*TB2z
        elif basis_name == "C":
            F2 = Py*B2z + B2y*Pz
            ThetaYF2 = TPy*B2z + TB2y*Pz
            ThetaZF2 = Py*TB2z + B2y*TPz
        else:
            raise ValueError("basis_name must be None, 'A', 'B', or 'C'")

        K2bar = (
            F2*(1 + 2*L*x)
            + L*y*ThetaYF2
            + L*z*ThetaZF2
        )/(2*L**3*V)

    V0 = V.deriv_value((0, 0, 0))

    G0 = [
        [normal_deriv_value(K0, [row, column]) for column in range(3)]
        for row in range(3)
    ]

    G1 = [
        [twisted_deriv_value(K1bar, [row, column], 1) for column in range(3)]
        for row in range(3)
    ]

    G2 = [
        [twisted_deriv_value(K2bar, [row, column], 2) for column in range(3)]
        for row in range(3)
    ]

    M = [
        [normal_deriv_value(V, [row, column]) for column in range(3)]
        for row in range(3)
    ]

    e = [-(2*Xv+1), RF(1), RF(1)]
    o = [RF(0), RF(1), RF(-1)]
    radial = [Xv, RF(1), RF(1)]

    mu_e_1 = dot(G1, e)/dot(M, e)
    mu_o_1 = dot(G1, o)/dot(M, o)
    mu_r_1 = dot(G1, radial)/dot(M, radial)

    sigma1 = mu_o_1 - mu_e_1
    D1 = (mu_e_1 + mu_o_1)/2 + 2*mu_r_1

    mu_e_2 = dot(G2, e)/dot(M, e)
    mu_o_2 = dot(G2, o)/dot(M, o)
    mu_r_2 = dot(G2, radial)/dot(M, radial)

    sigma2 = mu_o_2 - mu_e_2
    D2 = (mu_e_2 + mu_o_2)/2 + 2*mu_r_2

    A0 = []
    A1 = []
    A2 = []

    for q in range(3):
        a0q = RF(0)
        a1q = RF(0)
        a2q = RF(0)

        for index in range(3):
            a0q += e[index]*(
                normal_deriv_value(K0, [index, 1, q])
                - normal_deriv_value(K0, [index, 2, q])
            )

            a1q += e[index]*(
                twisted_deriv_value(K1bar, [index, 1, q], 1)
                - twisted_deriv_value(K1bar, [index, 2, q], 1)
            )

            a2q += e[index]*(
                twisted_deriv_value(K2bar, [index, 1, q], 2)
                - twisted_deriv_value(K2bar, [index, 2, q], 2)
            )

        A0.append(a0q)
        A1.append(a1q)
        A2.append(a2q)

    K4_0 = RF(0)
    K4_1 = RF(0)
    K4_2 = RF(0)

    for first in range(3):
        for second in range(3):
            pattern = [first, second]

            K4_0 += e[first]*e[second]*(
                normal_deriv_value(K0, pattern + [1, 1])
                - 2*normal_deriv_value(K0, pattern + [1, 2])
                + normal_deriv_value(K0, pattern + [2, 2])
            )

            K4_1 += e[first]*e[second]*(
                twisted_deriv_value(K1bar, pattern + [1, 1], 1)
                - 2*twisted_deriv_value(K1bar, pattern + [1, 2], 1)
                + twisted_deriv_value(K1bar, pattern + [2, 2], 1)
            )

            K4_2 += e[first]*e[second]*(
                twisted_deriv_value(K2bar, pattern + [1, 1], 2)
                - 2*twisted_deriv_value(K2bar, pattern + [1, 2], 2)
                + twisted_deriv_value(K2bar, pattern + [2, 2], 2)
            )

    Ginv0 = mat_inv_3(G0)
    Ginv1 = matscale(matmul(matmul(Ginv0, G1), Ginv0), RF(-1))

    Ginv2 = matadd(
        matmul(matmul(matmul(matmul(Ginv0, G1), Ginv0), G1), Ginv0),
        matscale(matmul(matmul(Ginv0, G2), Ginv0), RF(-1)),
    )

    Curv0 = -K4_0 + dot(Ginv0, A0)

    Curv1 = (
        -K4_1
        + dot(Ginv1, A0)
        + dot(Ginv0, A1, A0)
        + dot(Ginv0, A0, A1)
    )

    Curv2 = (
        -K4_2
        + dot(Ginv2, A0)
        + dot(Ginv1, A1, A0)
        + dot(Ginv1, A0, A1)
        + dot(Ginv0, A2, A0)
        + dot(Ginv0, A0, A2)
        + dot(Ginv0, A1, A1)
    )

    ge0 = dot(G0, e)
    ge1 = dot(G1, e)
    ge2 = dot(G2, e)

    go0 = dot(G0, o)
    go1 = dot(G1, o)
    go2 = dot(G2, o)

    Den0 = ge0*go0
    Den1 = ge1*go0 + ge0*go1
    Den2 = ge2*go0 + ge1*go1 + ge0*go2

    Bcurv0 = Curv0/Den0
    Bcurv1 = Curv1/Den0 - Curv0*Den1/Den0**2
    Bcurv2 = (
        Curv2/Den0
        - Curv1*Den1/Den0**2
        + Curv0*(Den1**2/Den0**3 - Den2/Den0**2)
    )

    Scurv1 = -Bcurv1/V0
    Scurv2 = -Bcurv2/V0

    return {
        "Scurv1": Scurv1,
        "Scurv2": Scurv2,
        "sigma1": sigma1,
        "D1": D1,
        "sigma2": sigma2,
        "D2": D2,
    }


# ----------------------------------------------------------------------------
# Cache the affine closure columns at all sample points
# ----------------------------------------------------------------------------

all_X = []
for value in TRAIN_X + VALIDATE_X + CANDIDATE_INTERPOLATION_X:
    if value not in all_X:
        all_X.append(value)

print("Building full p^2 geometry at {} radial points...".format(len(all_X)), flush=True)

DATA = {}

for index, X_value in enumerate(all_X, start=1):
    base = component_data_at_X(X_value, None)
    colA = component_data_at_X(X_value, "A")
    colB = component_data_at_X(X_value, "B")
    colC = component_data_at_X(X_value, "C")

    # Fixed quadratic term from K1 x K1.
    Rquad = base["Scurv2"]

    # Linear K2 columns, with alpha_0 sigma_2 + beta_0 D_2 removed.
    def corrected_column(total):
        return (
            total["Scurv2"] - base["Scurv2"]
            - alpha0*total["sigma2"]
            - beta0*total["D2"]
        )

    DATA[X_value] = {
        "Rquad": Rquad,
        "rA": corrected_column(colA),
        "rB": corrected_column(colB),
        "rC": corrected_column(colC),
        "sigma1": base["sigma1"],
        "D1": base["D1"],
    }

    print("  completed {}/{} at X={}".format(index, len(all_X), X_value), flush=True)

print(flush=True)

# ----------------------------------------------------------------------------
# Solve the full natural-span closure problem
# ----------------------------------------------------------------------------
#
# At order p^2 we seek
#
#   Rquad + a*rA + b*rB + c*rC - alpha1*sigma1 - beta1*D1 = 0.
#
# Unknown vector:
#
#   u = (a,b,c,alpha1,beta1).
# ----------------------------------------------------------------------------

A_train = matrix(RF, [
    [
        DATA[Xv]["rA"],
        DATA[Xv]["rB"],
        DATA[Xv]["rC"],
        -DATA[Xv]["sigma1"],
        -DATA[Xv]["D1"],
    ]
    for Xv in TRAIN_X
])

b_train = vector(RF, [-DATA[Xv]["Rquad"] for Xv in TRAIN_X])

print("Natural factorized-span solve:", flush=True)
print("  training matrix determinant =", A_train.det(), flush=True)

if abs(A_train.det()) < ZERO_TOL:
    raise ArithmeticError(
        "The five-point interpolation matrix is numerically singular. "
        "Change TRAIN_X or increase precision."
    )

solution = A_train.solve_right(b_train)
a_sol, b_sol, c_sol, alpha1_sol, beta1_sol = solution

print("  solved coefficients:", flush=True)
print("    a       =", a_sol, flush=True)
print("    b       =", b_sol, flush=True)
print("    c       =", c_sol, flush=True)
print("    alpha_1 =", alpha1_sol, flush=True)
print("    beta_1  =", beta1_sol, flush=True)


def residual_for_solution(Xv, values):
    aa, bb, cc, al1, be1 = values
    row = DATA[Xv]
    return (
        row["Rquad"]
        + aa*row["rA"]
        + bb*row["rB"]
        + cc*row["rC"]
        - al1*row["sigma1"]
        - be1*row["D1"]
    )


validation_residuals = [
    residual_for_solution(Xv, solution)
    for Xv in VALIDATE_X
]

max_validation = max(abs(value) for value in validation_residuals)
span_preserves = bool(max_validation < ZERO_TOL)

print("  max validation residual =", max_validation, flush=True)
print("  natural factorized span restores closure?", span_preserves, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Test named candidate gluing prescriptions
# ----------------------------------------------------------------------------


def test_named_candidate(label, abc):
    aa, bb, cc = [RF(value) for value in abc]
    XA, XB = CANDIDATE_INTERPOLATION_X

    rhsA = (
        DATA[XA]["Rquad"]
        + aa*DATA[XA]["rA"]
        + bb*DATA[XA]["rB"]
        + cc*DATA[XA]["rC"]
    )

    rhsB = (
        DATA[XB]["Rquad"]
        + aa*DATA[XB]["rA"]
        + bb*DATA[XB]["rB"]
        + cc*DATA[XB]["rC"]
    )

    interp = matrix(RF, [
        [DATA[XA]["sigma1"], DATA[XA]["D1"]],
        [DATA[XB]["sigma1"], DATA[XB]["D1"]],
    ])

    rhs = vector(RF, [rhsA, rhsB])

    if abs(interp.det()) < ZERO_TOL:
        raise ArithmeticError("Candidate alpha/beta interpolation matrix singular.")

    shifts = interp.solve_right(rhs)
    al1, be1 = shifts

    values = vector(RF, [aa, bb, cc, al1, be1])
    residuals = [residual_for_solution(Xv, values) for Xv in VALIDATE_X]
    max_residual = max(abs(value) for value in residuals)
    preserves = bool(max_residual < ZERO_TOL)

    print("-"*78, flush=True)
    print(label, flush=True)
    print("  (a,b,c) =", (aa, bb, cc), flush=True)
    print("  fitted alpha_1 =", al1, flush=True)
    print("  fitted beta_1  =", be1, flush=True)
    print("  max validation residual =", max_residual, flush=True)
    print("  preserves order-p^2 two-channel closure?", preserves, flush=True)

    return {
        "label": label,
        "abc": (aa, bb, cc),
        "alpha1": al1,
        "beta1": be1,
        "max_residual": max_residual,
        "preserves": preserves,
    }


print("Named candidate tests:", flush=True)

candidate_results = []

candidate_results.append(test_named_candidate(
    "Candidate 1: pure primitive product P2(y)P2(z)",
    (1, 0, 0),
))

candidate_results.append(test_named_candidate(
    "Candidate 2: primitive product + GV double cover",
    (1, QQ(1)/8, 0),
))

candidate_results.append(test_named_candidate(
    "Candidate 3: uninserted product Z2(y)Z2(z)",
    (1, QQ(1)/64, QQ(1)/8),
))

# ----------------------------------------------------------------------------
# Compact interpretation
# ----------------------------------------------------------------------------

print("="*78, flush=True)
print("RESULT", flush=True)
print(flush=True)

if span_preserves:
    print(
        "A unique combination in the natural factorized degree-two span",
        flush=True,
    )
    print(
        "restores the two-channel closure at all independent validation points.",
        flush=True,
    )
    print(
        "This is a candidate ratio for the missing (2) and (1,1) contact",
        flush=True,
    )
    print(
        "contributions, but it must still be compared with the exact relative",
        flush=True,
    )
    print("Schoen degeneration formula before being called geometric.", flush=True)
else:
    print(
        "No combination in the natural factorized span",
        flush=True,
    )
    print(
        "{P2*P2, B(q^2)*B(q^2), sym(P2*B(q^2))}",
        flush=True,
    )
    print(
        "restores the order-p^2 two-channel closure.",
        flush=True,
    )
    print(
        "Therefore any restoration by the full Schoen gluing would require",
        flush=True,
    )
    print(
        "genuinely non-factorized relative-contact/Jacobi data, or additional",
        flush=True,
    )
    print("spectral channels at degree two.", flush=True)

print(flush=True)
print("Specific candidate verdicts:", flush=True)
for result in candidate_results:
    print(
        "  {} -> {}".format(result["label"], result["preserves"]),
        flush=True,
    )

print(flush=True)
print("="*78, flush=True)
print("SUCCESS", flush=True)
print(
    "The full p^2 curvature calculation and factorized-gluing feasibility test",
    flush=True,
)
print("completed.", flush=True)
print("Copy the complete compact output back into the chat.", flush=True)
print("="*78, flush=True)
