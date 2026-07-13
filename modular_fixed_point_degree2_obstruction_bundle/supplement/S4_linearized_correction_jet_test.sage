from sage.all import *
from math import factorial, comb

# ============================================================================
# LINEARIZED CORRECTION-JET CLOSURE TEST
# ============================================================================
#
# PURPOSE
# -------
# Test whether a first correction to the one-variable HSS/eta series preserves
# the two-channel spectral--curvature closure at the square-lattice point.
#
# We perturb the six normalized finite jets by
#
#       C_m(eps) = C0 * (u_m + eps*j_m),       m=0,...,5,
#
# where u_m are the eta-product/Ramanujan jets and j_0,...,j_5 are completely
# formal correction jets.  The calculation is exact to O(eps^2).
#
# We ALSO allow the closure coefficients themselves to shift:
#
#       alpha(eps) = alpha_0 + eps*alpha_1,
#       beta(eps)  = beta_0  + eps*beta_1.
#
# The script solves alpha_1 and beta_1 from two radial points, then checks the
# residual as an identity in the generic radial variable X.
#
# OUTPUT
# ------
# 1. Whether an arbitrary correction jet preserves closure (normally: no).
# 2. Whether the first corrected residual lies in <delta,E6>.
# 3. The exact rank and nullity of the linear constraints on (j0,...,j5) at
#    the square-lattice fixed point.
# 4. Sanity checks showing that pure rescaling and the formal E4-tangent are
#    in the preserving kernel.
#
# INTERPRETATION
# --------------
# This is a UNIVERSAL FINITE-JET TEST.  It does not insert a genuine E8/Jacobi
# correction because no such correction is present in the supplied source.
# Once actual correction jets are known, paste them into test_fixed_correction()
# at the bottom of the cell.
# ============================================================================

ORDER = 4
PRINT_CONSTRAINT_ROWS = False
PRINT_KERNEL_BASIS = False
PRINT_SHIFT_COEFFICIENTS = False

# Interpolation points used to solve alpha_1 and beta_1.
X_A = QQ(5)/4
X_B = QQ(3)/2

print("="*78)
print("LINEARIZED CORRECTION-JET CLOSURE TEST")
print("C_m(eps)=C0*(u_m+eps*j_m), m=0,...,5")
print("alpha and beta are allowed first-order shifts")
print("="*78)
print()

# ----------------------------------------------------------------------------
# Exact coefficient field
# ----------------------------------------------------------------------------

names = [
    "X", "L", "E4", "delta", "E6", "C0",
    "j0", "j1", "j2", "j3", "j4", "j5",
]

R = PolynomialRing(QQ, names)
K = R.fraction_field()

gens = R.gens()
(
    rX, rL, rE4, rdelta, rE6, rC0,
    rj0, rj1, rj2, rj3, rj4, rj5,
) = gens

(
    X, L, E4, delta, E6, C0,
    j0, j1, j2, j3, j4, j5,
) = [K(g) for g in gens]

j = [j0, j1, j2, j3, j4, j5]
E2 = (delta + 6)/L

print("Coefficient-ring variables:")
print(" ", names)
print("E2 represented as (delta+6)/L")
print()

# Dual-number coefficient ring: exact modulo eps^2.
PS = PowerSeriesRing(K, "eps", default_prec=2)
eps = PS.gen()

# ----------------------------------------------------------------------------
# Rational differentiation and Ramanujan jets
# ----------------------------------------------------------------------------

def dK(expr, rvar):
    """Differentiate a K-expression with respect to a polynomial generator."""
    expr = K(expr)
    num = R(expr.numerator())
    den = R(expr.denominator())
    return K(num.derivative(rvar)*den - num*den.derivative(rvar)) / K(den)**2


theta_delta = L*(E2**2 - E4)/12
theta_E4 = (E2*E4 - E6)/3
theta_E6 = (E2*E6 - E4**2)/2


def theta(expr):
    """Ramanujan theta=q*d/dq in variables delta,E4,E6."""
    return (
        dK(expr, rdelta)*theta_delta
        + dK(expr, rE4)*theta_E4
        + dK(expr, rE6)*theta_E6
    )


Alog = (1 - E2)/6

u = [K(1)]
for _ in range(5):
    u.append(K(theta(u[-1]) + Alog*u[-1]))

# Perturbed absolute jets, exact to first order in eps.
Ceps = [
    PS(C0*u[m]) + eps*PS(C0*j[m])
    for m in range(6)
]

print("Built eta-product jets u_0,...,u_5 and formal corrections j_0,...,j_5.")
print()

# ============================================================================
# Truncated Taylor-jet algebra in (x,y,z), with dual-number coefficients
# ============================================================================

def total_degree(index):
    return int(index[0]) + int(index[1]) + int(index[2])


class Jet:
    def __init__(self, coeffs=None):
        self.c = {}
        if coeffs:
            for key, value in coeffs.items():
                kk = tuple(int(a) for a in key)
                vv = PS(value)
                if total_degree(kk) <= ORDER and vv != 0:
                    self.c[kk] = vv

    @staticmethod
    def const(value):
        return Jet({(0, 0, 0): PS(value)})

    @staticmethod
    def var(index, base):
        key = [0, 0, 0]
        key[index] = 1
        return Jet({(0, 0, 0): PS(base), tuple(key): PS(1)})

    def scale(self, value):
        value = PS(value)
        return Jet({key: value*entry for key, entry in self.c.items()})

    def __add__(self, other):
        other = tojet(other)
        out = dict(self.c)
        for key, value in other.c.items():
            out[key] = out.get(key, PS(0)) + value
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
                    out[key] = out.get(key, PS(0)) + va*vb
        return Jet(out)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, Jet):
            return self.scale(PS(1)/PS(other))
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
        a0 = self.c.get((0, 0, 0), PS(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.inv()")

        relative = (self - Jet.const(a0)).scale(PS(1)/a0)
        out = Jet.const(0)
        term = Jet.const(1)
        sign = PS(1)

        for _ in range(ORDER + 1):
            out = out + term.scale(sign)
            term = term*relative
            sign = -sign

        return out.scale(PS(1)/a0)

    def log(self):
        """Formal log; the irrelevant constant log(a0) is omitted."""
        a0 = self.c.get((0, 0, 0), PS(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.log()")

        relative = (self - Jet.const(a0)).scale(PS(1)/a0)
        out = Jet.const(0)
        term = Jet.const(1)

        for power in range(1, ORDER + 1):
            term = term*relative
            out = out + term.scale(PS(QQ((-1)**(power+1))/power))

        return out

    def deriv_value(self, counts):
        counts = tuple(int(a) for a in counts)
        coefficient = self.c.get(counts, PS(0))
        multiplier = ZZ(1)
        for count in counts:
            multiplier *= factorial(count)
        return coefficient*multiplier


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


def mat_neg(A):
    return [[-A[row][column] for column in range(3)] for row in range(3)]


def dot(A, vector, other=None):
    if other is None:
        other = vector

    return sum(
        vector[row]*A[row][column]*other[column]
        for row in range(3)
        for column in range(3)
    )

# ----------------------------------------------------------------------------
# Derivative helpers and perturbed B jets
# ----------------------------------------------------------------------------

def normal_deriv_value(jet, indices):
    return jet.deriv_value((
        indices.count(0),
        indices.count(1),
        indices.count(2),
    ))


def pfree_deriv_value(jet, indices):
    """D_x=partial_x-L; D_y=partial_y; D_z=partial_z."""
    a = indices.count(0)
    b = indices.count(1)
    c = indices.count(2)

    total = PS(0)
    for power in range(a + 1):
        total += (
            comb(a, power)
            *(-L)**(a-power)
            *jet.deriv_value((power, b, c))
        )
    return total


def Bjet_y():
    out = Jet.const(0)
    dy = Jet({(0, 1, 0): PS(1)})

    for m in range(ORDER + 1):
        out += (
            (-L)**m*Ceps[m]/factorial(m)
        )*(dy**m)

    return out


def Bjet_z():
    out = Jet.const(0)
    dz = Jet({(0, 0, 1): PS(1)})

    for m in range(ORDER + 1):
        out += (
            (-L)**m*Ceps[m]/factorial(m)
        )*(dz**m)

    return out


def Thetajet_y():
    out = Jet.const(0)
    dy = Jet({(0, 1, 0): PS(1)})

    for m in range(ORDER + 1):
        out += (
            (-L)**m*Ceps[m+1]/factorial(m)
        )*(dy**m)

    return out


def Thetajet_z():
    out = Jet.const(0)
    dz = Jet({(0, 0, 1): PS(1)})

    for m in range(ORDER + 1):
        out += (
            (-L)**m*Ceps[m+1]/factorial(m)
        )*(dz**m)

    return out

# ----------------------------------------------------------------------------
# Series and fraction helpers
# ----------------------------------------------------------------------------

def coeff0(series):
    return K(PS(series)[0])


def coeff1(series):
    return K(PS(series)[1])


def primitive_numerator(expr):
    expr = K(expr)
    num = R(expr.numerator())

    try:
        return num.primitive_part()
    except Exception:
        content = num.content()
        return R(num/content) if content != 0 else num


def substitute_fraction(expr, substitutions):
    """Substitute polynomial generators into a reduced K-expression."""
    expr = K(expr)
    num = R(expr.numerator()).subs(substitutions)
    den = R(expr.denominator()).subs(substitutions)
    return K(num)/K(den)


def substitute_X(expr, value):
    return substitute_fraction(expr, {rX: QQ(value)})

# ============================================================================
# Build perturbed channels at generic X
# ============================================================================

def channels_at_generic_X():
    x = Jet.var(0, X)
    y = Jet.var(1, K(1))
    z = Jet.var(2, K(1))

    V = 9*x*y*z + QQ(3)/2*y*y*z + QQ(3)/2*y*z*z
    K0 = -V.log()

    By = Bjet_y()
    Bz = Bjet_z()
    Ty = Thetajet_y()
    Tz = Thetajet_z()

    K1bar = (
        By*Bz*(1 + L*x)
        + L*y*Ty*Bz
        + L*z*By*Tz
    )/(2*L**3*V)

    V0 = V.deriv_value((0, 0, 0))

    G0 = [
        [normal_deriv_value(K0, [row, column]) for column in range(3)]
        for row in range(3)
    ]

    G1 = [
        [pfree_deriv_value(K1bar, [row, column]) for column in range(3)]
        for row in range(3)
    ]

    M = [
        [normal_deriv_value(V, [row, column]) for column in range(3)]
        for row in range(3)
    ]

    e = [-(2*X+1), K(1), K(1)]
    o = [K(0), K(1), K(-1)]
    radial = [X, K(1), K(1)]

    mu_e = dot(G1, e)/dot(M, e)
    mu_o = dot(G1, o)/dot(M, o)
    mu_r = dot(G1, radial)/dot(M, radial)

    sigma = mu_o - mu_e
    Dchannel = (mu_e + mu_o)/2 + 2*mu_r

    A0 = []
    A1 = []

    for q in range(3):
        a0q = PS(0)
        a1q = PS(0)

        for index in range(3):
            a0q += e[index]*(
                normal_deriv_value(K0, [index, 1, q])
                - normal_deriv_value(K0, [index, 2, q])
            )

            a1q += e[index]*(
                pfree_deriv_value(K1bar, [index, 1, q])
                - pfree_deriv_value(K1bar, [index, 2, q])
            )

        A0.append(a0q)
        A1.append(a1q)

    K4_0 = PS(0)
    K4_1 = PS(0)

    for first in range(3):
        for second in range(3):
            K4_0 += e[first]*e[second]*(
                normal_deriv_value(K0, [first, second, 1, 1])
                - 2*normal_deriv_value(K0, [first, second, 1, 2])
                + normal_deriv_value(K0, [first, second, 2, 2])
            )

            K4_1 += e[first]*e[second]*(
                pfree_deriv_value(K1bar, [first, second, 1, 1])
                - 2*pfree_deriv_value(K1bar, [first, second, 1, 2])
                + pfree_deriv_value(K1bar, [first, second, 2, 2])
            )

    Ginv0 = mat_inv_3(G0)
    Ginv1 = mat_neg(matmul(matmul(Ginv0, G1), Ginv0))

    Curv0 = -K4_0 + dot(Ginv0, A0)
    Curv1 = (
        -K4_1
        + dot(Ginv1, A0)
        + dot(Ginv0, A1, A0)
        + dot(Ginv0, A0, A1)
    )

    Den0 = dot(G0, e)*dot(G0, o)
    Den1 = dot(G1, e)*dot(G0, o) + dot(G0, e)*dot(G1, o)

    B0 = Curv0/Den0
    B1 = Curv1/Den0 - Curv0*Den1/Den0**2

    Scurv = -B1/V0

    return {
        "Scurv": PS(Scurv),
        "sigma": PS(sigma),
        "D": PS(Dchannel),
        "B0": PS(B0),
    }


print("Building the perturbed metric, spectral, and curvature channels...")
channels = channels_at_generic_X()
print("Channels built exactly modulo eps^2.")
print()

# Baseline sanity check.
B0_base = coeff0(channels["B0"])
print("Baseline B0=-1/3 sanity check:", bool(primitive_numerator(B0_base + QQ(1)/3) == 0))
assert primitive_numerator(B0_base + QQ(1)/3) == 0
print()

# ----------------------------------------------------------------------------
# Base and first-order channel pieces
# ----------------------------------------------------------------------------

S0 = coeff0(channels["Scurv"])
S1 = coeff1(channels["Scurv"])

sigma0 = coeff0(channels["sigma"])
sigma1 = coeff1(channels["sigma"])

D0 = coeff0(channels["D"])
D1 = coeff1(channels["D"])

alpha0 = (L**2*E4 + 12)/72
beta0 = -(
    7*L**4*E4**2 - 552*L**2*E4 - 13392
)/(288*(L**2*E4 + 36))

base_defect = K(S0 - alpha0*sigma0 - beta0*D0)
print("Baseline closure numerator lies in the fixed-point ideal?")
base_num = R(base_defect.numerator())
base_fixed = R(base_num.subs({rdelta: 0, rE6: 0}))
print("  baseline numerator vanishes at delta=E6=0:", bool(base_fixed == 0))
assert base_fixed == 0
print()

# Raw first-order defect before shifting alpha and beta.
raw1 = K(S1 - alpha0*sigma1 - beta0*D1)

# ----------------------------------------------------------------------------
# Solve alpha_1 and beta_1 from two generic radial points
# ----------------------------------------------------------------------------

sigA = substitute_X(sigma0, X_A)
sigB = substitute_X(sigma0, X_B)
DA = substitute_X(D0, X_A)
DB = substitute_X(D0, X_B)
rawA = substitute_X(raw1, X_A)
rawB = substitute_X(raw1, X_B)

interp_det = K(sigA*DB - sigB*DA)
interp_det_fixed = substitute_fraction(interp_det, {rdelta: 0, rE6: 0})

print("Interpolation system:")
print("  X_A =", X_A)
print("  X_B =", X_B)
print("  determinant nonzero generically?", bool(interp_det != 0))
print("  determinant nonzero on fixed locus?", bool(interp_det_fixed != 0))

assert interp_det != 0
assert interp_det_fixed != 0

alpha1 = K((rawA*DB - rawB*DA)/interp_det)
beta1 = K((sigA*rawB - sigB*rawA)/interp_det)

if PRINT_SHIFT_COEFFICIENTS:
    print()
    print("alpha_1 =")
    print(alpha1)
    print()
    print("beta_1 =")
    print(beta1)

# Corrected first-order residual as a function of generic X.
residual1 = K(raw1 - alpha1*sigma0 - beta1*D0)

# It must vanish at the interpolation points by construction.
resA = substitute_X(residual1, X_A)
resB = substitute_X(residual1, X_B)

print("  residual vanishes at X_A?", bool(resA == 0))
print("  residual vanishes at X_B?", bool(resB == 0))
assert resA == 0
assert resB == 0
print()

# ----------------------------------------------------------------------------
# Strong fixed-point ideal-membership test
# ----------------------------------------------------------------------------

N1 = R(residual1.numerator())
N1_delta0 = R(N1.subs({rdelta: 0}))
N1_fixed = R(N1_delta0.subs({rE6: 0}))

membership = bool(N1_fixed == 0)

print("Strong corrected-defect test:")
print("  generic first-order residual identically zero?", bool(residual1 == 0))
print("  numerator term count:", len(N1.dict()))
print("  N1|delta=E6=0 is zero?", membership)

if membership:
    A_certificate, rem_delta = (N1 - N1_delta0).quo_rem(rdelta)
    B_certificate, rem_E6 = N1_delta0.quo_rem(rE6)

    certificate_ok = bool(
        rem_delta == 0
        and rem_E6 == 0
        and N1 == rdelta*A_certificate + rE6*B_certificate
    )

    print("  exact N1=delta*A+E6*B certificate?", certificate_ok)
    assert certificate_ok
else:
    print("  arbitrary correction jets do NOT preserve fixed-point closure.")

print()

# ----------------------------------------------------------------------------
# Reduce the fixed-locus obstruction and extract its linear constraints
# ----------------------------------------------------------------------------

residual_fixed = substitute_fraction(
    residual1,
    {rdelta: 0, rE6: 0, rC0: 1},
)

Pfixed_original = R(residual_fixed.numerator())

print("Fixed-locus corrected residual:")
print("  reduced numerator zero?", bool(Pfixed_original == 0))
print("  reduced numerator term count:", len(Pfixed_original.dict()))
print()

# Build QQ(L,E4)[X,j0,...,j5].
BaseR = PolynomialRing(QQ, names=("LL", "EE4"))
LL, EE4 = BaseR.gens()
BaseF = BaseR.fraction_field()

P = PolynomialRing(
    BaseF,
    names=("XX", "jj0", "jj1", "jj2", "jj3", "jj4", "jj5"),
)
XX, jj0, jj1, jj2, jj3, jj4, jj5 = P.gens()
jj = [jj0, jj1, jj2, jj3, jj4, jj5]
PF = P.fraction_field()

phi_values = [
    PF(XX),
    PF(BaseF(LL)),
    PF(BaseF(EE4)),
    PF(0),
    PF(0),
    PF(1),
    PF(jj0), PF(jj1), PF(jj2),
    PF(jj3), PF(jj4), PF(jj5),
]

phi = R.hom(phi_values, PF)
Pfixed_image = phi(Pfixed_original)

if Pfixed_image.denominator() != 1:
    raise ArithmeticError("Expected a polynomial after the fixed-locus map.")

Pfixed = P(Pfixed_image.numerator())

# Verify linearity in the correction jets and collect rows by X degree.
max_x_degree = Pfixed.degree(XX) if Pfixed != 0 else -1

constant_by_x = {
    degree: BaseF(0)
    for degree in range(max_x_degree + 1)
}

rows_by_x = {
    degree: [BaseF(0) for _ in range(6)]
    for degree in range(max_x_degree + 1)
}

nonlinear_terms = []

for exponent_tuple, coefficient in Pfixed.dict().items():
    x_degree = int(exponent_tuple[0])
    correction_exponents = exponent_tuple[1:]
    correction_degree = sum(int(value) for value in correction_exponents)

    if correction_degree == 0:
        constant_by_x[x_degree] += BaseF(coefficient)

    elif correction_degree == 1:
        positions = [
            index
            for index, exponent in enumerate(correction_exponents)
            if exponent == 1
        ]

        if len(positions) != 1:
            nonlinear_terms.append((exponent_tuple, coefficient))
        else:
            rows_by_x[x_degree][positions[0]] += BaseF(coefficient)

    else:
        nonlinear_terms.append((exponent_tuple, coefficient))

print("Constraint extraction:")
print("  degree in X:", max_x_degree)
print("  nonlinear correction terms found?", bool(len(nonlinear_terms) > 0))

assert len(nonlinear_terms) == 0
assert all(value == 0 for value in constant_by_x.values())

# Keep only nonzero coefficient rows.
active_degrees = []
active_rows = []

for degree in range(max_x_degree + 1):
    row = rows_by_x[degree]
    if any(entry != 0 for entry in row):
        active_degrees.append(degree)
        active_rows.append(row)

M = Matrix(BaseF, active_rows)
rank = M.rank()
nullity = 6 - rank

print("  active X-coefficient equations:", len(active_rows))
print("  active X degrees:", active_degrees)
print("  exact constraint-matrix rank:", rank)
print("  preserving-kernel dimension:", nullity)
print()

if PRINT_CONSTRAINT_ROWS:
    print("Constraint rows, columns=(j0,j1,j2,j3,j4,j5):")
    for degree, row in zip(active_degrees, active_rows):
        print("  [X^{}]:".format(degree))
        print("   ", row)
    print()

kernel = M.right_kernel()

if PRINT_KERNEL_BASIS:
    print("Kernel basis over QQ(L,E4):")
    for vector in kernel.basis():
        print(vector)
    print()

# ----------------------------------------------------------------------------
# Sanity directions that MUST preserve closure
# ----------------------------------------------------------------------------

# Map a K-expression at delta=E6=0,C0=1 into BaseF.
def to_BaseF_fixed(expr):
    fixed = substitute_fraction(
        expr,
        {rdelta: 0, rE6: 0, rC0: 1, rX: 0},
    )

    # The u_m do not depend on X or correction jets.
    num = R(fixed.numerator())
    den = R(fixed.denominator())

    # Direct substitution into BaseF.
    num_value = num(
        0, BaseF(LL), BaseF(EE4), 0, 0, 1,
        0, 0, 0, 0, 0, 0,
    )
    den_value = den(
        0, BaseF(LL), BaseF(EE4), 0, 0, 1,
        0, 0, 0, 0, 0, 0,
    )

    return BaseF(num_value)/BaseF(den_value)


u_fixed = vector(BaseF, [to_BaseF_fixed(u[m]) for m in range(6)])

# Formal E4 tangent of the fixed normalized jets.
# Differentiate first in K, then specialize.
u_E4_tangent = vector(
    BaseF,
    [to_BaseF_fixed(dK(u[m], rE4)) for m in range(6)],
)

scale_test = M*u_fixed
E4_test = M*u_E4_tangent

print("Mandatory kernel sanity checks:")
print("  pure B-rescaling direction lies in kernel?", bool(scale_test == 0))
print("  formal E4-tangent direction lies in kernel?", bool(E4_test == 0))

assert scale_test == 0
assert E4_test == 0

independent_sanity_directions = Matrix(
    BaseF,
    [list(u_fixed), list(u_E4_tangent)],
).rank()

print("  these two sanity directions are independent?", bool(independent_sanity_directions == 2))
assert independent_sanity_directions == 2
print()

# ----------------------------------------------------------------------------
# Reusable concrete-correction tester
# ----------------------------------------------------------------------------

def test_fixed_correction(jet_values, label="user correction"):
    """
    Test six normalized correction jets against the exact fixed-locus
    constraint matrix.

    INPUT
    -----
    jet_values : list/tuple of six expressions in LL and EE4, coercible to
                 BaseF=QQ(LL,EE4).

    EXAMPLE
    -------
    test_fixed_correction(
        [1, 0, 0, 0, 0, 0],
        label="change only C0",
    )

    For an actual HSS correction, supply

        [j0(LL,EE4), ..., j5(LL,EE4)].
    """
    if len(jet_values) != 6:
        raise ValueError("Provide exactly six normalized correction jets.")

    vector_value = vector(BaseF, [BaseF(value) for value in jet_values])
    obstruction = M*vector_value
    preserves = bool(obstruction == 0)

    print("-"*78)
    print("Concrete correction test:", label)
    print("  preserves two-channel closure at first order?", preserves)

    if not preserves:
        nonzero = [
            (active_degrees[index], value)
            for index, value in enumerate(obstruction)
            if value != 0
        ]

        print("  number of nonzero coefficient obstructions:", len(nonzero))
        print("  first nonzero obstruction(s):")

        for degree, value in nonzero[:5]:
            print("    X^{} : {}".format(degree, value))

    return preserves, obstruction


# Two built-in tests.
test_fixed_correction(list(u_fixed), label="pure rescaling B -> (1+eps)B")
test_fixed_correction(list(u_E4_tangent), label="formal E4 tangent")

print()
print("="*78)
print("SUCCESS")
print()

if membership:
    print("An arbitrary formal correction jet preserves fixed-point ideal closure.")
else:
    print("An arbitrary formal correction jet does NOT preserve closure.")
    print("The preserving corrections form the exact kernel reported above.")

print()
print("Next use:")
print("  Insert six genuine normalized HSS/Jacobi correction jets into")
print("  test_fixed_correction([j0,...,j5], label='first HSS correction').")
print()
print("Please copy the complete printed output back into the chat.")
print("="*78)
