from sage.all import *
from math import factorial, comb
import mpmath as mp

# ============================================================
# FORMAL FINITE-JET OBSTRUCTION EXTRACTOR
#
# Goal:
#   Extract the jet-level relation responsible for the exact-looking
#   HSS closure on v=(x,1,1).
#
# We treat
#
#   L = 2*pi
#   C_m = (theta^m B)(e^{-2*pi})
#
# as formal variables, and compute the determinant obstruction:
#
#   det [[F(x1), S(x1), D(x1)],
#        [F(x2), S(x2), D(x2)],
#        [F(x3), S(x3), D(x3)]]
#
# where
#
#   F = Scurv_bar
#   S = sigma_bar = (mu_o - mu_e)/p
#   D = Dchan_bar = (Delta1/(2V))/p
#
# If the obstruction has a common factor across several x3 values,
# and that factor vanishes on the true HSS jets, that factor is the
# candidate relation we need to prove from the HSS/Picard-Fuchs/modular data.
# ============================================================

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

ORDER = 4
JMAX = 5  # need C0,...,C5

X1 = QQ(6)/5
X2 = QQ(2)

# Start with two points. Add more after first success.
OBS_POINTS = [QQ(5)/4, QQ(7)/5, QQ(3)/2]

DO_GCD = True
DO_FACTOR_GCD = True
DO_FACTOR_OBSTRUCTIONS = False   # can be very verbose

# If something stalls, first set:
#   OBS_POINTS = [QQ(5)/4]
#   DO_GCD = False
#   DO_FACTOR_GCD = False

# HSS B coefficients through U^50.
B_COEFFS = [
    ZZ(9),
    ZZ(36),
    ZZ(126),
    ZZ(360),
    ZZ(945),
    ZZ(2268),
    ZZ(5166),
    ZZ(11160),
    ZZ(23220),
    ZZ(46620),
    ZZ(90972),
    ZZ(172872),
    ZZ(321237),
    ZZ(584640),
    ZZ(1044810),
    ZZ(1835856),
    ZZ(3177153),
    ZZ(5421132),
    ZZ(9131220),
    ZZ(15195600),
    ZZ(25006653),
    ZZ(40722840),
    ZZ(65670768),
    ZZ(104930280),
    ZZ(166214205),
    ZZ(261141300),
    ZZ(407118726),
    ZZ(630048384),
    ZZ(968272605),
    ZZ(1478208420),
    ZZ(2242463580),
    ZZ(3381344280),
    ZZ(5069259342),
    ZZ(7557818940),
    ZZ(11208455370),
    ZZ(16538048640),
    ZZ(24282822798),
    ZZ(35487134928),
    ZZ(51626878470),
    ZZ(74779896240),
    ZZ(107861179482),
    ZZ(154945739844),
    ZZ(221711362038),
    ZZ(316042958880),
    ZZ(448856366490),
    ZZ(635216766732),
    ZZ(895854679650),
    ZZ(1259213600736),
    ZZ(1764210946995),
    ZZ(2463949037340),
    ZZ(3430694064888),
]

# ------------------------------------------------------------
# Polynomial ring for formal jet constants
# ------------------------------------------------------------

print("Setting up polynomial ring...")

names = ["L"] + ["C{}".format(i) for i in range(JMAX + 1)]
R = PolynomialRing(QQ, names)
K = R.fraction_field()

L = K(R.gen(0))
C = [K(R.gen(1+i)) for i in range(JMAX + 1)]

print("Ring variables:", names)
print()

# ============================================================
# Generic local Taylor-jet algebra over K
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
        return Jet({(0,0,0): K(a)})

    @staticmethod
    def var(idx, base):
        k = [0,0,0]
        k[idx] = 1
        return Jet({
            (0,0,0): K(base),
            tuple(k): K(1)
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
        a0 = self.c.get((0,0,0), K(0))
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
        Formal log. Constant log(a0) is irrelevant for derivatives,
        so we omit it. All derivatives of order >=1 are correct.
        """
        a0 = self.c.get((0,0,0), K(0))
        if a0 == 0:
            raise ZeroDivisionError("zero constant term in Jet.log()")

        r = (self - Jet.const(a0)).scale(K(1)/a0)

        out = Jet.const(0)
        term = Jet.const(1)

        for k in range(1, ORDER + 1):
            term = term * r
            out = out + term.scale(K((-1)**(k+1)) / K(k))

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
# Jet derivative helpers
# ============================================================

def pfree_deriv_value(jet, inds):
    """
    D_x = partial_x - L, D_y=partial_y, D_z=partial_z.
    """
    a = inds.count(0)
    b = inds.count(1)
    c = inds.count(2)

    total = K(0)

    for j in range(a + 1):
        total += K(comb(a, j)) * ((-L)**(a-j)) * jet.deriv_value((j,b,c))

    return total


def normal_deriv_value(jet, inds):
    return jet.deriv_value((inds.count(0), inds.count(1), inds.count(2)))


def Bjet_y():
    out = Jet.const(0)
    dy = Jet({(0,1,0): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m] / K(factorial(m))) * (dy**m)

    return out


def Bjet_z():
    out = Jet.const(0)
    dz = Jet({(0,0,1): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m] / K(factorial(m))) * (dz**m)

    return out


def Thetajet_y():
    out = Jet.const(0)
    dy = Jet({(0,1,0): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m+1] / K(factorial(m))) * (dy**m)

    return out


def Thetajet_z():
    out = Jet.const(0)
    dz = Jet({(0,0,1): K(1)})

    for m in range(ORDER + 1):
        out = out + ((-L)**m * C[m+1] / K(factorial(m))) * (dz**m)

    return out


# ============================================================
# Formal channels at a fixed x
# ============================================================

def channels_at_x(x0):
    """
    Return formal K-elements:
      F = Scurv_bar
      S = sigma_bar
      D = Dchan_bar
      B0 sanity value
      denominator sanity residual
    at v=(x0,1,1).
    """
    x0 = QQ(x0)
    s = QQ(1)

    x = Jet.var(0, K(x0))
    y = Jet.var(1, K(s))
    z = Jet.var(2, K(s))

    V = 9*x*y*z + K(QQ(3)/2)*y*y*z + K(QQ(3)/2)*y*z*z
    K0 = -V.log()

    By = Bjet_y()
    Bz = Bjet_z()
    Ty = Thetajet_y()
    Tz = Thetajet_z()

    K1bar = (
        By*Bz*(1 + L*x)
        + L*y*Ty*Bz
        + L*z*By*Tz
    ) / (2 * L**3 * V)

    V0 = V.deriv_value((0,0,0))

    G0 = [
        [normal_deriv_value(K0, [i,j]) for j in range(3)]
        for i in range(3)
    ]

    G1bar = [
        [pfree_deriv_value(K1bar, [i,j]) for j in range(3)]
        for i in range(3)
    ]

    M = [
        [normal_deriv_value(V, [i,j]) for j in range(3)]
        for i in range(3)
    ]

    e = [-(2*K(x0) + 1), K(1), K(1)]
    o = [K(0), K(1), K(-1)]
    r = [K(x0), K(1), K(1)]

    mu_e = dot(G1bar, e) / dot(M, e)
    mu_o = dot(G1bar, o) / dot(M, o)
    mu_r = dot(G1bar, r) / dot(M, r)

    sigma = mu_o - mu_e
    Dchan = (mu_e + mu_o)/2 + 2*mu_r

    A0 = []
    A1 = []

    for q in range(3):
        A0q = K(0)
        A1q = K(0)

        for i in range(3):
            A0q += e[i] * (
                normal_deriv_value(K0, [i,1,q])
                - normal_deriv_value(K0, [i,2,q])
            )

            A1q += e[i] * (
                pfree_deriv_value(K1bar, [i,1,q])
                - pfree_deriv_value(K1bar, [i,2,q])
            )

        A0.append(A0q)
        A1.append(A1q)

    K4_0 = K(0)
    K4_1 = K(0)

    for i in range(3):
        for j in range(3):
            K4_0 += e[i]*e[j] * (
                normal_deriv_value(K0, [i,j,1,1])
                - 2*normal_deriv_value(K0, [i,j,1,2])
                + normal_deriv_value(K0, [i,j,2,2])
            )

            K4_1 += e[i]*e[j] * (
                pfree_deriv_value(K1bar, [i,j,1,1])
                - 2*pfree_deriv_value(K1bar, [i,j,1,2])
                + pfree_deriv_value(K1bar, [i,j,2,2])
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


def primitive_numerator(expr):
    """
    Return primitive numerator in R for a K-element.
    """
    expr = K(expr)
    num = R(expr.numerator())

    try:
        return num.primitive_part()
    except Exception:
        cont = num.content()
        if cont != 0:
            return R(num / cont)
        return num


def obstruction_polynomial(x3):
    """
    Obstruction determinant numerator for x1,x2,x3.
    """
    print("  building channels at x1 =", X1)
    ch1 = channels_at_x(X1)

    print("  building channels at x2 =", X2)
    ch2 = channels_at_x(X2)

    print("  building channels at x3 =", x3)
    ch3 = channels_at_x(x3)

    print("  sanity B0 at x1: B0+1/3 numerator zero?",
          bool(primitive_numerator(ch1["B0"] + K(QQ(1)/3)) == 0))

    print("  sanity den at x1 numerator zero?",
          bool(primitive_numerator(ch1["den_check"]) == 0))

    Mdet = matrix(K, [
        [ch1["F"], ch1["S"], ch1["D"]],
        [ch2["F"], ch2["S"], ch2["D"]],
        [ch3["F"], ch3["S"], ch3["D"]],
    ]).det()

    num = primitive_numerator(Mdet)
    return num


# ============================================================
# Actual HSS jet evaluation helpers
# ============================================================

def rational_to_mpf(q):
    q = QQ(q)
    return mp.mpf(str(q.numerator())) / mp.mpf(str(q.denominator()))


def compute_actual_hss_jets():
    mp.mp.dps = 100
    L0 = 2 * mp.pi
    q = mp.e**(-L0)

    Cvals = []

    for m in range(JMAX + 1):
        total = mp.mpf(0)

        for n, b in enumerate(B_COEFFS):
            total += mp.mpf(str(b)) * (mp.mpf(n)**m) * (q**n)

        Cvals.append(total)

    return L0, Cvals


def eval_poly_mpmath(poly, L0, Cvals):
    """
    Evaluate poly in R at actual numerical L,C.
    Also return a scale = sum abs(terms), so abs(val)/scale is meaningful.
    """
    mp.mp.dps = 100

    vals = [L0] + Cvals

    total = mp.mpf(0)
    scale = mp.mpf(0)

    for exp, coeff in poly.dict().items():
        coeff_mp = rational_to_mpf(coeff)
        term = coeff_mp

        for i, ei in enumerate(exp):
            if ei:
                term *= vals[i] ** int(ei)

        total += term
        scale += abs(term)

    rel = abs(total) / max(scale, mp.mpf("1e-100"))

    return total, scale, rel


# ============================================================
# Main run
# ============================================================

print("="*70)
print("FORMAL FINITE-JET OBSTRUCTION EXTRACTOR")
print("="*70)
print("OBS_POINTS:", OBS_POINTS)
print()

obs_polys = []

for x3 in OBS_POINTS:
    print("="*70)
    print("Computing obstruction for x3 =", x3)
    print("="*70)

    P = obstruction_polynomial(x3)
    obs_polys.append(P)

    print()
    print("Obstruction polynomial summary for x3 =", x3)
    print("  zero?       ", bool(P == 0))
    print("  degree:     ", P.degree())
    print("  term count: ", len(P.dict()))

    for name, gen in zip(names, R.gens()):
        print("  degree in {}: {}".format(name, P.degree(gen)))

    if DO_FACTOR_OBSTRUCTIONS:
        print()
        print("Factoring obstruction polynomial...")
        try:
            print(P.factor())
        except Exception as err:
            print("factor failed:", err)

    print()

# ------------------------------------------------------------
# GCD across obstruction polynomials
# ------------------------------------------------------------

if DO_GCD and len(obs_polys) >= 2:
    print("="*70)
    print("Computing gcd across obstruction polynomials...")
    print("="*70)

    Gpoly = obs_polys[0]

    for P in obs_polys[1:]:
        print("  gcd step...")
        Gpoly = Gpoly.gcd(P)
        try:
            Gpoly = Gpoly.primitive_part()
        except Exception:
            pass

        print("    current gcd degree:", Gpoly.degree())
        print("    current gcd term count:", len(Gpoly.dict()))

    print()
    print("GCD summary:")
    print("  zero?       ", bool(Gpoly == 0))
    print("  degree:     ", Gpoly.degree())
    print("  term count: ", len(Gpoly.dict()))

    for name, gen in zip(names, R.gens()):
        print("  degree in {}: {}".format(name, Gpoly.degree(gen)))

    # Evaluate on actual HSS jets.
    print()
    print("Evaluating GCD on actual 51-term HSS jets...")
    L0, Cvals = compute_actual_hss_jets()
    val, scale, rel = eval_poly_mpmath(Gpoly, L0, Cvals)

    print("  value =", mp.nstr(val, 30))
    print("  scale =", mp.nstr(scale, 30))
    print("  relative |value|/scale =", mp.nstr(rel, 30))

    if DO_FACTOR_GCD:
        print()
        print("Factoring gcd...")
        try:
            fac = Gpoly.factor()
            print(fac)

            print()
            print("Evaluating individual factors on actual HSS jets:")
            for idx, item in enumerate(fac):
                factor_poly, mult = item
                valf, scalef, relf = eval_poly_mpmath(factor_poly, L0, Cvals)

                print("factor", idx, "multiplicity", mult)
                print("  degree:", factor_poly.degree())
                print("  terms:", len(factor_poly.dict()))
                print("  value:", mp.nstr(valf, 30))
                print("  relative:", mp.nstr(relf, 30))

        except Exception as err:
            print("factor failed:", err)

else:
    Gpoly = None
    print("Skipped gcd step.")

# ------------------------------------------------------------
# Evaluate each obstruction on true HSS jets too
# ------------------------------------------------------------

print()
print("="*70)
print("Evaluating each obstruction on actual 51-term HSS jets")
print("="*70)

L0, Cvals = compute_actual_hss_jets()

for x3, P in zip(OBS_POINTS, obs_polys):
    val, scale, rel = eval_poly_mpmath(P, L0, Cvals)

    print("x3 =", x3)
    print("  value =", mp.nstr(val, 30))
    print("  scale =", mp.nstr(scale, 30))
    print("  relative |value|/scale =", mp.nstr(rel, 30))

print()
print("="*70)
print("Interpretation")
print("="*70)
print("- If obstruction polynomials are nonzero, the closure is not formal.")
print("- A nontrivial gcd suggests a common jet-level relation.")
print("- If the gcd/factor evaluates near zero on actual HSS jets, that is the candidate relation.")
print("- If gcd=1, the relation may be more complicated or the obstruction may not have a simple common polynomial factor.")
print("Done.")

# ============================================================
# ELIMINATE C5 FROM THE OBSTRUCTION POLYNOMIALS
#
# Run this AFTER the formal finite-jet obstruction extractor.
#
# We have obstruction polynomials P_i(L,C0,...,C5).
# Each was degree 1 in C5.
#
# If P_i = A_i*C5 + B_i,
# then eliminating C5 between P_i and P_j gives:
#
#   Q_ij = A_j*B_i - A_i*B_j.
#
# The true HSS jets should make Q_ij vanish.
# ============================================================

print("="*70)
print("C5 ELIMINATION STEP")
print("="*70)

try:
    obs_polys
except NameError:
    raise RuntimeError("obs_polys is not defined. Run the obstruction extractor first.")

try:
    R
    names
except NameError:
    raise RuntimeError("R/names not found. Run this in the same session as the obstruction extractor.")

# Identify C5 generator.
C5_gen = R.gen(names.index("C5"))

def coeff_in_var(poly, var, deg):
    """
    Coefficient of var^deg in poly.
    """
    return R(poly).coefficient({var: deg})

def split_linear(poly, var):
    """
    Return A,B such that poly = A*var + B.
    Assumes degree in var <= 1.
    """
    P = R(poly)
    d = P.degree(var)
    if d > 1:
        raise ValueError("Polynomial is not linear in {}".format(var))
    A = coeff_in_var(P, var, 1)
    B = coeff_in_var(P, var, 0)
    return R(A), R(B)

AB = []

print("Splitting obstruction polynomials as A*C5+B...")
for idx, P in enumerate(obs_polys):
    A, B = split_linear(P, C5_gen)
    AB.append((A, B))

    print("P{}: deg_C5={} | terms A={} | terms B={}".format(
        idx,
        P.degree(C5_gen),
        len(A.dict()),
        len(B.dict())
    ))

print()
print("Building pairwise C5-elimination polynomials Q_ij...")

elim_C5 = []

for i in range(len(AB)):
    for j in range(i+1, len(AB)):
        Ai, Bi = AB[i]
        Aj, Bj = AB[j]

        Q = R(Aj*Bi - Ai*Bj)

        try:
            Q = Q.primitive_part()
        except Exception:
            cont = Q.content()
            if cont != 0:
                Q = R(Q/cont)

        elim_C5.append(((i,j), Q))

        print()
        print("Q_{}{} summary:".format(i,j))
        print("  zero?      ", bool(Q == 0))
        print("  degree:    ", Q.degree())
        print("  term count:", len(Q.dict()))
        for nm, gen in zip(names, R.gens()):
            print("  degree in {}: {}".format(nm, Q.degree(gen)))

# Evaluate on actual HSS jets.
print()
print("="*70)
print("Evaluating C5-elimination polynomials on actual HSS jets")
print("="*70)

L0, Cvals = compute_actual_hss_jets()

for (i,j), Q in elim_C5:
    val, scale, rel = eval_poly_mpmath(Q, L0, Cvals)

    print("Q_{}{}:".format(i,j))
    print("  value =", mp.nstr(val, 30))
    print("  scale =", mp.nstr(scale, 30))
    print("  relative =", mp.nstr(rel, 30))

# GCD among elimination polynomials.
print()
print("="*70)
print("GCD among C5-elimination polynomials")
print("="*70)

if len(elim_C5) >= 2:
    GQ = elim_C5[0][1]
    for _, Q in elim_C5[1:]:
        print("gcd step...")
        GQ = GQ.gcd(Q)
        try:
            GQ = GQ.primitive_part()
        except Exception:
            pass
        print("  current degree:", GQ.degree())
        print("  current terms:", len(GQ.dict()))

    print()
    print("GQ summary:")
    print("  zero?      ", bool(GQ == 0))
    print("  degree:    ", GQ.degree())
    print("  term count:", len(GQ.dict()))
    for nm, gen in zip(names, R.gens()):
        print("  degree in {}: {}".format(nm, GQ.degree(gen)))

    val, scale, rel = eval_poly_mpmath(GQ, L0, Cvals)
    print()
    print("GQ actual-HSS evaluation:")
    print("  value =", mp.nstr(val, 30))
    print("  scale =", mp.nstr(scale, 30))
    print("  relative =", mp.nstr(rel, 30))

    print()
    print("Trying to factor GQ...")
    try:
        fac = GQ.factor()
        print(fac)

        print()
        print("Evaluating individual GQ factors:")
        for idx, item in enumerate(fac):
            factor_poly, mult = item
            valf, scalef, relf = eval_poly_mpmath(factor_poly, L0, Cvals)
            print("factor", idx, "multiplicity", mult)
            print("  degree:", factor_poly.degree())
            print("  terms:", len(factor_poly.dict()))
            print("  value:", mp.nstr(valf, 30))
            print("  relative:", mp.nstr(relf, 30))

    except Exception as err:
        print("factor failed:", err)

else:
    print("Need at least two elimination polynomials.")

print()
print("Done.")

# ============================================================
# IMPOSE THE LOG-DERIVATIVE RELATION AND SOLVE FOR C5
#
# Run this in the SAME Sage session after:
#   - obstruction extractor
#   - C5 elimination
#
# We found the key factor:
#
#   R1 = -L*C0 + 6*L*C1 + 6*C0.
#
# Equivalently:
#
#   C1 = (L-6)*C0/(6L).
#
# This code substitutes that relation into each original obstruction
# polynomial P_i and solves the resulting linear equation for C5.
# ============================================================

print("="*70)
print("IMPOSE R1 AND SOLVE ORIGINAL OBSTRUCTIONS FOR C5")
print("="*70)

try:
    obs_polys
    R
    names
    compute_actual_hss_jets
except NameError:
    raise RuntimeError("Run this in the same session as the obstruction extractor.")

# New reduced ring: variables are L, C0, C2, C3, C4, C5.
red_names = ["L", "C0", "C2", "C3", "C4", "C5"]
Rred = PolynomialRing(QQ, red_names)
Kred = Rred.fraction_field()

Lr  = Kred(Rred.gen(0))
C0r = Kred(Rred.gen(1))
C2r = Kred(Rred.gen(2))
C3r = Kred(Rred.gen(3))
C4r = Kred(Rred.gen(4))
C5r = Kred(Rred.gen(5))

# Substitute C1 = (L-6)*C0/(6L).
C1_expr = ((Lr - 6) * C0r) / (6 * Lr)

# Map old variables [L,C0,C1,C2,C3,C4,C5] into reduced field.
phi_vals = [
    Lr,
    C0r,
    C1_expr,
    C2r,
    C3r,
    C4r,
    C5r,
]

phi_rel = R.hom(phi_vals, Kred)

def primitive_poly_red(expr):
    """
    Convert reduced fraction expression to a primitive numerator polynomial.
    """
    expr = Kred(expr)
    num = Rred(expr.numerator())
    try:
        return num.primitive_part()
    except Exception:
        cont = num.content()
        if cont != 0:
            return Rred(num / cont)
        return num

def split_linear_red(poly, var):
    """
    Split poly = A*var + B.
    """
    P = Rred(poly)
    d = P.degree(var)
    if d > 1:
        raise ValueError("Polynomial not linear in requested variable.")
    A = P.coefficient({var: 1})
    B = P.coefficient({var: 0})
    return Rred(A), Rred(B)

def rational_to_mpf(q):
    q = QQ(q)
    return mp.mpf(str(q.numerator())) / mp.mpf(str(q.denominator()))

def eval_red_poly_mpmath(poly, L0, Cvals):
    """
    Evaluate a polynomial in Rred at:
      L, C0, C2, C3, C4, C5.
    """
    vals = [
        L0,
        Cvals[0],
        Cvals[2],
        Cvals[3],
        Cvals[4],
        Cvals[5],
    ]

    total = mp.mpf(0)
    scale = mp.mpf(0)

    for exp, coeff in Rred(poly).dict().items():
        term = rational_to_mpf(coeff)
        for i, ei in enumerate(exp):
            if ei:
                term *= vals[i] ** int(ei)
        total += term
        scale += abs(term)

    rel = abs(total) / max(scale, mp.mpf("1e-100"))
    return total, scale, rel

def eval_red_frac_mpmath(expr, L0, Cvals):
    """
    Evaluate fraction expression in Kred.
    """
    expr = Kred(expr)
    num = Rred(expr.numerator())
    den = Rred(expr.denominator())

    nval, _, _ = eval_red_poly_mpmath(num, L0, Cvals)
    dval, _, _ = eval_red_poly_mpmath(den, L0, Cvals)

    return nval / dval

# Actual HSS jets from 51 terms.
L0, Cvals = compute_actual_hss_jets()

print("Actual C5 =", mp.nstr(Cvals[5], 50))
print()

C5_formulas = []

for idx, P in enumerate(obs_polys):
    print("="*70)
    print("Processing obstruction P{}".format(idx))
    print("="*70)

    # Substitute R1 relation.
    P_rel_expr = phi_rel(P)

    # Clear denominator.
    P_rel = primitive_poly_red(P_rel_expr)

    print("After imposing R1:")
    print("  zero?      ", bool(P_rel == 0))
    print("  degree:    ", P_rel.degree())
    print("  term count:", len(P_rel.dict()))
    for nm, gen in zip(red_names, Rred.gens()):
        print("  degree in {}: {}".format(nm, P_rel.degree(gen)))

    # Evaluate reduced obstruction at actual HSS jets.
    val, scale, rel = eval_red_poly_mpmath(P_rel, L0, Cvals)
    print()
    print("  actual-HSS evaluation:")
    print("    value =", mp.nstr(val, 30))
    print("    scale =", mp.nstr(scale, 30))
    print("    relative =", mp.nstr(rel, 30))

    # Split as A*C5+B.
    A, B = split_linear_red(P_rel, Rred.gen(5))

    print()
    print("  Split P_rel = A*C5+B:")
    print("    A terms:", len(A.dict()), "degree:", A.degree())
    print("    B terms:", len(B.dict()), "degree:", B.degree())

    C5_pred = -Kred(B) / Kred(A)
    C5_formulas.append(C5_pred)

    pred_val = eval_red_frac_mpmath(C5_pred, L0, Cvals)
    err = pred_val - Cvals[5]
    relerr = abs(err) / max(abs(Cvals[5]), mp.mpf("1e-100"))

    print()
    print("  C5 predicted from P{}:".format(idx))
    print("    pred =", mp.nstr(pred_val, 50))
    print("    actual C5 =", mp.nstr(Cvals[5], 50))
    print("    abs error =", mp.nstr(err, 30))
    print("    rel error =", mp.nstr(relerr, 30))

    print()

# Compare C5 formulas pairwise.
print("="*70)
print("PAIRWISE COMPARISON OF C5 FORMULAS")
print("="*70)

for i in range(len(C5_formulas)):
    for j in range(i+1, len(C5_formulas)):
        diff = Kred(C5_formulas[i] - C5_formulas[j])
        num = primitive_poly_red(diff)

        print("formula {} - formula {}:".format(i, j))
        print("  numerator zero?", bool(num == 0))
        print("  numerator degree:", num.degree())
        print("  numerator terms:", len(num.dict()))

        val, scale, rel = eval_red_poly_mpmath(num, L0, Cvals)
        print("  actual-HSS numerator relative:", mp.nstr(rel, 30))
        print()

# Optional: factor the first C5 relation numerator.
print("="*70)
print("FACTOR FIRST REDUCED OBSTRUCTION")
print("="*70)

try:
    P0_rel = primitive_poly_red(phi_rel(obs_polys[0]))
    fac = P0_rel.factor()
    print(fac)
except Exception as err:
    print("factor failed:")
    print(err)

print()
print("Interpretation:")
print("- If all C5 formulas agree symbolically, then R1 makes the obstruction system rank one.")
print("- If actual C5 satisfies the formula to high precision, then the original closure follows from R1 plus this C5 relation.")
print("- The next proof target would be to derive R1 and the C5 relation from the Picard-Fuchs/modular equation for B at q=e^{-2*pi}.")
print("Done.")

# ============================================================
# CONDITIONAL FORMAL PROOF CHECK
#
# Run this in the SAME Sage session after the previous obstruction code.
#
# We impose:
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
# Then it checks whether the obstruction polynomials vanish formally.
# ============================================================

print("="*70)
print("CONDITIONAL FORMAL PROOF CHECK")
print("="*70)

try:
    obs_polys
    R
    names
except NameError:
    raise RuntimeError("Run this in the same session as the obstruction extractor.")

# Reduced ring with independent variables L,C0,C2,C3,C4.
proof_names = ["L", "C0", "C2", "C3", "C4"]
Rproof = PolynomialRing(QQ, proof_names)
Kproof = Rproof.fraction_field()

Lp  = Kproof(Rproof.gen(0))
C0p = Kproof(Rproof.gen(1))
C2p = Kproof(Rproof.gen(2))
C3p = Kproof(Rproof.gen(3))
C4p = Kproof(Rproof.gen(4))

# Relation R1:
#   -L C0 + 6 L C1 + 6 C0 = 0
# gives:
C1p = ((Lp - 6) * C0p) / (6 * Lp)

# Relation R5 solved for C5.
#
# R5 =
# -54L^2C2^2 + 12L^2C0C3 + 108L^2C2C3
# -45L^2C0C4 + 54L^2C0C5
# +LC0^2 -72LC0C2 +216LC0C3 -270LC0C4 -8C0^2 = 0.
#
C5p = (
    54*Lp**2*C2p**2
    - 12*Lp**2*C0p*C3p
    - 108*Lp**2*C2p*C3p
    + 45*Lp**2*C0p*C4p
    - Lp*C0p**2
    + 72*Lp*C0p*C2p
    - 216*Lp*C0p*C3p
    + 270*Lp*C0p*C4p
    + 8*C0p**2
) / (54*Lp**2*C0p)

# Map old variables [L,C0,C1,C2,C3,C4,C5] into proof field.
phi_vals = [
    Lp,
    C0p,
    C1p,
    C2p,
    C3p,
    C4p,
    C5p,
]

phi_proof = R.hom(phi_vals, Kproof)

def proof_num(expr):
    expr = Kproof(expr)
    num = Rproof(expr.numerator())
    try:
        return num.primitive_part()
    except Exception:
        cont = num.content()
        if cont != 0:
            return Rproof(num / cont)
        return num

all_zero = True

for idx, P in enumerate(obs_polys):
    print("-"*70)
    print("Checking obstruction P{}".format(idx))

    P_sub = phi_proof(P)
    num = proof_num(P_sub)

    is_zero = bool(num == 0)
    all_zero = all_zero and is_zero

    print("  numerator zero?", is_zero)

    if not is_zero:
        print("  degree:", num.degree())
        print("  term count:", len(num.dict()))
        for nm, gen in zip(proof_names, Rproof.gens()):
            print("  degree in {}: {}".format(nm, num.degree(gen)))

print()
print("="*70)
print("Result")
print("="*70)
print("All obstruction polynomials vanish under R1 and R5?", all_zero)

if all_zero:
    print()
    print("SUCCESS:")
    print("For the tested obstruction polynomials, the closure follows formally")
    print("from the two jet identities R1=0 and R5=0.")
else:
    print()
    print("Not all obstructions vanished. More relations are needed.")

print()
print("Done.")
