from sage.all import *

# ============================================================================
# MODULAR q-PATH TRANSVERSALITY TEST AT tau = i
# ============================================================================
#
# This script studies the first derivative of the CLEARED closure-defect
# numerator along the actual imaginary modular path
#
#       tau(y) = i*y,     y = 1 at the square lattice.
#
# The preceding computations established
#
#       N = delta*A + E6*B,
#       delta = L*E2 - 6,
#
# and therefore, at delta=E6=0,
#
#       dN/dy = delta'(1)*A_delta(X) + E6'(1)*B_E6(X).
#
# Here L=2*pi is held fixed, exactly as in the jet algebra used to construct
# N.  The E4 derivative contributes nothing at first order because
#
#       N(X,L,E4,0,0) == 0
#
# identically as a function of E4.
#
# Along tau=i*y, the Ramanujan equations give
#
#       dE2/dy = (pi/6)*(E4 - E2^2),
#       dE6/dy = pi*(E4^2 - E2*E6).
#
# At y=1:
#
#       E2(i)=3/pi,   E6(i)=0,
#
# so, writing Y=L^2*E4(i),
#
#       delta'(1) = (Y-36)/12,
#       E6'(1)    = Y^2/(2*L^3).
#
# The script proves rigorously that the resulting cubic first variation has
# all coefficients strictly negative at the physical square-lattice values.
# Consequently it is nonzero (indeed negative) for every X>0.
#
# It also includes, as a clearly labelled comparison, the derivative obtained
# from the completed combination
#
#       delta_hat(y) = 2*pi*y*E2(i*y) - 6,
#
# whose derivative at y=1 is (Y+36)/12.  This comparison is NOT needed for the
# main fixed-L jet-path theorem, but it checks that the same transversality
# persists in the modular-covariant normalization.
# ============================================================================

BITS = 256
BISECTION_STEPS = 70
PRINT_FULL_POLYNOMIALS = False

print("="*78)
print("MODULAR q-PATH TRANSVERSALITY TEST")
print("Path: tau(y)=i*y through y=1")
print("Main coordinate: delta=L*E2-6 with L=2*pi fixed")
print("="*78)
print()

# ----------------------------------------------------------------------------
# Exact symbolic setup in QQ(L,Y)[X], with Y=L^2*E4(i)
# ----------------------------------------------------------------------------

R = PolynomialRing(QQ, names=("L", "Y"))
rL, rY = R.gens()
K = R.fraction_field()
L = K(rL)
Y = K(rY)

PX = PolynomialRing(K, "X")
X = PX.gen()

# The cubic P from the preceding rigidity calculation, rewritten using
# Y = L^2*E4.  Recall:
#
#       A_delta(X) / C0^2 = 3*P(X).

p3 = -6912*L**3*(Y + 36)

p2 = -72*L**2*(
    96*L*Y
    + 3456*L
    + 7*Y**2
    - 648*Y
    - 16848
)

p1 = -24*L*(
    96*L**2*Y
    + 3456*L**2
    + 14*L*Y**2
    - 1296*L*Y
    - 33696*L
    - 63*Y**2
    - 1368*Y
    - 14256
)

p0 = (
    -256*L**3*Y
    - 9216*L**3
    - 56*L**2*Y**2
    + 5184*L**2*Y
    + 134784*L**2
    + 504*L*Y**2
    + 10944*L*Y
    + 114048*L
    + 13*Y**3
    + 2916*Y**2
    + 88560*Y
    + 1135296
)

P = PX(p3*X**3 + p2*X**2 + p1*X + p0)

Abar = PX(3*P)  # A_delta / C0^2

Bbar = PX(
    -48*L**3*(Y + 36)
    *(Y + 30*L*X + 10*L + 18)
)  # B_E6 / C0^2

# Actual fixed-L jet-path derivatives.
ddelta_fixed = K((Y - 36)/12)
dE6_dy = K(Y**2/(2*L**3))

# Optional modular-completed comparison.
ddelta_completed = K((Y + 36)/12)

T_fixed = PX(ddelta_fixed*Abar + dE6_dy*Bbar)
T_completed = PX(ddelta_completed*Abar + dE6_dy*Bbar)

print("Exact symbolic checks:")
print("  deg_X T_fixed     =", T_fixed.degree())
print("  deg_X T_completed =", T_completed.degree())
print("  T_fixed is nonzero?", bool(T_fixed != 0))
print("  T_completed is nonzero?", bool(T_completed != 0))

assert T_fixed.degree() == 3
assert T_completed.degree() == 3
assert T_fixed != 0
assert T_completed != 0

# Since ddelta_completed - ddelta_fixed = 6, the two comparison cubics differ
# by exactly 6*Abar.
comparison_identity = bool(T_completed - T_fixed == 6*Abar)
print("  T_completed - T_fixed = 6*Abar ?", comparison_identity)
assert comparison_identity

# Exact leading coefficients.
lead_fixed_expected = K(-1728*L**3*(Y - 36)*(Y + 36))
lead_completed_expected = K(-1728*L**3*(Y + 36)**2)

print("  fixed-path leading factor correct?",
      bool(T_fixed[3] == lead_fixed_expected))
print("  completed-path leading factor correct?",
      bool(T_completed[3] == lead_completed_expected))

assert T_fixed[3] == lead_fixed_expected
assert T_completed[3] == lead_completed_expected

if PRINT_FULL_POLYNOMIALS:
    print()
    print("T_fixed(X) / C0^2 =")
    print(T_fixed)
    print()
    print("T_completed(X) / C0^2 =")
    print(T_completed)

print()

# ----------------------------------------------------------------------------
# Rigorous physical specialization
# ----------------------------------------------------------------------------

RB = RealBallField(BITS)
pi_b = RB.pi()
gamma_quarter_b = RB.gamma(QQ(1)/4)

L_b = 2*pi_b
E4_b = 3*gamma_quarter_b**8/(2*pi_b)**6
Y_b = L_b**2*E4_b

print("Rigorous square-lattice constants:")
print("  L=2*pi     =", L_b)
print("  E4(i)      =", E4_b)
print("  Y=L^2*E4   =", Y_b)
print()


def certified_sign(ball):
    """Return +1 or -1 when a real ball has certified sign; otherwise 0."""
    ball = RB(ball)
    if ball.lower() > 0:
        return 1
    if ball.upper() < 0:
        return -1
    return 0


def eval_K_at_physical(value):
    """Evaluate an element of QQ(L,Y) rigorously at the physical point."""
    value = K(value)
    num = value.numerator()
    den = value.denominator()
    return RB(num(L_b, Y_b))/RB(den(L_b, Y_b))


def physical_coefficients(poly):
    """Return descending real-ball coefficients of a polynomial in X."""
    return [
        eval_K_at_physical(poly[k])
        for k in range(poly.degree(), -1, -1)
    ]


def eval_ball_poly(coeffs_desc, x_value):
    """Horner evaluation of descending real-ball coefficients."""
    xb = RB(x_value)
    total = RB(0)
    for coeff in coeffs_desc:
        total = total*xb + coeff
    return total


def cubic_discriminant(coeffs_desc):
    """Discriminant of a cubic a*x^3+b*x^2+c*x+d."""
    a, b, c, d = coeffs_desc
    return (
        18*a*b*c*d
        - 4*b**3*d
        + b**2*c**2
        - 4*a*c**3
        - 27*a**2*d**2
    )


def isolate_root(coeffs_desc, left, right, steps):
    """Certified bisection for an exact rational bracket."""
    left = QQ(left)
    right = QQ(right)

    s_left = certified_sign(eval_ball_poly(coeffs_desc, left))
    s_right = certified_sign(eval_ball_poly(coeffs_desc, right))

    if s_left == 0 or s_right == 0 or s_left == s_right:
        raise ValueError(
            "Initial root bracket is not certified with opposite signs."
        )

    completed = 0

    for _ in range(steps):
        mid = (left + right)/2
        s_mid = certified_sign(eval_ball_poly(coeffs_desc, mid))

        if s_mid == 0:
            break

        if s_mid == s_left:
            left = mid
            s_left = s_mid
        else:
            right = mid
            s_right = s_mid

        completed += 1

    assert certified_sign(eval_ball_poly(coeffs_desc, left)) != 0
    assert certified_sign(eval_ball_poly(coeffs_desc, right)) != 0
    assert (
        certified_sign(eval_ball_poly(coeffs_desc, left))
        != certified_sign(eval_ball_poly(coeffs_desc, right))
    )

    return left, right, completed


# Derivative values, evaluated independently from the compact formulas.
ddelta_fixed_b = (Y_b - 36)/12
ddelta_completed_b = (Y_b + 36)/12
dE6_dy_b = Y_b**2/(2*L_b**3)
dE4_dy_b = -2*E4_b

print("Ramanujan path derivatives at y=1:")
print("  d(delta)/dy, fixed L       =", ddelta_fixed_b)
print("  d(delta_hat)/dy, completed =", ddelta_completed_b)
print("  dE4/dy                     =", dE4_dy_b)
print("  dE6/dy                     =", dE6_dy_b)

assert ddelta_fixed_b.lower() > 0
assert ddelta_completed_b.lower() > 0
assert dE4_dy_b.upper() < 0
assert dE6_dy_b.lower() > 0

print()


def analyze_path(label, poly, root_left, root_right):
    """Run the rigorous sign, discriminant, and root-isolation tests."""
    print("-"*78)
    print(label)
    print("-"*78)

    coeffs = physical_coefficients(poly)
    signs = [certified_sign(c) for c in coeffs]

    print("Certified coefficient signs [X^3,X^2,X,1]:", signs)
    print("Expected signs:                              [-1,-1,-1,-1]")

    assert signs == [-1, -1, -1, -1]

    print("All coefficients are strictly negative.")
    print("Therefore the first variation is strictly negative for every X>=0.")

    discr = cubic_discriminant(coeffs)
    discr_sign = certified_sign(discr)

    print("Cubic discriminant =", discr)
    print("Discriminant certified negative?", bool(discr_sign == -1))

    assert discr_sign == -1

    print("Hence the cubic has exactly one real root.")
    print("Since there is no nonnegative root, that real root is negative.")

    left, right, completed = isolate_root(
        coeffs,
        root_left,
        root_right,
        BISECTION_STEPS
    )

    RF = RealField(100)
    midpoint = (left + right)/2
    width = right - left

    print("Certified negative-root isolation:")
    print("  bisection steps =", completed)
    print("  rational left   =", left)
    print("  rational right  =", right)
    print("  decimal interval = [")
    print("   ", RF(left))
    print("   ", RF(right))
    print("  ]")
    print("  midpoint         =", RF(midpoint))
    print("  interval width   =", RF(width))
    print()

    return {
        "coefficients": coeffs,
        "signs": signs,
        "discriminant": discr,
        "root_left": left,
        "root_right": right,
    }


fixed_result = analyze_path(
    "MAIN TEST: fixed-L jet coordinate delta=L*E2-6",
    T_fixed,
    -QQ(663)/1000,
    -QQ(331)/500
)

completed_result = analyze_path(
    "COMPARISON: completed coordinate delta_hat=2*pi*y*E2-6",
    T_completed,
    -QQ(501)/1000,
    -QQ(1)/2
)

print("="*78)
print("SUCCESS")
print()
print("For the actual fixed-L modular q-path tau=i*y:")
print("  dN/dy at y=1 is a nonzero cubic in X.")
print("  Every coefficient is rigorously negative at the physical point.")
print("  Therefore dN/dy < 0 for every X>0 (after dividing by C0^2>0).")
print("  The square-lattice closure is transverse along this modular path.")
print()
print("The optional completed-coordinate comparison has the same property.")
print()
print("Please copy the complete printed output back into the chat.")
print("="*78)
