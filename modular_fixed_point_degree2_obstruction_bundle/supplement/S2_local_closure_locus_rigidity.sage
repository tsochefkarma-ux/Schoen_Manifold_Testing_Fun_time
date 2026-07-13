from sage.all import *

# ============================================================================
# LOCAL CLOSURE-LOCUS AND FIRST-ORDER RIGIDITY TEST
# ============================================================================
#
# Input from the preceding fixed-point ideal computation:
#
#   N = delta*A + E6*B,
#   delta = L*E2 - 6,
#
# where the first-order coefficients at delta=E6=0 are
#
#   A_delta(X) = (dN/ddelta)|fixed,
#   B_E6(X)    = (dN/dE6)|fixed.
#
# This cell does four things:
#
#   1. Reconstructs A_delta(X) and B_E6(X) exactly.
#   2. Proves first-order FUNCTIONAL rigidity of the closure identity in X.
#      The X^3 and X^1 coefficient equations have a nonzero 2x2 Jacobian.
#   3. Specializes rigorously to L=2*pi and E4=E4(i), using real-ball
#      arithmetic, and proves that A_delta has exactly one positive root.
#   4. Proves that B_E6 never vanishes for X>0 and isolates the unique
#      positive root of A_delta in a certified rational interval.
#
# IMPORTANT INTERPRETATION
# ------------------------
# A_delta does vanish at one positive X.  This does NOT spoil rigidity of the
# closure as an identity in X: B_E6 is nonzero there, and the polynomial
# coefficient map has rank two at the fixed point.
# ============================================================================

BITS = 256
BISECTION_STEPS = 70
PRINT_FULL_POLYNOMIALS = False

print("="*78)
print("LOCAL CLOSURE-LOCUS AND FIRST-ORDER RIGIDITY TEST")
print("Physical point: L=2*pi, E4=E4(i)")
print("="*78)
print()

# ----------------------------------------------------------------------------
# Exact symbolic first-order coefficients
# ----------------------------------------------------------------------------

S = PolynomialRing(QQ, names=("L", "E4", "C0"))
rL, rE4, rC0 = S.gens()
F = S.fraction_field()
L, E4, C0 = [F(g) for g in S.gens()]

PX = PolynomialRing(F, "X")
X = PX.gen()

# The parenthesized cubic printed by 05_fixed_point_defect_membership.sage.
P = PX(
    -504*X**2*L**6*E4**2
    -6912*X**3*L**5*E4
    -336*X*L**6*E4**2
    +13*L**6*E4**3
    -6912*X**2*L**5*E4
    +1512*X*L**5*E4**2
    -56*L**6*E4**2
    +46656*X**2*L**4*E4
    -2304*X*L**5*E4
    +504*L**5*E4**2
    -248832*X**3*L**3
    +31104*X*L**4*E4
    -256*L**5*E4
    +2916*L**4*E4**2
    -248832*X**2*L**3
    +32832*X*L**3*E4
    +5184*L**4*E4
    +1213056*X**2*L**2
    -82944*X*L**3
    +10944*L**3*E4
    +808704*X*L**2
    -9216*L**3
    +88560*L**2*E4
    +342144*X*L
    +134784*L**2
    +114048*L
    +1135296
)

A_delta = PX(3*C0**2*P)
B_E6 = PX(
    -48*C0**2*L**3*(L**2*E4 + 36)
    *(L**2*E4 + 30*X*L + 10*L + 18)
)

print("Exact polynomial degrees:")
print("  deg_X A_delta =", A_delta.degree())
print("  deg_X B_E6    =", B_E6.degree())
assert A_delta.degree() == 3
assert B_E6.degree() == 1
print()

# Exact leading/coefficient checks.
a3_expected = -20736*C0**2*L**3*(L**2*E4 + 36)
b1_expected = -1440*C0**2*L**4*(L**2*E4 + 36)

print("Exact coefficient checks:")
print("  coeff_X^3(A_delta) has expected factorization?",
      bool(A_delta[3] == a3_expected))
print("  coeff_X^1(B_E6) has expected factorization?",
      bool(B_E6[1] == b1_expected))
assert A_delta[3] == a3_expected
assert B_E6[1] == b1_expected
print()

# For the coefficient equations [X^3 N = 0, X^1 N = 0], the first-order
# Jacobian with columns (delta,E6) is triangular:
#
#       [ A_delta[3]      0       ]
#       [ A_delta[1]   B_E6[1]    ]
#
# Its determinant is A_delta[3]*B_E6[1].

jacobian_det = F(A_delta[3]*B_E6[1])
jacobian_expected = F(
    29859840*C0**4*L**7*(L**2*E4 + 36)**2
)

print("First-order functional-rigidity certificate:")
print("  Jacobian determinant has expected factorization?",
      bool(jacobian_det == jacobian_expected))
print("  Jacobian determinant is nonzero in QQ(L,E4,C0)?",
      bool(jacobian_det != 0))
assert jacobian_det == jacobian_expected
assert jacobian_det != 0
print("  determinant =")
print("   ", jacobian_expected)
print()

# A generic gcd check is supplementary.  The rank proof above is simpler and
# does not rely on this check.
generic_gcd = A_delta.gcd(B_E6)
print("Generic gcd over QQ(L,E4,C0)[X]:", generic_gcd)
print()

if PRINT_FULL_POLYNOMIALS:
    print("A_delta(X) =")
    print(A_delta)
    print()
    print("B_E6(X) =")
    print(B_E6)
    print()

# ----------------------------------------------------------------------------
# Rigorous physical specialization using Arb real balls
# ----------------------------------------------------------------------------

RB = RealBallField(BITS)
pi_b = RB.pi()
gamma_quarter_b = RB.gamma(QQ(1)/4)
L_b = 2*pi_b
E4_b = 3*gamma_quarter_b**8/(2*pi_b)**6
Y_b = L_b**2*E4_b       # Y = L^2 E4(i)

print("Rigorous physical constants (real-ball enclosures):")
print("  L       =", L_b)
print("  E4(i)   =", E4_b)
print("  Y=L^2E4 =", Y_b)
print()


def certified_sign(ball):
    """Return +1/-1 if the real ball has certified sign, otherwise 0."""
    ball = RB(ball)
    if ball.lower() > 0:
        return 1
    if ball.upper() < 0:
        return -1
    return 0


# Write P(X)=p3*X^3+p2*X^2+p1*X+p0 using Y=L^2 E4.
# A_delta=3*C0^2*P, so for C0 != 0 the roots and signs are those of P.
p3 = -6912*L_b**3*(Y_b + 36)

p2 = -72*L_b**2*(
    96*L_b*Y_b + 3456*L_b
    + 7*Y_b**2 - 648*Y_b - 16848
)

p1 = -24*L_b*(
    96*L_b**2*Y_b + 3456*L_b**2
    + 14*L_b*Y_b**2 - 1296*L_b*Y_b - 33696*L_b
    - 63*Y_b**2 - 1368*Y_b - 14256
)

p0 = (
    -256*L_b**3*Y_b
    -9216*L_b**3
    -56*L_b**2*Y_b**2
    +5184*L_b**2*Y_b
    +134784*L_b**2
    +504*L_b*Y_b**2
    +10944*L_b*Y_b
    +114048*L_b
    +13*Y_b**3
    +2916*Y_b**2
    +88560*Y_b
    +1135296
)

p_coeffs = [p3, p2, p1, p0]
p_signs = [certified_sign(c) for c in p_coeffs]

print("Certified coefficient signs of P, descending in X:")
print("  [p3,p2,p1,p0] signs =", p_signs)
print("  expected              = [-1, -1, +1, +1]")
assert p_signs == [-1, -1, 1, 1]
print()

# There is exactly one sign change.  By Descartes' rule of signs, P has
# exactly one positive real root, counted with multiplicity.
sign_changes = sum(
    1 for a, b in zip(p_signs, p_signs[1:]) if a != b
)
print("Descartes sign changes:", sign_changes)
assert sign_changes == 1
print("Conclusion: P has exactly one positive real root.")
print()


def P_ball(x_value):
    """Rigorous real-ball evaluation of P at an exact rational x_value."""
    xb = RB(x_value)
    return ((p3*xb + p2)*xb + p1)*xb + p0


# A broad rational bracket with strongly certified opposite signs.
left = QQ(79)/125          # 0.632
right = QQ(633)/1000       # 0.633
s_left = certified_sign(P_ball(left))
s_right = certified_sign(P_ball(right))

print("Initial certified root bracket:")
print("  left  =", left, " sign =", s_left)
print("  right =", right, " sign =", s_right)
assert s_left == 1
assert s_right == -1
print()

# Certified bisection.  Since there is exactly one positive root and the
# endpoint signs are opposite, every retained interval contains that root.
completed_steps = 0
for _ in range(BISECTION_STEPS):
    mid = (left + right)/2
    s_mid = certified_sign(P_ball(mid))

    if s_mid == 0:
        print("Bisection stopped because the midpoint ball contains zero.")
        print("Increase BITS if a narrower interval is needed.")
        break

    if s_mid == s_left:
        left = mid
        s_left = s_mid
    else:
        right = mid
        s_right = s_mid

    completed_steps += 1

assert certified_sign(P_ball(left)) == 1
assert certified_sign(P_ball(right)) == -1

RF = RealField(100)
root_mid = (left + right)/2
root_width = right - left

print("Certified positive-root isolation:")
print("  bisection steps completed =", completed_steps)
print("  rational left endpoint    =", left)
print("  rational right endpoint   =", right)
print("  decimal interval          = [")
print("   ", RF(left))
print("   ", RF(right))
print("  ]")
print("  midpoint approximation    =", RF(root_mid))
print("  interval width            =", RF(root_width))
print()

# ----------------------------------------------------------------------------
# E6-direction nonvanishing on X>0
# ----------------------------------------------------------------------------
#
# B_E6(X) = -48*C0^2*L^3*(Y+36)*(Y+30*L*X+10*L+18).
# For L>0, Y>0, C0 != 0 and X>0, every factor except -48 is positive.
# Thus B_E6(X)<0 throughout X>0.

B_root_b = -(Y_b + 10*L_b + 18)/(30*L_b)
print("E6-linear coefficient analysis:")
print("  unique zero of B_E6(X) =", B_root_b)
print("  zero is certified negative?", bool(B_root_b.upper() < 0))
assert B_root_b.upper() < 0
print("  therefore B_E6(X) is nonzero (indeed negative) for every X>0.")
print()

# Evaluate the B inner factor throughout the isolated A-root interval.
# It is affine increasing in X, so its minimum is at the left endpoint.
B_inner_at_left = Y_b + 30*L_b*RB(left) + 10*L_b + 18
print("At the positive root of A_delta:")
print("  lower bound for B inner factor =", B_inner_at_left.lower())
print("  B inner factor certified positive?", bool(B_inner_at_left.lower() > 0))
assert B_inner_at_left.lower() > 0
print("  hence the E6-linear response remains nonzero there.")
print()

# Physical Jacobian determinant, omitting only the positive C0^4 factor.
jacobian_physical_reduced = (
    29859840*L_b**7*(Y_b + 36)**2
)
print("Physical rank-two Jacobian check:")
print("  determinant / C0^4 =", jacobian_physical_reduced)
print("  certified positive?", bool(jacobian_physical_reduced.lower() > 0))
assert jacobian_physical_reduced.lower() > 0
print()

print("="*78)
print("SUCCESS")
print("1. A_delta has exactly one positive root, rigorously isolated above.")
print("2. B_E6 is nonzero for all X>0, including at that root.")
print("3. The X^3/X^1 coefficient Jacobian has rank two at the fixed point.")
print("Therefore the closure identity as a FUNCTION OF X is infinitesimally")
print("rigid in the two modular directions delta and E6.")
print("Please copy the complete printed output back into the chat.")
print("="*78)
