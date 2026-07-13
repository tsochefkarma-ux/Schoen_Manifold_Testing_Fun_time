from sage.all import *

# ============================================================================
# ROUTE A, STAGE 1:
# FINITE q-EXPANSION OF THE OBERDIECK--PIXTON WDVV SERIES M_2
# ============================================================================
#
# Published input:
#
#   M_1(tau,z) = Theta_E8(tau,z) / Delta(tau)^(1/2),
#
# and, for a root b with <b,b> = -2 in E8(-1),
#
#   4 <b,b> M_2
#     = sum_{a,c} eta^{ac}
#       [ D_b D_a M_1 * D_b D_c M_1
#         - D_a M_1 * D_b^2 D_c M_1 ].
#
# We use the divisor equation and the orthogonal splitting
#
#   H^2(R) = <W,F>  (+)  E8(-1),
#   W^2=F^2=0, W.F=1.
#
# On base degree 1:
#
#   D_F M_1 = M_1,
#   D_W M_1 = D_q M_1.
#
# After restricting the elliptic variable to z = s*b, with
# zeta = exp(2*pi*i*s), the WDVV right-hand side becomes
#
#   H = 2 (D_b D_q M_1)(D_b M_1)
#       - (D_q M_1)(D_b^2 M_1)
#       - M_1(D_b^2 D_q M_1),
#
#   E = - <D_b grad M_1, D_b grad M_1>_Euclidean
#       + <grad M_1, D_b^2 grad M_1>_Euclidean,
#
#   M_2 = (H+E)/(-8).
#
# The code computes exact truncated q-series with Laurent-polynomial
# dependence on zeta.  It does NOT yet assert that z=0 or zeta=-1 is the
# HSS restriction that produced B(q)=9*prod(1-q^n)^(-4).  A compatibility
# gate at the end checks these two simplest specializations and refuses to
# identify them when their leading series do not match.
#
# Runtime:
#   QPREC=3 is intentionally modest and should be quick on CoCalc.
#   Increase to 4 or 5 only after the first successful run.
# ============================================================================

QPREC = 3
PRINT_LAURENT_LEVELS = True

print("="*78, flush=True)
print("ROUTE A: OBERDIECK--PIXTON WDVV RECONSTRUCTION OF M_2", flush=True)
print("q precision through q^{} after removing the leading q powers".format(QPREC), flush=True)
print("="*78, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# E8 lattice enumeration
# ----------------------------------------------------------------------------
# We use the standard model
#
#   E8 = {x in Z^8 : sum x_i even}
#        union
#        {x in (Z+1/2)^8 : sum x_i even}.
#
# Write y=2x.  Then y is either all even or all odd, sum(y)=0 mod 4,
# and
#
#   q exponent = ||x||^2/2 = sum(y_i^2)/8.
#
# The chosen root is b=(1,-1,0,...,0), so
#
#   <x,b> = x_1-x_2 = (y_1-y_2)/2.
# ----------------------------------------------------------------------------


def enumerate_E8_vectors(max_q_degree):
    max_sq = 8*ZZ(max_q_degree)
    bound = floor(sqrt(max_sq))
    output = []

    for parity in [0, 1]:
        choices = [
            y for y in range(-bound, bound + 1)
            if y % 2 == parity
        ]

        current = [0]*8

        def recurse(position, square_sum):
            if position == 8:
                if sum(current) % 4 != 0:
                    return
                if square_sum % 8 != 0:
                    return

                q_degree = ZZ(square_sum // 8)
                if q_degree > max_q_degree:
                    return

                y_tuple = tuple(ZZ(v) for v in current)
                lam = vector(QQ, [QQ(v)/2 for v in y_tuple])
                elliptic_power = ZZ((y_tuple[0] - y_tuple[1]) // 2)

                output.append((q_degree, elliptic_power, lam))
                return

            for value in choices:
                next_square = square_sum + value*value
                if next_square <= max_sq:
                    current[position] = value
                    recurse(position + 1, next_square)

        recurse(0, 0)

    return output


print("Enumerating E8 lattice vectors...", flush=True)
E8_vectors = enumerate_E8_vectors(QPREC)
print("  vectors found:", len(E8_vectors), flush=True)

counts_by_degree = [0]*(QPREC + 1)
for degree, elliptic_power, lam in E8_vectors:
    counts_by_degree[degree] += 1

expected_theta_counts = [1] + [240*sigma(n, 3) for n in range(1, QPREC + 1)]

print("  theta coefficients at z=0:", counts_by_degree, flush=True)
print("  expected E4 coefficients:  ", expected_theta_counts, flush=True)
print("  E8 theta enumeration valid?", counts_by_degree == expected_theta_counts, flush=True)

assert counts_by_degree == expected_theta_counts
print(flush=True)

# ----------------------------------------------------------------------------
# Reciprocal Delta^(1/2)
# ----------------------------------------------------------------------------
#
#   Delta^(1/2) = q^(1/2) prod_{n>=1}(1-q^n)^12,
#
# so
#
#   M_1 = q^(-1/2) * Theta_E8 * prod(1-q^n)^(-12).
#
# We store only the integer-power coefficient series after q^(-1/2).
# ----------------------------------------------------------------------------

PS = PowerSeriesRing(QQ, "q", default_prec=QPREC + 1)
q = PS.gen()

P12 = PS(1)
for n in range(1, QPREC + 1):
    P12 *= (1 - q**n)**(-12)
P12 = P12.add_bigoh(QPREC + 1)
P12_coeffs = [QQ(P12[n]) for n in range(QPREC + 1)]

print("prod(1-q^n)^(-12) coefficients:", P12_coeffs, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Sparse q/zeta series helpers
# ----------------------------------------------------------------------------
# A scalar series is a dictionary
#
#   (n,m) -> coefficient
#
# representing
#
#   q^leading_power * sum coefficient*q^n*zeta^m.
#
# A vector series has the same keys and QQ^8 vector values.
# ----------------------------------------------------------------------------

ZERO8 = vector(QQ, [0]*8)


def scalar_add_term(series, key, value):
    value = QQ(value)
    if value == 0:
        return
    new_value = series.get(key, QQ(0)) + value
    if new_value == 0:
        series.pop(key, None)
    else:
        series[key] = new_value


def vector_add_term(series, key, value):
    value = vector(QQ, value)
    if value == ZERO8:
        return
    new_value = series.get(key, ZERO8) + value
    if new_value == ZERO8:
        series.pop(key, None)
    else:
        series[key] = new_value


def scalar_scale(series, factor):
    factor = QQ(factor)
    if factor == 0:
        return {}
    return {
        key: factor*value
        for key, value in series.items()
        if factor*value != 0
    }


def scalar_sum(*series_list):
    output = {}
    for series in series_list:
        for key, value in series.items():
            scalar_add_term(output, key, value)
    return output


def scalar_product(left, right, max_degree):
    output = {}
    for (n1, m1), c1 in left.items():
        for (n2, m2), c2 in right.items():
            n = n1 + n2
            if n <= max_degree:
                scalar_add_term(output, (n, m1 + m2), c1*c2)
    return output


def vector_dot_product(left, right, max_degree):
    output = {}
    for (n1, m1), v1 in left.items():
        for (n2, m2), v2 in right.items():
            n = n1 + n2
            if n <= max_degree:
                scalar_add_term(output, (n, m1 + m2), v1.dot_product(v2))
    return output


def Dq_M1_part(series):
    """
    Apply D_q=q*d/dq to a series carrying the common leading q^(-1/2).
    A coefficient at stored degree n is multiplied by n-1/2.
    """
    output = {}
    for (n, m), value in series.items():
        scalar_add_term(output, (n, m), (QQ(n) - QQ(1)/2)*value)
    return output


def evaluate_zeta(series, zeta_value, max_degree):
    """Evaluate zeta at +1 or -1 and return q coefficients."""
    if zeta_value not in [1, -1]:
        raise ValueError("This compact evaluator currently supports zeta=+1 or -1.")

    coefficients = [QQ(0)]*(max_degree + 1)
    for (n, m), value in series.items():
        sign = 1 if (zeta_value == 1 or m % 2 == 0) else -1
        coefficients[n] += sign*value
    return coefficients


def is_zeta_symmetric(series):
    keys = set(series.keys())
    for n, m in keys:
        if series.get((n, m), 0) != series.get((n, -m), 0):
            return False
    return True


def level_dictionary(series, degree):
    return {
        m: series[(degree, m)]
        for n, m in sorted(series)
        if n == degree and series[(degree, m)] != 0
    }

# ----------------------------------------------------------------------------
# Build M_1 and the elliptic moment series needed by WDVV
# ----------------------------------------------------------------------------
# Stored scalar series:
#   f       : M_1
#   f_b     : D_b M_1
#   f_bb    : D_b^2 M_1
#
# Stored vector series:
#   grad    : grad_z M_1
#   grad_b  : D_b grad_z M_1
#   grad_bb : D_b^2 grad_z M_1
#
# All carry the common leading factor q^(-1/2).
# ----------------------------------------------------------------------------

f = {}
f_b = {}
f_bb = {}
grad = {}
grad_b = {}
grad_bb = {}

print("Building M_1 and derivative moments...", flush=True)

for theta_degree, elliptic_power, lam in E8_vectors:
    for partition_degree, partition_coefficient in enumerate(P12_coeffs):
        total_degree = theta_degree + partition_degree
        if total_degree > QPREC:
            continue

        key = (total_degree, elliptic_power)
        c = partition_coefficient
        m = QQ(elliptic_power)

        scalar_add_term(f, key, c)
        scalar_add_term(f_b, key, c*m)
        scalar_add_term(f_bb, key, c*m**2)

        vector_add_term(grad, key, c*lam)
        vector_add_term(grad_b, key, c*m*lam)
        vector_add_term(grad_bb, key, c*m**2*lam)

print("  scalar M_1 terms:", len(f), flush=True)
print("  vector gradient terms:", len(grad), flush=True)
print("  M_1 zeta symmetry?", is_zeta_symmetric(f), flush=True)
assert is_zeta_symmetric(f)
print(flush=True)

# ----------------------------------------------------------------------------
# WDVV reconstruction of M_2
# ----------------------------------------------------------------------------

print("Constructing the WDVV right-hand side...", flush=True)

Dq_f = Dq_M1_part(f)
Dq_f_b = Dq_M1_part(f_b)
Dq_f_bb = Dq_M1_part(f_bb)

H_term = scalar_sum(
    scalar_scale(scalar_product(Dq_f_b, f_b, QPREC), 2),
    scalar_scale(scalar_product(Dq_f, f_bb, QPREC), -1),
    scalar_scale(scalar_product(f, Dq_f_bb, QPREC), -1),
)

E_term = scalar_sum(
    scalar_scale(vector_dot_product(grad_b, grad_b, QPREC), -1),
    vector_dot_product(grad, grad_bb, QPREC),
)

WDVV_rhs = scalar_sum(H_term, E_term)

# Since 4<b,b> = 4*(-2) = -8.
M2 = scalar_scale(WDVV_rhs, QQ(-1)/8)

print("  H-plane contribution terms:", len(H_term), flush=True)
print("  E8 contribution terms:     ", len(E_term), flush=True)
print("  reconstructed M_2 terms:   ", len(M2), flush=True)
print("  M_2 zeta symmetry?", is_zeta_symmetric(M2), flush=True)
assert is_zeta_symmetric(M2)
print(flush=True)

# ----------------------------------------------------------------------------
# Compact output
# ----------------------------------------------------------------------------
# M_1 has leading q^(-1/2), so the printed coefficients are for q^(1/2) M_1.
# M_2 has leading q^(-1), so the printed coefficients are for q M_2.
# ----------------------------------------------------------------------------

M1_at_1 = evaluate_zeta(f, 1, QPREC)
M1_at_minus1 = evaluate_zeta(f, -1, QPREC)
M2_at_1 = evaluate_zeta(M2, 1, QPREC)
M2_at_minus1 = evaluate_zeta(M2, -1, QPREC)

print("Specialized coefficient lists:", flush=True)
print("  q^(1/2) M1 | zeta=+1:", M1_at_1, flush=True)
print("  q^(1/2) M1 | zeta=-1:", M1_at_minus1, flush=True)
print("  q M2       | zeta=+1:", M2_at_1, flush=True)
print("  q M2       | zeta=-1:", M2_at_minus1, flush=True)
print(flush=True)

if PRINT_LAURENT_LEVELS:
    print("Low q-level Laurent data for q*M2:", flush=True)
    for degree in range(QPREC + 1):
        print("  q^{} : {}".format(degree, level_dictionary(M2, degree)), flush=True)
    print(flush=True)

# ----------------------------------------------------------------------------
# Compatibility gate with the current leading HSS series
# ----------------------------------------------------------------------------
# The current curvature note uses
#
#   B_HSS(q) = 9*prod(1-q^n)^(-4)
#            = 9 + 36q + 126q^2 + 360q^3 + ...
#
# If a proposed specialization of M_1 is to be the SAME one-variable series,
# then after one overall normalization its coefficients must match B_HSS.
# We test only the two simplest specializations zeta=+1 and zeta=-1.
# ----------------------------------------------------------------------------

P4 = PS(1)
for n in range(1, QPREC + 1):
    P4 *= (1 - q**n)**(-4)
P4 = P4.add_bigoh(QPREC + 1)
B_HSS = [QQ(9*P4[n]) for n in range(QPREC + 1)]


def normalized_to_constant_nine(coefficients):
    if coefficients[0] == 0:
        return None
    scale = QQ(9)/coefficients[0]
    return [scale*c for c in coefficients]


candidate_plus = normalized_to_constant_nine(M1_at_1)
candidate_minus = normalized_to_constant_nine(M1_at_minus1)

plus_matches = bool(candidate_plus == B_HSS)
minus_matches = bool(candidate_minus == B_HSS)

print("Compatibility gate against B_HSS=9*prod(1-q^n)^(-4):", flush=True)
print("  target coefficients:             ", B_HSS, flush=True)
print("  normalized M1 at zeta=+1:        ", candidate_plus, flush=True)
print("  normalized M1 at zeta=-1:        ", candidate_minus, flush=True)
print("  zeta=+1 matches current B_HSS?   ", plus_matches, flush=True)
print("  zeta=-1 matches current B_HSS?   ", minus_matches, flush=True)
print(flush=True)

# Keep useful objects available in the Sage session.
ROUTE_A_DATA = {
    "QPREC": QPREC,
    "E8_vectors": E8_vectors,
    "P12_coeffs": P12_coeffs,
    "M1": f,
    "M2": M2,
    "M1_zeta_plus": M1_at_1,
    "M1_zeta_minus": M1_at_minus1,
    "M2_zeta_plus": M2_at_1,
    "M2_zeta_minus": M2_at_minus1,
    "B_HSS": B_HSS,
    "simple_specialization_matches": {
        "+1": plus_matches,
        "-1": minus_matches,
    },
}

print("="*78, flush=True)
print("SUCCESS", flush=True)
print("The finite-q WDVV reconstruction of M_2 completed exactly.", flush=True)

if plus_matches or minus_matches:
    print("A simple zeta specialization matches the current leading B_HSS series.", flush=True)
    print("The corresponding M_2 series is a candidate correction for the next adapter.", flush=True)
else:
    print("Neither zeta=+1 nor zeta=-1 reproduces the current leading B_HSS series.", flush=True)
    print("Therefore the exact HSS E8-variable embedding must be supplied or derived", flush=True)
    print("before M_2 can be converted into the six correction jets j_0,...,j_5.", flush=True)

print("Copy the complete compact output back into the chat.", flush=True)
print("="*78, flush=True)
