from sage.all import *

# ============================================================================
# ROUTE A, STAGES 2--4:
# HSS GAMMA SPECIALIZATION, DEGREE-TWO GATE, AND CM JET SHAPE
# ============================================================================
#
# This cell performs four logically separate checks.
#
# (1) HSS gamma specialization
#
#       q^(3/2) M_1(3t,t*gamma) = 9 prod_{n>=1}(1-q^n)^(-4),
#       gamma=(1,1,1,1,1,1,1,-1).
#
# (2) WDVV stationary M_2 specialization
#
#     It reconstructs the Oberdieck--Pixton stationary series M_2 and
#     specializes q^3 M_2(3t,t*gamma).
#
#     IMPORTANT: M_2 has a point insertion.  It is not automatically the
#     uninserted HSS genus-zero coefficient Z_{0;2}.
#
# (3) Actual rational-elliptic-surface degree-two prepotential coefficient
#
#       Z_{0;2} = (E_2/72) * Z_{0;1}^2,
#
#     together with the multiple-cover-subtracted series
#
#       Ztilde_{0;2}(q) = Z_{0;2}(q) - (1/8) Z_{0;1}(q^2).
#
# (4) Exact square-lattice finite-jet shape of Z_{0;2}
#
#     If C0=B(e^{-2*pi}) and
#
#       j_m = theta^m Z_{0;2}(e^{-2*pi}) / C0,
#
#     then the code constructs j_m/C0 in QQ(L,E4), L=2*pi.
#     Since kernel membership is homogeneous, the vector j/C0 is enough for
#     the linearized correction-kernel test from cell 08.
#
# The final optional adapter tests this exact one-sided rational-surface jet
# shape when the successful cell 08 is still loaded in the same Sage session.
# It does NOT claim that this is already the full Schoen p^2 bi-Jacobi term.
# ============================================================================

OUT_PREC = 4
PRINT_CM_JET_EXPRESSIONS = False

print("="*78, flush=True)
print("HSS GAMMA SPECIALIZATION AND DEGREE-TWO COMPATIBILITY GATE", flush=True)
print("q precision through q^{}".format(OUT_PREC), flush=True)
print("="*78, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Small exact series helpers
# ----------------------------------------------------------------------------


def scalar_add(array, degree, value):
    if 0 <= degree < len(array) and value != 0:
        array[degree] += QQ(value)


def vector_add(array, degree, value):
    if 0 <= degree < len(array):
        array[degree] += vector(QQ, value)


def convolution(left, right, max_degree):
    output = [QQ(0)]*(max_degree + 1)
    for i, a in enumerate(left):
        if a == 0:
            continue
        for j in range(max_degree + 1 - i):
            b = right[j]
            if b != 0:
                output[i+j] += a*b
    return output


def vector_dot_convolution(left, right, max_degree):
    output = [QQ(0)]*(max_degree + 1)
    for i, a in enumerate(left):
        for j in range(max_degree + 1 - i):
            output[i+j] += a.dot_product(right[j])
    return output


# ----------------------------------------------------------------------------
# Safe lattice cutoff
# ----------------------------------------------------------------------------
#
# For a vector lambda with d=||lambda||^2/2, the shifted HSS exponent is
#
#       e = 3d + <lambda,gamma> + 3r,
#
# where r>=0 comes from prod(1-q_tau^n)^(-12).
# Since ||gamma||^2=8,
#
#       <lambda,gamma> >= -4 sqrt(d).
#
# Hence any contribution with e<=N must satisfy
#
#       3d - 4 sqrt(d) <= N.
#
# The following cutoff is therefore conservative.
# ----------------------------------------------------------------------------

safe_root = (QQ(4) + sqrt(QQ(16) + QQ(12)*OUT_PREC))/QQ(6)
THETA_MAX_DEGREE = ZZ(ceil(safe_root**2)) + 1
PARTITION_MAX_DEGREE = OUT_PREC//3 + 1

print("Safe theta-degree cutoff:", THETA_MAX_DEGREE, flush=True)
print("Partition-degree cutoff:", PARTITION_MAX_DEGREE, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# E8 lattice enumeration
# ----------------------------------------------------------------------------


def enumerate_E8_vectors(max_degree):
    max_square = 8*ZZ(max_degree)
    coordinate_bound = floor(sqrt(max_square))
    output = []

    for parity in [0, 1]:
        choices = [
            value for value in range(-coordinate_bound, coordinate_bound + 1)
            if value % 2 == parity
        ]

        current = [0]*8

        def recurse(position, square_sum):
            if position == 8:
                if sum(current) % 4 != 0:
                    return
                if square_sum % 8 != 0:
                    return

                degree = ZZ(square_sum // 8)
                if degree > max_degree:
                    return

                lam = vector(QQ, [QQ(value)/2 for value in current])
                output.append((degree, lam))
                return

            for value in choices:
                next_square = square_sum + value*value
                if next_square <= max_square:
                    current[position] = value
                    recurse(position + 1, next_square)

        recurse(0, 0)

    return output


print("Enumerating E8 lattice vectors...", flush=True)
A12_E8_vectors = enumerate_E8_vectors(THETA_MAX_DEGREE)
print("  vectors found:", len(A12_E8_vectors), flush=True)

counts = [0]*(THETA_MAX_DEGREE + 1)
for degree, lam in A12_E8_vectors:
    counts[degree] += 1

expected_counts = [1] + [
    240*sigma(n, 3)
    for n in range(1, THETA_MAX_DEGREE + 1)
]

theta_ok = bool(counts == expected_counts)
print("  E8 theta-count check valid?", theta_ok, flush=True)
assert theta_ok
print(flush=True)

# ----------------------------------------------------------------------------
# eta denominator coefficients
# ----------------------------------------------------------------------------

PS = PowerSeriesRing(QQ, "qA12", default_prec=OUT_PREC + 1)
qA12 = PS.gen()

P12 = PS(1)
for n in range(1, PARTITION_MAX_DEGREE + 1):
    P12 *= (1 - qA12**n)**(-12)
P12 = P12.add_bigoh(PARTITION_MAX_DEGREE + 1)
P12_coeffs = [QQ(P12[n]) for n in range(PARTITION_MAX_DEGREE + 1)]

# ----------------------------------------------------------------------------
# HSS vectors gamma and the WDVV root b
# ----------------------------------------------------------------------------

gamma = vector(QQ, [1, 1, 1, 1, 1, 1, 1, -1])
b_root = vector(QQ, [1, -1, 0, 0, 0, 0, 0, 0])

assert gamma.dot_product(gamma) == 8
assert b_root.dot_product(b_root) == 2

ZERO8 = vector(QQ, [0]*8)

# Shifted one-variable series.  Each array index is the final q exponent after
# multiplying M_1 by q^(3/2).
A12_f = [QQ(0)]*(OUT_PREC + 1)
A12_fb = [QQ(0)]*(OUT_PREC + 1)
A12_fbb = [QQ(0)]*(OUT_PREC + 1)
A12_Dqf = [QQ(0)]*(OUT_PREC + 1)
A12_Dqfb = [QQ(0)]*(OUT_PREC + 1)
A12_Dqfbb = [QQ(0)]*(OUT_PREC + 1)

A12_grad = [vector(QQ, [0]*8) for _ in range(OUT_PREC + 1)]
A12_gradb = [vector(QQ, [0]*8) for _ in range(OUT_PREC + 1)]
A12_gradbb = [vector(QQ, [0]*8) for _ in range(OUT_PREC + 1)]

print("Building the gamma-specialized derivative moments...", flush=True)

for theta_degree, lam in A12_E8_vectors:
    gamma_power = lam.dot_product(gamma)
    b_power = lam.dot_product(b_root)

    if gamma_power.denominator() != 1 or b_power.denominator() != 1:
        raise ArithmeticError("Expected integral E8 pairings with gamma and b.")

    gamma_power = ZZ(gamma_power)
    b_power = ZZ(b_power)

    for partition_degree, partition_coefficient in enumerate(P12_coeffs):
        exponent = 3*(theta_degree + partition_degree) + gamma_power

        if exponent < 0 or exponent > OUT_PREC:
            continue

        q_tau_weight = QQ(theta_degree + partition_degree) - QQ(1)/2
        coefficient = QQ(partition_coefficient)

        scalar_add(A12_f, exponent, coefficient)
        scalar_add(A12_fb, exponent, coefficient*b_power)
        scalar_add(A12_fbb, exponent, coefficient*b_power**2)

        scalar_add(A12_Dqf, exponent, coefficient*q_tau_weight)
        scalar_add(A12_Dqfb, exponent, coefficient*q_tau_weight*b_power)
        scalar_add(A12_Dqfbb, exponent, coefficient*q_tau_weight*b_power**2)

        vector_add(A12_grad, exponent, coefficient*lam)
        vector_add(A12_gradb, exponent, coefficient*b_power*lam)
        vector_add(A12_gradbb, exponent, coefficient*b_power**2*lam)

print("  derivative moments built.", flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Stage 1: exact HSS embedding check
# ----------------------------------------------------------------------------

P4 = PS(1)
for n in range(1, OUT_PREC + 1):
    P4 *= (1 - qA12**n)**(-4)
P4 = P4.add_bigoh(OUT_PREC + 1)

B_HSS_series = PS(9)*P4
B_HSS_coeffs = [QQ(B_HSS_series[n]) for n in range(OUT_PREC + 1)]

embedding_ok = bool(A12_f == B_HSS_coeffs)

print("Stage 1: HSS gamma embedding", flush=True)
print("  q^(3/2) M1(3t,t*gamma):", A12_f, flush=True)
print("  9*prod(1-q^n)^(-4):   ", B_HSS_coeffs, flush=True)
print("  exact coefficient match?", embedding_ok, flush=True)
assert embedding_ok
print(flush=True)

# ----------------------------------------------------------------------------
# Stage 2: stationary WDVV M2 under the same specialization
# ----------------------------------------------------------------------------

H1 = convolution(A12_Dqfb, A12_fb, OUT_PREC)
H2 = convolution(A12_Dqf, A12_fbb, OUT_PREC)
H3 = convolution(A12_f, A12_Dqfbb, OUT_PREC)

E1 = vector_dot_convolution(A12_gradb, A12_gradb, OUT_PREC)
E2moment = vector_dot_convolution(A12_grad, A12_gradbb, OUT_PREC)

A12_stationary_M2 = [
    -QQ(1)/8 * (
        2*H1[n] - H2[n] - H3[n] - E1[n] + E2moment[n]
    )
    for n in range(OUT_PREC + 1)
]

print("Stage 2: specialized stationary WDVV M2", flush=True)
print("  q^3 M2(3t,t*gamma):", A12_stationary_M2, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Stage 3: actual uninserted degree-two genus-zero coefficient
# ----------------------------------------------------------------------------
#
# HST give
#
#       Z_{0;2} = (E2/72) * B^2.
#
# This contains multiple covers.  For n=2,
#
#       Ztilde_{0;2} = Z_{0;2} - (1/8) B(q^2).
# ----------------------------------------------------------------------------

E2_series = PS(1)
for n in range(1, OUT_PREC + 1):
    E2_series += -24*sigma(n, 1)*qA12**n
E2_series = E2_series.add_bigoh(OUT_PREC + 1)

Z02_series = (E2_series * B_HSS_series**2 / 72).add_bigoh(OUT_PREC + 1)
Z02_coeffs = [QQ(Z02_series[n]) for n in range(OUT_PREC + 1)]

B_q2_series = PS(0)
for n in range(OUT_PREC + 1):
    if 2*n <= OUT_PREC:
        B_q2_series += B_HSS_coeffs[n]*qA12**(2*n)
B_q2_series = B_q2_series.add_bigoh(OUT_PREC + 1)

Z02_primitive_series = (Z02_series - B_q2_series/8).add_bigoh(OUT_PREC + 1)
Z02_primitive_coeffs = [
    QQ(Z02_primitive_series[n])
    for n in range(OUT_PREC + 1)
]

stationary_equals_uninserted = bool(A12_stationary_M2 == Z02_coeffs)
stationary_equals_primitive = bool(A12_stationary_M2 == Z02_primitive_coeffs)

print("Stage 3: uninserted rational-surface degree-two series", flush=True)
print("  Z_0;2=(E2/72)B^2:          ", Z02_coeffs, flush=True)
print("  Ztilde_0;2=Z_0;2-B(q^2)/8:", Z02_primitive_coeffs, flush=True)
print("  stationary M2 equals Z_0;2?", stationary_equals_uninserted, flush=True)
print("  stationary M2 equals Ztilde?", stationary_equals_primitive, flush=True)
print(flush=True)

if stationary_equals_uninserted or stationary_equals_primitive:
    raise AssertionError(
        "Unexpected identification: re-check conventions before proceeding."
    )

# ----------------------------------------------------------------------------
# Stage 4: exact CM jet shape of Z_{0;2}
# ----------------------------------------------------------------------------

A12_RJ = PolynomialRing(QQ, ["A12_L", "A12_E2", "A12_E4", "A12_E6"])
A12_KJ = A12_RJ.fraction_field()
A12_L, A12_E2, A12_E4, A12_E6 = A12_KJ.gens()


def A12_dK(expression, variable):
    expression = A12_KJ(expression)
    numerator = A12_RJ(expression.numerator())
    denominator = A12_RJ(expression.denominator())

    return A12_KJ(
        numerator.derivative(variable)*denominator
        - numerator*denominator.derivative(variable)
    ) / A12_KJ(denominator)**2


def A12_theta(expression):
    return (
        A12_dK(expression, A12_RJ.gen(1))
        * (A12_E2**2 - A12_E4)/12
        + A12_dK(expression, A12_RJ.gen(2))
        * (A12_E2*A12_E4 - A12_E6)/3
        + A12_dK(expression, A12_RJ.gen(3))
        * (A12_E2*A12_E6 - A12_E4**2)/2
    )


# theta log B and theta log Z02.
A12_logB = (1 - A12_E2)/6
A12_logZ02 = A12_theta(A12_E2)/A12_E2 + 2*A12_logB

# v_m = theta^m Z02 / Z02.
A12_v = [A12_KJ(1)]
for m in range(5):
    A12_v.append(A12_KJ(A12_theta(A12_v[-1]) + A12_logZ02*A12_v[-1]))


def A12_cm_sub(expression):
    expression = A12_KJ(expression)
    hom = A12_RJ.hom(
        [A12_L, 6/A12_L, A12_E4, 0],
        A12_KJ,
    )

    numerator = A12_RJ(expression.numerator())
    denominator = A12_RJ(expression.denominator())

    return A12_KJ(hom(numerator))/A12_KJ(hom(denominator))


A12_v_cm = [A12_cm_sub(value) for value in A12_v]

# At tau=i,
#
#   Z02/B = (E2/72)B = C0/(12L).
#
# Therefore
#
#   j_m/C0 = v_m/(12L).
A12_degree2_shape = [
    A12_KJ(value/(12*A12_L))
    for value in A12_v_cm
]

print("Stage 4: exact CM finite-jet shape", flush=True)
print("  constructed six values (j_m/C0), m=0,...,5.", flush=True)

if PRINT_CM_JET_EXPRESSIONS:
    for m, value in enumerate(A12_degree2_shape):
        print("  (j_{}/C0) = {}".format(m, value), flush=True)
else:
    print("  exact expressions retained in A12_degree2_shape.", flush=True)

# Numerical values for compact inspection.
A12_RB = RealBallField(160)
A12_L_num = 2*A12_RB.pi()
A12_E4_num = 3*A12_RB.gamma(QQ(1)/4)**8/(2*A12_RB.pi())**6

A12_eval_ring = PolynomialRing(QQ, ["A12_L_eval", "A12_E4_eval"])
A12_eval_field = A12_eval_ring.fraction_field()
A12_L_eval, A12_E4_eval = A12_eval_field.gens()


def A12_to_two_variable(expression):
    expression = A12_KJ(expression)

    # After CM substitution the expression depends only on A12_L and A12_E4.
    source_num = A12_RJ(expression.numerator())
    source_den = A12_RJ(expression.denominator())

    map_to_two = A12_RJ.hom(
        [A12_L_eval, 0, A12_E4_eval, 0],
        A12_eval_field,
    )

    return A12_eval_field(map_to_two(source_num))/A12_eval_field(map_to_two(source_den))


def A12_numeric(expression):
    two = A12_to_two_variable(expression)
    num = two.numerator()(A12_L_num, A12_E4_num)
    den = two.denominator()(A12_L_num, A12_E4_num)
    return A12_RB(num)/A12_RB(den)


A12_shape_numeric = [A12_numeric(value) for value in A12_degree2_shape]
print("  numerical (j_m/C0) values:", flush=True)
for m, value in enumerate(A12_shape_numeric):
    print("    m={} : {}".format(m, value), flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Optional exact adapter to the cell-08 correction kernel
# ----------------------------------------------------------------------------

A12_kernel_test_available = all(
    name in globals()
    for name in ["BaseR", "BaseF", "M"]
)

A12_kernel_preserves = None
A12_kernel_obstruction = None

print("Optional cell-08 adapter:", flush=True)
print("  correction-kernel data available?", A12_kernel_test_available, flush=True)

if A12_kernel_test_available:
    base_generators = list(BaseR.gens())

    if len(base_generators) != 2:
        raise RuntimeError("Expected BaseR=QQ(L,E4) from cell 08.")

    base_L, base_E4 = [BaseF(generator) for generator in base_generators]

    # Map the CM expressions to BaseF by substituting the two surviving names.
    adapter_ring = PolynomialRing(QQ, ["adapter_L", "adapter_E4"])
    adapter_field = adapter_ring.fraction_field()
    adapter_L, adapter_E4 = adapter_field.gens()

    mapped_shape = []

    for expression in A12_degree2_shape:
        two = A12_to_two_variable(expression)

        # Convert through strings only between two rational-function fields
        # with the same ordered generators; this avoids global-name clashes.
        two_num = two.numerator()
        two_den = two.denominator()

        phi_adapter = A12_eval_ring.hom(
            [adapter_L, adapter_E4],
            adapter_field,
        )

        adapted = adapter_field(phi_adapter(two_num))/adapter_field(phi_adapter(two_den))

        phi_base = adapter_ring.hom(
            [base_L, base_E4],
            BaseF,
        )

        mapped = BaseF(phi_base(adapted.numerator()))/BaseF(phi_base(adapted.denominator()))
        mapped_shape.append(mapped)

    A12_kernel_obstruction = M*vector(BaseF, mapped_shape)
    A12_kernel_preserves = bool(A12_kernel_obstruction == 0)

    print("  one-sided Z_0;2 jet shape preserves closure?", A12_kernel_preserves, flush=True)
    print("  nonzero obstruction entries:", sum(1 for x in A12_kernel_obstruction if x != 0), flush=True)
else:
    print("  Run this cell after cell 08 to obtain the exact kernel verdict.", flush=True)

print(flush=True)

# Keep compact reusable data.
ROUTE_A_STAGE2_DATA = {
    "OUT_PREC": OUT_PREC,
    "gamma": gamma,
    "B_coefficients": B_HSS_coeffs,
    "stationary_M2_coefficients": A12_stationary_M2,
    "Z02_coefficients": Z02_coeffs,
    "Z02_primitive_coefficients": Z02_primitive_coeffs,
    "degree2_shape_j_over_C0": A12_degree2_shape,
    "degree2_shape_numeric": A12_shape_numeric,
    "kernel_test_available": A12_kernel_test_available,
    "kernel_preserves": A12_kernel_preserves,
    "kernel_obstruction": A12_kernel_obstruction,
}

print("="*78, flush=True)
print("SUCCESS", flush=True)
print("1. The HSS gamma specialization reproduces B(q) exactly.", flush=True)
print("2. The stationary WDVV M2 is not the uninserted Z_0;2 series.", flush=True)
print("3. The actual rational-surface Z_0;2 and its primitive subtraction were built.", flush=True)
print("4. Six exact square-lattice degree-two jet-shape values were constructed.", flush=True)
print(flush=True)
print("Important: this is a one-sided rational-elliptic-surface correction shape.", flush=True)
print("The full Schoen p^2 term still requires the E8 x E8 bi-Jacobi/gluing lift.", flush=True)
print("Copy the complete compact output back into the chat.", flush=True)
print("="*78, flush=True)
