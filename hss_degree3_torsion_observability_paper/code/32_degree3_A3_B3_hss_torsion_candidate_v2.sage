from sage.all import *
import os
import sys

# =============================================================================
# CELL 32: DEGREE-THREE A3/B3 HSS TORSION-CANDIDATE SERIES
# =============================================================================
#
# Purpose
# -------
# Construct the six weight-16, index-3 one-sided holomorphic basis series
# along the canonical order-three shift candidate
#
#     tau = 3 t,
#     z_r = t*gamma + (r/3)*gamma,       r = 0,1,2,
#     gamma = (1,1,1,1,1,1,1,-1).
#
# The character tags are
#
#     one    <-> r=0,
#     omega  <-> r=1,
#     omega2 <-> r=2.
#
# IMPORTANT SCOPE
# ---------------
# This is a mathematically canonical gamma/3 torsion-point candidate.  The
# identification of these three evaluations with the geometric BKOS torsion
# labels must still be checked against torsion recombination and independent
# low-degree data.  The cell does not assume that identification as a theorem.
#
# Sakai generators used
# ---------------------
#
#   A1 = Theta_E8,
#
#   A2 = 8/9 [ A1(2tau,2z)
#              + 1/16 sum_{k=0}^1 A1((tau+k)/2,z) ],
#
#   A3 = 27/28 A1(3tau,3z)
#              + 1/81 sum_{k=0}^2 A1((tau+k)/3,z),
#
#   B2 = 32/5 [ e1 A1(2tau,2z)
#              + e3/16 A1(tau/2,z)
#              + e2/16 A1((tau+1)/2,z) ],
#
#   B3 = 81/80 [ h0(tau)^2 A1(3tau,3z)
#              - 1/3^5 sum_{k=0}^2 h0((tau+k)/3)^2
#                                      A1((tau+k)/3,z) ].
#
# It then constructs
#
#   E6^2 A3, E4^3 A3, E4 E6 B3,
#   E4^2 A1 A2, E6 A1 B2, E4 A1^3,
#
# and computes the exact q-coefficient rank of the two shared sectors
# (omega,omega) and (omega,omega^2).
#
# Coefficients are exact in Q(zeta_12); no floating-point arithmetic is used.
# =============================================================================

print("="*79, flush=True)
print("DEGREE-THREE A3/B3 HSS TORSION-CANDIDATE SERIES", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# Sector q cutoff.  A command-line integer overrides the default.
SECTOR_Q_MAX = ZZ(sys.argv[1]) if len(sys.argv) > 1 else ZZ(24)
if SECTOR_Q_MAX < 8:
    raise ValueError("Use a sector q cutoff of at least 8")

# Compute one-sided series farther than the final sector cutoff because the
# Laurent generators can begin at negative q powers.
GEN_MIN = ZZ(-8)
GEN_MAX = ZZ(SECTOR_Q_MAX + 16)
CHECK_MAX = min(ZZ(20), GEN_MAX)

K12 = CyclotomicField(12)
zeta12 = K12.gen()
omega = zeta12**4
omega2 = omega**2

CHARACTERS = ["one", "omega", "omega2"]
CHARACTER_R = {"one": ZZ(0), "omega": ZZ(1), "omega2": ZZ(2)}

BASIS_LABELS = [
    "E6^2*A3",
    "E4^3*A3",
    "E4*E6*B3",
    "E4^2*A1*A2",
    "E6*A1*B2",
    "E4*A1^3",
]

print("sector q cutoff: q^{}".format(SECTOR_Q_MAX), flush=True)
print("one-sided working range: q^{} through q^{}".format(GEN_MIN, GEN_MAX), flush=True)
print("coefficient field: {}".format(K12), flush=True)
print("candidate shift: z_r=t*gamma+(r/3)*gamma", flush=True)

# -----------------------------------------------------------------------------
# Sparse exact Laurent-series helpers.
# -----------------------------------------------------------------------------

def add_coeff(target, exponent, coefficient):
    exponent = ZZ(exponent)
    coefficient = K12(coefficient)
    if coefficient == 0:
        return
    target[exponent] = target.get(exponent, K12(0)) + coefficient
    if target[exponent] == 0:
        del target[exponent]


def series_scale(series, scalar):
    scalar = K12(scalar)
    if scalar == 0:
        return {}
    return {e: scalar*c for e, c in series.items() if scalar*c != 0}


def series_add(*series_list):
    output = {}
    for series in series_list:
        for exponent, coefficient in series.items():
            add_coeff(output, exponent, coefficient)
    return output


def series_convolve(left, right, min_exp, max_exp):
    min_exp = ZZ(min_exp)
    max_exp = ZZ(max_exp)
    output = {}
    for e1, c1 in left.items():
        if c1 == 0:
            continue
        for e2, c2 in right.items():
            if c2 == 0:
                continue
            exponent = e1 + e2
            if min_exp <= exponent <= max_exp:
                add_coeff(output, exponent, c1*c2)
    return output


def series_power(series, power, min_exp, max_exp):
    output = {ZZ(0): K12(1)}
    for _ in range(ZZ(power)):
        output = series_convolve(output, series, min_exp, max_exp)
    return output


def coefficient_list(series, q_min, q_max):
    return [series.get(ZZ(e), K12(0)) for e in range(ZZ(q_min), ZZ(q_max) + 1)]


def series_min(series):
    return min(series.keys()) if series else None


def series_max(series):
    return max(series.keys()) if series else None


def equal_on_range(left, right, q_min, q_max):
    return all(
        left.get(ZZ(e), K12(0)) == right.get(ZZ(e), K12(0))
        for e in range(ZZ(q_min), ZZ(q_max) + 1)
    )

# -----------------------------------------------------------------------------
# Level-one modular forms at tau=3t.
# -----------------------------------------------------------------------------

def divisor_sigma(n, power):
    n = ZZ(n)
    if n <= 0:
        return ZZ(0)
    return sum(ZZ(d)**ZZ(power) for d in divisors(n))


E4 = {ZZ(0): K12(1)}
E6 = {ZZ(0): K12(1)}
for exponent in range(3, GEN_MAX + 1, 3):
    n = ZZ(exponent // 3)
    add_coeff(E4, exponent, 240*divisor_sigma(n, 3))
    add_coeff(E6, exponent, -504*divisor_sigma(n, 5))

# -----------------------------------------------------------------------------
# Compressed E8 lattice records (d,g,count).
#
# d=(lambda,lambda)/2 and g=(lambda,gamma).  Since ||gamma||^2=8,
# |g| <= 4 sqrt(d).  The slowest transformed A1 term has exponent d+g.
# -----------------------------------------------------------------------------

def degree_bound(target):
    d = ZZ(1)
    target = RR(target + 10)
    while RR(d) - 4*sqrt(RR(d)) <= target:
        d += 1
    return d + 3


D_BOUND = degree_bound(GEN_MAX)
MAX_NORM4 = ZZ(8*D_BOUND)
COORD_BOUND = ZZ(ceil(sqrt(RR(2*D_BOUND)))) + 2
GAMMA_SIGNS = [1, 1, 1, 1, 1, 1, 1, -1]

print("\nE8 lattice dynamic-programming bounds:", flush=True)
print("  d <= {}".format(D_BOUND), flush=True)
print("  coordinate bound = {}".format(COORD_BOUND), flush=True)


def coordinate_options(sign, half_coset):
    options = []
    if not half_coset:
        for n in range(-COORD_BOUND, COORD_BOUND + 1):
            norm4 = ZZ(4*n*n)
            if norm4 <= MAX_NORM4:
                options.append((norm4, ZZ(2*sign*n), ZZ(n % 2)))
    else:
        for n in range(-COORD_BOUND - 1, COORD_BOUND + 1):
            odd = ZZ(2*n + 1)
            norm4 = odd*odd
            if norm4 <= MAX_NORM4:
                options.append((norm4, ZZ(sign*odd), ZZ(n % 2)))
    return options


def build_coset(half_coset):
    states = {(ZZ(0), ZZ(0), ZZ(0)): ZZ(1)}
    for coordinate, sign in enumerate(GAMMA_SIGNS):
        options = coordinate_options(sign, half_coset)
        new_states = {}
        for (norm4, gamma2, parity), count in states.items():
            for dn, dg, dp in options:
                next_norm4 = norm4 + dn
                if next_norm4 > MAX_NORM4:
                    continue
                key = (next_norm4, gamma2 + dg, (parity + dp) % 2)
                new_states[key] = new_states.get(key, ZZ(0)) + count
        states = new_states
        print(
            "  {} coset coordinate {}/8: {} states".format(
                "half" if half_coset else "integer",
                coordinate + 1,
                len(states),
            ),
            flush=True,
        )
    return states


integer_states = build_coset(False)
half_states = build_coset(True)


def records_from_states(states_list):
    combined = {}
    for states in states_list:
        for (norm4, gamma2, parity), count in states.items():
            if parity != 0:
                continue
            if norm4 % 8 != 0 or gamma2 % 2 != 0:
                continue
            d = ZZ(norm4 // 8)
            g = ZZ(gamma2 // 2)
            key = (d, g)
            combined[key] = combined.get(key, ZZ(0)) + count
    return [(d, g, K12(c)) for (d, g), c in sorted(combined.items())]


RECORDS = records_from_states([integer_states, half_states])
print("  retained (d,g) records: {}".format(len(RECORDS)), flush=True)

# -----------------------------------------------------------------------------
# Theta-series constructors.
# -----------------------------------------------------------------------------

def theta_hss(r):
    output = {}
    for d, g, count in RECORDS:
        exponent = 3*d + g
        if GEN_MIN <= exponent <= GEN_MAX:
            add_coeff(output, exponent, count*omega**(r*g))
    return output


def theta_zero(scale_d, q_min=0, q_max=GEN_MAX):
    output = {}
    for d, g, count in RECORDS:
        exponent = ZZ(scale_d)*d
        if q_min <= exponent <= q_max:
            add_coeff(output, exponent, count)
    return output


def theta_A3_direct(use_gamma):
    output = {}
    for d, g, count in RECORDS:
        exponent = 9*d + (3*g if use_gamma else 0)
        if GEN_MIN <= exponent <= GEN_MAX:
            # 3*(r/3)*gamma is an E8 lattice shift, hence phase 1.
            add_coeff(output, exponent, count)
    return output


def theta_A3_divided(k, r, use_gamma):
    output = {}
    for d, g, count in RECORDS:
        exponent = d + (g if use_gamma else 0)
        if GEN_MIN <= exponent <= GEN_MAX:
            phase_power = k*d + (r*g if use_gamma else 0)
            add_coeff(output, exponent, count*omega**phase_power)
    return output

# -----------------------------------------------------------------------------
# A2/B2 theta constants and transformed theta series, using exponent units 1/8.
# -----------------------------------------------------------------------------

UNIT_MIN = ZZ(8*GEN_MIN)
UNIT_MAX = ZZ(8*(GEN_MAX + 8))


def theta_constant_units(kind, max_units):
    output = {}
    if kind in ("theta3", "theta4"):
        bound = ZZ(ceil(sqrt(RR(max_units)/12))) + 3
        for n in range(-bound, bound + 1):
            exponent = ZZ(12*n*n)
            if exponent > max_units:
                continue
            coefficient = K12(-1 if kind == "theta4" and n % 2 else 1)
            add_coeff(output, exponent, coefficient)
        return output
    if kind == "theta2":
        bound = ZZ(ceil((sqrt(RR(max_units)/3) + 1)/2)) + 3
        for n in range(-bound, bound + 1):
            odd = ZZ(2*n + 1)
            exponent = ZZ(3*odd*odd)
            if exponent <= max_units:
                add_coeff(output, exponent, 1)
        return output
    raise ValueError("unknown theta constant {}".format(kind))


theta2 = theta_constant_units("theta2", UNIT_MAX + 64)
theta3 = theta_constant_units("theta3", UNIT_MAX + 64)
theta4 = theta_constant_units("theta4", UNIT_MAX + 64)

theta2_4 = series_power(theta2, 4, 0, UNIT_MAX + 64)
theta3_4 = series_power(theta3, 4, 0, UNIT_MAX + 64)
theta4_4 = series_power(theta4, 4, 0, UNIT_MAX + 64)

e1_units = series_scale(series_add(theta3_4, theta4_4), QQ(1)/12)
e2_units = series_scale(series_add(theta2_4, series_scale(theta4_4, -1)), QQ(1)/12)
e3_units = series_scale(series_add(series_scale(theta2_4, -1), series_scale(theta3_4, -1)), QQ(1)/12)


def collapse_units(series, label):
    bad = [(e, c) for e, c in sorted(series.items()) if e % 8 != 0 and c != 0]
    if bad:
        raise ArithmeticError(
            "{} retained nonintegral q powers; first terms {}".format(label, bad[:8])
        )
    output = {}
    for exponent, coefficient in series.items():
        if exponent % 8 == 0:
            q_exp = ZZ(exponent // 8)
            if GEN_MIN <= q_exp <= GEN_MAX:
                add_coeff(output, q_exp, coefficient)
    return output


def sakai_A2_B2(r, use_gamma):
    theta_6_2 = {}
    theta_half_0 = {}
    theta_half_1 = {}

    for d, g, count in RECORDS:
        gg = g if use_gamma else ZZ(0)
        phase_r = r*g if use_gamma else ZZ(0)

        exponent_direct = ZZ(8*(6*d + 2*gg))
        if UNIT_MIN <= exponent_direct <= UNIT_MAX:
            add_coeff(theta_6_2, exponent_direct, count*omega**(2*phase_r))

        exponent_half = ZZ(12*d + 8*gg)
        if UNIT_MIN <= exponent_half <= UNIT_MAX:
            base = count*omega**phase_r
            add_coeff(theta_half_0, exponent_half, base)
            add_coeff(theta_half_1, exponent_half, base*((-1)**d))

    A2_units = series_scale(
        series_add(
            theta_6_2,
            series_scale(theta_half_0, QQ(1)/16),
            series_scale(theta_half_1, QQ(1)/16),
        ),
        QQ(8)/9,
    )

    B2_inner = series_add(
        series_convolve(e1_units, theta_6_2, UNIT_MIN, UNIT_MAX),
        series_scale(
            series_convolve(e3_units, theta_half_0, UNIT_MIN, UNIT_MAX),
            QQ(1)/16,
        ),
        series_scale(
            series_convolve(e2_units, theta_half_1, UNIT_MIN, UNIT_MAX),
            QQ(1)/16,
        ),
    )
    B2_units = series_scale(B2_inner, QQ(32)/5)

    return (
        collapse_units(A2_units, "A2 r={} gamma={}".format(r, use_gamma)),
        collapse_units(B2_units, "B2 r={} gamma={}".format(r, use_gamma)),
    )

# -----------------------------------------------------------------------------
# h0 series entering B3.
# -----------------------------------------------------------------------------

H0_MAX = ZZ(GEN_MAX + 8)


def h0_at_3t():
    output = {}
    bound = ZZ(ceil(sqrt(RR(H0_MAX)))) + 3

    # theta3(6t) theta3(18t)
    for n in range(-bound, bound + 1):
        for m in range(-bound, bound + 1):
            exponent = ZZ(3*n*n + 9*m*m)
            if exponent <= H0_MAX:
                add_coeff(output, exponent, 1)

    # theta2(6t) theta2(18t)
    for n in range(-bound, bound + 1):
        a = ZZ(2*n + 1)
        for m in range(-bound, bound + 1):
            b = ZZ(2*m + 1)
            numerator = ZZ(3*(a*a + 3*b*b))
            if numerator % 4 != 0:
                raise ArithmeticError("unexpected h0(3t) fractional exponent")
            exponent = numerator // 4
            if exponent <= H0_MAX:
                add_coeff(output, exponent, 1)
    return output


def h0_divided(k):
    output = {}
    bound = ZZ(ceil(sqrt(RR(4*H0_MAX)))) + 3

    # theta3(2(t+k/3)) theta3(6(t+k/3))
    for n in range(-bound, bound + 1):
        for m in range(-bound, bound + 1):
            exponent = ZZ(n*n + 3*m*m)
            if exponent <= H0_MAX:
                add_coeff(output, exponent, omega**(k*n*n))

    # theta2(2(t+k/3)) theta2(6(t+k/3))
    for n in range(-bound, bound + 1):
        a = ZZ(2*n + 1)
        for m in range(-bound, bound + 1):
            b = ZZ(2*m + 1)
            numerator = a*a + 3*b*b
            if numerator % 4 != 0:
                raise ArithmeticError("unexpected h0 divided fractional exponent")
            exponent = ZZ(numerator // 4)
            if exponent <= H0_MAX:
                # Exact phase exp(2*pi*i*k*(a^2+3b^2)/12).
                add_coeff(output, exponent, zeta12**(k*numerator))
    return output


h0_3 = h0_at_3t()
h0_div = [h0_divided(k) for k in range(3)]
h0_3_sq = series_convolve(h0_3, h0_3, 0, H0_MAX)
h0_div_sq = [series_convolve(h, h, 0, H0_MAX) for h in h0_div]

# -----------------------------------------------------------------------------
# A3 and B3.
# -----------------------------------------------------------------------------

def sakai_A3_B3(r, use_gamma):
    direct = theta_A3_direct(use_gamma)
    divided = [theta_A3_divided(k, r, use_gamma) for k in range(3)]

    # Sakai's normalized index-three Hecke transform:
    #   A3 = 27/28 * (A1(3*tau,3*z) + 3^(-4) sum_k A1((tau+k)/3,z)).
    # Hence the coefficient of each divided term is (27/28)/81 = 1/84.
    A3 = series_add(
        series_scale(direct, QQ(27)/28),
        series_scale(series_add(*divided), QQ(1)/84),
    )

    first_B = series_convolve(h0_3_sq, direct, GEN_MIN, GEN_MAX)
    divided_B_terms = [
        series_convolve(h0_div_sq[k], divided[k], GEN_MIN, GEN_MAX)
        for k in range(3)
    ]
    B3 = series_scale(
        series_add(
            first_B,
            series_scale(series_add(*divided_B_terms), -QQ(1)/243),
        ),
        QQ(81)/80,
    )
    return A3, B3

# -----------------------------------------------------------------------------
# Normalization checks at z=0.
# -----------------------------------------------------------------------------

print("\nPART I. SAKAI NORMALIZATION CHECKS", flush=True)
print("-"*79, flush=True)

A1_zero = theta_zero(3)
A2_zero, B2_zero = sakai_A2_B2(ZZ(0), False)
A3_zero, B3_zero = sakai_A3_B3(ZZ(0), False)

checks = {
    "A1(3t,0)=E4(3t)": equal_on_range(A1_zero, E4, 0, CHECK_MAX),
    "A2(3t,0)=E4(3t)": equal_on_range(A2_zero, E4, 0, CHECK_MAX),
    "A3(3t,0)=E4(3t)": equal_on_range(A3_zero, E4, 0, CHECK_MAX),
    "B2(3t,0)=E6(3t)": equal_on_range(B2_zero, E6, 0, CHECK_MAX),
    "B3(3t,0)=E6(3t)": equal_on_range(B3_zero, E6, 0, CHECK_MAX),
}
for label, passed in checks.items():
    print("  {:28s} {}".format(label, passed), flush=True)
if not all(checks.values()):
    raise ArithmeticError("A Sakai normalization check failed")

# -----------------------------------------------------------------------------
# Build the character-specialized generators and six basis series.
# -----------------------------------------------------------------------------

print("\nPART II. CHARACTER-SPECIALIZED ONE-SIDED BASIS", flush=True)
print("-"*79, flush=True)

forms = {}
basis = {}

for character in CHARACTERS:
    r = CHARACTER_R[character]
    A1 = theta_hss(r)
    A2, B2 = sakai_A2_B2(r, True)
    A3, B3 = sakai_A3_B3(r, True)

    forms[character] = {
        "A1": A1,
        "A2": A2,
        "B2": B2,
        "A3": A3,
        "B3": B3,
    }

    E6_sq = series_convolve(E6, E6, 0, GEN_MAX)
    E4_sq = series_convolve(E4, E4, 0, GEN_MAX)
    E4_cube = series_convolve(E4_sq, E4, 0, GEN_MAX)
    E4E6 = series_convolve(E4, E6, 0, GEN_MAX)

    h1 = series_convolve(E6_sq, A3, GEN_MIN, GEN_MAX)
    h2 = series_convolve(E4_cube, A3, GEN_MIN, GEN_MAX)
    h3 = series_convolve(E4E6, B3, GEN_MIN, GEN_MAX)
    h4 = series_convolve(
        E4_sq,
        series_convolve(A1, A2, GEN_MIN, GEN_MAX),
        GEN_MIN,
        GEN_MAX,
    )
    h5 = series_convolve(
        E6,
        series_convolve(A1, B2, GEN_MIN, GEN_MAX),
        GEN_MIN,
        GEN_MAX,
    )
    h6 = series_convolve(
        E4,
        series_power(A1, 3, GEN_MIN, GEN_MAX),
        GEN_MIN,
        GEN_MAX,
    )

    basis[character] = {
        BASiS_LABEL: series
        for BASiS_LABEL, series in zip(BASIS_LABELS, [h1, h2, h3, h4, h5, h6])
    }

    print("  {}:".format(character), flush=True)
    for label in BASIS_LABELS:
        series = basis[character][label]
        print(
            "    {:16s} support [{},{}]".format(
                label,
                series_min(series),
                series_max(series),
            ),
            flush=True,
        )

# -----------------------------------------------------------------------------
# Shared symmetric 21-column sector matrix.
# -----------------------------------------------------------------------------

print("\nPART III. SHARED TWO-SECTOR q-SERIES RANK", flush=True)
print("-"*79, flush=True)

pairs = [(i, j) for i in range(6) for j in range(i, 6)]
assert len(pairs) == 21


def symmetric_basis_product(left_character, right_character, i, j):
    left_i = basis[left_character][BASIS_LABELS[i]]
    right_j = basis[right_character][BASIS_LABELS[j]]
    first = series_convolve(left_i, right_j, 2*GEN_MIN, SECTOR_Q_MAX)
    if i == j:
        return first
    left_j = basis[left_character][BASIS_LABELS[j]]
    right_i = basis[right_character][BASIS_LABELS[i]]
    second = series_convolve(left_j, right_i, 2*GEN_MIN, SECTOR_Q_MAX)
    return series_add(first, second)


sector_specs = [
    ("(omega,omega)", "omega", "omega"),
    ("(omega,omega2)", "omega", "omega2"),
]

sector_columns = {}
sector_min = ZZ(0)
first_min = True
for sector_name, left_character, right_character in sector_specs:
    columns = [
        symmetric_basis_product(left_character, right_character, i, j)
        for i, j in pairs
    ]
    sector_columns[sector_name] = columns
    mins = [series_min(column) for column in columns if column]
    if mins:
        this_min = min(mins)
        sector_min = this_min if first_min else min(sector_min, this_min)
        first_min = False

rows_by_sector = {}
for sector_name, _, _ in sector_specs:
    rows = []
    for q_power in range(sector_min, SECTOR_Q_MAX + 1):
        rows.append([
            column.get(ZZ(q_power), K12(0))
            for column in sector_columns[sector_name]
        ])
    rows_by_sector[sector_name] = rows

rank_by_sector = {}
for sector_name, _, _ in sector_specs:
    M = matrix(K12, rows_by_sector[sector_name])
    rank_by_sector[sector_name] = ZZ(M.rank())
    print("  {:18s} rank = {}".format(sector_name, rank_by_sector[sector_name]), flush=True)

combined_rows = []
row_labels = []
for sector_name, _, _ in sector_specs:
    for q_power, row in zip(range(sector_min, SECTOR_Q_MAX + 1), rows_by_sector[sector_name]):
        combined_rows.append(row)
        row_labels.append((sector_name, ZZ(q_power)))

combined_matrix = matrix(K12, combined_rows)
combined_rank = ZZ(combined_matrix.rank())
combined_nullity = ZZ(21 - combined_rank)

print("  q-power range used: [{} , {}]".format(sector_min, SECTOR_Q_MAX), flush=True)
print("  combined rows: {}".format(combined_matrix.nrows()), flush=True)
print("  shared unknowns: 21", flush=True)
print("  combined exact rank: {}".format(combined_rank), flush=True)
print("  remaining evaluation nullity: {}".format(combined_nullity), flush=True)

# Incremental rank profile: record only points at which the rank rises.
rank_profile = []
previous_rank = ZZ(-1)
for cutoff in range(sector_min, SECTOR_Q_MAX + 1):
    incremental_rows = []
    for sector_name, _, _ in sector_specs:
        start = ZZ(cutoff - sector_min + 1)
        incremental_rows.extend(rows_by_sector[sector_name][:start])
    rank_value = ZZ(matrix(K12, incremental_rows).rank())
    if rank_value != previous_rank:
        rank_profile.append((ZZ(cutoff), rank_value))
        previous_rank = rank_value

print("  incremental rank jumps:", flush=True)
for cutoff, rank_value in rank_profile:
    print("    through q^{:>3}: rank {}".format(cutoff, rank_value), flush=True)

# -----------------------------------------------------------------------------
# Save a shifted coefficient-list data file plus the full Laurent result.
# -----------------------------------------------------------------------------

all_mins = [
    series_min(basis[character][label])
    for character in CHARACTERS
    for label in BASIS_LABELS
]
GLOBAL_MIN = min(all_mins)
SHIFT = ZZ(max(0, -GLOBAL_MIN))
SHIFTED_MAX = ZZ(GEN_MAX + SHIFT)

shifted_series = {}
for character in CHARACTERS:
    shifted_series[character] = {}
    for label in BASIS_LABELS:
        source_series = basis[character][label]
        shifted_series[character][label] = [
            source_series.get(ZZ(n - SHIFT), K12(0))
            for n in range(SHIFTED_MAX + 1)
        ]

DATA_PATH = os.path.join(RESULTS_DIR, "degree3_one_sided_torsion_basis.sobj")
FULL_PATH = os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.sobj")
SUMMARY_PATH = os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.txt")

compatibility_data = {
    "q_order": SHIFTED_MAX,
    "q_shift": SHIFT,
    "original_q_min": GLOBAL_MIN,
    "original_q_max": GEN_MAX,
    "characters": CHARACTERS,
    "basis_labels": BASIS_LABELS,
    "series": shifted_series,
    "coefficient_ring": str(K12),
    "coefficient_field_conductor": ZZ(12),
    "torsion_candidate": "z_r=t*gamma+(r/3)*gamma",
    "notes": (
        "Every one-sided series was multiplied by the common q^q_shift so "
        "coefficient lists begin at q^0.  This common shift does not change "
        "the shared-sector evaluation rank."
    ),
}
save(compatibility_data, DATA_PATH)

full_result = {
    "sector_q_max": SECTOR_Q_MAX,
    "one_sided_range": (GEN_MIN, GEN_MAX),
    "coefficient_field": K12,
    "zeta12": zeta12,
    "omega": omega,
    "gamma": vector(ZZ, [1,1,1,1,1,1,1,-1]),
    "torsion_candidate": "z_r=t*gamma+(r/3)*gamma",
    "normalization_checks": checks,
    "forms": forms,
    "basis": basis,
    "basis_labels": BASIS_LABELS,
    "pairs": pairs,
    "sector_columns": sector_columns,
    "sector_min": sector_min,
    "rank_by_sector": rank_by_sector,
    "combined_rank": combined_rank,
    "combined_nullity": combined_nullity,
    "rank_profile": rank_profile,
    "row_labels": row_labels,
}
save(full_result, FULL_PATH)

summary_lines = [
    "DEGREE-THREE A3/B3 HSS TORSION-CANDIDATE SERIES",
    "candidate: z_r=t*gamma+(r/3)*gamma, r=0,1,2",
    "coefficient field: Q(zeta_12)",
    "sector q range: [{} , {}]".format(sector_min, SECTOR_Q_MAX),
    "",
    "Sakai normalization checks:",
]
summary_lines.extend("  {}: {}".format(label, passed) for label, passed in checks.items())
summary_lines.extend([
    "",
    "sector ranks:",
    "  (omega,omega): {}".format(rank_by_sector["(omega,omega)"]),
    "  (omega,omega2): {}".format(rank_by_sector["(omega,omega2)"]),
    "combined shared-sector rank: {} / 21".format(combined_rank),
    "remaining evaluation nullity: {}".format(combined_nullity),
    "",
    "rank jumps:",
])
summary_lines.extend("  through q^{}: {}".format(qp, rv) for qp, rv in rank_profile)
summary_lines.extend([
    "",
    "scope: this is the canonical gamma/3 torsion-point candidate;",
    "the BKOS torsion-label identification remains to be independently checked.",
])

with open(SUMMARY_PATH, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

print("\nSaved:", flush=True)
print("  {}".format(DATA_PATH), flush=True)
print("  {}".format(FULL_PATH), flush=True)
print("  {}".format(SUMMARY_PATH), flush=True)
print("\nCELL 32 COMPLETE", flush=True)
