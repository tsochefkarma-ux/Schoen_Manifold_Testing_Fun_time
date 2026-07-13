# Sage symbols are already loaded by cell 08.
# Do not re-import sage.all here: it can overwrite user ring names such as R.
from itertools import combinations

# ============================================================================
# EXPLICIT LINEARIZED CORRECTION-KERNEL STRUCTURE
# ============================================================================
#
# FOLLOW-UP CELL
# --------------
# Run 08_linearized_correction_jet_test.sage first in the SAME Sage session,
# then paste this entire cell underneath it.
#
# The preceding cell must have defined:
#
#   M, BaseR, BaseF, LL, EE4,
#   active_degrees, active_rows,
#   kernel, u_fixed, u_E4_tangent,
#   alpha1, beta1, phi, P, PF, jj,
#   u, theta, Alog, dK,
#   rL, rE4, rdelta, rE6,
#   K and the fixed-locus map phi,
#   rC0, rX, to_BaseF_fixed.
#
# PURPOSE
# -------
# 1. Print two exact independent obstruction functionals Lambda_1,Lambda_2.
# 2. Print the reduced row-echelon form of the constraint matrix.
# 3. Build a four-vector kernel basis whose first two vectors are:
#       - overall eta-series rescaling;
#       - the formal E4 tangent.
# 4. Extract two complementary nontrivial kernel directions.
# 5. Print the induced first-order shifts alpha_1 and beta_1 for every basis
#    vector.
# 6. Test natural candidate directions:
#       - delta tangent;
#       - E6 tangent;
#       - L tangent;
#       - q-translation / theta-shift;
#       - jet-index weighting.
#
# All calculations are exact over QQ(L,E4).
# ============================================================================

PRINT_RREF = False
PRINT_COMPACT_OBSTRUCTIONS = True
PRINT_FULL_BASIS = False
PRINT_COMPLEMENTARY_BASIS = True
PRINT_SHIFT_FORMS = False
PRINT_CANDIDATE_COORDINATES = True

print("="*78)
print("EXPLICIT LINEARIZED CORRECTION-KERNEL STRUCTURE")
print("Exact coefficient field: QQ(L,E4)")
print("="*78)
print()

# ----------------------------------------------------------------------------
# Confirm that the previous cell is loaded
# ----------------------------------------------------------------------------

required_names = [
    "M", "BaseR", "BaseF", "LL", "EE4",
    "active_degrees", "active_rows",
    "u_fixed", "u_E4_tangent", "alpha1", "beta1", "K",
    "phi", "P", "PF", "jj", "u", "theta", "Alog", "dK",
    "rL", "rE4", "rdelta", "rE6", "rC0", "rX",
    "to_BaseF_fixed",
]

missing = [name for name in required_names if name not in globals()]

if missing:
    raise RuntimeError(
        "Run 08_linearized_correction_jet_test.sage first in the same "
        "Sage session. Missing names: {}".format(missing)
    )

# Earlier versions of this follow-up cell used `from sage.all import *`, which
# can overwrite the polynomial-ring variable R with Sage's external R-language
# interface.  Recover the exact source ring directly from the homomorphism phi
# created by cell 08.  This also repairs inherited helpers such as dK() and
# theta(), whose function bodies resolve the global name R when called.
R = phi.domain()

assert K == R.fraction_field()
assert M.base_ring() == BaseF
assert M.ncols() == 6
assert M.rank() == 2
preserving_kernel = M.right_kernel()
assert preserving_kernel.dimension() == 4

j_names = ["j0", "j1", "j2", "j3", "j4", "j5"]

print("Loaded preceding constraint system successfully.")
print("  constraint-matrix shape:", M.nrows(), "x", M.ncols())
print("  rank:", M.rank())
print("  kernel dimension:", preserving_kernel.dimension())
print()

# ----------------------------------------------------------------------------
# Exact simplification and formatting helpers
# ----------------------------------------------------------------------------

def exact_polynomial_quotient(a, b):
    """Return a/b in BaseR, asserting exact polynomial division."""
    a = BaseR(a)
    b = BaseR(b)
    q, rem = a.quo_rem(b)
    if rem != 0:
        raise ArithmeticError("non-exact polynomial division")
    return BaseR(q)


def primitive_polynomial_vector(vector_entries):
    """
    Clear all rational-function denominators and divide out the polynomial gcd.

    The result is a primitive vector over QQ[L,E4], defined up to a nonzero
    rational scalar. It spans the same one-dimensional direction over BaseF.
    """
    entries = [BaseF(value) for value in vector_entries]

    common_denominator = BaseR(1)
    for value in entries:
        denominator = BaseR(value.denominator())
        common_denominator = common_denominator.lcm(denominator)

    polynomials = []
    for value in entries:
        numerator = BaseR(value.numerator())
        denominator = BaseR(value.denominator())
        multiplier = exact_polynomial_quotient(
            common_denominator,
            denominator,
        )
        polynomials.append(BaseR(numerator*multiplier))

    nonzero = [value for value in polynomials if value != 0]
    if not nonzero:
        return vector(BaseR, [0 for _ in polynomials])

    common_gcd = nonzero[0]
    for value in nonzero[1:]:
        common_gcd = common_gcd.gcd(value)

    if common_gcd != 0 and not common_gcd.is_unit():
        polynomials = [
            exact_polynomial_quotient(value, common_gcd)
            if value != 0 else BaseR(0)
            for value in polynomials
        ]

    # Normalize the sign using the first nonzero leading coefficient.
    for value in polynomials:
        if value != 0:
            if value.leading_coefficient() < 0:
                polynomials = [-entry for entry in polynomials]
            break

    return vector(BaseR, polynomials)


def factored_fraction(value):
    """Readable exact factorization of an element of QQ(L,E4)."""
    value = BaseF(value)
    if value == 0:
        return "0"

    numerator = BaseR(value.numerator())
    denominator = BaseR(value.denominator())

    numerator_text = str(factor(numerator))
    denominator_text = str(factor(denominator))

    if denominator == 1:
        return numerator_text

    return "({})/({})".format(numerator_text, denominator_text)


def format_linear_functional(row):
    """Format row dot (j0,...,j5) as an exact linear expression."""
    pieces = []

    for name, coefficient in zip(j_names, row):
        coefficient = BaseF(coefficient)
        if coefficient == 0:
            continue
        pieces.append("({})*{}".format(
            factored_fraction(coefficient),
            name,
        ))

    return " + ".join(pieces) if pieces else "0"


def safe_factor_text(value):
    """Return a readable factorization, treating zero safely."""
    if value == 0:
        return "0"

    try:
        return str(factor(value))
    except (ArithmeticError, TypeError, NotImplementedError):
        return str(value)


def print_labelled_vector(label, vector_value, primitive=True):
    """Print a six-vector with labelled jet components."""
    print(label)

    if primitive:
        displayed = primitive_polynomial_vector(vector_value)
        ring_label = "primitive polynomial representative"
    else:
        displayed = vector(BaseF, vector_value)
        ring_label = "rational-function representative"

    print("  ({})".format(ring_label))
    for name, value in zip(j_names, displayed):
        print("    {} = {}".format(name, safe_factor_text(value)))
    print()

    return displayed


def independent_rank(vectors):
    if len(vectors) == 0:
        return 0
    return Matrix(BaseF, [list(vector(BaseF, v)) for v in vectors]).rank()


def proportional(v, w):
    """Exact proportionality test for two nonzero vectors over BaseF."""
    v = vector(BaseF, v)
    w = vector(BaseF, w)

    if v == 0 or w == 0:
        return bool(v == 0 and w == 0)

    return bool(Matrix(BaseF, [list(v), list(w)]).rank() == 1)

# ----------------------------------------------------------------------------
# Reduced row-echelon form and compact independent obstruction equations
# ----------------------------------------------------------------------------

rref = M.echelon_form()
rref_rows = [
    vector(BaseF, rref.row(index))
    for index in range(rref.nrows())
    if rref.row(index) != 0
]

assert len(rref_rows) == 2
assert Matrix(BaseF, rref_rows).rank() == 2

print("Reduced row-echelon data:")
print("  pivot columns:", list(rref.pivots()))
print()

if PRINT_RREF:
    for index, row in enumerate(rref_rows, start=1):
        print("  Lambda_RREF_{} = 0".format(index))
        print("   ", format_linear_functional(row))
        print()

# Among the four original X-coefficient equations, choose the independent pair
# with the shortest primitive factored representation.

def row_complexity(row):
    primitive = primitive_polynomial_vector(row)
    return sum(len(str(factor(entry))) for entry in primitive if entry != 0)

independent_pairs = []

for first, second in combinations(range(len(active_rows)), 2):
    pair_matrix = Matrix(
        BaseF,
        [active_rows[first], active_rows[second]],
    )

    if pair_matrix.rank() == 2:
        score = (
            row_complexity(active_rows[first])
            + row_complexity(active_rows[second])
        )
        independent_pairs.append((score, first, second))

assert independent_pairs
independent_pairs.sort(key=lambda item: item[0])
_, best_first, best_second = independent_pairs[0]

compact_rows = [
    vector(BaseF, active_rows[best_first]),
    vector(BaseF, active_rows[best_second]),
]

compact_degrees = [
    active_degrees[best_first],
    active_degrees[best_second],
]

assert Matrix(BaseF, compact_rows).rank() == 2
compact_kernel = Matrix(BaseF, compact_rows).right_kernel()
assert compact_kernel.dimension() == preserving_kernel.dimension()
assert all(M*vector(BaseF, value) == 0 for value in compact_kernel.basis())

print("Compact independent obstruction equations:")
print("  selected original X degrees:", compact_degrees)
print()

if PRINT_COMPACT_OBSTRUCTIONS:
    for number, (degree, row) in enumerate(
        zip(compact_degrees, compact_rows),
        start=1,
    ):
        primitive = primitive_polynomial_vector(row)
        print("  Lambda_{} from [X^{}] = 0".format(number, degree))
        print("   ", format_linear_functional(
            [BaseF(entry) for entry in primitive]
        ))
        print()

# ----------------------------------------------------------------------------
# Construct a basis beginning with the two known geometric directions
# ----------------------------------------------------------------------------

scale_direction = vector(BaseF, u_fixed)
E4_direction = vector(BaseF, u_E4_tangent)

assert M*scale_direction == 0
assert M*E4_direction == 0
assert independent_rank([scale_direction, E4_direction]) == 2

chosen_basis = [scale_direction, E4_direction]
extra_directions = []

for candidate in preserving_kernel.basis():
    candidate = vector(BaseF, candidate)
    old_rank = independent_rank(chosen_basis)
    new_rank = independent_rank(chosen_basis + [candidate])

    if new_rank > old_rank:
        chosen_basis.append(candidate)
        extra_directions.append(candidate)

    if len(extra_directions) == 2:
        break

assert len(chosen_basis) == 4
assert len(extra_directions) == 2
assert independent_rank(chosen_basis) == 4
assert all(M*vector_value == 0 for vector_value in chosen_basis)

print("Kernel-basis construction:")
print("  first direction: overall eta-series rescaling")
print("  second direction: formal E4 tangent")
print("  two exact complementary directions extracted")
print("  resulting basis rank:", independent_rank(chosen_basis))
print()

if PRINT_FULL_BASIS:
    displayed_basis = []

    displayed_basis.append(print_labelled_vector(
        "K1: overall rescaling direction",
        chosen_basis[0],
    ))

    displayed_basis.append(print_labelled_vector(
        "K2: formal E4-tangent direction",
        chosen_basis[1],
    ))

    displayed_basis.append(print_labelled_vector(
        "K3: first complementary preserving direction",
        chosen_basis[2],
    ))

    displayed_basis.append(print_labelled_vector(
        "K4: second complementary preserving direction",
        chosen_basis[3],
    ))

elif PRINT_COMPLEMENTARY_BASIS:
    print("Known directions K1 and K2 retained but not expanded.")
    print("  K1 = overall eta-series rescaling")
    print("  K2 = formal E4 tangent")
    print()

    print_labelled_vector(
        "K3: first complementary preserving direction",
        chosen_basis[2],
    )

    print_labelled_vector(
        "K4: second complementary preserving direction",
        chosen_basis[3],
    )

# ----------------------------------------------------------------------------
# Map alpha_1 and beta_1 to the fixed-locus correction-jet ring
# ----------------------------------------------------------------------------

def map_K_expression_to_PF_fixed(expr):
    """
    Map a K-expression directly to PF using the homomorphism phi created
    by cell 08.  That map already performs

        delta -> 0,  E6 -> 0,  C0 -> 1,
        X -> XX,     j_m -> jj_m,
        L -> LL,     E4 -> EE4.

    We deliberately avoid the global name R, because `from sage.all import *`
    may bind R to Sage's external R-language interface rather than to the
    polynomial ring created in cell 08.
    """
    expr = K(expr)

    numerator = phi(expr.numerator())
    denominator = phi(expr.denominator())

    return PF(numerator)/PF(denominator)


alpha1_form = map_K_expression_to_PF_fixed(alpha1)
beta1_form = map_K_expression_to_PF_fixed(beta1)


def evaluate_PF_on_jet_vector(expr, vector_value):
    """Evaluate a PF expression on one correction-jet vector."""
    vector_value = vector(BaseF, vector_value)

    evaluation = P.hom(
        [BaseF(0)] + [BaseF(value) for value in vector_value],
        BaseF,
    )

    numerator = evaluation(P(PF(expr).numerator()))
    denominator = evaluation(P(PF(expr).denominator()))

    return BaseF(numerator)/BaseF(denominator)


# Sanity: alpha1_form and beta1_form must vanish at zero correction.
zero_vector = vector(BaseF, [0]*6)
assert evaluate_PF_on_jet_vector(alpha1_form, zero_vector) == 0
assert evaluate_PF_on_jet_vector(beta1_form, zero_vector) == 0

if PRINT_SHIFT_FORMS:
    print("Fixed-locus alpha_1 linear form:")
    print(alpha1_form)
    print()
    print("Fixed-locus beta_1 linear form:")
    print(beta1_form)
    print()

print("Induced closure-coefficient shifts on the chosen kernel basis:")

for index, vector_value in enumerate(chosen_basis, start=1):
    alpha_shift = evaluate_PF_on_jet_vector(alpha1_form, vector_value)
    beta_shift = evaluate_PF_on_jet_vector(beta1_form, vector_value)

    print("  K{}:".format(index))
    print("    alpha_1 =", factored_fraction(alpha_shift))
    print("    beta_1  =", factored_fraction(beta_shift))

print()

# ----------------------------------------------------------------------------
# Coordinates in the chosen kernel basis
# ----------------------------------------------------------------------------

basis_matrix = Matrix(
    BaseF,
    6,
    4,
    lambda row, column: chosen_basis[column][row],
)

assert basis_matrix.rank() == 4


def kernel_coordinates(vector_value):
    """Coordinates of a kernel vector in K1,K2,K3,K4."""
    vector_value = vector(BaseF, vector_value)

    if M*vector_value != 0:
        raise ValueError("vector is not in the preserving kernel")

    coordinates = basis_matrix.solve_right(vector_value)
    assert basis_matrix*coordinates == vector_value
    return vector(BaseF, coordinates)

# ----------------------------------------------------------------------------
# Natural candidate directions
# ----------------------------------------------------------------------------

# Extend the normalized jet recurrence once more for q-translation.
u_extended = list(u)
u_extended.append(theta(u_extended[-1]) + Alog*u_extended[-1])
assert len(u_extended) == 7

# Formal tangents of the fixed normalized jets.
delta_tangent = vector(
    BaseF,
    [to_BaseF_fixed(dK(u[m], rdelta)) for m in range(6)],
)

E6_tangent = vector(
    BaseF,
    [to_BaseF_fixed(dK(u[m], rE6)) for m in range(6)],
)

L_tangent = vector(
    BaseF,
    [to_BaseF_fixed(dK(u[m], rL)) for m in range(6)],
)

# Infinitesimal q-translation: B(q) -> B(q*exp(eps)), hence
# theta^m B -> theta^m B + eps*theta^(m+1)B.
q_translation = vector(
    BaseF,
    [to_BaseF_fixed(u_extended[m+1]) for m in range(6)],
)

# Simple grading/index-weight direction.
index_weight = vector(
    BaseF,
    [BaseF(m)*u_fixed[m] for m in range(6)],
)

candidate_directions = [
    ("formal delta tangent", delta_tangent),
    ("formal E6 tangent", E6_tangent),
    ("formal L tangent with geometry held fixed", L_tangent),
    ("q-translation / theta-shift", q_translation),
    ("jet-index weighting m*u_m", index_weight),
]

print("Natural candidate-direction tests:")
print()

for label, candidate in candidate_directions:
    obstruction = M*candidate
    preserves = bool(obstruction == 0)

    print("  {}".format(label))
    print("    lies in preserving kernel?", preserves)

    if preserves:
        coordinates = kernel_coordinates(candidate)

        print("    coordinates in (K1,K2,K3,K4):")
        print("     ", [factored_fraction(value) for value in coordinates])

        alpha_shift = evaluate_PF_on_jet_vector(alpha1_form, candidate)
        beta_shift = evaluate_PF_on_jet_vector(beta1_form, candidate)

        print("    alpha_1 =", factored_fraction(alpha_shift))
        print("    beta_1  =", factored_fraction(beta_shift))

        matches = []
        for basis_index, basis_vector in enumerate(chosen_basis, start=1):
            if proportional(candidate, basis_vector):
                matches.append("K{}".format(basis_index))

        if matches:
            print("    proportional to:", ", ".join(matches))

    else:
        nonzero_entries = [
            (active_degrees[row_index], obstruction[row_index])
            for row_index in range(len(obstruction))
            if obstruction[row_index] != 0
        ]

        print("    nonzero obstruction coefficients:")
        for degree, value in nonzero_entries:
            print("      [X^{}] {}".format(
                degree,
                factored_fraction(value),
            ))

    print()

# ----------------------------------------------------------------------------
# Quotient information: verify that K3,K4 are genuinely new modulo K1,K2
# ----------------------------------------------------------------------------

known_span_rank = independent_rank([scale_direction, E4_direction])
full_span_rank = independent_rank(chosen_basis)

assert known_span_rank == 2
assert full_span_rank == 4

print("Quotient-kernel check:")
print("  dim span(rescaling,E4 tangent) =", known_span_rank)
print("  dim full preserving kernel     =", full_span_rank)
print("  quotient dimension             =", full_span_rank-known_span_rank)
print()

# Explicitly verify that neither extra direction lies in the known span.
known_matrix = Matrix(
    BaseF,
    6,
    2,
    lambda row, column: chosen_basis[column][row],
)

for number, extra in enumerate(extra_directions, start=3):
    augmented = known_matrix.augment(matrix(BaseF, 6, 1, list(extra)))
    outside = bool(augmented.rank() == 3)
    print("  K{} is outside known two-dimensional span? {}".format(
        number,
        outside,
    ))
    assert outside

print()
print("="*78)
print("SUCCESS")
print()
print("Two exact obstruction functionals and a four-vector preserving basis")
print("have been extracted over QQ(L,E4).")
print()
print("K1 and K2 are the known rescaling and E4-tangent directions.")
print("K3 and K4 are exact complementary preserving directions.")
print()
print("The output above reports whether natural modular/coordinate candidates")
print("identify either of the two additional directions.")
print()
print("Please copy the complete printed output back into the chat.")
print("="*78)
