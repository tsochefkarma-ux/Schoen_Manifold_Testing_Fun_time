from sage.all import *
import os

# =============================================================================
# CELL 34 v3: CERTIFIED-CUTOFF ONE-SIDED RELATION AND 15D SECTOR QUOTIENT
# =============================================================================
#
# Corrects the original Cell 34, which accidentally treated Cell 32's internal
# convolution buffer as certified q-series data.  This version:
#
#   * uses only q <= sector_q_max for the one-sided rank computation;
#   * labels coefficients above sector_q_max as an unsafe working buffer;
#   * constructs the induced 15 x 21 symmetric-square quotient;
#   * reconstructs the six sector matrices directly through the certified
#     sector cutoff;
#   * verifies factorization by kernel annihilation and an exact right section
#     of the quotient map, rather than by multiplying truncated quotient series;
#   * compares the quotient kernel with the exact parked-sector kernel from
#     Cell 33.
#
# Scope: exact for the certified finite q-range and the canonical gamma/3
# torsion candidate.  It is not yet an all-orders Jacobi-form identity.
# =============================================================================

print("="*79, flush=True)
print("DEGREE-THREE CERTIFIED-CUTOFF 15D QUOTIENT DIAGNOSTIC", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

CELL32_PATH = os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.sobj")
CELL33_PATH = os.path.join(RESULTS_DIR, "degree3_evaluation_kernel_and_sector_lift.sobj")
if not os.path.exists(CELL32_PATH):
    raise OSError("Missing Cell-32 result: {}".format(CELL32_PATH))
if not os.path.exists(CELL33_PATH):
    raise OSError("Missing Cell-33 result: {}".format(CELL33_PATH))

data32 = load(CELL32_PATH)
data33 = load(CELL33_PATH)

K = data32["coefficient_field"]
basis = data32["basis"]
basis_labels = list(data32["basis_labels"])
pairs6 = [tuple(pair) for pair in data32["pairs"]]
characters = ["one", "omega", "omega2"]
q_min, q_buffer_max = [ZZ(x) for x in data32["one_sided_range"]]
q_cert = ZZ(data32["sector_q_max"])
sector_min = ZZ(data32["sector_min"])

if len(basis_labels) != 6 or len(pairs6) != 21:
    raise ArithmeticError("Unexpected Cell-32 dimensions")
if q_cert > q_buffer_max:
    raise ArithmeticError("Certified cutoff exceeds stored working range")

print("coefficient field: {}".format(K), flush=True)
print("certified q cutoff: q^{}".format(q_cert), flush=True)
print("stored working buffer: q^{} through q^{}".format(q_min, q_buffer_max), flush=True)
if q_buffer_max > q_cert:
    print("unsafe upper buffer ignored: q^{} through q^{}".format(q_cert + 1, q_buffer_max), flush=True)

# -----------------------------------------------------------------------------
# Helpers.
# -----------------------------------------------------------------------------

def add_coeff(target, exponent, coefficient):
    exponent = ZZ(exponent)
    coefficient = K(coefficient)
    if coefficient == 0:
        return
    target[exponent] = target.get(exponent, K(0)) + coefficient
    if target[exponent] == 0:
        del target[exponent]


def series_add(*series_list):
    output = {}
    for series in series_list:
        for exponent, coefficient in series.items():
            add_coeff(output, exponent, coefficient)
    return output


def series_convolve(left, right, min_exp, max_exp):
    output = {}
    min_exp = ZZ(min_exp)
    max_exp = ZZ(max_exp)
    for e1, c1 in left.items():
        if c1 == 0:
            continue
        for e2, c2 in right.items():
            if c2 == 0:
                continue
            exponent = ZZ(e1 + e2)
            if min_exp <= exponent <= max_exp:
                add_coeff(output, exponent, c1*c2)
    return output


def matrix_for_character(character, top_power):
    return matrix(K, [
        [basis[character][label].get(ZZ(exponent), K(0)) for label in basis_labels]
        for exponent in range(q_min, ZZ(top_power) + 1)
    ])


def stack_matrices(matrices, ncols):
    rows = []
    for M in matrices:
        rows.extend(M.rows())
    return matrix(K, rows) if rows else matrix(K, 0, ncols)


def primitive_normalize(vector_in):
    values = list(vector_in)
    rationals = []
    rational = True
    for value in values:
        try:
            rationals.append(QQ(value))
        except (TypeError, ValueError):
            rational = False
            break
    if rational:
        den = lcm([x.denominator() for x in rationals])
        integers = [ZZ(den*x) for x in rationals]
        nz = [abs(x) for x in integers if x != 0]
        common = gcd(nz) if nz else ZZ(1)
        integers = [x // common for x in integers]
        first = next((x for x in integers if x != 0), ZZ(1))
        if first < 0:
            integers = [-x for x in integers]
        return vector(K, integers), vector(ZZ, integers), True
    first_index = next(i for i, x in enumerate(values) if x != 0)
    normalized = vector(K, [x/values[first_index] for x in values])
    return normalized, normalized, False


def compact(value, limit=500):
    text = str(value)
    return text if len(text) <= limit else "<exact expression stored; {} chars>".format(len(text))

# -----------------------------------------------------------------------------
# Part I. Certified one-sided rank and unique relation.
# -----------------------------------------------------------------------------

print("\nPART I. CERTIFIED ONE-SIDED NULL RELATION", flush=True)
print("-"*79, flush=True)

M_char = {ch: matrix_for_character(ch, q_cert) for ch in characters}
M_combined = stack_matrices([M_char[ch] for ch in characters], 6)
char_ranks = {ch: ZZ(M_char[ch].rank()) for ch in characters}
combined_rank = ZZ(M_combined.rank())
combined_kernel = M_combined.right_kernel()
combined_nullity = ZZ(combined_kernel.dimension())

for ch in characters:
    print("  {:7s} rank through certified q^{} = {} / 6".format(ch, q_cert, char_ranks[ch]), flush=True)
print("  combined rank through certified q^{} = {} / 6".format(q_cert, combined_rank), flush=True)
print("  common nullity = {}".format(combined_nullity), flush=True)

# Diagnostic only: show what happens if the unsafe buffer is included.
buffer_rank = None
first_unsafe_jump = None
if q_buffer_max > q_cert:
    M_buffer = stack_matrices([matrix_for_character(ch, q_buffer_max) for ch in characters], 6)
    buffer_rank = ZZ(M_buffer.rank())
    previous = combined_rank
    for cutoff in range(q_cert + 1, q_buffer_max + 1):
        r = ZZ(stack_matrices([matrix_for_character(ch, cutoff) for ch in characters], 6).rank())
        if r > previous:
            first_unsafe_jump = (ZZ(cutoff), r)
            break
    print("  diagnostic rank through unsafe buffer q^{} = {} / 6".format(q_buffer_max, buffer_rank), flush=True)
    if first_unsafe_jump is not None:
        print("  first rank jump outside certified range: q^{} -> rank {}".format(*first_unsafe_jump), flush=True)

if combined_rank != 5 or combined_nullity != 1:
    raise ArithmeticError(
        "Certified range does not have the expected rank-five one-sided image"
    )

null_raw = combined_kernel.basis()[0]
null_vector, null_display, rational_null = primitive_normalize(null_raw)

relation_checks = {}
for ch in characters:
    zero = vector(K, M_char[ch].nrows(), [0]*M_char[ch].nrows())
    relation_checks[ch] = bool(M_char[ch]*null_vector == zero)
    print("  {:7s} relation verified through q^{}? {}".format(ch, q_cert, relation_checks[ch]), flush=True)
if not all(relation_checks.values()):
    raise ArithmeticError("Certified relation failed a character check")

print("\n  normalized relation coefficients:", flush=True)
for label, coefficient in zip(basis_labels, null_display):
    print("    {:16s} {}".format(label, compact(coefficient)), flush=True)

relation_terms = [
    "({})*{}".format(compact(coefficient), label)
    for label, coefficient in zip(basis_labels, null_display)
    if coefficient != 0
]
relation_text = " + ".join(relation_terms) + " = 0"
print("  relation: {}".format(relation_text), flush=True)

# Rank jumps only within the certified range.
rank_profile = []
previous_rank = ZZ(-1)
for cutoff in range(q_min, q_cert + 1):
    rank_value = ZZ(stack_matrices([matrix_for_character(ch, cutoff) for ch in characters], 6).rank())
    if rank_value != previous_rank:
        rank_profile.append((ZZ(cutoff), rank_value))
        previous_rank = rank_value

# -----------------------------------------------------------------------------
# Part II. Construct the five-dimensional one-sided quotient and Sym^2 map.
# -----------------------------------------------------------------------------

print("\nPART II. FIVE-DIMENSIONAL CERTIFIED QUOTIENT", flush=True)
print("-"*79, flush=True)

pivot_indices = []
current_rank = ZZ(0)
for column_index in range(6):
    trial = pivot_indices + [column_index]
    trial_rank = ZZ(M_combined.matrix_from_columns(trial).rank())
    if trial_rank > current_rank:
        pivot_indices.append(column_index)
        current_rank = trial_rank
    if current_rank == 5:
        break

if len(pivot_indices) != 5:
    raise ArithmeticError("Could not select five independent certified columns")
dropped_index = [i for i in range(6) if i not in pivot_indices][0]
if null_vector[dropped_index] == 0:
    raise ArithmeticError("Dropped basis coefficient vanishes in relation")

replacement = vector(K, [
    -null_vector[index]/null_vector[dropped_index]
    for index in pivot_indices
])
pivot_labels = [basis_labels[index] for index in pivot_indices]
dropped_label = basis_labels[dropped_index]

print("  retained basis: {}".format(pivot_labels), flush=True)
print("  eliminated basis: {}".format(dropped_label), flush=True)

# h_original = T * g_quotient on the certified coefficient range.
T = matrix(K, 6, 5, 0)
for a, original_index in enumerate(pivot_indices):
    T[original_index, a] = K(1)
for a, coefficient in enumerate(replacement):
    T[dropped_index, a] = coefficient

embedding_checks = {}
for ch in characters:
    M6 = M_char[ch]
    M5 = M6.matrix_from_columns(pivot_indices)
    embedding_checks[ch] = bool(M6 == M5*T.transpose())
    print("  {:7s} certified embedding exact? {}".format(ch, embedding_checks[ch]), flush=True)
if not all(embedding_checks.values()):
    raise ArithmeticError("Certified one-sided quotient embedding failed")

pairs5 = [(i, j) for i in range(5) for j in range(i, 5)]
quotient_columns = []
for i, j in pairs6:
    C = matrix(K, 6, 6, 0)
    C[i, j] = K(1)
    C[j, i] = K(1)
    Ceff = T.transpose()*C*T
    quotient_columns.append(vector(K, [Ceff[a, b] for a, b in pairs5]))

Q = matrix(K, 15, 21, lambda row, col: quotient_columns[col][row])
Q_rank = ZZ(Q.rank())
Q_kernel = Q.right_kernel()
Q_kernel_matrix = Q_kernel.basis_matrix().echelon_form()
Q_nullity = ZZ(Q_kernel.dimension())
print("  symmetric quotient map: {} x {}".format(Q.nrows(), Q.ncols()), flush=True)
print("  quotient rank/nullity: {}/{}".format(Q_rank, Q_nullity), flush=True)
if Q_rank != 15 or Q_nullity != 6:
    raise ArithmeticError("Unexpected symmetric quotient dimensions")

# Exact right section S with Q*S=I_15.
pivot_columns_Q = list(Q.pivots())
if len(pivot_columns_Q) != 15:
    raise ArithmeticError("Could not find 15 pivot columns of Q")
Qp = Q.matrix_from_columns(pivot_columns_Q)
Qp_inv = Qp.inverse()
S = matrix(K, 21, 15, 0)
for local_row, original_row in enumerate(pivot_columns_Q):
    for col in range(15):
        S[original_row, col] = Qp_inv[local_row, col]
if Q*S != identity_matrix(K, 15):
    raise ArithmeticError("Failed to construct an exact right section of Q")

# -----------------------------------------------------------------------------
# Part III. Rebuild certified sector matrices and verify quotient factorization.
# -----------------------------------------------------------------------------

print("\nPART III. CERTIFIED SECTOR FACTORIZATION", flush=True)
print("-"*79, flush=True)

sector_specs = [
    ("(one,one)", "one", "one"),
    ("(one,omega)", "one", "omega"),
    ("(one,omega2)", "one", "omega2"),
    ("(omega,omega)", "omega", "omega"),
    ("(omega,omega2)", "omega", "omega2"),
    ("(omega2,omega2)", "omega2", "omega2"),
]


def sector_matrix(left, right):
    columns = []
    for i, j in pairs6:
        first = series_convolve(
            basis[left][basis_labels[i]],
            basis[right][basis_labels[j]],
            sector_min,
            q_cert,
        )
        if i == j:
            columns.append(first)
        else:
            second = series_convolve(
                basis[left][basis_labels[j]],
                basis[right][basis_labels[i]],
                sector_min,
                q_cert,
            )
            columns.append(series_add(first, second))
    return matrix(K, [
        [column.get(ZZ(exponent), K(0)) for column in columns]
        for exponent in range(sector_min, q_cert + 1)
    ])

sector_matrices21 = {}
sector_matrices15 = {}
sector_factorization = {}
sector_kernel_annihilation = {}
sector_rank21 = {}
sector_rank15 = {}

zero_kernel_target = None
for name, left, right in sector_specs:
    M21 = sector_matrix(left, right)
    annihilation = M21*Q_kernel_matrix.transpose()
    if zero_kernel_target is None or zero_kernel_target.nrows() != M21.nrows():
        zero_kernel_target = zero_matrix(K, M21.nrows(), Q_kernel_matrix.nrows())
    kills_kernel = bool(annihilation == zero_kernel_target)
    M15 = M21*S
    factors = bool(M21 == M15*Q)

    sector_matrices21[name] = M21
    sector_matrices15[name] = M15
    sector_kernel_annihilation[name] = kills_kernel
    sector_factorization[name] = factors
    sector_rank21[name] = ZZ(M21.rank())
    sector_rank15[name] = ZZ(M15.rank())

    print("  {:18s} kills ker(Q)? {}  factors? {}  ranks {}/{}".format(
        name, kills_kernel, factors, sector_rank21[name], sector_rank15[name]
    ), flush=True)

all_factor = all(sector_factorization.values())
print("  all six sector maps factor through Q? {}".format(all_factor), flush=True)

parked21 = stack_matrices([
    sector_matrices21["(omega,omega)"],
    sector_matrices21["(omega,omega2)"],
], 21)
parked15 = stack_matrices([
    sector_matrices15["(omega,omega)"],
    sector_matrices15["(omega,omega2)"],
], 15)
parked_rank21 = ZZ(parked21.rank())
parked_rank15 = ZZ(parked15.rank())
print("  parked rank in 21D coordinates: {} / 21".format(parked_rank21), flush=True)
print("  parked rank in 15D quotient: {} / 15".format(parked_rank15), flush=True)

# -----------------------------------------------------------------------------
# Part IV. Compare with the independently computed Cell-33 parked kernel.
# -----------------------------------------------------------------------------

print("\nPART IV. KERNEL IDENTIFICATION", flush=True)
print("-"*79, flush=True)

parked_kernel33 = data33["parked_kernel_matrix"]
if parked_kernel33.ncols() != 21:
    raise ArithmeticError("Unexpected Cell-33 kernel width")
union_rank = ZZ(matrix(K, list(parked_kernel33.rows()) + list(Q_kernel_matrix.rows())).rank())
kernels_equal = bool(
    parked_kernel33.nrows() == 6
    and Q_kernel_matrix.nrows() == 6
    and union_rank == 6
)
print("  Cell-33 parked kernel dimension: {}".format(parked_kernel33.nrows()), flush=True)
print("  quotient kernel dimension: {}".format(Q_kernel_matrix.nrows()), flush=True)
print("  rank of union of kernel bases: {}".format(union_rank), flush=True)
print("  kernels exactly equal through certified q^{}? {}".format(q_cert, kernels_equal), flush=True)

# -----------------------------------------------------------------------------
# Save results regardless of whether the factorization hypothesis passes.
# -----------------------------------------------------------------------------

RESULT_PATH = os.path.join(RESULTS_DIR, "degree3_certified_15d_quotient_v3.sobj")
SUMMARY_PATH = os.path.join(RESULTS_DIR, "degree3_certified_15d_quotient_v3.txt")

result = {
    "cell32_path": CELL32_PATH,
    "cell33_path": CELL33_PATH,
    "coefficient_field": K,
    "torsion_candidate": data32["torsion_candidate"],
    "certified_q_range": (q_min, q_cert),
    "working_buffer_range": (q_min, q_buffer_max),
    "unsafe_buffer_rank": buffer_rank,
    "first_unsafe_rank_jump": first_unsafe_jump,
    "basis_labels6": basis_labels,
    "one_sided_character_ranks": char_ranks,
    "one_sided_combined_rank": combined_rank,
    "one_sided_null_vector": null_vector,
    "one_sided_null_display": null_display,
    "one_sided_null_is_rational": rational_null,
    "relation_text": relation_text,
    "relation_checks": relation_checks,
    "rank_profile_certified": rank_profile,
    "pivot_indices": [ZZ(i) for i in pivot_indices],
    "pivot_labels5": pivot_labels,
    "dropped_index": ZZ(dropped_index),
    "dropped_label": dropped_label,
    "replacement_coefficients": replacement,
    "one_sided_embedding_T": T,
    "pairs6": pairs6,
    "pairs5": pairs5,
    "symmetric_quotient_Q": Q,
    "quotient_section_S": S,
    "quotient_rank": Q_rank,
    "quotient_nullity": Q_nullity,
    "sector_factorization": sector_factorization,
    "sector_kernel_annihilation": sector_kernel_annihilation,
    "sector_rank21": sector_rank21,
    "sector_rank15": sector_rank15,
    "parked_rank21": parked_rank21,
    "parked_rank15": parked_rank15,
    "kernels_equal": kernels_equal,
}
save(result, RESULT_PATH)

summary = [
    "DEGREE-THREE CERTIFIED-CUTOFF 15D QUOTIENT DIAGNOSTIC",
    "candidate: {}".format(data32["torsion_candidate"]),
    "certified q range: [{} , {}]".format(q_min, q_cert),
    "stored working buffer: [{} , {}]".format(q_min, q_buffer_max),
    "",
    "one-sided certified ranks:",
]
for ch in characters:
    summary.append("  {}: {} / 6".format(ch, char_ranks[ch]))
summary.extend([
    "  combined: {} / 6".format(combined_rank),
    "  common nullity: {}".format(combined_nullity),
    "  unsafe-buffer rank (diagnostic only): {}".format(buffer_rank),
    "  first unsafe rank jump: {}".format(first_unsafe_jump),
    "",
    "certified null relation:",
])
for label, coefficient in zip(basis_labels, null_display):
    summary.append("  {}: {}".format(label, compact(coefficient)))
summary.extend([
    "",
    "symmetric quotient:",
    "  rank: {} / 15".format(Q_rank),
    "  kernel dimension: {}".format(Q_nullity),
    "  all six sectors factor: {}".format(all_factor),
    "  parked rank in 21D: {} / 21".format(parked_rank21),
    "  parked rank in quotient: {} / 15".format(parked_rank15),
    "  quotient kernel equals Cell-33 parked kernel: {}".format(kernels_equal),
    "",
    "scope: exact through the certified finite q cutoff only;",
    "no all-orders one-sided identity is claimed by this cell.",
])
with open(SUMMARY_PATH, "w") as handle:
    handle.write("\n".join(summary) + "\n")

print("\nSaved:", flush=True)
print("  {}".format(RESULT_PATH), flush=True)
print("  {}".format(SUMMARY_PATH), flush=True)

if all_factor and parked_rank15 == 15 and kernels_equal:
    print("\nCERTIFIED FINITE-RANGE QUOTIENT: PASS", flush=True)
else:
    print("\nCERTIFIED FINITE-RANGE QUOTIENT: DID NOT FULLY PASS", flush=True)
    print("Inspect the saved summary before drawing a quotient conclusion.", flush=True)
