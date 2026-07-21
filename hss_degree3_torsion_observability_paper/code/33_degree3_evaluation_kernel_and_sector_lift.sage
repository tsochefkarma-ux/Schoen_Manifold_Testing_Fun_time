from sage.all import *
import os
import itertools

# =============================================================================
# CELL 33: DEGREE-THREE EVALUATION KERNEL AND SECTOR-LIFT CENSUS
# =============================================================================
#
# Loads the exact Cell-32 torsion-candidate series and:
#   1. reconstructs the 21-column shared symmetric evaluation matrices;
#   2. computes the exact six-dimensional kernel of the two parked sectors;
#   3. measures one-sided visibility at 1, omega, omega^2;
#   4. tests how the four remaining symmetric character sectors act on the
#      six invisible directions;
#   5. finds inclusion-minimal additional sector sets that lift the shared
#      coefficient rank as far as possible (and, if possible, to 21).
#
# All arithmetic is exact in the coefficient field stored by Cell 32.
# No floating-point arithmetic and no large formal-CM checkpoint are used.
# =============================================================================

print("="*79, flush=True)
print("DEGREE-THREE EVALUATION KERNEL AND SECTOR-LIFT CENSUS", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

INPUT_PATH = os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.sobj")
if not os.path.exists(INPUT_PATH):
    raise OSError("Missing Cell-32 result: {}".format(INPUT_PATH))

print("Loading {}".format(INPUT_PATH), flush=True)
data = load(INPUT_PATH)

K = data["coefficient_field"]
basis = data["basis"]
basis_labels = list(data["basis_labels"])
pairs = [tuple(pair) for pair in data["pairs"]]
q_max = ZZ(data["sector_q_max"])
characters = ["one", "omega", "omega2"]

if len(basis_labels) != 6 or len(pairs) != 21:
    raise ArithmeticError("Unexpected Cell-32 basis dimensions")

pair_labels = ["c{}{}".format(i+1, j+1) for i, j in pairs]

# -----------------------------------------------------------------------------
# Sparse Laurent-series helpers.
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


def support_min(series):
    return min(series.keys()) if series else None


all_one_sided_mins = [
    support_min(basis[ch][label])
    for ch in characters
    for label in basis_labels
    if basis[ch][label]
]
one_sided_min = min(all_one_sided_mins)
sector_min = ZZ(2*one_sided_min)


def one_sided_matrix(character):
    rows = []
    for exponent in range(one_sided_min, q_max + 1):
        rows.append([
            basis[character][label].get(ZZ(exponent), K(0))
            for label in basis_labels
        ])
    return matrix(K, rows)


def symmetric_sector_columns(left_character, right_character):
    columns = []
    for i, j in pairs:
        first = series_convolve(
            basis[left_character][basis_labels[i]],
            basis[right_character][basis_labels[j]],
            sector_min,
            q_max,
        )
        if i == j:
            columns.append(first)
        else:
            second = series_convolve(
                basis[left_character][basis_labels[j]],
                basis[right_character][basis_labels[i]],
                sector_min,
                q_max,
            )
            columns.append(series_add(first, second))
    return columns


def sector_matrix(left_character, right_character):
    columns = symmetric_sector_columns(left_character, right_character)
    rows = []
    for exponent in range(sector_min, q_max + 1):
        rows.append([
            column.get(ZZ(exponent), K(0))
            for column in columns
        ])
    return matrix(K, rows)


def stack_matrices(matrices):
    nonempty = [M for M in matrices if M.nrows() > 0]
    if not nonempty:
        return matrix(K, 0, 21)
    rows = []
    for M in nonempty:
        rows.extend(M.rows())
    return matrix(K, rows)


# -----------------------------------------------------------------------------
# Part I. One-sided visibility.
# -----------------------------------------------------------------------------

print("\nPART I. ONE-SIDED CHARACTER VISIBILITY", flush=True)
print("-"*79, flush=True)

one_sided_matrices = {ch: one_sided_matrix(ch) for ch in characters}
one_sided_ranks = {}
for ch in characters:
    one_sided_ranks[ch] = ZZ(one_sided_matrices[ch].rank())
    print("  {:7s} rank = {} / 6".format(ch, one_sided_ranks[ch]), flush=True)

combined_one_sided = stack_matrices([one_sided_matrices[ch] for ch in characters])
combined_one_sided_rank = ZZ(combined_one_sided.rank())
print("  all three characters combined rank = {} / 6".format(combined_one_sided_rank), flush=True)

one_sided_kernel_dimensions = {
    ch: ZZ(6 - one_sided_ranks[ch]) for ch in characters
}

# -----------------------------------------------------------------------------
# Part II. All six symmetric character sectors.
# -----------------------------------------------------------------------------

print("\nPART II. SYMMETRIC CHARACTER-SECTOR RANKS", flush=True)
print("-"*79, flush=True)

sector_specs = [
    ("(one,one)", "one", "one"),
    ("(one,omega)", "one", "omega"),
    ("(one,omega2)", "one", "omega2"),
    ("(omega,omega)", "omega", "omega"),
    ("(omega,omega2)", "omega", "omega2"),
    ("(omega2,omega2)", "omega2", "omega2"),
]

sector_matrices = {}
sector_ranks = {}
for name, left, right in sector_specs:
    M = sector_matrix(left, right)
    sector_matrices[name] = M
    sector_ranks[name] = ZZ(M.rank())
    print("  {:18s} rank = {} / 21".format(name, sector_ranks[name]), flush=True)

# -----------------------------------------------------------------------------
# Part III. Exact kernel of the two parked sectors.
# -----------------------------------------------------------------------------

print("\nPART III. EXACT PARKED-SECTOR KERNEL", flush=True)
print("-"*79, flush=True)

parked_names = ["(omega,omega)", "(omega,omega2)"]
parked_matrix = stack_matrices([sector_matrices[name] for name in parked_names])
parked_rank = ZZ(parked_matrix.rank())
parked_kernel = parked_matrix.right_kernel()
parked_kernel_matrix = parked_kernel.basis_matrix().echelon_form()
parked_nullity = ZZ(parked_kernel_matrix.nrows())

print("  parked combined rank = {} / 21".format(parked_rank), flush=True)
print("  exact kernel dimension = {}".format(parked_nullity), flush=True)
if parked_rank != ZZ(data["combined_rank"]):
    raise ArithmeticError("Reconstructed parked rank disagrees with Cell 32")
if parked_nullity != ZZ(data["combined_nullity"]):
    raise ArithmeticError("Reconstructed parked nullity disagrees with Cell 32")

kernel_descriptions = []
for row_index, row in enumerate(parked_kernel_matrix.rows()):
    support = [j for j, coefficient in enumerate(row) if coefficient != 0]
    terms = {pair_labels[j]: row[j] for j in support}
    kernel_descriptions.append({
        "index": ZZ(row_index),
        "support_indices": [ZZ(j) for j in support],
        "support_labels": [pair_labels[j] for j in support],
        "terms": terms,
    })
    print(
        "  kernel vector {:>2}: support size {:>2}; {}".format(
            row_index + 1,
            len(support),
            ", ".join(pair_labels[j] for j in support),
        ),
        flush=True,
    )

# Columns of kernel_embedding span the invisible coefficient directions.
kernel_embedding = parked_kernel_matrix.transpose()

# -----------------------------------------------------------------------------
# Part IV. How the remaining sectors lift the six-dimensional kernel.
# -----------------------------------------------------------------------------

print("\nPART IV. LIFT OF THE SIX INVISIBLE DIRECTIONS", flush=True)
print("-"*79, flush=True)

additional_names = [
    "(one,one)",
    "(one,omega)",
    "(one,omega2)",
    "(omega2,omega2)",
]

lift_ranks = {}
augmented_ranks = {}
for name in additional_names:
    restricted = sector_matrices[name] * kernel_embedding
    lift_ranks[name] = ZZ(restricted.rank())
    augmented = stack_matrices([parked_matrix, sector_matrices[name]])
    augmented_ranks[name] = ZZ(augmented.rank())
    print(
        "  {:18s} lifts {} / 6; total rank {} / 21".format(
            name,
            lift_ranks[name],
            augmented_ranks[name],
        ),
        flush=True,
    )

# Search every subset of the four additional sectors.
subset_records = []
for subset_size in range(len(additional_names) + 1):
    for subset in itertools.combinations(additional_names, subset_size):
        matrices = [parked_matrix] + [sector_matrices[name] for name in subset]
        total_rank = ZZ(stack_matrices(matrices).rank())
        restricted_matrices = [sector_matrices[name] * kernel_embedding for name in subset]
        lift_rank = ZZ(stack_matrices(restricted_matrices).rank()) if subset else ZZ(0)
        subset_records.append({
            "subset": tuple(subset),
            "size": ZZ(subset_size),
            "lift_rank": lift_rank,
            "total_rank": total_rank,
            "remaining_nullity": ZZ(21 - total_rank),
        })

max_total_rank = max(record["total_rank"] for record in subset_records)
full_rank_records = [record for record in subset_records if record["total_rank"] == 21]

if full_rank_records:
    minimum_full_size = min(record["size"] for record in full_rank_records)
    minimal_full_rank_sets = [
        record for record in full_rank_records if record["size"] == minimum_full_size
    ]
else:
    minimum_full_size = None
    minimal_full_rank_sets = []

maximal_records = [record for record in subset_records if record["total_rank"] == max_total_rank]
minimum_max_size = min(record["size"] for record in maximal_records)
minimal_max_rank_sets = [
    record for record in maximal_records if record["size"] == minimum_max_size
]

print("\n  maximum rank obtainable from all six sectors: {} / 21".format(max_total_rank), flush=True)
if minimal_full_rank_sets:
    print("  minimum additional sectors needed for full rank: {}".format(minimum_full_size), flush=True)
    for record in minimal_full_rank_sets:
        print("    {}".format(list(record["subset"])), flush=True)
else:
    print("  full rank 21 is not reached by these six candidate sectors", flush=True)
    print("  smallest additional sets attaining rank {}:".format(max_total_rank), flush=True)
    for record in minimal_max_rank_sets:
        print(
            "    {} (remaining nullity {})".format(
                list(record["subset"]),
                record["remaining_nullity"],
            ),
            flush=True,
        )

# -----------------------------------------------------------------------------
# Part V. Rank stabilization on the current six-dimensional kernel.
# -----------------------------------------------------------------------------

print("\nPART V. q-CUTOFF STABILIZATION OF KERNEL LIFT", flush=True)
print("-"*79, flush=True)

# For each additional sector, record only cutoffs at which its restriction to
# the parked kernel gains rank.  This distinguishes genuine invisibility from
# an insufficient q truncation.
lift_profiles = {}
for name in additional_names:
    M = sector_matrices[name]
    profile = []
    previous = ZZ(-1)
    for row_count in range(1, M.nrows() + 1):
        value = ZZ((M[:row_count, :] * kernel_embedding).rank())
        if value != previous:
            exponent = ZZ(sector_min + row_count - 1)
            profile.append((exponent, value))
            previous = value
    lift_profiles[name] = profile
    print("  {}:".format(name), flush=True)
    for exponent, value in profile:
        print("    through q^{:>3}: lift rank {}".format(exponent, value), flush=True)

# -----------------------------------------------------------------------------
# Save exact data and compact summary.
# -----------------------------------------------------------------------------

RESULT_PATH = os.path.join(RESULTS_DIR, "degree3_evaluation_kernel_and_sector_lift.sobj")
SUMMARY_PATH = os.path.join(RESULTS_DIR, "degree3_evaluation_kernel_and_sector_lift.txt")

result = {
    "input_path": INPUT_PATH,
    "coefficient_field": K,
    "q_range": (sector_min, q_max),
    "basis_labels": basis_labels,
    "pairs": pairs,
    "pair_labels": pair_labels,
    "one_sided_ranks": one_sided_ranks,
    "one_sided_kernel_dimensions": one_sided_kernel_dimensions,
    "combined_one_sided_rank": combined_one_sided_rank,
    "sector_ranks": sector_ranks,
    "parked_sector_names": parked_names,
    "parked_rank": parked_rank,
    "parked_nullity": parked_nullity,
    "parked_kernel_matrix": parked_kernel_matrix,
    "kernel_descriptions": kernel_descriptions,
    "additional_sector_names": additional_names,
    "lift_ranks": lift_ranks,
    "augmented_ranks": augmented_ranks,
    "subset_records": subset_records,
    "max_total_rank": max_total_rank,
    "minimal_full_rank_sets": minimal_full_rank_sets,
    "minimal_max_rank_sets": minimal_max_rank_sets,
    "lift_profiles": lift_profiles,
}
save(result, RESULT_PATH)

summary_lines = [
    "DEGREE-THREE EVALUATION KERNEL AND SECTOR-LIFT CENSUS",
    "coefficient field: {}".format(K),
    "q range: [{} , {}]".format(sector_min, q_max),
    "",
    "one-sided ranks:",
]
for ch in characters:
    summary_lines.append("  {}: {} / 6".format(ch, one_sided_ranks[ch]))
summary_lines.extend([
    "  all three combined: {} / 6".format(combined_one_sided_rank),
    "",
    "symmetric sector ranks:",
])
for name, _, _ in sector_specs:
    summary_lines.append("  {}: {} / 21".format(name, sector_ranks[name]))
summary_lines.extend([
    "",
    "parked sectors:",
    "  rank: {} / 21".format(parked_rank),
    "  exact kernel dimension: {}".format(parked_nullity),
    "",
    "additional-sector lift:",
])
for name in additional_names:
    summary_lines.append(
        "  {}: lifts {} / 6; total rank {} / 21".format(
            name,
            lift_ranks[name],
            augmented_ranks[name],
        )
    )
summary_lines.extend([
    "",
    "maximum rank from all six candidate sectors: {} / 21".format(max_total_rank),
])
if minimal_full_rank_sets:
    summary_lines.append(
        "minimum additional sectors for full rank: {}".format(minimum_full_size)
    )
    for record in minimal_full_rank_sets:
        summary_lines.append("  {}".format(list(record["subset"])))
else:
    summary_lines.append("full rank 21 not reached")
    summary_lines.append("smallest sets attaining maximum rank:")
    for record in minimal_max_rank_sets:
        summary_lines.append(
            "  {} (remaining nullity {})".format(
                list(record["subset"]),
                record["remaining_nullity"],
            )
        )
summary_lines.extend([
    "",
    "kernel basis supports:",
])
for description in kernel_descriptions:
    summary_lines.append(
        "  k{}: {}".format(
            description["index"] + 1,
            ", ".join(description["support_labels"]),
        )
    )
summary_lines.extend([
    "",
    "scope: exact for the Cell-32 canonical gamma/3 torsion candidate;",
    "geometric identification of the character labels remains a separate check.",
])

with open(SUMMARY_PATH, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

print("\nSaved:", flush=True)
print("  {}".format(RESULT_PATH), flush=True)
print("  {}".format(SUMMARY_PATH), flush=True)
print("\nCELL 33 COMPLETE", flush=True)
