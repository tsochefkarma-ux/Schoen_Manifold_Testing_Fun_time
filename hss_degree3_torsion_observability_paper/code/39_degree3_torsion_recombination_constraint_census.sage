# Cell 39 — degree-three torsion recombination / constraint-budget census
#
# Purpose:
#   Work entirely in the audited 15-dimensional observable quotient.
#   Rebuild the six symmetric character-evaluation sector matrices from the
#   exact Cell-32 q-series, then measure:
#     * the rank of the trivial-character (torsion-summed) series;
#     * the rank of every character sector and every torsion-class Fourier
#       projector;
#     * the smallest character-sector sets reaching full rank 15;
#     * the absolute ceiling after adding the rank-two scalar transport data.
#
# This is exact cyclotomic linear algebra.  It does not assume a new closure.

from sage.all import *
import os
import itertools

ROOT = os.getcwd()
RESULTS_DIR = os.path.join(ROOT, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

CELL32_PATH = os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.sobj")
CELL34_PATH = os.path.join(RESULTS_DIR, "degree3_certified_15d_quotient_v3.sobj")
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_torsion_recombination_constraint_census_cell39.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_torsion_recombination_constraint_census_cell39.txt")

for path in [CELL32_PATH, CELL34_PATH]:
    if not os.path.exists(path + ".sobj") and not os.path.exists(path):
        raise IOError("Missing required input: {}".format(path))

data32 = load(CELL32_PATH)
data34 = load(CELL34_PATH)

K = data32["coefficient_field"]
omega = K(data32["omega"])
characters = ["one", "omega", "omega2"]
char_index = {name: ZZ(i) for i, name in enumerate(characters)}

certified_q_max = ZZ(data32["sector_q_max"])
if "certified_q_max" in data34:
    certified_q_max = min(certified_q_max, ZZ(data34["certified_q_max"]))

basis_all = data32["basis"]
labels5 = list(data34["pivot_labels5"])
if len(labels5) != 5:
    raise ArithmeticError("Expected five reduced one-sided basis labels")

pairs5 = [(i, j) for i in range(5) for j in range(i, 5)]
pair_labels = ["{}*{}".format(labels5[i], labels5[j]) for i, j in pairs5]
if len(pairs5) != 15:
    raise ArithmeticError("Expected fifteen symmetric quotient coordinates")


def sget(series, exponent):
    return K(series.get(ZZ(exponent), K(0)))


def series_min(series):
    nz = [ZZ(e) for e, value in series.items() if value != 0]
    return min(nz) if nz else ZZ(0)


def convolve(A, B, exp_min, exp_max):
    out = {}
    keysA = [(ZZ(e), K(v)) for e, v in A.items() if v != 0]
    keysB = [(ZZ(e), K(v)) for e, v in B.items() if v != 0]
    for ea, va in keysA:
        for eb, vb in keysB:
            e = ea + eb
            if e < exp_min or e > exp_max:
                continue
            out[e] = out.get(e, K(0)) + va * vb
    return {e: v for e, v in out.items() if v != 0}


def add_series(A, B, scale=K(1)):
    out = dict(A)
    for e, v in B.items():
        out[e] = out.get(e, K(0)) + scale * K(v)
        if out[e] == 0:
            del out[e]
    return out


# A safe common lower bound for all two-sided products.
one_sided_min = min(
    series_min(basis_all[ch][label])
    for ch in characters for label in labels5
)
sector_min = ZZ(2 * one_sided_min)


def symmetric_column(left, right, i, j):
    A = basis_all[left][labels5[i]]
    B = basis_all[right][labels5[j]]
    first = convolve(A, B, sector_min, certified_q_max)
    if i == j:
        return first
    C = basis_all[left][labels5[j]]
    D = basis_all[right][labels5[i]]
    second = convolve(C, D, sector_min, certified_q_max)
    return add_series(first, second)


def sector_matrix(left, right):
    columns = [symmetric_column(left, right, i, j) for i, j in pairs5]
    rows = []
    for exponent in range(sector_min, certified_q_max + 1):
        rows.append([sget(column, exponent) for column in columns])
    return matrix(K, rows)


def vstack(matrices):
    matrices = [M for M in matrices if M is not None and M.nrows() > 0]
    if not matrices:
        return matrix(K, 0, 15)
    rows = []
    for M in matrices:
        rows.extend([list(row) for row in M.rows()])
    return matrix(K, rows)


def linear_combination(weighted_matrices):
    result = matrix(K, certified_q_max - sector_min + 1, 15, 0)
    for weight, M in weighted_matrices:
        result += K(weight) * M
    return result


print("=" * 79)
print("CELL 39: DEGREE-THREE TORSION-RECOMBINATION CONSTRAINT CENSUS")
print("=" * 79)
print("  coefficient field      : {}".format(K))
print("  certified q range      : [{} , {}]".format(sector_min, certified_q_max))
print("  reduced one-sided basis: {}".format(labels5))
print("  observable coefficients: 15")

# -----------------------------------------------------------------------------
# Part I. All ordered character-evaluation matrices.
# -----------------------------------------------------------------------------
ordered = {}
for left in characters:
    for right in characters:
        ordered[(left, right)] = sector_matrix(left, right)

# Exact symmetry audit E_(r,s)=E_(s,r) for the symmetric potential.
symmetry_checks = {}
for left in characters:
    for right in characters:
        symmetry_checks[(left, right)] = bool(
            ordered[(left, right)] == ordered[(right, left)]
        )
if not all(symmetry_checks.values()):
    raise ArithmeticError("Symmetric sector matrices failed exchange audit")

sector_reps = [
    ("one", "one"),
    ("one", "omega"),
    ("one", "omega2"),
    ("omega", "omega"),
    ("omega", "omega2"),
    ("omega2", "omega2"),
]
sector_names = {
    pair: "({},{})".format(pair[0], pair[1]) for pair in sector_reps
}
sector_ranks = {sector_names[pair]: ZZ(ordered[pair].rank()) for pair in sector_reps}

print("\nPART I. CHARACTER-EVALUATION SECTOR RANKS")
print("-" * 79)
for pair in sector_reps:
    print("  {:18s}: {} / 15".format(sector_names[pair], sector_ranks[sector_names[pair]]))

# The ordinary torsion-summed series is evaluation at trivial characters.
torsion_summed_matrix = ordered[("one", "one")]
torsion_summed_rank = ZZ(torsion_summed_matrix.rank())
scalar_transport_ceiling = ZZ(2)
torsion_plus_scalar_upper_bound = min(ZZ(15), torsion_summed_rank + scalar_transport_ceiling)
torsion_plus_scalar_remaining_lower_bound = ZZ(15) - torsion_plus_scalar_upper_bound

print("\nPART II. PUBLISHED TORSION-SUMMED CONSTRAINT BUDGET")
print("-" * 79)
print("  torsion-summed evaluation        : (one,one)")
print("  exact q-series map rank          : {} / 15".format(torsion_summed_rank))
print("  scalar-transport absolute ceiling: +{}".format(scalar_transport_ceiling))
print("  best possible combined rank      : <= {} / 15".format(torsion_plus_scalar_upper_bound))
print("  unavoidable remaining dimension  : >= {}".format(torsion_plus_scalar_remaining_lower_bound))

# -----------------------------------------------------------------------------
# Part III. Minimal sets of character-evaluation sectors reaching full rank.
# -----------------------------------------------------------------------------
subset_records = []
for size in range(1, len(sector_reps) + 1):
    for subset in itertools.combinations(sector_reps, size):
        rank = ZZ(vstack([ordered[pair] for pair in subset]).rank())
        subset_records.append({
            "subset": [sector_names[pair] for pair in subset],
            "rank": rank,
            "size": ZZ(size),
        })

full_rank_records = [record for record in subset_records if record["rank"] == 15]
minimal_character_sector_count = min(
    [record["size"] for record in full_rank_records], default=None
)
minimal_character_sector_sets = [
    record["subset"] for record in full_rank_records
    if record["size"] == minimal_character_sector_count
]

# Minimal additions when the torsion-summed sector is already known.
base_pair = ("one", "one")
additional_reps = [pair for pair in sector_reps if pair != base_pair]
base_augmented = []
for size in range(0, len(additional_reps) + 1):
    for subset in itertools.combinations(additional_reps, size):
        matrices = [ordered[base_pair]] + [ordered[pair] for pair in subset]
        rank = ZZ(vstack(matrices).rank())
        base_augmented.append({
            "additional": [sector_names[pair] for pair in subset],
            "rank": rank,
            "additional_count": ZZ(size),
        })

base_full = [record for record in base_augmented if record["rank"] == 15]
minimal_additional_count = min(
    [record["additional_count"] for record in base_full], default=None
)
minimal_additional_sets = [
    record["additional"] for record in base_full
    if record["additional_count"] == minimal_additional_count
]

print("\nPART III. MINIMAL CHARACTER DATA FOR FULL OBSERVABILITY")
print("-" * 79)
print("  minimal number of character sectors for rank 15: {}".format(minimal_character_sector_count))
for subset in minimal_character_sector_sets:
    print("    {}".format(subset))
print("  with (one,one) already known, extra sectors needed: {}".format(minimal_additional_count))
for subset in minimal_additional_sets:
    print("    {}".format(subset))

# -----------------------------------------------------------------------------
# Part IV. Exact Z3 x Z3 Fourier projectors (torsion-class series).
# -----------------------------------------------------------------------------
# Evaluation convention:
#   E_(r,s) = sum_(a,b) N_(a,b) omega^(r a + s b).
# Hence the unnormalised inverse Fourier projector is
#   P_(a,b) = sum_(r,s) omega^(-r a - s b) E_(r,s) = 9 N_(a,b).
projectors = {}
projector_ranks = {}
for a in range(3):
    for b in range(3):
        weighted = []
        for r, left in enumerate(characters):
            for s, right in enumerate(characters):
                weight = omega ** ZZ(-(r * a + s * b))
                weighted.append((weight, ordered[(left, right)]))
        Pab = linear_combination(weighted)
        projectors[(ZZ(a), ZZ(b))] = Pab
        projector_ranks[(ZZ(a), ZZ(b))] = ZZ(Pab.rank())

print("\nPART IV. TORSION-CLASS FOURIER-PROJECTOR RANKS")
print("-" * 79)
for a in range(3):
    for b in range(3):
        print("  class ({},{}) : {} / 15".format(a, b, projector_ranks[(ZZ(a), ZZ(b))]))

# Orbit representatives under exchange (a,b)<->(b,a).
class_reps = [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
class_subset_records = []
for size in range(1, len(class_reps) + 1):
    for subset in itertools.combinations(class_reps, size):
        rank = ZZ(vstack([projectors[(ZZ(a),ZZ(b))] for a,b in subset]).rank())
        class_subset_records.append({
            "subset": ["N_{}_{}".format(a,b) for a,b in subset],
            "rank": rank,
            "size": ZZ(size),
        })
class_full = [record for record in class_subset_records if record["rank"] == 15]
minimal_class_count = min([record["size"] for record in class_full], default=None)
minimal_class_sets = [
    record["subset"] for record in class_full if record["size"] == minimal_class_count
]

print("  minimal torsion-class series needed for rank 15: {}".format(minimal_class_count))
for subset in minimal_class_sets:
    print("    {}".format(subset))

# -----------------------------------------------------------------------------
# Part V. Natural orbit sums.
# -----------------------------------------------------------------------------
diagonal_character_sum = linear_combination([
    (1, ordered[("one","one")]),
    (1, ordered[("omega","omega")]),
    (1, ordered[("omega2","omega2")]),
])
offdiagonal_character_sum = linear_combination([
    (1, ordered[("one","omega")]),
    (1, ordered[("omega","one")]),
    (1, ordered[("one","omega2")]),
    (1, ordered[("omega2","one")]),
    (1, ordered[("omega","omega2")]),
    (1, ordered[("omega2","omega")]),
])
all_character_sum = diagonal_character_sum + offdiagonal_character_sum

orbit_sum_ranks = {
    "diagonal_character_sum": ZZ(diagonal_character_sum.rank()),
    "offdiagonal_character_sum": ZZ(offdiagonal_character_sum.rank()),
    "all_character_sum": ZZ(all_character_sum.rank()),
    "diagonal_plus_offdiagonal": ZZ(vstack([diagonal_character_sum, offdiagonal_character_sum]).rank()),
}

print("\nPART V. NATURAL CHARACTER-ORBIT SUMS")
print("-" * 79)
for name, rank in orbit_sum_ranks.items():
    print("  {:30s}: {} / 15".format(name, rank))

# -----------------------------------------------------------------------------
# Save.
# -----------------------------------------------------------------------------
verdict = {
    "observable_dimension": ZZ(15),
    "torsion_summed_rank": torsion_summed_rank,
    "scalar_transport_rank_ceiling": scalar_transport_ceiling,
    "torsion_plus_scalar_rank_upper_bound": torsion_plus_scalar_upper_bound,
    "remaining_dimension_lower_bound": torsion_plus_scalar_remaining_lower_bound,
    "minimal_character_sector_count_for_full_rank": minimal_character_sector_count,
    "minimal_additional_character_count_given_torsion_sum": minimal_additional_count,
    "minimal_torsion_class_count_for_full_rank": minimal_class_count,
}

result = {
    "schema_version": ZZ(1),
    "scope": "exact degree-three torsion recombination rank census in the audited 15D HSS quotient",
    "coefficient_field": K,
    "certified_q_range": (sector_min, certified_q_max),
    "basis_labels5": labels5,
    "pair_labels15": pair_labels,
    "sector_ranks": sector_ranks,
    "symmetry_checks": symmetry_checks,
    "torsion_summed_rank": torsion_summed_rank,
    "subset_records": subset_records,
    "minimal_character_sector_count": minimal_character_sector_count,
    "minimal_character_sector_sets": minimal_character_sector_sets,
    "base_augmented_records": base_augmented,
    "minimal_additional_count_given_torsion_sum": minimal_additional_count,
    "minimal_additional_sets_given_torsion_sum": minimal_additional_sets,
    "projector_ranks": projector_ranks,
    "class_subset_records": class_subset_records,
    "minimal_class_count": minimal_class_count,
    "minimal_class_sets": minimal_class_sets,
    "orbit_sum_ranks": orbit_sum_ranks,
    "verdict": verdict,
}
save(result, OUT_SOBJ)

summary_lines = [
    "CELL 39: DEGREE-THREE TORSION-RECOMBINATION CONSTRAINT CENSUS",
    "certified q range: [{} , {}]".format(sector_min, certified_q_max),
    "observable quotient dimension: 15",
    "torsion-summed sector (one,one) rank: {} / 15".format(torsion_summed_rank),
    "scalar transport absolute rank ceiling: 2",
    "torsion sum + scalar transport best possible rank: <= {} / 15".format(torsion_plus_scalar_upper_bound),
    "remaining dimension after those inputs: >= {}".format(torsion_plus_scalar_remaining_lower_bound),
    "minimal character sectors for full rank: {}".format(minimal_character_sector_count),
    "minimal extra sectors given (one,one): {}".format(minimal_additional_count),
    "minimal torsion-class series for full rank: {}".format(minimal_class_count),
    "natural orbit-sum ranks: {}".format(orbit_sum_ranks),
    "conclusion: quantify external enumerative data before proposing a new transport law",
]
with open(OUT_TXT, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

print("\nPART VI. CONSTRAINT-BUDGET VERDICT")
print("-" * 79)
print("  torsion-summed data alone fixes at most {} of 15 directions".format(torsion_summed_rank))
print("  adding the present scalar transport can fix at most {} of 15".format(torsion_plus_scalar_upper_bound))
print("  at least {} directions require other independent input".format(torsion_plus_scalar_remaining_lower_bound))
print("\nSaved:")
print("  {}".format(OUT_SOBJ))
print("  {}".format(OUT_TXT))
