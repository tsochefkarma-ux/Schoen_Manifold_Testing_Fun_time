# Cell 40 — exact degree-three reconstruction atlas
#
# Purpose:
#   Work in the audited 15-dimensional observable quotient and determine:
#     * which single character sector, when combined with the ordinary
#       torsion-summed sector (one,one), gives full rank 15;
#     * the earliest certified q cutoff at which full rank is reached;
#     * an exact set of 9 torsion-summed coefficients plus 6 coefficients
#       from the additional sector forming an invertible 15x15 system;
#     * the exact inverse reconstruction matrix c = R y;
#     * which pairs of torsion-class Fourier projectors give full rank.
#
# This cell does not insert enumerative values.  It produces the exact linear
# maps that will reconstruct the 15 observable coefficients once the selected
# q-series data are supplied.

from sage.all import *
import os
import itertools

ROOT = os.getcwd()
RESULTS_DIR = os.path.join(ROOT, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

CELL32_PATH = os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.sobj")
CELL34_PATH = os.path.join(RESULTS_DIR, "degree3_certified_15d_quotient_v3.sobj")
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_exact_reconstruction_atlas_cell40.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_exact_reconstruction_atlas_cell40.txt")

for path in [CELL32_PATH, CELL34_PATH]:
    if not os.path.exists(path) and not os.path.exists(path + ".sobj"):
        raise IOError("Missing required input: {}".format(path))

D32 = load(CELL32_PATH)
D34 = load(CELL34_PATH)

K = D32["coefficient_field"]
omega = K(D32["omega"])
characters = ["one", "omega", "omega2"]

certified_q_max = ZZ(D32["sector_q_max"])
if "certified_q_max" in D34:
    certified_q_max = min(certified_q_max, ZZ(D34["certified_q_max"]))

basis_all = D32["basis"]
labels5 = list(D34["pivot_labels5"])
if len(labels5) != 5:
    raise ArithmeticError("Expected five reduced one-sided basis labels")

pairs5 = [(i, j) for i in range(5) for j in range(i, 5)]
pair_labels15 = ["{}*{}".format(labels5[i], labels5[j]) for i, j in pairs5]


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


one_sided_min = min(
    series_min(basis_all[ch][label])
    for ch in characters for label in labels5
)
sector_min = ZZ(2 * one_sided_min)
exponents = list(range(sector_min, certified_q_max + 1))


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
    return matrix(K, [
        [sget(column, exponent) for column in columns]
        for exponent in exponents
    ])


def vstack(matrices):
    matrices = [M for M in matrices if M is not None and M.nrows() > 0]
    if not matrices:
        return matrix(K, 0, 15)
    rows = []
    for M in matrices:
        rows.extend([list(row) for row in M.rows()])
    return matrix(K, rows)


def rows_through(M, qcut):
    if qcut < sector_min:
        return matrix(K, 0, M.ncols())
    count = min(M.nrows(), ZZ(qcut - sector_min + 1))
    return M.matrix_from_rows(range(count))


def greedy_independent_rows(M, row_labels, initial_rows=None, target_rank=None):
    selected_rows = [] if initial_rows is None else [list(row) for row in initial_rows]
    selected_labels = []
    current = matrix(K, selected_rows) if selected_rows else matrix(K, 0, M.ncols())
    current_rank = ZZ(current.rank())
    if target_rank is None:
        target_rank = ZZ(M.ncols())
    for idx in range(M.nrows()):
        candidate = vstack([current, matrix(K, [list(M.row(idx))])])
        new_rank = ZZ(candidate.rank())
        if new_rank > current_rank:
            selected_rows.append(list(M.row(idx)))
            selected_labels.append(row_labels[idx])
            current = candidate
            current_rank = new_rank
            if current_rank >= target_rank:
                break
    return current, selected_labels, current_rank


def earliest_full_cutoff(M1, M2, target=15):
    for qcut in exponents:
        if vstack([rows_through(M1, qcut), rows_through(M2, qcut)]).rank() >= target:
            return ZZ(qcut)
    return None


def linear_combination(weighted_matrices):
    result = matrix(K, len(exponents), 15, 0)
    for weight, M in weighted_matrices:
        result += K(weight) * M
    return result


print("=" * 79)
print("CELL 40: DEGREE-THREE EXACT RECONSTRUCTION ATLAS")
print("=" * 79)
print("  coefficient field       : {}".format(K))
print("  certified q range       : [{} , {}]".format(sector_min, certified_q_max))
print("  observable dimension    : 15")
print("  reduced one-sided basis : {}".format(labels5))

ordered = {(left, right): sector_matrix(left, right)
           for left in characters for right in characters}

# Exchange audit.
for left in characters:
    for right in characters:
        if ordered[(left, right)] != ordered[(right, left)]:
            raise ArithmeticError("Exchange audit failed for ({},{})".format(left, right))

base_pair = ("one", "one")
base_name = "(one,one)"
base = ordered[base_pair]
base_rank = ZZ(base.rank())

# Select an early exact basis for the nine-dimensional torsion-summed rowspace.
base_row_labels = [(base_name, ZZ(e)) for e in exponents]
base_basis_matrix, base_selected_labels, base_basis_rank = greedy_independent_rows(
    base, base_row_labels, target_rank=base_rank
)
if base_basis_rank != base_rank:
    raise ArithmeticError("Failed to extract a basis of the torsion-summed rowspace")

print("\nPART I. TORSION-SUMMED BASE")
print("-" * 79)
print("  rank                         : {} / 15".format(base_rank))
print("  independent coefficients used: {}".format(len(base_selected_labels)))
print("  selected q exponents         : {}".format([e for _, e in base_selected_labels]))

sector_reps = [
    ("one", "omega"),
    ("one", "omega2"),
    ("omega", "omega"),
    ("omega", "omega2"),
    ("omega2", "omega2"),
]

sector_atlas = {}
full_rank_sector_names = []

print("\nPART II. ONE-EXTRA-SECTOR RECONSTRUCTION MAPS")
print("-" * 79)
for pair in sector_reps:
    name = "({},{})".format(pair[0], pair[1])
    M = ordered[pair]
    combined_rank = ZZ(vstack([base, M]).rank())
    increment = combined_rank - base_rank
    qcut = earliest_full_cutoff(base, M, target=15) if combined_rank == 15 else None

    extra_row_labels = [(name, ZZ(e)) for e in exponents]
    combined_basis, extra_selected_labels, final_rank = greedy_independent_rows(
        M,
        extra_row_labels,
        initial_rows=[list(row) for row in base_basis_matrix.rows()],
        target_rank=15,
    )

    reconstruction = None
    inverse_verified = False
    if final_rank == 15:
        if combined_basis.nrows() != 15 or combined_basis.ncols() != 15:
            raise ArithmeticError("Expected a square 15x15 reconstruction matrix")
        reconstruction = combined_basis.inverse()
        inverse_verified = bool(
            combined_basis * reconstruction == identity_matrix(K, 15)
            and reconstruction * combined_basis == identity_matrix(K, 15)
        )
        if not inverse_verified:
            raise ArithmeticError("Exact inverse verification failed for {}".format(name))
        full_rank_sector_names.append(name)

    sector_atlas[name] = {
        "pair": pair,
        "sector_rank": ZZ(M.rank()),
        "combined_rank_with_torsion_sum": combined_rank,
        "rank_increment": increment,
        "earliest_full_q_cutoff": qcut,
        "base_selected_rows": list(base_selected_labels),
        "extra_selected_rows": list(extra_selected_labels),
        "extra_coefficient_count": ZZ(len(extra_selected_labels)),
        "system_matrix": combined_basis if final_rank == 15 else None,
        "reconstruction_inverse": reconstruction,
        "inverse_verified": inverse_verified,
    }

    print("  {:18s}: sector rank {:2d}, combined {:2d}, increment {:2d}".format(
        name, ZZ(M.rank()), combined_rank, increment
    ))
    if final_rank == 15:
        print("    full reconstruction       : YES")
        print("    earliest full q cutoff    : q^{}".format(qcut))
        print("    extra coefficients needed : {}".format(len(extra_selected_labels)))
        print("    extra q exponents         : {}".format([e for _, e in extra_selected_labels]))
        print("    exact inverse verified    : {}".format(inverse_verified))
    else:
        print("    full reconstruction       : NO")

# Choose the cheapest viable sector by earliest cutoff, then lexical name.
viable = [
    (record["earliest_full_q_cutoff"], name)
    for name, record in sector_atlas.items()
    if record["combined_rank_with_torsion_sum"] == 15
]
preferred_sector = None
if viable:
    preferred_sector = sorted(viable, key=lambda item: (item[0], item[1]))[0][1]

print("\nPART III. PREFERRED CHARACTER-SECTOR COMPLETION")
print("-" * 79)
print("  viable single extra sectors: {}".format(full_rank_sector_names))
print("  preferred by earliest q depth: {}".format(preferred_sector))
if preferred_sector is not None:
    rec = sector_atlas[preferred_sector]
    print("  torsion-summed coefficients : {}".format(rec["base_selected_rows"]))
    print("  extra-sector coefficients   : {}".format(rec["extra_selected_rows"]))
    print("  reconstruction convention   : c = reconstruction_inverse * y")

# -----------------------------------------------------------------------------
# Exact inverse-Fourier torsion-class projectors and full-rank pairs.
# -----------------------------------------------------------------------------
projectors = {}
for a in range(3):
    for b in range(3):
        weighted = []
        for r, left in enumerate(characters):
            for s, right in enumerate(characters):
                weighted.append((omega ** ZZ(-(r*a + s*b)), ordered[(left, right)]))
        projectors[(ZZ(a), ZZ(b))] = linear_combination(weighted)

class_reps = [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
class_pair_atlas = []
for p1, p2 in itertools.combinations(class_reps, 2):
    M1 = projectors[(ZZ(p1[0]), ZZ(p1[1]))]
    M2 = projectors[(ZZ(p2[0]), ZZ(p2[1]))]
    rank = ZZ(vstack([M1, M2]).rank())
    if rank == 15:
        qcut = earliest_full_cutoff(M1, M2, target=15)
        class_pair_atlas.append({
            "classes": ("N_{}_{}".format(*p1), "N_{}_{}".format(*p2)),
            "earliest_full_q_cutoff": qcut,
            "rank": rank,
        })

class_pair_atlas.sort(key=lambda r: (r["earliest_full_q_cutoff"], r["classes"]))

print("\nPART IV. TWO TORSION-CLASS SERIES SUFFICE")
print("-" * 79)
print("  full-rank class pairs: {}".format(len(class_pair_atlas)))
for record in class_pair_atlas[:20]:
    print("  {} + {} : rank 15 by q^{}".format(
        record["classes"][0], record["classes"][1], record["earliest_full_q_cutoff"]
    ))
if len(class_pair_atlas) > 20:
    print("  ... {} additional full-rank pairs stored in result object".format(
        len(class_pair_atlas) - 20
    ))

# -----------------------------------------------------------------------------
# Save exact atlas.
# -----------------------------------------------------------------------------
result = {
    "schema_version": ZZ(1),
    "scope": "exact reconstruction atlas in the audited 15D degree-three HSS torsion quotient",
    "coefficient_field": K,
    "certified_q_range": (sector_min, certified_q_max),
    "basis_labels5": labels5,
    "pair_labels15": pair_labels15,
    "torsion_summed_rank": base_rank,
    "torsion_summed_selected_rows": list(base_selected_labels),
    "torsion_summed_basis_matrix": base_basis_matrix,
    "sector_atlas": sector_atlas,
    "viable_single_extra_sectors": full_rank_sector_names,
    "preferred_extra_sector": preferred_sector,
    "torsion_class_full_rank_pairs": class_pair_atlas,
}
save(result, OUT_SOBJ)

summary = [
    "CELL 40: DEGREE-THREE EXACT RECONSTRUCTION ATLAS",
    "certified q range: [{} , {}]".format(sector_min, certified_q_max),
    "observable dimension: 15",
    "torsion-summed rank: {} / 15".format(base_rank),
    "missing directions after torsion sum: {}".format(15 - base_rank),
    "viable single extra character sectors: {}".format(full_rank_sector_names),
    "preferred extra sector by earliest q depth: {}".format(preferred_sector),
    "two torsion-class series can reconstruct: {}".format(bool(class_pair_atlas)),
    "interpretation: one independent additional sector supplies the six directions invisible to the torsion sum",
    "caution: the atlas proves sufficiency of data, not that the additional sector is already known independently",
]
with open(OUT_TXT, "w") as handle:
    handle.write("\n".join(summary) + "\n")

print("\nPART V. VERDICT")
print("-" * 79)
print("  torsion sum fixes                 : {} / 15".format(base_rank))
print("  one viable extra sector fixes     : remaining {} directions".format(15 - base_rank))
print("  exact reconstruction maps saved   : {}".format(bool(full_rank_sector_names)))
print("  this is a data-sufficiency theorem: not yet a derivation from transport")
print("\nSaved:")
print("  {}".format(OUT_SOBJ))
print("  {}".format(OUT_TXT))
