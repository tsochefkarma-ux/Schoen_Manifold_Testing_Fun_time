from sage.all import *
from itertools import combinations
import glob
import os

# ============================================================================
# CELL 28F: AGGREGATE SPECIALIZE-FIRST WITNESSES
# ============================================================================
# Uses exact good-reduction witnesses produced by Cell 28E.
#
# Certification rules:
#   * a nonzero m x m coefficient minor at one good specialization proves
#     generic row rank m; for m <= 3 this makes the affine system surjective
#     and hence consistent, with geometry-family dimension 3-m;
#   * a nonzero 4 x 4 augmented minor proves inconsistency because the
#     coefficient matrix has only three geometry columns;
#   * beta_-3 = rho3 alpha_-3 is an exact characteristic-zero relation already
#     certified in Cell 27, so when both cubic rows occur we remove beta_-3
#     before applying the rules above.
# ============================================================================

INPUT_GLOB = "results/28e_specialize_first/point_*.sobj"
OUTPUT_TEXT = "results/degree2_inverse_selection_specialize_first_summary.txt"
OUTPUT_SOBJ = "results/degree2_inverse_selection_specialize_first.sobj"

POLAR_NAMES = [
    "alpha_U_-3", "beta_U_-3", "alpha_U_-2",
    "beta_U_-2", "alpha_U_-1",
]
CONSTANT_NAMES = ["alpha_const", "beta_const"]

paths = sorted(glob.glob(INPUT_GLOB))
if not paths:
    raise RuntimeError("No specialize-first witness files found")

records = [load(path) for path in paths]
print("="*78)
print("SPECIALIZE-FIRST INVERSE-SELECTION CERTIFICATE")
print("="*78)
print("witness files: {}".format(len(records)))
for path, record in zip(paths, records):
    print("  {}: p={}, point={}".format(
        os.path.basename(path), record["prime"], record["parameter_values"]
    ))
print()


def canonical_rows(names):
    names = list(names)
    # Exact formal relation, including intercept and all geometry slopes:
    # beta_-3 = rho3 alpha_-3.
    if "alpha_U_-3" in names and "beta_U_-3" in names:
        names.remove("beta_U_-3")
    return tuple(names)


def matrix_ranks_at(record, names):
    p = Integer(record["prime"])
    Fp = GF(p)
    rows = {
        name: [Fp(value) for value in values]
        for name, values in record["rows_mod_p"].items()
    }
    if not names:
        return 0, 0, None
    M = matrix(Fp, [rows[name][1:4] for name in names])
    A = matrix(Fp, [rows[name][1:4] + [-rows[name][0]] for name in names])
    solution = None
    if len(names) == 3 and M.rank() == 3:
        solution = tuple(M.solve_right(vector(Fp, [-rows[name][0] for name in names])))
    return Integer(M.rank()), Integer(A.rank()), solution


def classify(names):
    effective = canonical_rows(names)
    m = len(effective)
    ranks = []
    solutions = []
    for index, record in enumerate(records):
        rank_m, rank_a, solution = matrix_ranks_at(record, effective)
        ranks.append((rank_m, rank_a, index))
        if solution is not None:
            solutions.append((index, tuple(Integer(v) for v in solution)))
    max_m = max(item[0] for item in ranks)
    max_a = max(item[1] for item in ranks)

    if m <= 3 and max_m == m:
        status = "CERTIFIED CONSISTENT"
        dimension = 3 - m
        witness = next(item[2] for item in ranks if item[0] == m)
    elif m >= 4 and max_a >= 4:
        status = "CERTIFIED INCONSISTENT"
        dimension = None
        witness = next(item[2] for item in ranks if item[1] >= 4)
    else:
        status = "UNRESOLVED BY CURRENT WITNESSES"
        dimension = None
        witness = None

    return {
        "original_rows": tuple(names),
        "effective_rows": effective,
        "effective_count": m,
        "max_coefficient_rank": max_m,
        "max_augmented_rank": max_a,
        "status": status,
        "geometry_dimension": dimension,
        "witness_index": witness,
        "modular_unique_solutions": solutions,
    }


hierarchy = {
    "NO CUBIC POLE": ("alpha_U_-3", "beta_U_-3"),
    "NO CUBIC OR QUADRATIC POLES": (
        "alpha_U_-3", "beta_U_-3", "alpha_U_-2", "beta_U_-2"
    ),
    "POLE-FREE TRANSPORT": tuple(POLAR_NAMES),
    "PURE UNIVERSAL TRANSPORT": tuple(POLAR_NAMES + CONSTANT_NAMES),
}

print("HIERARCHY")
print("-"*78)
hierarchy_results = {}
for label, names in hierarchy.items():
    result = classify(names)
    hierarchy_results[label] = result
    print(label)
    print("  requested rows : {}".format(list(names)))
    print("  independent rows after exact cubic relation: {}".format(
        list(result["effective_rows"])
    ))
    print("  max ranks M/(M|r): {}/{}".format(
        result["max_coefficient_rank"], result["max_augmented_rank"]
    ))
    print("  status: {}".format(result["status"]))
    if result["geometry_dimension"] is not None:
        print("  geometry-family dimension: {}".format(result["geometry_dimension"]))
    if result["witness_index"] is not None:
        print("  certified by witness point {:02d}".format(result["witness_index"]))
    if result["modular_unique_solutions"]:
        print("  unique geometry residues at witness primes:")
        for index, solution in result["modular_unique_solutions"]:
            print("    point {:02d}: {}".format(index, solution))
    print()

print("EXHAUSTIVE POLAR SUBSET CLASSIFICATION")
print("-"*78)
subset_results = {}
for size in range(6):
    for subset in combinations(POLAR_NAMES, size):
        subset_results[tuple(subset)] = classify(subset)

for subset in sorted(subset_results, key=lambda x: (len(x), x)):
    result = subset_results[subset]
    print("  {:58s} {}".format(str(subset), result["status"]))

consistent_subsets = {
    subset for subset, result in subset_results.items()
    if result["status"] == "CERTIFIED CONSISTENT"
}
maximal_consistent = []
for subset in consistent_subsets:
    if not any(set(subset) < set(other) for other in consistent_subsets):
        maximal_consistent.append(subset)
maximal_consistent.sort(key=lambda x: (-len(x), x))

print()
print("INCLUSION-MAXIMAL CERTIFIED COMPATIBLE ZERO SETS")
print("-"*78)
if maximal_consistent:
    for subset in maximal_consistent:
        result = subset_results[subset]
        print("  {}  [dimension {}]".format(
            list(subset), result["geometry_dimension"]
        ))
else:
    print("  none certified")

print()
print("NAMED-GEOMETRY NONVANISHING WITNESSES")
print("-"*78)
named_labels = sorted(records[0]["named_values_mod_p"])
named_results = {}
for label in named_labels:
    nonzero_rows = set()
    witnesses = {}
    for name in POLAR_NAMES + CONSTANT_NAMES:
        for index, record in enumerate(records):
            value = Integer(record["named_values_mod_p"][label][name])
            if value != 0:
                nonzero_rows.add(name)
                witnesses[name] = index
                break
    named_results[label] = {
        "certified_nonzero_rows": sorted(nonzero_rows),
        "witness_indices": witnesses,
        "pole_free_disproved": any(name in nonzero_rows for name in POLAR_NAMES),
    }
    print("{}".format(label))
    print("  certified nonzero rows: {}".format(sorted(nonzero_rows)))
    print("  pole-free disproved? {}".format(
        named_results[label]["pole_free_disproved"]
    ))

unresolved = [
    subset for subset, result in subset_results.items()
    if result["status"].startswith("UNRESOLVED")
]

summary_lines = [
    "DEGREE-TWO SPECIALIZE-FIRST INVERSE-SELECTION CERTIFICATE",
    "witness count: {}".format(len(records)),
    "method: specialize parameters before geometry; no formal connection map loaded",
    "exact dependency used: beta_U_-3 = rho3 alpha_U_-3",
    "",
]
for label, result in hierarchy_results.items():
    summary_lines.append("{}: {}".format(label, result["status"]))
    summary_lines.append("  effective rows: {}".format(list(result["effective_rows"])))
    summary_lines.append("  max ranks: {}/{}".format(
        result["max_coefficient_rank"], result["max_augmented_rank"]
    ))
    if result["geometry_dimension"] is not None:
        summary_lines.append("  geometry dimension: {}".format(
            result["geometry_dimension"]
        ))
summary_lines.extend([
    "",
    "maximal certified compatible polar zero sets:",
])
for subset in maximal_consistent:
    summary_lines.append("  {} (dimension {})".format(
        list(subset), subset_results[subset]["geometry_dimension"]
    ))
summary_lines.extend([
    "",
    "unresolved polar subsets: {}".format(len(unresolved)),
])

os.makedirs("results", exist_ok=True)
with open(OUTPUT_TEXT, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

certificate = {
    "method": "specialize-before-geometry",
    "witness_paths": paths,
    "hierarchy_results": hierarchy_results,
    "subset_results": subset_results,
    "maximal_consistent_subsets": maximal_consistent,
    "named_results": named_results,
    "unresolved_subsets": unresolved,
}
save(certificate, OUTPUT_SOBJ)

print()
print("FINAL SUMMARY")
print("-"*78)
for line in summary_lines:
    print(line)
print()
print("wrote {}".format(OUTPUT_TEXT))
print("wrote {}".format(OUTPUT_SOBJ))
