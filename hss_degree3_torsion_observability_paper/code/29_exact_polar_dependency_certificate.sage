from sage.all import *
from itertools import combinations
import glob
import os

# ============================================================================
# CELL 29: EXACT POLAR-DEPENDENCY AND MINIMALITY CERTIFICATE
# ============================================================================
#
# This is a lightweight corollary of Cells 27 and 28E.
#
# Exact inputs already established upstream:
#   (i)  beta_-3 = rho3 alpha_-3 in every affine column;
#   (ii) every defect column has F_-4 = 0;
#   (iii) the geometry slopes have F_-3 = 0, while the intercept has a
#         nonzero U^-3 coefficient;
#   (iv) the specialize-first files are exact good-reduction witnesses.
#
# The U^-4 and U^-3 Laurent equations then determine all generic linear
# dependencies among the five polar rows.  The eight finite-field witnesses
# certify that the two slope basis rows and the intercept residual are truly
# independent.  Consequently all 32 polar zero subsets can be classified
# exactly without loading the 467 MB formal checkpoint.
# ============================================================================

INPUT_GLOB = "results/28e_specialize_first/point_*.sobj"
OUTPUT_TEXT = "results/degree2_exact_polar_dependency_summary.txt"
OUTPUT_SOBJ = "results/degree2_exact_polar_dependency_certificate.sobj"

POLAR_NAMES = [
    "alpha_U_-3", "beta_U_-3", "alpha_U_-2",
    "beta_U_-2", "alpha_U_-1",
]
CONSTANT_NAMES = ["alpha_const", "beta_const"]

paths = sorted(glob.glob(INPUT_GLOB))
if not paths:
    raise RuntimeError(
        "No specialize-first witnesses found at {}. Run Cell 28E/V4 first."
        .format(INPUT_GLOB)
    )
records = [load(path) for path in paths]

print("="*78)
print("CELL 29: EXACT POLAR-DEPENDENCY CERTIFICATE")
print("="*78)
print("witness files: {}".format(len(records)))
print()

# ----------------------------------------------------------------------------
# PART I. Derive the dependency coefficients from the finite Laurent channels.
# ----------------------------------------------------------------------------

R = PolynomialRing(QQ, names=("L", "Y", "C0"))
L, Y, C0 = R.gens()
K = R.fraction_field()
L, Y, C0 = map(K, (L, Y, C0))

s_m2 = -C0^2*(Y + 108)/(972*L^3)
s_m1 = -C0^2*(Y + 84)/(1944*L^2)
s_0  = -C0^2/(243*L)
d_m2 =  C0^2*(Y + 36)/(648*L^3)
d_m1 =  C0^2*(Y + 36)/(1944*L^2)

rho = K(-s_m2/d_m2)

# U^-4 = 0 after beta_-3 = rho alpha_-3:
#   (s_-1 + rho d_-1) alpha_-3 + s_-2 alpha_-2
#       + d_-2 beta_-2 = 0.
b2_from_a3 = K(-(s_m1 + rho*d_m1)/d_m2)
b2_from_a2 = K(-s_m2/d_m2)

# U^-3 = 0 for each geometry slope:
#   s_0 alpha_-3 + s_-1 alpha_-2 + d_-1 beta_-2
#       + s_-2 alpha_-1 = 0.
a1_from_a3 = K(-(s_0 + d_m1*b2_from_a3)/s_m2)
a1_from_a2 = K(-(s_m1 + d_m1*b2_from_a2)/s_m2)

rho_expected = K(2*(Y + 108)/(3*(Y + 36)))
b2_a3_expected = K(L/9)
b2_a2_expected = rho_expected
a1_a3_expected = K(L^2*(Y - 36)/(18*(Y + 108)))
a1_a2_expected = K(-L*(Y + 36)/(6*(Y + 108)))

symbolic_checks = {
    "rho": rho == rho_expected,
    "beta_-2 coefficient of alpha_-3": b2_from_a3 == b2_a3_expected,
    "beta_-2 coefficient of alpha_-2": b2_from_a2 == b2_a2_expected,
    "alpha_-1 coefficient of alpha_-3": a1_from_a3 == a1_a3_expected,
    "alpha_-1 coefficient of alpha_-2": a1_from_a2 == a1_a2_expected,
}
if not all(symbolic_checks.values()):
    raise ArithmeticError("Laurent dependency simplification failed")

Delta = K(Y^2 + 72*Y - 2160)
selected_pair_determinant = K(-L^2*Delta/(18*(Y + 36)*(Y + 108)))

print("PART I. EXACT LAURENT DEPENDENCIES")
print("-"*78)
print("  rho3 = 2*(Y+108)/(3*(Y+36))")
print("  beta_-3 = rho3*alpha_-3")
print()
print("  beta_-2 = (L/9)*alpha_-3 + rho3*alpha_-2")
print("  equivalently:")
print("    9*(Y+36)*beta_-2")
print("      - L*(Y+36)*alpha_-3 - 6*(Y+108)*alpha_-2 = 0")
print()
print("  on the three geometry-slope columns:")
print("  alpha_-1 = A*alpha_-3 + B*alpha_-2")
print("    A = L^2*(Y-36)/(18*(Y+108))")
print("    B = -L*(Y+36)/(6*(Y+108))")
print()
print("  det[beta_-2, alpha_-1] in the slope basis =")
print("    -L^2*(Y^2+72*Y-2160)/(18*(Y+36)*(Y+108))")
print("  symbolic derivation passed? {}".format(all(symbolic_checks.values())))
print()

# ----------------------------------------------------------------------------
# PART II. Verify the exact relations at every specialize-first witness.
# ----------------------------------------------------------------------------

witness_checks = []
rank2_witness = None
rank3_augmented_witness = None
pure_universal_witness = None
residual_nonzero_witness = None

for record_index, record in enumerate(records):
    p = Integer(record["prime"])
    F = GF(p)
    parameter_values = record["parameter_values"]
    Lp = F(parameter_values[0])
    E4p = F(parameter_values[1])
    Yp = Lp^2*E4p

    denominators = [Lp, Yp + 36, Yp + 108]
    if any(value == 0 for value in denominators):
        raise ArithmeticError("bad reduction in stored witness {}".format(record_index))

    rhop = F(2)*(Yp + 108)/(F(3)*(Yp + 36))
    Ap = Lp^2*(Yp - 36)/(F(18)*(Yp + 108))
    Bp = -Lp*(Yp + 36)/(F(6)*(Yp + 108))

    rows = {
        name: vector(F, [F(Integer(value)) for value in values])
        for name, values in record["rows_mod_p"].items()
    }

    cubic_ok = rows["beta_U_-3"] == rhop*rows["alpha_U_-3"]
    quadratic_residual = (
        rows["beta_U_-2"]
        - (Lp/F(9))*rows["alpha_U_-3"]
        - rhop*rows["alpha_U_-2"]
    )
    quadratic_ok = quadratic_residual.is_zero()

    simple_residual = (
        rows["alpha_U_-1"]
        - Ap*rows["alpha_U_-3"]
        - Bp*rows["alpha_U_-2"]
    )
    simple_slopes_ok = all(simple_residual[j] == 0 for j in range(1, 4))

    effective = [
        rows["alpha_U_-3"], rows["alpha_U_-2"],
        rows["beta_U_-2"], rows["alpha_U_-1"],
    ]
    Mpolar = matrix(F, [list(row[1:4]) for row in effective])
    Apolar = matrix(F, [list(row[1:4]) + [-row[0]] for row in effective])

    if rank2_witness is None and Mpolar.rank() == 2:
        rank2_witness = record_index
    if rank3_augmented_witness is None and Apolar.rank() == 3:
        rank3_augmented_witness = record_index
    if residual_nonzero_witness is None and simple_residual[0] != 0:
        residual_nonzero_witness = record_index

    pure_names = POLAR_NAMES + CONSTANT_NAMES
    Mpure = matrix(F, [list(rows[name][1:4]) for name in pure_names])
    Apure = matrix(F, [list(rows[name][1:4]) + [-rows[name][0]]
                       for name in pure_names])
    if pure_universal_witness is None and Apure.rank() >= 4:
        pure_universal_witness = record_index

    checks = {
        "cubic": cubic_ok,
        "quadratic": quadratic_ok,
        "simple_slopes": simple_slopes_ok,
        "simple_residual_intercept": Integer(simple_residual[0]),
        "polar_coefficient_rank": Integer(Mpolar.rank()),
        "polar_augmented_rank": Integer(Apolar.rank()),
        "pure_coefficient_rank": Integer(Mpure.rank()),
        "pure_augmented_rank": Integer(Apure.rank()),
    }
    witness_checks.append(checks)

    if not (cubic_ok and quadratic_ok and simple_slopes_ok):
        raise ArithmeticError(
            "structural relation failed at witness {}".format(record_index)
        )

if rank2_witness is None:
    raise RuntimeError("No witness proves the generic polar coefficient rank is 2")
if rank3_augmented_witness is None or residual_nonzero_witness is None:
    raise RuntimeError("No witness proves the simple-pole intercept residual is nonzero")
if pure_universal_witness is None:
    raise RuntimeError("No witness certifies pure universal transport is inconsistent")

print("PART II. GOOD-REDUCTION WITNESSES")
print("-"*78)
print("  all witnesses satisfy the cubic relation: True")
print("  all witnesses satisfy the U^-4 relation: True")
print("  all witnesses satisfy the slope U^-3 relation: True")
print("  coefficient-rank-2 witness: point {:02d}".format(rank2_witness))
print("  augmented-rank-3 witness: point {:02d}".format(rank3_augmented_witness))
print("  nonzero intercept-residual witness: point {:02d}".format(
    residual_nonzero_witness
))
print("  pure-universal rank-4 witness: point {:02d}".format(
    pure_universal_witness
))
print()

# ----------------------------------------------------------------------------
# PART III. Exact abstract polar matroid.
# ----------------------------------------------------------------------------
#
# The slope rows live in a two-dimensional basis (alpha_-3, alpha_-2).
# The augmented rows need one extra coordinate for the nonzero intercept-only
# U^-3 residual.  Good-reduction witnesses prove these basis directions are
# genuinely independent over the characteristic-zero formal field.
# ----------------------------------------------------------------------------

S = PolynomialRing(QQ, names=("l", "y"))
l, y = S.gens()
Q = S.fraction_field()
l, y = map(Q, (l, y))

r = Q(2*(y + 108)/(3*(y + 36)))
Acoef = Q(l^2*(y - 36)/(18*(y + 108)))
Bcoef = Q(-l*(y + 36)/(6*(y + 108)))

# Coefficient/slope row coordinates in basis (alpha_-3, alpha_-2).
coefficient_model = {
    "alpha_U_-3": vector(Q, [1, 0]),
    "beta_U_-3":  vector(Q, [r, 0]),
    "alpha_U_-2": vector(Q, [0, 1]),
    "beta_U_-2":  vector(Q, [l/9, r]),
    "alpha_U_-1": vector(Q, [Acoef, Bcoef]),
}

# Augmented row coordinates in basis (intercept residual, alpha_-3, alpha_-2).
# The residual coordinate of alpha_-1 is normalized to 1; only its nonvanishing
# matters for rank and consistency.
augmented_model = {
    "alpha_U_-3": vector(Q, [0, 1, 0]),
    "beta_U_-3":  vector(Q, [0, r, 0]),
    "alpha_U_-2": vector(Q, [0, 0, 1]),
    "beta_U_-2":  vector(Q, [0, l/9, r]),
    "alpha_U_-1": vector(Q, [1, Acoef, Bcoef]),
}


def model_rank(rows, model, width):
    if not rows:
        return Integer(0)
    return Integer(matrix(Q, [list(model[name]) for name in rows]).rank())


def classify_polar_subset(subset):
    rank_m = model_rank(subset, coefficient_model, 2)
    rank_a = model_rank(subset, augmented_model, 3)
    consistent = (rank_m == rank_a)
    return {
        "rows": tuple(subset),
        "coefficient_rank": rank_m,
        "augmented_rank": rank_a,
        "consistent": bool(consistent),
        "geometry_dimension": Integer(3-rank_m) if consistent else None,
    }

subset_results = {}
for size in range(len(POLAR_NAMES) + 1):
    for subset in combinations(POLAR_NAMES, size):
        subset_results[tuple(subset)] = classify_polar_subset(subset)

consistent_subsets = {
    subset for subset, result in subset_results.items() if result["consistent"]
}
maximal_consistent = []
for subset in consistent_subsets:
    if not any(set(subset) < set(other) for other in consistent_subsets):
        maximal_consistent.append(subset)
maximal_consistent.sort(key=lambda item: (-len(item), item))

unresolved = []

print("PART III. COMPLETE EXACT POLAR-SUBSET CLASSIFICATION")
print("-"*78)
for subset in sorted(subset_results, key=lambda item: (len(item), item)):
    result = subset_results[subset]
    status = "CONSISTENT" if result["consistent"] else "INCONSISTENT"
    suffix = (
        "dimension {}".format(result["geometry_dimension"])
        if result["consistent"] else "no geometry"
    )
    print("  {:58s} {:12s} ranks {}/{}; {}".format(
        str(subset), status,
        result["coefficient_rank"], result["augmented_rank"], suffix
    ))

print()
print("  inclusion-maximal compatible polar zero sets:")
for subset in maximal_consistent:
    print("    {} (dimension {})".format(
        list(subset), subset_results[subset]["geometry_dimension"]
    ))
print("  unresolved polar subsets: 0")
print()

# ----------------------------------------------------------------------------
# PART IV. Hierarchy and pure-universal conclusion.
# ----------------------------------------------------------------------------

hierarchy = {
    "NO CUBIC POLE": ("alpha_U_-3", "beta_U_-3"),
    "NO CUBIC OR QUADRATIC POLES": (
        "alpha_U_-3", "beta_U_-3", "alpha_U_-2", "beta_U_-2"
    ),
    "POLE-FREE TRANSPORT": tuple(POLAR_NAMES),
}

hierarchy_results = {
    label: classify_polar_subset(rows) for label, rows in hierarchy.items()
}

pure_universal_result = {
    "rows": tuple(POLAR_NAMES + CONSTANT_NAMES),
    "coefficient_rank_lower_bound": 3,
    "augmented_rank_lower_bound": 4,
    "consistent": False,
    "witness_index": pure_universal_witness,
}

print("PART IV. FINAL MINIMALITY HIERARCHY")
print("-"*78)
for label in ["NO CUBIC POLE", "NO CUBIC OR QUADRATIC POLES",
              "POLE-FREE TRANSPORT"]:
    result = hierarchy_results[label]
    status = "CONSISTENT" if result["consistent"] else "INCONSISTENT"
    print("  {}: {}".format(label, status))
    print("    ranks M/(M|r): {}/{}".format(
        result["coefficient_rank"], result["augmented_rank"]
    ))
    if result["consistent"]:
        print("    geometry-family dimension: {}".format(
            result["geometry_dimension"]
        ))
print("  PURE UNIVERSAL TRANSPORT: INCONSISTENT")
print("    exact good-reduction witness has augmented rank 4")
print()

summary_lines = [
    "DEGREE-TWO EXACT POLAR-DEPENDENCY CERTIFICATE",
    "witness count: {}".format(len(records)),
    "method: Laurent structural relations plus exact good-reduction witnesses",
    "",
    "EXACT RELATIONS",
    "  beta_-3 = rho3 alpha_-3",
    "  rho3 = 2(Y+108)/(3(Y+36))",
    "  beta_-2 = (L/9) alpha_-3 + rho3 alpha_-2",
    "  slope(alpha_-1) = A slope(alpha_-3) + B slope(alpha_-2)",
    "  A = L^2(Y-36)/(18(Y+108))",
    "  B = -L(Y+36)/(6(Y+108))",
    "  alpha_-1 has a nonzero intercept-only residual at U^-3",
    "",
]
for label in ["NO CUBIC POLE", "NO CUBIC OR QUADRATIC POLES",
              "POLE-FREE TRANSPORT"]:
    result = hierarchy_results[label]
    status = "CERTIFIED CONSISTENT" if result["consistent"] else "CERTIFIED INCONSISTENT"
    summary_lines.append("{}: {}".format(label, status))
    summary_lines.append("  exact ranks: {}/{}".format(
        result["coefficient_rank"], result["augmented_rank"]
    ))
    if result["consistent"]:
        summary_lines.append("  geometry dimension: {}".format(
            result["geometry_dimension"]
        ))
summary_lines.extend([
    "PURE UNIVERSAL TRANSPORT: CERTIFIED INCONSISTENT",
    "  witness augmented rank: 4",
    "",
    "maximal compatible polar zero sets:",
])
for subset in maximal_consistent:
    summary_lines.append("  {} (dimension {})".format(
        list(subset), subset_results[subset]["geometry_dimension"]
    ))
summary_lines.extend([
    "",
    "unresolved polar subsets: 0",
])

os.makedirs("results", exist_ok=True)
with open(OUTPUT_TEXT, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

certificate = {
    "method": "exact-Laurent-relations-plus-good-reduction",
    "witness_paths": paths,
    "symbolic_checks": symbolic_checks,
    "dependency_coefficients": {
        "rho3": rho,
        "beta_m2_from_alpha_m3": b2_from_a3,
        "beta_m2_from_alpha_m2": b2_from_a2,
        "alpha_m1_from_alpha_m3_on_slopes": a1_from_a3,
        "alpha_m1_from_alpha_m2_on_slopes": a1_from_a2,
        "polar_pair_discriminant": selected_pair_determinant,
    },
    "witness_checks": witness_checks,
    "rank2_witness": rank2_witness,
    "rank3_augmented_witness": rank3_augmented_witness,
    "residual_nonzero_witness": residual_nonzero_witness,
    "pure_universal_witness": pure_universal_witness,
    "hierarchy_results": hierarchy_results,
    "pure_universal_result": pure_universal_result,
    "subset_results": subset_results,
    "maximal_consistent_subsets": maximal_consistent,
    "unresolved_subsets": unresolved,
}
save(certificate, OUTPUT_SOBJ)

print("FINAL SUMMARY")
print("-"*78)
for line in summary_lines:
    print(line)
print()
print("wrote {}".format(OUTPUT_TEXT))
print("wrote {}".format(OUTPUT_SOBJ))
