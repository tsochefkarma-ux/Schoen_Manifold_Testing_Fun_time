from sage.all import *
import os

# =============================================================================
# CELL 30: DEGREE-THREE QUASI-JACOBI BOOTSTRAP DIMENSION CENSUS
# =============================================================================
#
# Purpose
# -------
# Before constructing the p^3 HSS curvature response, determine the finite
# holomorphic ambiguity space which transport would have to reconstruct.
#
# For the genus-g, base-degree-k Schoen potential, the one-sided E8 Jacobi
# numerator has weight
#
#     ell = 2*g - 2 + 6*k
#
# and E8 index k.  At (g,k)=(0,3), this gives weight 16 and index 3.
#
# The holomorphic W(E8)-invariant index-three module is free over
# M_* = QQ[E4,E6] on
#
#     A3, B3, A1*A2, A1*B2, A1^3
#
# of weights 4, 6, 8, 10, 12.  This cell enumerates all modular multiples at
# total weight 16, verifies the equivalent Sakai-monomial enumeration, builds
# the symmetric bi-Jacobi basis, and writes a configurable reconstruction
# ledger for the two parked torsion-character sectors.
#
# This is a counting/certificate cell only.  It does not yet construct A3 or
# B3 q-expansions and does not assume a degree-three transport module.
# =============================================================================

print("="*79, flush=True)
print("DEGREE-THREE QUASI-JACOBI BOOTSTRAP DIMENSION CENSUS", flush=True)
print("="*79, flush=True)

# ----------------------------------------------------------------------------
# User/configuration block
# ----------------------------------------------------------------------------

GENUS = ZZ(0)
BASE_DEGREE = ZZ(3)
TARGET_WEIGHT = ZZ(2*GENUS - 2 + 6*BASE_DEGREE)
TARGET_INDEX = BASE_DEGREE

# The two unresolved torsion-character sectors named in the research ledger.
TORSION_SECTORS = ["(omega,omega)", "(omega,omega^2)"]

# Possible prior information can be entered here once established exactly.
# These numbers are deliberately zero by default: this cell reports what must
# still be supplied by HAE, torsion sums, transport, and independent checks.
KNOWN_ONE_SIDED_LINEAR_CONSTRAINTS_PER_SECTOR = ZZ(0)
KNOWN_SYMMETRIC_BI_LINEAR_CONSTRAINTS_PER_SECTOR = ZZ(0)
CROSS_SECTOR_LINEAR_CONSTRAINTS = ZZ(0)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

print("genus g              = {}".format(GENUS), flush=True)
print("base degree k        = {}".format(BASE_DEGREE), flush=True)
print("target Jacobi weight = {}".format(TARGET_WEIGHT), flush=True)
print("target E8 index      = {}".format(TARGET_INDEX), flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# Modular-form monomials at level one
# ----------------------------------------------------------------------------

def modular_monomials(weight):
    """Return labels E4^a E6^b of exact modular weight `weight`."""
    weight = ZZ(weight)
    if weight < 0 or weight % 2:
        return []
    output = []
    for a in range(weight//4 + 1):
        for b in range(weight//6 + 1):
            if 4*a + 6*b == weight:
                factors = []
                if a:
                    factors.append("E4" if a == 1 else "E4^{}".format(a))
                if b:
                    factors.append("E6" if b == 1 else "E6^{}".format(b))
                output.append("*".join(factors) if factors else "1")
    return output

# Exact index-three holomorphic module generators over M_*.
module_generators = [
    {"label": "A3",       "weight": ZZ(4),  "index": ZZ(3)},
    {"label": "B3",       "weight": ZZ(6),  "index": ZZ(3)},
    {"label": "A1*A2",    "weight": ZZ(8),  "index": ZZ(3)},
    {"label": "A1*B2",    "weight": ZZ(10), "index": ZZ(3)},
    {"label": "A1^3",     "weight": ZZ(12), "index": ZZ(3)},
]

one_sided_basis = []
for generator in module_generators:
    residual_weight = TARGET_WEIGHT - generator["weight"]
    for modular_factor in modular_monomials(residual_weight):
        if modular_factor == "1":
            label = generator["label"]
        else:
            label = modular_factor + "*" + generator["label"]
        one_sided_basis.append({
            "label": label,
            "module_generator": generator["label"],
            "modular_factor": modular_factor,
            "weight": TARGET_WEIGHT,
            "index": TARGET_INDEX,
        })

print("PART I. EXACT ONE-SIDED HOLOMORPHIC AMBIGUITY", flush=True)
print("-"*79, flush=True)
print("module generators over M_*:", flush=True)
for generator in module_generators:
    print("  {:8s} weight {:2d}, index {}".format(
        generator["label"], int(generator["weight"]), int(generator["index"])
    ), flush=True)

print("\nweight-{} index-{} basis:".format(TARGET_WEIGHT, TARGET_INDEX), flush=True)
for position, term in enumerate(one_sided_basis):
    print("  h{} = {}".format(position + 1, term["label"]), flush=True)

ONE_SIDED_DIMENSION = ZZ(len(one_sided_basis))
print("\none-sided holomorphic ambiguity dimension = {}".format(
    ONE_SIDED_DIMENSION
), flush=True)

# ----------------------------------------------------------------------------
# Independent cross-check by raw Sakai-monomial enumeration
# ----------------------------------------------------------------------------

raw_generators = [
    ("E4", ZZ(4), ZZ(0)),
    ("E6", ZZ(6), ZZ(0)),
    ("A1", ZZ(4), ZZ(1)),
    ("A2", ZZ(4), ZZ(2)),
    ("A3", ZZ(4), ZZ(3)),
    ("B2", ZZ(6), ZZ(2)),
    ("B3", ZZ(6), ZZ(3)),
]


def raw_monomial_label(exponents):
    factors = []
    for exponent, (name, _, _) in zip(exponents, raw_generators):
        if exponent == 0:
            continue
        factors.append(name if exponent == 1 else "{}^{}".format(name, exponent))
    return "*".join(factors) if factors else "1"


def enumerate_raw_monomials(target_weight, target_index):
    # Bounds are tiny because all generator weights are positive.
    bounds = [target_weight // weight + 1 for _, weight, _ in raw_generators]
    output = []

    def recurse(position, remaining_weight, remaining_index, exponents):
        if position == len(raw_generators):
            if remaining_weight == 0 and remaining_index == 0:
                output.append(tuple(exponents))
            return

        _, weight, index = raw_generators[position]
        maximum = min(
            remaining_weight // weight,
            remaining_index // index if index > 0 else remaining_weight // weight,
        )
        for exponent in range(maximum + 1):
            recurse(
                position + 1,
                remaining_weight - exponent*weight,
                remaining_index - exponent*index,
                exponents + [ZZ(exponent)],
            )

    recurse(0, ZZ(target_weight), ZZ(target_index), [])
    return output

raw_exponents = enumerate_raw_monomials(TARGET_WEIGHT, TARGET_INDEX)
raw_labels = sorted(raw_monomial_label(exponents) for exponents in raw_exponents)
module_labels = sorted(term["label"] for term in one_sided_basis)

print("\nPART II. INDEPENDENT SAKAI-MONOMIAL CROSS-CHECK", flush=True)
print("-"*79, flush=True)
for label in raw_labels:
    print("  {}".format(label), flush=True)
print("raw monomial count = {}".format(len(raw_labels)), flush=True)
print("matches module enumeration? {}".format(raw_labels == module_labels), flush=True)
assert raw_labels == module_labels

# ----------------------------------------------------------------------------
# Symmetric E8 x E8 bi-Jacobi ambiguity
# ----------------------------------------------------------------------------

symmetric_bi_basis = []
ordered_bi_basis = []
for i, left in enumerate(one_sided_basis):
    for j, right in enumerate(one_sided_basis):
        ordered_bi_basis.append((left["label"], right["label"]))
        if i <= j:
            if i == j:
                label = "{}(y)*{}(z)".format(left["label"], right["label"])
            else:
                label = (
                    "{}(y)*{}(z) + {}(y)*{}(z)".format(
                        left["label"], right["label"],
                        right["label"], left["label"],
                    )
                )
            symmetric_bi_basis.append(label)

ORDERED_BI_DIMENSION = ZZ(len(ordered_bi_basis))
SYMMETRIC_BI_DIMENSION = ZZ(len(symmetric_bi_basis))

print("\nPART III. BI-JACOBI AMBIGUITY COUNTS", flush=True)
print("-"*79, flush=True)
print("ordered tensor-product dimension  = {}".format(ORDERED_BI_DIMENSION), flush=True)
print("symmetric tensor-square dimension = {}".format(SYMMETRIC_BI_DIMENSION), flush=True)
print("expected symmetric dimension n(n+1)/2 = {}".format(
    ONE_SIDED_DIMENSION*(ONE_SIDED_DIMENSION + 1)//2
), flush=True)
assert SYMMETRIC_BI_DIMENSION == ONE_SIDED_DIMENSION*(ONE_SIDED_DIMENSION + 1)//2

print("\nfirst symmetric basis elements:", flush=True)
for position, label in enumerate(symmetric_bi_basis[:min(10, len(symmetric_bi_basis))]):
    print("  H{} = {}".format(position + 1, label), flush=True)
if len(symmetric_bi_basis) > 10:
    print("  ... {} further elements".format(len(symmetric_bi_basis) - 10), flush=True)

# ----------------------------------------------------------------------------
# Reconstruction thresholds for the parked torsion sectors
# ----------------------------------------------------------------------------

sector_count = ZZ(len(TORSION_SECTORS))
unknowns_one_sided = sector_count*ONE_SIDED_DIMENSION
unknowns_symmetric_bi = sector_count*SYMMETRIC_BI_DIMENSION

remaining_one_sided = max(
    ZZ(0),
    unknowns_one_sided
    - sector_count*KNOWN_ONE_SIDED_LINEAR_CONSTRAINTS_PER_SECTOR
    - CROSS_SECTOR_LINEAR_CONSTRAINTS,
)
remaining_symmetric_bi = max(
    ZZ(0),
    unknowns_symmetric_bi
    - sector_count*KNOWN_SYMMETRIC_BI_LINEAR_CONSTRAINTS_PER_SECTOR
    - CROSS_SECTOR_LINEAR_CONSTRAINTS,
)

print("\nPART IV. DEGREE-THREE BOOTSTRAP RANK THRESHOLDS", flush=True)
print("-"*79, flush=True)
print("torsion sectors: {}".format(", ".join(TORSION_SECTORS)), flush=True)
print("\nScenario A: each parked character is a one-sided index-three series", flush=True)
print("  total unknown coefficients              = {}".format(unknowns_one_sided), flush=True)
print("  currently entered prior constraints     = {}".format(
    unknowns_one_sided - remaining_one_sided
), flush=True)
print("  independent new equations needed        = {}".format(remaining_one_sided), flush=True)
print("  square reconstruction target rank        = {}".format(remaining_one_sided), flush=True)

print("\nScenario B: each parked character is a symmetric bi-Jacobi ambiguity", flush=True)
print("  total unknown coefficients              = {}".format(unknowns_symmetric_bi), flush=True)
print("  currently entered prior constraints     = {}".format(
    unknowns_symmetric_bi - remaining_symmetric_bi
), flush=True)
print("  independent new equations needed        = {}".format(remaining_symmetric_bi), flush=True)
print("  square reconstruction target rank        = {}".format(remaining_symmetric_bi), flush=True)

# Degree-two exact calibration, included as a structural design constraint.
degree2_calibration = {
    "bare_transport_selective_constraints": ZZ(0),
    "no_cubic_codimension": ZZ(1),
    "no_cubic_or_quadratic_codimension": ZZ(2),
    "pole_free_consistent": False,
    "pure_universal_consistent": False,
}

print("\nPART V. STRUCTURAL LEDGER FOR THE NEXT COMPUTATION", flush=True)
print("-"*79, flush=True)
print("Degree-two calibration carried forward:", flush=True)
print("  bare closure is not selective", flush=True)
print("  deepest-pole cancellation has codimension 1", flush=True)
print("  cancellation through quadratic order has codimension 2", flush=True)
print("  complete pole-freeness is inconsistent", flush=True)
print("  pure universal transport is inconsistent", flush=True)

constraint_sources = [
    "torsion-character sum/recombination",
    "exchange and character symmetries",
    "holomorphic-anomaly particular solution",
    "allowed Laurent support at U=0",
    "exact residue dependencies",
    "universal positive-power tail",
    "restricted master-coordinate count",
    "independent low-q or enumerative checkpoints",
]

print("\nConstraint sources to measure at degree three:", flush=True)
for source in constraint_sources:
    print("  - {}".format(source), flush=True)

print("\nDECISION RULE", flush=True)
print("-"*79, flush=True)
print("For a chosen sector model with N unknown coefficients:", flush=True)
print("  transport/HAE rank < N  : underdetermined family", flush=True)
print("  transport/HAE rank = N  : unique candidate series", flush=True)
print("  transport/HAE rank > N  : overdetermined; surplus equations are predictions", flush=True)

# ----------------------------------------------------------------------------
# Save machine-readable and human-readable outputs
# ----------------------------------------------------------------------------

result = {
    "genus": GENUS,
    "base_degree": BASE_DEGREE,
    "target_weight": TARGET_WEIGHT,
    "target_index": TARGET_INDEX,
    "module_generators": module_generators,
    "one_sided_basis": one_sided_basis,
    "one_sided_dimension": ONE_SIDED_DIMENSION,
    "raw_monomial_labels": raw_labels,
    "symmetric_bi_basis": symmetric_bi_basis,
    "symmetric_bi_dimension": SYMMETRIC_BI_DIMENSION,
    "ordered_bi_dimension": ORDERED_BI_DIMENSION,
    "torsion_sectors": TORSION_SECTORS,
    "unknowns_two_one_sided_sectors": unknowns_one_sided,
    "unknowns_two_symmetric_bi_sectors": unknowns_symmetric_bi,
    "remaining_rank_one_sided_scenario": remaining_one_sided,
    "remaining_rank_symmetric_bi_scenario": remaining_symmetric_bi,
    "degree2_calibration": degree2_calibration,
    "constraint_sources": constraint_sources,
}

sobj_path = os.path.join(RESULTS_DIR, "degree3_bootstrap_dimension_census.sobj")
save(result, sobj_path)

summary_lines = [
    "DEGREE-THREE QUASI-JACOBI BOOTSTRAP DIMENSION CENSUS",
    "genus/base degree: ({},{})".format(GENUS, BASE_DEGREE),
    "target one-sided weight/index: {}/{}".format(TARGET_WEIGHT, TARGET_INDEX),
    "",
    "one-sided holomorphic basis (dimension {}):".format(ONE_SIDED_DIMENSION),
]
summary_lines.extend("  {}".format(term["label"]) for term in one_sided_basis)
summary_lines.extend([
    "",
    "ordered E8 x E8 ambiguity dimension: {}".format(ORDERED_BI_DIMENSION),
    "symmetric E8 x E8 ambiguity dimension: {}".format(SYMMETRIC_BI_DIMENSION),
    "",
    "parked torsion sectors: {}".format(", ".join(TORSION_SECTORS)),
    "two one-sided sectors: {} unknowns".format(unknowns_one_sided),
    "two symmetric bi-sectors: {} unknowns".format(unknowns_symmetric_bi),
    "",
    "next decision: determine which sector model is geometrically correct,",
    "then compute the rank of genuinely selective degree-three conditions.",
])

summary_path = os.path.join(RESULTS_DIR, "degree3_bootstrap_dimension_census.txt")
with open(summary_path, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

print("\nSaved:", flush=True)
print("  {}".format(sobj_path), flush=True)
print("  {}".format(summary_path), flush=True)
print("\nCELL 30 COMPLETE", flush=True)
