from sage.all import *
import os

# =============================================================================
# CELL 31: SHARED SYMMETRIC BI-JACOBI TORSION-SECTOR MODEL
# =============================================================================
#
# Purpose
# -------
# Correct the degree-three ambiguity ledger after Cell 30.
#
# The Schoen degree-k potential is one E8 x E8 bi-quasi-Jacobi object, symmetric
# under interchange of the two rational elliptic surfaces.  Therefore the
# parked character sectors (omega,omega) and (omega,omega^2) are two
# evaluations/projections of ONE symmetric coefficient matrix, not two
# independent 21-parameter forms.
#
# For the six one-sided weight-16, index-3 basis elements h_i, write
#
#   H(y,z) = sum_i c_ii h_i(y)h_i(z)
#          + sum_{i<j} c_ij [h_i(y)h_j(z)+h_j(y)h_i(z)].
#
# There are 6*7/2 = 21 shared coefficients c_ij.
#
# This cell:
#   1. constructs the shared 21-element symmetric basis;
#   2. writes the exact character-evaluation formulas;
#   3. creates a machine-readable template for future one-sided q-series data;
#   4. if that data already exists, constructs the sector coefficient matrix
#      and computes its exact rank over QQ.
#
# It is a small structural cell and does not construct A3 or B3 itself.
# =============================================================================

print("="*79, flush=True)
print("DEGREE-THREE SHARED BI-JACOBI TORSION-SECTOR MODEL", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

CENSUS_PATH = os.path.join(RESULTS_DIR, "degree3_bootstrap_dimension_census.sobj")
DATA_PATH = os.path.join(RESULTS_DIR, "degree3_one_sided_torsion_basis.sobj")
TEMPLATE_PATH = os.path.join(RESULTS_DIR, "degree3_one_sided_torsion_basis_template.sobj")
MODEL_PATH = os.path.join(RESULTS_DIR, "degree3_shared_bijacobi_model.sobj")
SUMMARY_PATH = os.path.join(RESULTS_DIR, "degree3_shared_bijacobi_model.txt")

DEFAULT_BASIS = [
    "E6^2*A3",
    "E4^3*A3",
    "E4*E6*B3",
    "E4^2*A1*A2",
    "E6*A1*B2",
    "E4*A1^3",
]

# -----------------------------------------------------------------------------
# Load the Cell 30 basis when available.
# -----------------------------------------------------------------------------

if os.path.exists(CENSUS_PATH):
    census = load(CENSUS_PATH)
    one_sided_basis = [entry["label"] for entry in census["one_sided_basis"]]
    source = CENSUS_PATH
else:
    census = None
    one_sided_basis = list(DEFAULT_BASIS)
    source = "built-in Cell 30 basis"

if sorted(one_sided_basis) != sorted(DEFAULT_BASIS):
    raise RuntimeError(
        "Unexpected Cell 30 basis. Expected the six weight-16 index-3 terms."
    )

# Preserve the canonical order used in the human ledger.
one_sided_basis = list(DEFAULT_BASIS)
n = len(one_sided_basis)
assert n == 6

print("basis source: {}".format(source), flush=True)
print("one-sided dimension: {}".format(n), flush=True)
for i, label in enumerate(one_sided_basis):
    print("  h{} = {}".format(i + 1, label), flush=True)

# -----------------------------------------------------------------------------
# Shared symmetric tensor basis.
# -----------------------------------------------------------------------------

pairs = []
coefficient_names = []
symmetric_basis = []
for i in range(n):
    for j in range(i, n):
        name = "c{}_{}".format(i + 1, j + 1)
        coefficient_names.append(name)
        pairs.append((i, j))
        if i == j:
            formula = "h{0}(y)*h{0}(z)".format(i + 1)
        else:
            formula = (
                "h{0}(y)*h{1}(z) + h{1}(y)*h{0}(z)".format(
                    i + 1, j + 1
                )
            )
        symmetric_basis.append({
            "coefficient": name,
            "pair": (ZZ(i), ZZ(j)),
            "left_label": one_sided_basis[i],
            "right_label": one_sided_basis[j],
            "formula": formula,
        })

shared_dimension = ZZ(len(symmetric_basis))
assert shared_dimension == n*(n + 1)//2 == 21

print("\nPART I. SHARED SYMMETRIC AMBIGUITY", flush=True)
print("-"*79, flush=True)
print("shared coefficient count = {}".format(shared_dimension), flush=True)
print("independent-two-sector count from Cell 30 = 42", flush=True)
print("correct shared-potential count = 21", flush=True)
print("overcount removed = 21", flush=True)

for k, entry in enumerate(symmetric_basis):
    print("  H{:02d}: {:5s}  {}".format(
        k + 1, entry["coefficient"], entry["formula"]
    ), flush=True)

# -----------------------------------------------------------------------------
# Exact evaluation formula for a character pair.
# -----------------------------------------------------------------------------

# Each sector is represented by two one-sided specialization vectors
# u=(u_1,...,u_6), v=(v_1,...,v_6).  The coefficient of c_ij is
# u_i v_i for i=j and u_i v_j + u_j v_i for i<j.

def sector_row(left_vector, right_vector):
    if len(left_vector) != n or len(right_vector) != n:
        raise ValueError("sector vectors must have length {}".format(n))
    row = []
    for i, j in pairs:
        if i == j:
            row.append(left_vector[i]*right_vector[j])
        else:
            row.append(
                left_vector[i]*right_vector[j]
                + left_vector[j]*right_vector[i]
            )
    return row

print("\nPART II. CHARACTER-SECTOR MAP", flush=True)
print("-"*79, flush=True)
print("For any one-sided specialization vectors h^(chi_1), h^(chi_2):", flush=True)
print("  Z_(chi_1,chi_2) = E(chi_1,chi_2) * c", flush=True)
print("where c is the same 21-component vector in every sector.", flush=True)
print("", flush=True)
print("Parked sectors:", flush=True)
print("  (omega,omega)   : E(omega,omega) * c", flush=True)
print("  (omega,omega^2) : E(omega,omega^2) * c", flush=True)
print("They are not assigned separate coefficient vectors.", flush=True)

# Symbolic formula strings, suitable for a paper or downstream code generator.
sector_formula_terms = {}
for sector_name, left_tag, right_tag in [
    ("(omega,omega)", "w", "w"),
    ("(omega,omega^2)", "w", "w2"),
]:
    terms = []
    for name, (i, j) in zip(coefficient_names, pairs):
        if i == j:
            factor = "h{0}_{1}*h{0}_{2}".format(i + 1, left_tag, right_tag)
            # The formatter above cannot place two tags; replace explicitly.
            factor = "h{}_{}*h{}_{}".format(i + 1, left_tag, j + 1, right_tag)
        else:
            factor = "h{}_{}*h{}_{} + h{}_{}*h{}_{}".format(
                i + 1, left_tag, j + 1, right_tag,
                j + 1, left_tag, i + 1, right_tag,
            )
        terms.append("{}*({})".format(name, factor))
    sector_formula_terms[sector_name] = terms

# -----------------------------------------------------------------------------
# Template for the next q-series cell.
# -----------------------------------------------------------------------------

# Expected data schema:
# {
#   "q_order": N,
#   "characters": ["one", "omega", "omega2"],
#   "basis_labels": [... six labels ...],
#   "series": {
#       "one":   {label: [a_0,...,a_N]},
#       "omega": {label: [a_0,...,a_N]},
#       "omega2":{label: [a_0,...,a_N]},
#   },
#   "coefficient_ring": "QQ"  # initially
# }

template = {
    "q_order": None,
    "characters": ["one", "omega", "omega2"],
    "basis_labels": one_sided_basis,
    "series": {
        character: {label: None for label in one_sided_basis}
        for character in ["one", "omega", "omega2"]
    },
    "coefficient_ring": "QQ",
    "notes": (
        "Replace each None by a coefficient list [q^0,...,q^N]. "
        "All six basis series must use the same N."
    ),
}

if not os.path.exists(TEMPLATE_PATH):
    save(template, TEMPLATE_PATH)
    print("\nWrote q-series input template:", flush=True)
    print("  {}".format(TEMPLATE_PATH), flush=True)
else:
    print("\nExisting q-series template retained:", flush=True)
    print("  {}".format(TEMPLATE_PATH), flush=True)

# -----------------------------------------------------------------------------
# Optional data mode: build the exact q-coefficient evaluation matrix.
# -----------------------------------------------------------------------------

def convolution_truncated(a, b, N):
    output = [QQ(0)]*(N + 1)
    for r in range(N + 1):
        total = QQ(0)
        for k in range(r + 1):
            total += a[k]*b[r-k]
        output[r] = total
    return output


def symmetric_product_series(left_series, right_series, i, j, N):
    a_i = left_series[one_sided_basis[i]]
    b_j = right_series[one_sided_basis[j]]
    first = convolution_truncated(a_i, b_j, N)
    if i == j:
        return first
    a_j = left_series[one_sided_basis[j]]
    b_i = right_series[one_sided_basis[i]]
    second = convolution_truncated(a_j, b_i, N)
    return [first[r] + second[r] for r in range(N + 1)]

rank_report = None
if os.path.exists(DATA_PATH):
    print("\nPART III. EXACT q-SERIES EVALUATION RANK", flush=True)
    print("-"*79, flush=True)
    data = load(DATA_PATH)
    N = ZZ(data["q_order"])
    if data["basis_labels"] != one_sided_basis:
        raise RuntimeError("q-series data basis order does not match Cell 31")

    series = data["series"]
    for character in ["one", "omega", "omega2"]:
        if character not in series:
            raise RuntimeError("missing character {}".format(character))
        for label in one_sided_basis:
            coefficients = series[character][label]
            if coefficients is None or len(coefficients) != N + 1:
                raise RuntimeError(
                    "bad coefficient list for {} / {}".format(character, label)
                )

    sector_specs = [
        ("(omega,omega)", "omega", "omega"),
        ("(omega,omega2)", "omega", "omega2"),
    ]

    rows = []
    row_labels = []
    sector_basis_series = {}
    for sector_name, left_character, right_character in sector_specs:
        basis_columns = []
        for i, j in pairs:
            basis_columns.append(
                symmetric_product_series(
                    series[left_character],
                    series[right_character],
                    i, j, N,
                )
            )
        sector_basis_series[sector_name] = basis_columns
        for q_power in range(N + 1):
            rows.append([column[q_power] for column in basis_columns])
            row_labels.append((sector_name, ZZ(q_power)))

    evaluation_matrix = matrix(QQ, rows)
    rank_value = ZZ(evaluation_matrix.rank())
    nullity = shared_dimension - rank_value
    rank_report = {
        "q_order": N,
        "row_labels": row_labels,
        "matrix": evaluation_matrix,
        "rank": rank_value,
        "nullity": nullity,
        "sector_basis_series": sector_basis_series,
    }

    print("q-order used = {}".format(N), flush=True)
    print("combined sector equations = {}".format(evaluation_matrix.nrows()), flush=True)
    print("shared unknowns = {}".format(shared_dimension), flush=True)
    print("exact evaluation rank = {}".format(rank_value), flush=True)
    print("remaining ambiguity dimension = {}".format(nullity), flush=True)
else:
    print("\nPART III. q-SERIES DATA NOT YET PRESENT", flush=True)
    print("-"*79, flush=True)
    print("Expected input:", flush=True)
    print("  {}".format(DATA_PATH), flush=True)
    print("Cell 32 must construct the six one-sided torsion-specialized", flush=True)
    print("basis series before the sector evaluation rank can be measured.", flush=True)

# -----------------------------------------------------------------------------
# Save model and summary.
# -----------------------------------------------------------------------------

model = {
    "one_sided_basis": one_sided_basis,
    "one_sided_dimension": ZZ(n),
    "symmetric_basis": symmetric_basis,
    "shared_dimension": shared_dimension,
    "coefficient_names": coefficient_names,
    "pairs": pairs,
    "parked_sectors": ["(omega,omega)", "(omega,omega2)"],
    "sector_formula_terms": sector_formula_terms,
    "independent_sector_overcount": ZZ(42),
    "correct_shared_unknown_count": shared_dimension,
    "q_series_data_path": DATA_PATH,
    "q_series_template_path": TEMPLATE_PATH,
    "rank_report": rank_report,
}
save(model, MODEL_PATH)

summary_lines = [
    "DEGREE-THREE SHARED BI-JACOBI TORSION-SECTOR MODEL",
    "one-sided basis dimension: 6",
    "symmetric shared coefficient dimension: 21",
    "independent two-sector count 42: rejected for the shared Schoen potential",
    "parked sectors are two evaluations of one coefficient vector c",
    "",
    "sector maps:",
    "  Z_(omega,omega)   = E_(omega,omega) c",
    "  Z_(omega,omega^2) = E_(omega,omega^2) c",
    "",
]
if rank_report is None:
    summary_lines.extend([
        "q-series evaluation rank: pending",
        "next input: results/degree3_one_sided_torsion_basis.sobj",
        "next construction: A3 and B3 HSS torsion specializations",
    ])
else:
    summary_lines.extend([
        "q-order: {}".format(rank_report["q_order"]),
        "combined exact evaluation rank: {}".format(rank_report["rank"]),
        "remaining shared ambiguity: {}".format(rank_report["nullity"]),
    ])

with open(SUMMARY_PATH, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

print("\nSaved:", flush=True)
print("  {}".format(MODEL_PATH), flush=True)
print("  {}".format(SUMMARY_PATH), flush=True)
print("\nCELL 31 COMPLETE", flush=True)
