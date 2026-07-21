# -*- coding: utf-8 -*-
# =============================================================================
# CELL 36: DEGREE-TWO FULL-SPACE BASELINE AND DEGREE-THREE SELECTIVITY LEDGER
# =============================================================================
#
# Purpose
# -------
# 1. Prove that the degree-two one-sided holomorphic ambiguity at genus zero
#    has dimension two, and that {P2, B(q^2)} is an exact independent HSS
#    basis.  Hence {FA,FB,FC} is the complete symmetric degree-two ambiguity.
#
# 2. Record the exact degree-two zero-selectivity baseline:
#       dim(connection)=dim(defect Laurent support)=rank(transport)=9,
#    so bare closure imposes no equation on the complete three-dimensional
#    degree-two ambiguity space.
#
# 3. Load the audited all-orders degree-three torsion relation and construct
#    the exact fifteen-coordinate observable basis.
#
# 4. Write a machine-readable matrix interface for the future degree-three
#    geometry workers.  Once matrices T (admissible connection response) and
#    G (fifteen ambiguity response columns) are supplied, the genuinely
#    selective rank is
#
#       s = rank([T | G]) - rank(T).
#
#    This cell intentionally does not invent a degree-three connection module.
# =============================================================================

from sage.all import *
import os

print("="*79, flush=True)
print("CELL 36: FULL-SPACE BASELINE AND DEGREE-THREE SELECTIVITY LEDGER", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

OUT_SOBJ = os.path.join(
    RESULTS_DIR,
    "degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj",
)
OUT_TXT = os.path.join(
    RESULTS_DIR,
    "degree2_full_space_baseline_and_degree3_selectivity_ledger.txt",
)
TEMPLATE_PATH = os.path.join(
    RESULTS_DIR,
    "degree3_transport_selectivity_input_template.sobj",
)
MATRIX_INPUT_PATH = os.path.join(
    RESULTS_DIR,
    "degree3_transport_selectivity_matrices.sobj",
)

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

def modular_monomials(weight):
    """Return labels for level-one modular monomials of exact weight."""
    weight = ZZ(weight)
    if weight < 0 or weight % 2:
        return []
    labels = []
    for a in range(weight // 4 + 1):
        for b in range(weight // 6 + 1):
            if 4*a + 6*b != weight:
                continue
            factors = []
            if a:
                factors.append("E4" if a == 1 else "E4^{}".format(a))
            if b:
                factors.append("E6" if b == 1 else "E6^{}".format(b))
            labels.append("*".join(factors) if factors else "1")
    return labels


def stack_columns(left, right):
    """Horizontal concatenation over the common exact base ring."""
    if left.nrows() != right.nrows():
        raise ValueError("T and G must have the same number of rows")
    K = left.base_ring()
    if right.base_ring() != K:
        right = matrix(K, right)
    return left.augment(right)


def selectivity_report(T, G, label):
    """Compute the exact quotient rank rank([T|G])-rank(T)."""
    if not hasattr(T, "rank") or not hasattr(G, "rank"):
        raise TypeError("{} matrices must be exact Sage matrices".format(label))
    if T.nrows() != G.nrows():
        raise ValueError("{}: T and G row counts differ".format(label))
    if G.ncols() != 15:
        raise ValueError("{}: G must have exactly 15 ambiguity columns".format(label))
    rank_T = ZZ(T.rank())
    rank_union = ZZ(stack_columns(T, G).rank())
    selective_rank = ZZ(rank_union - rank_T)
    if selective_rank < 0 or selective_rank > 15:
        raise ArithmeticError("{}: impossible selective rank".format(label))
    return {
        "label": label,
        "defect_coordinate_count": ZZ(T.nrows()),
        "connection_column_count": ZZ(T.ncols()),
        "ambiguity_column_count": ZZ(G.ncols()),
        "rank_T": rank_T,
        "rank_TG": rank_union,
        "selective_rank": selective_rank,
        "remaining_ambiguity_dimension": ZZ(15 - selective_rank),
        "unique_observable_reconstruction": bool(selective_rank == 15),
    }

# -----------------------------------------------------------------------------
# Part I. Complete degree-two one-sided ambiguity.
# -----------------------------------------------------------------------------

print("\nPART I. COMPLETE DEGREE-TWO ONE-SIDED AMBIGUITY", flush=True)
print("-"*79, flush=True)

# Genus zero, base degree two: weight = -2 + 6*2 = 10, index 2.
weight2 = ZZ(10)
index2 = ZZ(2)
index2_module = [
    ("A2", ZZ(4)),
    ("B2", ZZ(6)),
    ("A1^2", ZZ(8)),
]

one_sided_degree2 = []
for generator, generator_weight in index2_module:
    for modular_factor in modular_monomials(weight2 - generator_weight):
        label = generator if modular_factor == "1" else modular_factor + "*" + generator
        one_sided_degree2.append(label)

expected_degree2 = ["E6*A2", "E4*B2"]
if one_sided_degree2 != expected_degree2:
    raise ArithmeticError(
        "Unexpected weight-10 index-2 basis: {}".format(one_sided_degree2)
    )

print("  target weight/index: {}/{}".format(weight2, index2), flush=True)
print("  module enumeration:", flush=True)
for label in one_sided_degree2:
    print("    {}".format(label), flush=True)
print("  A1^2 contribution absent because M_2(SL2Z)=0", flush=True)
print("  exact one-sided dimension: 2", flush=True)

# -----------------------------------------------------------------------------
# Part II. Exact HSS q-series independence of P2 and B(q^2).
# -----------------------------------------------------------------------------

print("\nPART II. EXACT HSS BASIS {P2, B(q^2)}", flush=True)
print("-"*79, flush=True)

PREC = ZZ(12)
PS = PowerSeriesRing(QQ, "q36", default_prec=PREC + 1)
q = PS.gen()

B = PS(9)
for n in range(1, PREC + 1):
    B *= (1 - q**n)**(-4)
B = B.add_bigoh(PREC + 1)

E2 = PS(1)
for n in range(1, PREC + 1):
    E2 += -24*sigma(n, 1)*q**n
E2 = E2.add_bigoh(PREC + 1)

Z2 = (E2*B**2/72).add_bigoh(PREC + 1)
Bq2 = B(q**2).add_bigoh(PREC + 1)
P2 = (Z2 - Bq2/8).add_bigoh(PREC + 1)

# q^0 and q^1 already give an exact nonzero minor.
HSS_minor = matrix(QQ, [
    [QQ(P2[0]), QQ(Bq2[0])],
    [QQ(P2[1]), QQ(Bq2[1])],
])
HSS_minor_det = QQ(HSS_minor.det())
HSS_rank = ZZ(HSS_minor.rank())
if HSS_rank != 2 or HSS_minor_det == 0:
    raise ArithmeticError("P2 and B(q^2) are not independent")

print("  P2(q) begins       : {}".format(P2), flush=True)
print("  B(q^2) begins      : {}".format(Bq2), flush=True)
print("  q^0/q^1 minor      : {}".format(HSS_minor), flush=True)
print("  minor determinant  : {}".format(HSS_minor_det), flush=True)
print("  HSS one-sided rank : {} / 2".format(HSS_rank), flush=True)

symmetric_degree2_basis = [
    "FA=P2(y)P2(z)",
    "FC=P2(y)B(q_z^2)+B(q_y^2)P2(z)",
    "FB=B(q_y^2)B(q_z^2)",
]
print("  complete symmetric square dimension: 3", flush=True)
for label in symmetric_degree2_basis:
    print("    {}".format(label), flush=True)

# -----------------------------------------------------------------------------
# Part III. Exact zero-selectivity baseline at degree two.
# -----------------------------------------------------------------------------

print("\nPART III. DEGREE-TWO BARE-TRANSPORT SELECTIVITY", flush=True)
print("-"*79, flush=True)

# Exact result from the transfer determinant certificate.
defect_support2 = list(range(-5, 4))
defect_dimension2 = ZZ(len(defect_support2))
connection_dimension2 = ZZ(9)
transport_rank2 = ZZ(9)
transport_quotient_dimension2 = ZZ(defect_dimension2 - transport_rank2)
bare_selective_rank2 = ZZ(0)

if defect_dimension2 != 9 or transport_rank2 != 9:
    raise ArithmeticError("Degree-two transport ledger is inconsistent")

print("  defect Laurent support: U^-5 through U^3", flush=True)
print("  defect-space dimension: {}".format(defect_dimension2), flush=True)
print("  connection dimension   : {}".format(connection_dimension2), flush=True)
print("  exact transport rank   : {}".format(transport_rank2), flush=True)
print("  quotient dimension     : {}".format(transport_quotient_dimension2), flush=True)
print("  complete ambiguity dim : 3", flush=True)
print("  bare selective rank    : {}".format(bare_selective_rank2), flush=True)
print("  conclusion             : bare closure absorbs the full admissible space", flush=True)

# Optional consistency check from Cell 29 when present.
cell29_path = os.path.join(RESULTS_DIR, "degree2_exact_polar_dependency_certificate.sobj")
cell29_loaded = False
cell29_unresolved = None
if os.path.exists(cell29_path):
    cell29 = load(cell29_path)
    cell29_unresolved = list(cell29.get("unresolved_subsets", []))
    if cell29_unresolved:
        raise ArithmeticError("Cell 29 still reports unresolved polar subsets")
    cell29_loaded = True
    print("  Cell-29 polar calibration loaded: unresolved subsets = 0", flush=True)
else:
    print("  Cell-29 result not loaded; zero-selectivity result is independent of it", flush=True)

# -----------------------------------------------------------------------------
# Part IV. Audited degree-three observable space.
# -----------------------------------------------------------------------------

print("\nPART IV. AUDITED DEGREE-THREE OBSERVABLE SPACE", flush=True)
print("-"*79, flush=True)

audit_path = os.path.join(RESULTS_DIR, "degree3_torsion_modularity_audit_35b.sobj")
if not os.path.exists(audit_path):
    raise IOError("Missing Cell-35b audit result: {}".format(audit_path))
audit = load(audit_path)
if not bool(audit.get("all_orders_pass", False)):
    raise ArithmeticError("Cell 35b did not pass")
if bool(audit.get("global_jacobi_identity", True)):
    raise ArithmeticError("The torsion relation was incorrectly marked global")

relation_labels = list(audit["relation_labels"])
relation_vector = vector(QQ, audit["relation_vector"])
if relation_vector != vector(QQ, [3, -3, 0, -6, -10, 16]):
    raise ArithmeticError("Unexpected all-orders torsion relation")

quotient_path = os.path.join(RESULTS_DIR, "degree3_certified_15d_quotient_v3.sobj")
quotient_loaded = False
if os.path.exists(quotient_path):
    quotient = load(quotient_path)
    if ZZ(quotient.get("quotient_rank", -1)) != 15:
        raise ArithmeticError("Cell 34 quotient rank is not 15")
    if ZZ(quotient.get("parked_rank15", -1)) != 15:
        raise ArithmeticError("Parked sectors are not faithful on the quotient")
    if not bool(quotient.get("kernels_equal", False)):
        raise ArithmeticError("Cell 34 quotient and parked kernels differ")
    pivot_labels = list(quotient["pivot_labels5"])
    quotient_loaded = True
else:
    # Canonical relation has nonzero coefficient on h6, so eliminate h6.
    pivot_labels = relation_labels[:5]

if len(pivot_labels) != 5:
    raise ArithmeticError("Expected five one-sided observable basis labels")

observable_pairs = []
observable_labels15 = []
for i in range(5):
    for j in range(i, 5):
        observable_pairs.append((ZZ(i), ZZ(j)))
        if i == j:
            observable_labels15.append("g{}(y)g{}(z)".format(i + 1, j + 1))
        else:
            observable_labels15.append(
                "g{}(y)g{}(z)+g{}(y)g{}(z)".format(
                    i + 1, j + 1, j + 1, i + 1
                )
            )

if len(observable_labels15) != 15:
    raise ArithmeticError("Observable symmetric basis does not have dimension 15")

print("  Cell-35b all-orders audit: PASS", flush=True)
print("  relation scope: three-torsion HSS specialization kernel", flush=True)
print("  global Jacobi identity: False", flush=True)
print("  one-sided global basis dimension: 6", flush=True)
print("  one-sided torsion image dimension: 5", flush=True)
print("  symmetric observable dimension: 15", flush=True)
print("  reduced one-sided basis:", flush=True)
for i, label in enumerate(pivot_labels):
    print("    g{} = {}".format(i + 1, label), flush=True)
print("  Cell-34 finite matrix quotient loaded? {}".format(quotient_loaded), flush=True)

# -----------------------------------------------------------------------------
# Part V. Degree-three selectivity matrix interface.
# -----------------------------------------------------------------------------

print("\nPART V. DEGREE-THREE SELECTIVITY INTERFACE", flush=True)
print("-"*79, flush=True)

ledger_formula = "s = rank([T|G]) - rank(T)"
print("  observable ambiguity columns: 15", flush=True)
print("  exact selective-rank formula: {}".format(ledger_formula), flush=True)
print("  s=0   : transport absorbs every observable ambiguity direction", flush=True)
print("  1<=s<15: transport leaves a family of dimension 15-s", flush=True)
print("  s=15  : unique observable reconstruction is possible", flush=True)
print("  no degree-three module is assumed by this cell", flush=True)

template = {
    "schema_version": ZZ(1),
    "purpose": "degree-three exact transport selectivity matrices",
    "observable_dimension": ZZ(15),
    "observable_basis_labels": observable_labels15,
    "observable_pairs": observable_pairs,
    "one_sided_basis_labels": pivot_labels,
    "defect_coordinate_labels": None,
    "admissible_connection_column_labels": None,
    "structured_connection_column_labels": None,
    "T_admissible": None,
    "T_structured": None,
    "G_ambiguity": None,
    "source_vector": None,
    "required_formula": ledger_formula,
    "notes": (
        "All matrices must use the same exact base field and defect-coordinate "
        "row order. G_ambiguity must have exactly 15 columns. Fill this file "
        "from specialize-first degree-three geometry workers."
    ),
}
if not os.path.exists(TEMPLATE_PATH):
    save(template, TEMPLATE_PATH)
    print("  wrote input template: {}".format(TEMPLATE_PATH), flush=True)
else:
    print("  retained existing input template: {}".format(TEMPLATE_PATH), flush=True)

matrix_reports = {}
if os.path.exists(MATRIX_INPUT_PATH):
    print("\nPART VI. OPTIONAL FILLED-MATRIX SELECTIVITY REPORT", flush=True)
    print("-"*79, flush=True)
    matrix_data = load(MATRIX_INPUT_PATH)
    G = matrix_data.get("G_ambiguity", None)
    T_adm = matrix_data.get("T_admissible", None)
    T_str = matrix_data.get("T_structured", None)
    if G is None or T_adm is None:
        raise KeyError("Filled matrix input must contain G_ambiguity and T_admissible")
    matrix_reports["admissible"] = selectivity_report(T_adm, G, "admissible")
    report = matrix_reports["admissible"]
    print("  admissible rank(T)     : {}".format(report["rank_T"]), flush=True)
    print("  admissible rank([T|G]) : {}".format(report["rank_TG"]), flush=True)
    print("  bare selective rank    : {} / 15".format(report["selective_rank"]), flush=True)
    print("  remaining dimension    : {}".format(report["remaining_ambiguity_dimension"]), flush=True)
    if T_str is not None:
        matrix_reports["structured"] = selectivity_report(T_str, G, "structured")
        report = matrix_reports["structured"]
        print("  structured selective rank: {} / 15".format(
            report["selective_rank"]
        ), flush=True)
        print("  structured remaining dim : {}".format(
            report["remaining_ambiguity_dimension"]
        ), flush=True)
else:
    print("  no filled degree-three matrix file found yet", flush=True)
    print("  next producer: specialize-first degree-three response workers", flush=True)

# -----------------------------------------------------------------------------
# Save certificate and compact summary.
# -----------------------------------------------------------------------------

certificate = {
    "degree2": {
        "weight": weight2,
        "index": index2,
        "one_sided_module_basis": one_sided_degree2,
        "one_sided_dimension": ZZ(2),
        "hss_basis": ["P2", "B(q^2)"],
        "hss_minor": HSS_minor,
        "hss_minor_determinant": HSS_minor_det,
        "hss_rank": HSS_rank,
        "symmetric_basis": symmetric_degree2_basis,
        "symmetric_dimension": ZZ(3),
        "defect_support": defect_support2,
        "defect_dimension": defect_dimension2,
        "connection_dimension": connection_dimension2,
        "transport_rank": transport_rank2,
        "transport_quotient_dimension": transport_quotient_dimension2,
        "bare_selective_rank": bare_selective_rank2,
        "cell29_loaded": cell29_loaded,
        "cell29_unresolved_subsets": cell29_unresolved,
    },
    "degree3": {
        "audit_path": audit_path,
        "audit_all_orders_pass": bool(audit["all_orders_pass"]),
        "relation_scope": audit["scope"],
        "global_jacobi_identity": bool(audit["global_jacobi_identity"]),
        "relation_labels": relation_labels,
        "relation_vector": relation_vector,
        "one_sided_global_dimension": ZZ(6),
        "one_sided_observable_dimension": ZZ(5),
        "observable_symmetric_dimension": ZZ(15),
        "one_sided_observable_basis": pivot_labels,
        "observable_pairs": observable_pairs,
        "observable_basis_labels": observable_labels15,
        "quotient_path": quotient_path if quotient_loaded else None,
        "quotient_loaded": quotient_loaded,
        "selectivity_formula": ledger_formula,
        "matrix_input_path": MATRIX_INPUT_PATH,
        "matrix_reports": matrix_reports,
    },
}
save(certificate, OUT_SOBJ)

summary_lines = [
    "CELL 36: FULL-SPACE BASELINE AND DEGREE-THREE SELECTIVITY LEDGER",
    "",
    "DEGREE TWO",
    "  one-sided weight/index: 10/2",
    "  complete one-sided module dimension: 2",
    "  exact HSS basis: P2, B(q^2)",
    "  q^0/q^1 independence minor determinant: {}".format(HSS_minor_det),
    "  complete symmetric ambiguity dimension: 3",
    "  exact transport rank: 9 / 9",
    "  bare selectivity: 0",
    "  conclusion: the full admissible degree-two ambiguity is absorbed",
    "",
    "DEGREE THREE",
    "  Cell-35b all-orders torsion audit: PASS",
    "  relation is global Jacobi identity: False",
    "  relation scope: three-torsion HSS specialization kernel",
    "  global one-sided dimension: 6",
    "  observable one-sided dimension: 5",
    "  observable symmetric dimension: 15",
    "  parked sectors faithful on quotient: {}".format(quotient_loaded),
    "",
    "SELECTIVITY FORMULA",
    "  s = rank([T|G]) - rank(T)",
    "  s=0: no inverse information",
    "  s=15: unique observable reconstruction possible",
    "  filled degree-three matrices present: {}".format(bool(matrix_reports)),
]
with open(OUT_TXT, "w") as handle:
    handle.write("\n".join(summary_lines) + "\n")

print("\nFINAL SUMMARY", flush=True)
print("-"*79, flush=True)
for line in summary_lines:
    print(line, flush=True)
print("\nSaved:", flush=True)
print("  {}".format(OUT_SOBJ), flush=True)
print("  {}".format(OUT_TXT), flush=True)
print("  {}".format(TEMPLATE_PATH), flush=True)
