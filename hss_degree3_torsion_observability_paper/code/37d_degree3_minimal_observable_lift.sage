# -*- coding: utf-8 -*-
# =============================================================================
# CELL 37D: DEGREE-THREE MINIMAL LOCAL-OBSERVABLE LIFT CENSUS
# =============================================================================
#
# Input
# -----
#   results/degree3_information_loss_cell37c.sobj
#
# Purpose
# -------
# Cell 37C found
#
#   rank(curvature primitives) = 14 / 15,
#   rank(final scalar defect)   =  2 / 15.
#
# This cell determines the smallest geometrically meaningful local observable
# that lifts the fourteen-dimensional curvature-primitives stack to full rank.
# It tests, in order:
#
#   * frame-metric components g_ee,g_eo,g_er,g_oo,g_or,g_rr;
#   * coordinate metric entries G00,G01,G02,G11,G12,G22;
#   * connection-vector entries A0,A1,A2;
#   * individual Kähler-jet entries;
#   * individual raw CM-tensor entries.
#
# It also reconstructs the unique numerical null direction of the curvature-
# primitives matrix and verifies that every reported lift responds nontrivially
# to that direction.  No transport law is assumed here: this is an information
# census used to choose the next physical observable family non-circularly.
# =============================================================================

from sage.all import *
from itertools import combinations
import os
import sys

print("="*79, flush=True)
print("CELL 37D: DEGREE-THREE MINIMAL LOCAL-OBSERVABLE LIFT CENSUS", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

INPUT_PATH = os.path.join(RESULTS_DIR, "degree3_information_loss_cell37c.sobj")
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_minimal_observable_lift_cell37d.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_minimal_observable_lift_cell37d.txt")
TEMPLATE_OUT = os.path.join(RESULTS_DIR, "degree3_selected_observable_stack_cell37d.sobj")

if not os.path.exists(INPUT_PATH):
    raise IOError("Missing Cell 37C result: {}".format(INPUT_PATH))

cert = load(INPUT_PATH)
PREC_BITS = ZZ(cert["precision_bits"])
CC = ComplexField(PREC_BITS)
RR = RealField(PREC_BITS)
responses = cert["responses"]
labels15 = list(cert["observable_labels15"])
sector_names = list(cert["sector_names"])
EXP_MIN, EXP_MAX = map(ZZ, cert["census_range"])

if len(labels15) != 15:
    raise ArithmeticError("Expected fifteen ambiguity labels")
if len(sector_names) != 2:
    raise ArithmeticError("Expected two parked sectors")

# Use the same three tolerance levels as Cell 37C.
tolerances = {
    "1e-30": RR(10)^(-30),
    "1e-40": RR(10)^(-40),
    "1e-50": RR(10)^(-50),
}

print("\nPART I. INPUT AUDIT", flush=True)
print("-"*79, flush=True)
print("  precision bits      : {}".format(PREC_BITS), flush=True)
print("  Laurent window      : [{} , {}]".format(EXP_MIN, EXP_MAX), flush=True)
print("  ambiguity columns   : {}".format(len(labels15)), flush=True)
print("  parked sectors      : {}".format(sector_names), flush=True)


def component_order(block):
    if "component_orders" in cert and block in cert["component_orders"]:
        return list(cert["component_orders"][block])
    names = set()
    for sector in sector_names:
        for label in labels15:
            names.update(responses[sector][label][block].keys())
    return sorted(names)


def component_columns(block, component, sectors=None, lower=None, upper=None):
    """Fifteen column vectors for one stored component."""
    if sectors is None:
        sectors = sector_names
    if lower is None:
        lower = EXP_MIN
    if upper is None:
        upper = EXP_MAX
    columns = []
    for label in labels15:
        vec = []
        for sector in sectors:
            data = responses[sector][label][block][component]
            if block == "raw_tensor":
                vec.append(CC(data))
            else:
                vec.extend(CC(data.get(ZZ(e), 0)) for e in range(ZZ(lower), ZZ(upper)+1))
        columns.append(vec)
    return columns


def stack_column_families(families):
    """Concatenate row blocks, preserving the fifteen common columns."""
    out = [[] for _ in labels15]
    for family in families:
        if len(family) != len(labels15):
            raise ArithmeticError("Column-family width mismatch")
        for j in range(len(labels15)):
            out[j].extend(family[j])
    return out


def numerical_rank_columns(columns, tolerance, return_pivots=False):
    """Modified Gram-Schmidt rank after independent column normalization."""
    basis = []
    pivots = []
    residual_norms = []
    for idx, raw in enumerate(columns):
        v = vector(CC, raw)
        norm0 = sqrt(sum((abs(entry)^2 for entry in v), RR(0)))
        if norm0 == 0:
            continue
        v = v/CC(norm0)
        for q in basis:
            inner = sum((q[j].conjugate()*v[j] for j in range(len(v))), CC(0))
            v -= inner*q
        norm = sqrt(sum((abs(entry)^2 for entry in v), RR(0)))
        if norm > tolerance:
            basis.append(v/CC(norm))
            pivots.append(ZZ(idx))
            residual_norms.append(RR(norm))
    if return_pivots:
        return ZZ(len(basis)), pivots, residual_norms
    return ZZ(len(basis))


def columns_to_matrix(columns):
    if not columns:
        return matrix(CC, 0, len(labels15))
    nrows = len(columns[0])
    return matrix(CC, nrows, len(columns), lambda i,j: CC(columns[j][i]))


def vector_norm(v):
    return RR(sqrt(sum((abs(x)^2 for x in v), RR(0))))


# -----------------------------------------------------------------------------
# Reconstruct the frame-metric candidates from the stored coordinate metric.
# We keep two orders of safety at both Laurent boundaries because the e/r frame
# vectors are linear in U and their bilinears can shift by two powers.
# -----------------------------------------------------------------------------
LS = LaurentSeriesRing(CC, names=("u",), default_prec=int(EXP_MAX+8))
u = LS.gen()
X = (u-1)/3
FRAME_MIN = EXP_MIN + 2
FRAME_MAX = EXP_MAX - 2


def dict_to_series(data):
    return LS({ZZ(e):CC(v) for e,v in data.items()})


def metric_matrix(sector, label):
    block = responses[sector][label]["metric"]
    G00 = dict_to_series(block["G00"])
    G01 = dict_to_series(block["G01"])
    G02 = dict_to_series(block["G02"])
    G11 = dict_to_series(block["G11"])
    G12 = dict_to_series(block["G12"])
    G22 = dict_to_series(block["G22"])
    return [[G00,G01,G02],[G01,G11,G12],[G02,G12,G22]]


def bilinear(G, left, right):
    return sum((left[i]*G[i][j]*right[j] for i in range(3) for j in range(3)), LS(0))


evec = [-(2*X+1), LS(1), LS(1)]
ovec = [LS(0), LS(1), LS(-1)]
rvec = [X, LS(1), LS(1)]
frame_pairs = {
    "g_ee": (evec,evec),
    "g_eo": (evec,ovec),
    "g_er": (evec,rvec),
    "g_oo": (ovec,ovec),
    "g_or": (ovec,rvec),
    "g_rr": (rvec,rvec),
}


def frame_component_columns(name, sectors=None):
    if sectors is None:
        sectors = sector_names
    left,right = frame_pairs[name]
    columns = []
    for label in labels15:
        vec = []
        for sector in sectors:
            value = bilinear(metric_matrix(sector,label),left,right)
            vec.extend(CC(value[e]) for e in range(FRAME_MIN,FRAME_MAX+1))
        columns.append(vec)
    return columns


# Base curvature-primitives stack.
primitive_components = component_order("curvature_primitives")
primitive_families = [
    component_columns("curvature_primitives", comp)
    for comp in primitive_components
]
base_columns = stack_column_families(primitive_families)

print("\nPART II. FOURTEEN-DIMENSIONAL CURVATURE-PRIMITIVE BASE", flush=True)
print("-"*79, flush=True)
print("  primitive components: {}".format(primitive_components), flush=True)
base_ranks = {
    key:numerical_rank_columns(base_columns,tol)
    for key,tol in tolerances.items()
}
print("  ranks by tolerance: {}".format(base_ranks), flush=True)
if any(rank != 14 for rank in base_ranks.values()):
    raise ArithmeticError("Curvature-primitives rank is not stably fourteen")

# Search for the smallest primitive component subset that already has rank 14.
minimal_primitive_subsets = []
minimal_primitive_size = None
for size in range(1,len(primitive_components)+1):
    found = []
    for subset in combinations(primitive_components,size):
        cols = stack_column_families([
            component_columns("curvature_primitives",comp) for comp in subset
        ])
        ranks = [numerical_rank_columns(cols,tol) for tol in tolerances.values()]
        if all(rank == 14 for rank in ranks):
            found.append(tuple(subset))
    if found:
        minimal_primitive_size = ZZ(size)
        minimal_primitive_subsets = found
        break
print("  smallest component count retaining rank 14: {}".format(minimal_primitive_size), flush=True)
for subset in minimal_primitive_subsets:
    print("    {}".format(subset), flush=True)


# -----------------------------------------------------------------------------
# Unique null direction of the primitive stack.
# -----------------------------------------------------------------------------
print("\nPART III. UNIQUE PRIMITIVE-STACK NULL DIRECTION", flush=True)
print("-"*79, flush=True)
Mbase = columns_to_matrix(base_columns)
rank14,pivot_cols,_ = numerical_rank_columns(base_columns,tolerances["1e-50"],True)
if rank14 != 14:
    raise ArithmeticError("Could not identify fourteen independent primitive columns")
dependent_cols = [j for j in range(15) if j not in pivot_cols]
if len(dependent_cols) != 1:
    raise ArithmeticError("Expected exactly one dependent ambiguity column")
dep = dependent_cols[0]
P = Mbase.matrix_from_columns([int(j) for j in pivot_cols])

# Select fourteen independent rows of P using the same normalized MGS routine.
row_vectors = [list(P.row(i)) for i in range(P.nrows())]
row_rank,row_pivots,_ = numerical_rank_columns(row_vectors,tolerances["1e-50"],True)
if row_rank != 14:
    raise ArithmeticError("Could not select fourteen independent primitive rows")
row_ids = [int(i) for i in row_pivots[:14]]
A = P.matrix_from_rows(row_ids)
b = vector(CC,[Mbase[i,dep] for i in row_ids])
coeffs = A.solve_right(b)
nullvec = vector(CC,15)
nullvec[dep] = CC(1)
for k,j in enumerate(pivot_cols):
    nullvec[int(j)] = -coeffs[k]
scale = max(abs(x) for x in nullvec)
nullvec = nullvec/CC(scale)
null_residual = vector_norm(Mbase*nullvec)
print("  dependent column: {}".format(labels15[dep]), flush=True)
print("  normalized null residual: {}".format(null_residual), flush=True)
print("  visible null-vector support:", flush=True)
for label,value in zip(labels15,nullvec):
    if abs(value) > RR(10)^(-25):
        print("    {:45s} {}".format(label,value), flush=True)


# -----------------------------------------------------------------------------
# Candidate local observables.
# -----------------------------------------------------------------------------
print("\nPART IV. SINGLE-OBSERVABLE LIFT CENSUS", flush=True)
print("-"*79, flush=True)

candidates = []

# Prefer frame-covariant metric data, especially the off-diagonal mixing modes.
for name in ["g_eo","g_er","g_or","g_ee","g_oo","g_rr"]:
    candidates.append({
        "family":"frame_metric",
        "name":name,
        "columns":frame_component_columns(name),
        "window":(FRAME_MIN,FRAME_MAX),
    })

# Coordinate metric entries.
for name in component_order("metric"):
    candidates.append({
        "family":"coordinate_metric",
        "name":name,
        "columns":component_columns("metric",name),
        "window":(EXP_MIN,EXP_MAX),
    })

# Connection-vector entries.
for name in component_order("connection_vector"):
    candidates.append({
        "family":"connection_vector",
        "name":name,
        "columns":component_columns("connection_vector",name),
        "window":(EXP_MIN,EXP_MAX),
    })

# Individual Kähler-jet and raw-tensor entries are included as fallback audits.
for name in component_order("kahler_jet"):
    candidates.append({
        "family":"kahler_jet",
        "name":name,
        "columns":component_columns("kahler_jet",name),
        "window":(EXP_MIN,EXP_MAX),
    })
for name in component_order("raw_tensor"):
    candidates.append({
        "family":"raw_tensor",
        "name":name,
        "columns":component_columns("raw_tensor",name),
        "window":None,
    })

lift_results = []
for candidate in candidates:
    augmented = stack_column_families([base_columns,candidate["columns"]])
    ranks = {
        key:numerical_rank_columns(augmented,tol)
        for key,tol in tolerances.items()
    }
    C = columns_to_matrix(candidate["columns"])
    null_response = vector_norm(C*nullvec)
    null_scale = max([vector_norm(C.column(j)) for j in range(C.ncols())] + [RR(1)])
    relative_null_response = null_response/null_scale
    stable_lift = all(rank == 15 for rank in ranks.values())
    result = {
        "family":candidate["family"],
        "name":candidate["name"],
        "window":candidate["window"],
        "ranks":ranks,
        "stable_full_rank_lift":bool(stable_lift),
        "null_response_norm":null_response,
        "relative_null_response":relative_null_response,
    }
    lift_results.append(result)
    if stable_lift:
        print("  LIFT {:18s} {:18s} ranks={} null-response={}".format(
            candidate["family"],candidate["name"],ranks,relative_null_response
        ), flush=True)

stable_lifts = [entry for entry in lift_results if entry["stable_full_rank_lift"]]
if not stable_lifts:
    print("  no single stored local observable lifts rank 14 to 15", flush=True)

# Preferred choice: first stable candidate in the physically motivated ordering,
# with null-response magnitude used only to break ties within one family.
family_priority = {
    "frame_metric":0,
    "coordinate_metric":1,
    "connection_vector":2,
    "kahler_jet":3,
    "raw_tensor":4,
}
preferred = None
if stable_lifts:
    stable_lifts_sorted = sorted(
        stable_lifts,
        key=lambda entry:(
            family_priority[entry["family"]],
            -RR(entry["relative_null_response"]),
            entry["name"],
        )
    )
    preferred = stable_lifts_sorted[0]


# -----------------------------------------------------------------------------
# Sectorwise and cutoff stability for the preferred lift.
# -----------------------------------------------------------------------------
print("\nPART V. PREFERRED MINIMAL STACK", flush=True)
print("-"*79, flush=True)
preferred_audit = None
if preferred is None:
    print("  no preferred single-component lift available", flush=True)
else:
    print("  selected family/component: {}/{}".format(preferred["family"],preferred["name"]), flush=True)
    # Recover its columns.
    match = next(c for c in candidates if c["family"] == preferred["family"] and c["name"] == preferred["name"])
    selected_columns = match["columns"]
    combined_columns = stack_column_families([base_columns,selected_columns])
    sector_ranks = {}
    for sector in sector_names:
        primitive_sector = stack_column_families([
            component_columns("curvature_primitives",comp,[sector])
            for comp in primitive_components
        ])
        if preferred["family"] == "frame_metric":
            candidate_sector = frame_component_columns(preferred["name"],[sector])
        else:
            stored_block = {
                "coordinate_metric":"metric",
                "connection_vector":"connection_vector",
                "kahler_jet":"kahler_jet",
                "raw_tensor":"raw_tensor",
            }[preferred["family"]]
            candidate_sector = component_columns(stored_block,preferred["name"],[sector])
        sector_ranks[sector] = {
            "base":numerical_rank_columns(primitive_sector,tolerances["1e-40"]),
            "augmented":numerical_rank_columns(
                stack_column_families([primitive_sector,candidate_sector]),
                tolerances["1e-40"],
            ),
        }
    preferred_audit = {
        "selection":preferred,
        "combined_ranks":{
            key:numerical_rank_columns(combined_columns,tol)
            for key,tol in tolerances.items()
        },
        "sector_ranks_1e40":sector_ranks,
    }
    print("  combined ranks: {}".format(preferred_audit["combined_ranks"]), flush=True)
    print("  sector ranks: {}".format(sector_ranks), flush=True)
    print("  full fifteen-dimensional information recovered? {}".format(
        all(rank == 15 for rank in preferred_audit["combined_ranks"].values())
    ), flush=True)


# -----------------------------------------------------------------------------
# Verdict and saved interface for the next physical-transport cell.
# -----------------------------------------------------------------------------
print("\nPART VI. VERDICT", flush=True)
print("-"*79, flush=True)
print("  curvature-primitives rank: 14 / 15", flush=True)
print("  stable single-component lifts found: {}".format(len(stable_lifts)), flush=True)
if preferred:
    print("  preferred lift: {}/{}".format(preferred["family"],preferred["name"]), flush=True)
    print("  next physical question: derive a transport/covariance law for this added local observable", flush=True)
else:
    print("  next physical question: search minimal two-component local stacks", flush=True)

certificate = {
    "schema_version":ZZ(1),
    "scope":"high-precision minimal local-observable lift at the square-lattice CM point",
    "source_cell37c":INPUT_PATH,
    "precision_bits":PREC_BITS,
    "census_range":(EXP_MIN,EXP_MAX),
    "frame_safe_range":(FRAME_MIN,FRAME_MAX),
    "observable_labels15":labels15,
    "sector_names":sector_names,
    "primitive_components":primitive_components,
    "primitive_ranks":base_ranks,
    "minimal_primitive_component_count":minimal_primitive_size,
    "minimal_primitive_subsets":minimal_primitive_subsets,
    "primitive_null_vector":list(nullvec),
    "primitive_null_residual":null_residual,
    "lift_results":lift_results,
    "stable_lifts":stable_lifts,
    "preferred_lift":preferred,
    "preferred_audit":preferred_audit,
}
save(certificate,OUT_SOBJ)

save({
    "schema_version":ZZ(1),
    "source_cell37c":INPUT_PATH,
    "observable_labels15":labels15,
    "sector_names":sector_names,
    "base_observable":{
        "block":"curvature_primitives",
        "components":primitive_components,
    },
    "added_observable":preferred,
    "full_rank_expected":bool(preferred is not None),
    "next_required_step":"derive the geometrically admissible transport law for the selected local observable",
},TEMPLATE_OUT)

summary_lines = [
    "CELL 37D: DEGREE-THREE MINIMAL LOCAL-OBSERVABLE LIFT CENSUS",
    "precision bits: {}".format(PREC_BITS),
    "Laurent window: [{} , {}]".format(EXP_MIN,EXP_MAX),
    "curvature-primitives rank: 14 / 15",
    "minimal primitive component count retaining rank 14: {}".format(minimal_primitive_size),
    "primitive null residual: {}".format(null_residual),
    "stable single-component lifts: {}".format(len(stable_lifts)),
]
lift_counts_by_family = {}
for entry in stable_lifts:
    lift_counts_by_family[entry["family"]] = lift_counts_by_family.get(entry["family"],ZZ(0)) + 1
summary_lines.append("stable lifts by family: {}".format(lift_counts_by_family))
for entry in stable_lifts:
    if entry["family"] in ["frame_metric","coordinate_metric","connection_vector"]:
        summary_lines.append("  {}/{} -> 15/15".format(entry["family"],entry["name"]))
if preferred:
    summary_lines.extend([
        "preferred lift: {}/{}".format(preferred["family"],preferred["name"]),
        "preferred combined ranks: {}".format(preferred_audit["combined_ranks"]),
        "conclusion: one additional local observable restores all fifteen ambiguity directions",
        "next: derive a physical transport/covariance equation for the selected observable",
    ])
else:
    summary_lines.extend([
        "preferred lift: none",
        "conclusion: no single stored component restores full rank",
        "next: search minimal two-component observable stacks",
    ])

with open(OUT_TXT,"w") as handle:
    handle.write("\n".join(summary_lines)+"\n")

print("\nFINAL SUMMARY", flush=True)
print("-"*79, flush=True)
for line in summary_lines:
    print(line, flush=True)
print("\nSaved:", flush=True)
print("  {}".format(OUT_SOBJ), flush=True)
print("  {}".format(OUT_TXT), flush=True)
print("  {}".format(TEMPLATE_OUT), flush=True)
