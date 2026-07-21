# Cell 41 — minimal external-data bridge for the degree-three reconstruction
#
# Purpose:
#   Read the exact Cell-40 reconstruction atlas and turn it into a concrete
#   fifteen-number data request:
#       9 coefficients from the torsion-summed sector (one,one),
#       6 coefficients from one preferred nontrivial sector.
#
#   The cell exports:
#     * the exact selected q exponents;
#     * the exact 15x15 system and inverse matrices;
#     * a fillable Sage input template;
#     * a TSV worksheet for literature/data matching;
#     * a compact reconstruction bridge object.
#
#   If a completed input object already exists, the cell reconstructs the
#   fifteen observable coefficients exactly and verifies the supplied data.
#
# This cell does not assume a BKOS/HSS variable identification.  It only
# packages the exact linear algebra already certified in Cell 40.

from sage.all import *
import os

ROOT = os.getcwd()
RESULTS_DIR = os.path.join(ROOT, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

ATLAS_PATH = os.path.join(RESULTS_DIR, "degree3_exact_reconstruction_atlas_cell40.sobj")
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_minimal_external_data_bridge_cell41.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_minimal_external_data_bridge_cell41.txt")
OUT_TSV = os.path.join(RESULTS_DIR, "degree3_minimal_external_data_request_cell41.tsv")
OUT_TEMPLATE = os.path.join(RESULTS_DIR, "degree3_external_data_input_cell41_template.sage")
INPUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_external_data_input_cell41.sobj")
OUT_RECON = os.path.join(RESULTS_DIR, "degree3_reconstructed_observable_coefficients_cell41.sobj")
OUT_RECON_TXT = os.path.join(RESULTS_DIR, "degree3_reconstructed_observable_coefficients_cell41.txt")

if not os.path.exists(ATLAS_PATH) and not os.path.exists(ATLAS_PATH + ".sobj"):
    raise IOError("Missing required Cell-40 atlas: {}".format(ATLAS_PATH))

A = load(ATLAS_PATH)
K = A["coefficient_field"]
preferred = A["preferred_extra_sector"]
if preferred is None:
    raise ArithmeticError("Cell 40 did not identify a viable preferred sector")

record = A["sector_atlas"][preferred]
if not record["inverse_verified"]:
    raise ArithmeticError("Cell-40 inverse for {} was not verified".format(preferred))

base_rows = list(record["base_selected_rows"])
extra_rows = list(record["extra_selected_rows"])
selected_rows = base_rows + extra_rows
M = record["system_matrix"]
R = record["reconstruction_inverse"]
pair_labels15 = list(A["pair_labels15"])

if len(base_rows) != 9:
    raise ArithmeticError("Expected nine torsion-summed rows; got {}".format(len(base_rows)))
if len(extra_rows) != 6:
    raise ArithmeticError("Expected six extra-sector rows; got {}".format(len(extra_rows)))
if M.nrows() != 15 or M.ncols() != 15:
    raise ArithmeticError("Expected a 15x15 system matrix")
if M * R != identity_matrix(K, 15) or R * M != identity_matrix(K, 15):
    raise ArithmeticError("Exact Cell-40 inverse verification failed")


def row_key(row):
    sector, exponent = row
    return "{}|q^{}".format(sector, ZZ(exponent))


def coerce_value(value):
    if value is None:
        return None
    if isinstance(value, str):
        # Exact strings such as "7/3", "-2", or expressions in the
        # cyclotomic generator accepted by the coefficient field.
        try:
            return K(value)
        except Exception:
            return K(SR(value))
    return K(value)

print("=" * 79)
print("CELL 41: DEGREE-THREE MINIMAL EXTERNAL-DATA BRIDGE")
print("=" * 79)
print("  observable dimension       : 15")
print("  torsion-summed contribution : 9 independent coefficients")
print("  preferred extra sector     : {}".format(preferred))
print("  extra-sector contribution  : 6 independent coefficients")
print("  earliest full q cutoff     : q^{}".format(record["earliest_full_q_cutoff"]))
print("  exact inverse verified     : {}".format(record["inverse_verified"]))

print("\nPART I. MINIMAL COEFFICIENT REQUEST")
print("-" * 79)
print("  torsion-summed rows:")
for idx, row in enumerate(base_rows, 1):
    print("    B{:02d}: {}".format(idx, row_key(row)))
print("  preferred extra-sector rows:")
for idx, row in enumerate(extra_rows, 1):
    print("    E{:02d}: {}".format(idx, row_key(row)))

# Export a tab-separated worksheet.
with open(OUT_TSV, "w") as handle:
    handle.write("slot\tsource\tsector\tq_exponent\texact_value\tnotes\n")
    for idx, (sector, exponent) in enumerate(base_rows, 1):
        handle.write("B{:02d}\ttorsion_sum\t{}\t{}\t\t\n".format(
            idx, sector, ZZ(exponent)
        ))
    for idx, (sector, exponent) in enumerate(extra_rows, 1):
        handle.write("E{:02d}\textra_sector\t{}\t{}\t\t\n".format(
            idx, sector, ZZ(exponent)
        ))

# Export a fillable Sage template.  It is intentionally not executed until
# every None is replaced by an exact rational/cyclotomic expression.
with open(OUT_TEMPLATE, "w") as handle:
    handle.write("# Fill every None with an exact Sage value, then run this file with sage.\n")
    handle.write("from sage.all import *\n")
    handle.write("values = {\n")
    for row in selected_rows:
        handle.write("    {!r}: None,\n".format(row_key(row)))
    handle.write("}\n")
    handle.write("if any(v is None for v in values.values()):\n")
    handle.write("    raise ValueError('Replace every None with an exact coefficient before saving')\n")
    handle.write("save({'preferred_sector': %r, 'values': values}, %r)\n" % (
        preferred, INPUT_SOBJ
    ))
    handle.write("print('Saved completed Cell-41 input to %s')\n" % INPUT_SOBJ)

bridge = {
    "schema_version": ZZ(1),
    "scope": "minimal exact external-data bridge for the audited 15D degree-three quotient",
    "coefficient_field": K,
    "preferred_extra_sector": preferred,
    "earliest_full_q_cutoff": record["earliest_full_q_cutoff"],
    "basis_labels15": pair_labels15,
    "torsion_summed_rows": base_rows,
    "extra_sector_rows": extra_rows,
    "selected_rows": selected_rows,
    "selected_row_keys": [row_key(row) for row in selected_rows],
    "system_matrix": M,
    "reconstruction_inverse": R,
    "inverse_verified": True,
    "viable_extra_sectors": list(A["viable_single_extra_sectors"]),
    "torsion_class_full_rank_pairs": list(A["torsion_class_full_rank_pairs"]),
    "input_path": INPUT_SOBJ,
}
save(bridge, OUT_SOBJ)

summary = [
    "CELL 41: DEGREE-THREE MINIMAL EXTERNAL-DATA BRIDGE",
    "observable dimension: 15",
    "torsion-summed coefficients required: 9",
    "preferred extra sector: {}".format(preferred),
    "extra-sector coefficients required: 6",
    "earliest full q cutoff: q^{}".format(record["earliest_full_q_cutoff"]),
    "exact inverse verified: True",
    "data worksheet: {}".format(OUT_TSV),
    "fillable Sage template: {}".format(OUT_TEMPLATE),
]
with open(OUT_TXT, "w") as handle:
    handle.write("\n".join(summary) + "\n")

print("\nPART II. OUTPUT BRIDGE")
print("-" * 79)
print("  worksheet              : {}".format(OUT_TSV))
print("  fillable Sage template : {}".format(OUT_TEMPLATE))
print("  exact bridge object    : {}".format(OUT_SOBJ))

# Optional exact reconstruction when completed data are present.
if os.path.exists(INPUT_SOBJ) or os.path.exists(INPUT_SOBJ + ".sobj"):
    D = load(INPUT_SOBJ)
    values = D.get("values", {})
    missing = [key for key in bridge["selected_row_keys"] if key not in values or values[key] is None]
    if missing:
        raise ValueError("Completed Cell-41 input is missing values: {}".format(missing))

    y = vector(K, [coerce_value(values[key]) for key in bridge["selected_row_keys"]])
    c = R * y
    residual = M * c - y
    if residual != vector(K, 15, 0):
        raise ArithmeticError("Exact reconstruction residual is nonzero")

    recon = {
        "schema_version": ZZ(1),
        "scope": "exact reconstructed observable coefficient vector",
        "preferred_extra_sector": preferred,
        "basis_labels15": pair_labels15,
        "selected_rows": selected_rows,
        "input_vector": y,
        "observable_coefficients": c,
        "exact_residual_zero": True,
    }
    save(recon, OUT_RECON)
    with open(OUT_RECON_TXT, "w") as handle:
        handle.write("CELL 41 EXACT RECONSTRUCTION\n")
        handle.write("preferred extra sector: {}\n".format(preferred))
        handle.write("exact residual zero: True\n")
        for label, value in zip(pair_labels15, c):
            handle.write("{} = {}\n".format(label, value))

    print("\nPART III. EXACT RECONSTRUCTION")
    print("-" * 79)
    print("  completed input found  : {}".format(INPUT_SOBJ))
    print("  exact residual zero    : True")
    print("  reconstructed vector   : {}".format(OUT_RECON))
    print("  human-readable output  : {}".format(OUT_RECON_TXT))
else:
    print("\nPART III. EXACT RECONSTRUCTION")
    print("-" * 79)
    print("  completed input not present yet")
    print("  next action: fill and run {}".format(OUT_TEMPLATE))

print("\nCELL 41 COMPLETE")
