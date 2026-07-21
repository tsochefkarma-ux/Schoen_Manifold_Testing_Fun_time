# -*- coding: utf-8 -*-
# =============================================================================
# CELL 37A: DEGREE-THREE TORSION-SERIES CM-JET EXTRACTION
# =============================================================================
#
# Purpose
# -------
# Convert the audited Cell-32 canonical three-torsion q-series into numerical
# CM jets at t=i (q0=exp(-2*pi)) through theta^6, using only coefficients in
# the certified q-range.  The output supplies the one-sided and symmetric
# jet tensors required by the degree-three HSS curvature-response worker.
#
# This is deliberately a numerical/high-precision bridge, not yet the final
# exact formal-CM certificate.  Exact inputs retained here are:
#   * the coefficientwise torsion relation from Cell 35b;
#   * exact cyclotomic q-series coefficients through the certified cutoff.
# =============================================================================

from sage.all import *
import os
import sys

print("="*79, flush=True)
print("CELL 37A: DEGREE-THREE TORSION-SERIES CM-JET EXTRACTION", flush=True)
print("="*79, flush=True)

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

PREC_BITS = ZZ(sys.argv[1]) if len(sys.argv) > 1 else ZZ(256)
THETA_MAX = ZZ(sys.argv[2]) if len(sys.argv) > 2 else ZZ(6)
if PREC_BITS < 128:
    raise ValueError("Use at least 128 bits of precision")
if THETA_MAX < 6:
    raise ValueError("The degree-three curvature worker requires theta jets through 6")

CELL32_PATH = os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.sobj")
AUDIT_PATH = os.path.join(RESULTS_DIR, "degree3_torsion_modularity_audit_35b.sobj")
CELL36_PATH = os.path.join(
    RESULTS_DIR,
    "degree2_full_space_baseline_and_degree3_selectivity_ledger.sobj",
)
OUT_SOBJ = os.path.join(RESULTS_DIR, "degree3_cm_jets_cell37a.sobj")
OUT_TXT = os.path.join(RESULTS_DIR, "degree3_cm_jets_cell37a.txt")

for path in [CELL32_PATH, AUDIT_PATH, CELL36_PATH]:
    if not os.path.exists(path):
        raise IOError("Missing required result file: {}".format(path))

cell32 = load(CELL32_PATH)
audit = load(AUDIT_PATH)
cell36 = load(CELL36_PATH)

if not bool(audit.get("all_orders_pass", False)):
    raise ArithmeticError("Cell 35b audit did not pass")
if bool(audit.get("global_jacobi_identity", True)):
    raise ArithmeticError("The torsion relation was incorrectly marked global")

CC = ComplexField(PREC_BITS)
RR = RealField(PREC_BITS)
I = CC.gen()
piR = RR.pi()
L = RR(2)*piR
q0 = exp(-L)

CERTIFIED_Q_MAX = ZZ(cell32["sector_q_max"])
working_min, working_max = [ZZ(v) for v in cell32["one_sided_range"]]
if CERTIFIED_Q_MAX < 32:
    raise ArithmeticError("Certified cutoff is too shallow for the audited relation")

characters = list(cell32["basis"].keys())
expected_characters = ["one", "omega", "omega2"]
if sorted(characters) != sorted(expected_characters):
    raise ArithmeticError("Unexpected character labels: {}".format(characters))
characters = expected_characters

basis_labels6 = list(cell32["basis_labels"])
relation_labels = list(audit["relation_labels"])
relation_vector = vector(QQ, audit["relation_vector"])
if basis_labels6 != relation_labels:
    raise ArithmeticError("Cell-32 and Cell-35b basis orders differ")
if relation_vector != vector(QQ, [3, -3, 0, -6, -10, 16]):
    raise ArithmeticError("Unexpected audited torsion relation")

basis_labels5 = list(cell36["degree3"]["one_sided_observable_basis"])
if len(basis_labels5) != 5:
    raise ArithmeticError("Cell 36 does not contain five observable one-sided labels")
if not all(label in basis_labels6 for label in basis_labels5):
    raise ArithmeticError("Reduced basis labels are not contained in Cell-32 basis")

K12 = cell32["coefficient_field"]
zeta12_CC = CC.zeta(12)


def embed_cyclotomic(value):
    """Use the standard zeta_12 -> exp(2*pi*i/12) embedding."""
    value = K12(value)
    coeffs = value.list()
    return sum((CC(coeffs[j]) * zeta12_CC**j for j in range(len(coeffs))), CC(0))


def theta_sum(series, order, q_max):
    """Certified truncated theta^order value at q0."""
    order = ZZ(order)
    q_max = ZZ(q_max)
    total = CC(0)
    for exponent, coefficient in series.items():
        exponent = ZZ(exponent)
        if exponent > q_max:
            continue
        if exponent < working_min:
            raise ArithmeticError("Series exponent below declared working minimum")
        factor = ZZ(1) if order == 0 else exponent**order
        total += CC(factor) * embed_cyclotomic(coefficient) * CC(q0)**exponent
    return total


def exact_relation_on_range(character, q_min, q_max):
    basis = cell32["basis"][character]
    for exponent in range(ZZ(q_min), ZZ(q_max) + 1):
        residual = K12(0)
        for coefficient, label in zip(relation_vector, basis_labels6):
            residual += K12(coefficient) * basis[label].get(ZZ(exponent), K12(0))
        if residual != 0:
            return False, ZZ(exponent), residual
    return True, None, K12(0)


print("\nPART I. INPUT AND RANGE AUDIT", flush=True)
print("-"*79, flush=True)
print("  precision bits            : {}".format(PREC_BITS), flush=True)
print("  CM point                  : t=i", flush=True)
print("  q0                        : exp(-2*pi)", flush=True)
print("  q0 numerical              : {}".format(q0), flush=True)
print("  certified q maximum       : {}".format(CERTIFIED_Q_MAX), flush=True)
print("  stored working range      : [{} , {}]".format(working_min, working_max), flush=True)
print("  unsafe upper buffer ignored: [{} , {}]".format(CERTIFIED_Q_MAX + 1, working_max), flush=True)
print("  theta derivatives         : 0 through {}".format(THETA_MAX), flush=True)

print("\nPART II. EXACT COEFFICIENTWISE RELATION CHECK", flush=True)
print("-"*79, flush=True)
relation_checks = {}
for character in characters:
    ok, bad_exp, bad_value = exact_relation_on_range(character, working_min, CERTIFIED_Q_MAX)
    relation_checks[character] = bool(ok)
    print("  {:7s}: {}".format(character, ok), flush=True)
    if not ok:
        raise ArithmeticError(
            "Relation failed for {} at q^{}: {}".format(character, bad_exp, bad_value)
        )

print("\nPART III. ONE-SIDED CM JETS", flush=True)
print("-"*79, flush=True)

jets = {}
relation_jet_residuals = {}
convergence = {}
cutoffs = sorted(set([max(working_min, CERTIFIED_Q_MAX - 16),
                      max(working_min, CERTIFIED_Q_MAX - 8),
                      CERTIFIED_Q_MAX]))

for character in characters:
    jets[character] = {}
    basis = cell32["basis"][character]
    for label in basis_labels5:
        jets[character][label] = [
            theta_sum(basis[label], order, CERTIFIED_Q_MAX)
            for order in range(THETA_MAX + 1)
        ]

    relation_jet_residuals[character] = []
    for order in range(THETA_MAX + 1):
        residual = CC(0)
        for coefficient, label in zip(relation_vector, basis_labels6):
            residual += CC(coefficient) * theta_sum(
                basis[label], order, CERTIFIED_Q_MAX
            )
        relation_jet_residuals[character].append(residual)

    convergence[character] = {}
    for label in basis_labels5:
        convergence[character][label] = {}
        for order in range(THETA_MAX + 1):
            values = [theta_sum(basis[label], order, cutoff) for cutoff in cutoffs]
            differences = [abs(values[j+1] - values[j]) for j in range(len(values)-1)]
            convergence[character][label][order] = {
                "cutoffs": cutoffs,
                "values": values,
                "successive_differences": differences,
            }

    max_residual = max(abs(value) for value in relation_jet_residuals[character])
    print("  {:7s}: max relation jet residual = {}".format(character, max_residual), flush=True)

# A practical numerical stability statistic.
max_last_step = RR(0)
for character in characters:
    for label in basis_labels5:
        for order in range(THETA_MAX + 1):
            diffs = convergence[character][label][order]["successive_differences"]
            if diffs:
                max_last_step = max(max_last_step, RR(diffs[-1]))
print("  maximum q^{} -> q^{} jet change: {}".format(cutoffs[-2], cutoffs[-1], max_last_step), flush=True)

print("\nPART IV. LEADING HSS CHANNEL CM DATA", flush=True)
print("-"*79, flush=True)

# Compute E2,E4,E6 and B=9 prod(1-q^n)^(-4) at q0.  The q0 tail is tiny;
# use a generous numerical cutoff independent of the Cell-32 series buffer.
NMOD = ZZ(max(120, CERTIFIED_Q_MAX + 40))


def sigma_power(n, power):
    return sum(ZZ(d)**ZZ(power) for d in divisors(ZZ(n)))

E2_i = RR(1)
E4_i = RR(1)
E6_i = RR(1)
for n in range(1, NMOD + 1):
    qn = q0**n
    E2_i += RR(-24*sigma_power(n, 1))*qn
    E4_i += RR(240*sigma_power(n, 3))*qn
    E6_i += RR(-504*sigma_power(n, 5))*qn

B0 = RR(9)
for n in range(1, NMOD + 1):
    B0 *= (RR(1) - q0**n)**(-4)

# Ramanujan recurrence for theta=q d/dq.
Rram = PolynomialRing(QQ, names=("e2", "e4", "e6", "bb"))
e2r, e4r, e6r, bbr = Rram.gens()
ram_derivatives = {
    e2r: (e2r**2 - e4r)/12,
    e4r: (e2r*e4r - e6r)/3,
    e6r: (e2r*e6r - e4r**2)/2,
    bbr: ((1-e2r)/6)*bbr,
}


def theta_ram(poly):
    poly = Rram(poly)
    return Rram(sum(
        poly.derivative(generator)*ram_derivatives[generator]
        for generator in Rram.gens()
    ))


def eval_ram(poly):
    return RR(Rram(poly)(E2_i, E4_i, E6_i, B0))

B_jets = []
current = Rram(bbr)
for order in range(THETA_MAX + 2):
    B_jets.append(eval_ram(current))
    current = theta_ram(current)

Ycm = L**2 * E4_i
alpha0 = (Ycm + 12)/72
beta0 = -(7*Ycm**2 - 552*Ycm - 13392)/(288*(Ycm + 36))

print("  E2(i) versus 6/L residual: {}".format(E2_i - 6/L), flush=True)
print("  E6(i) residual           : {}".format(E6_i), flush=True)
print("  E4(i)                    : {}".format(E4_i), flush=True)
print("  B(q0)                    : {}".format(B0), flush=True)
print("  alpha0                   : {}".format(alpha0), flush=True)
print("  beta0                    : {}".format(beta0), flush=True)

print("\nPART V. FIFTEEN SYMMETRIC CM-JET TENSORS", flush=True)
print("-"*79, flush=True)

pairs = []
labels15 = []
for i in range(5):
    for j in range(i, 5):
        pairs.append((ZZ(i), ZZ(j)))
        if i == j:
            labels15.append("g{}(y)g{}(z)".format(i+1, j+1))
        else:
            labels15.append(
                "g{}(y)g{}(z)+g{}(y)g{}(z)".format(i+1, j+1, j+1, i+1)
            )

sector_pairs = [
    ("omega", "omega"),
    ("omega", "omega2"),
]
sector_tensors = {}
for left_character, right_character in sector_pairs:
    sector_name = "({},{})".format(left_character, right_character)
    sector_tensors[sector_name] = {}
    for pair_index, (i, j) in enumerate(pairs):
        left_i = jets[left_character][basis_labels5[i]]
        right_j = jets[right_character][basis_labels5[j]]
        if i == j:
            tensor = {
                (ZZ(a), ZZ(b)): left_i[a]*right_j[b]
                for a in range(THETA_MAX + 1)
                for b in range(THETA_MAX + 1)
            }
        else:
            left_j = jets[left_character][basis_labels5[j]]
            right_i = jets[right_character][basis_labels5[i]]
            tensor = {
                (ZZ(a), ZZ(b)): left_i[a]*right_j[b] + left_j[a]*right_i[b]
                for a in range(THETA_MAX + 1)
                for b in range(THETA_MAX + 1)
            }
        sector_tensors[sector_name][labels15[pair_index]] = tensor
    print("  {}: {} tensors, each {} x {}".format(
        sector_name, len(labels15), THETA_MAX + 1, THETA_MAX + 1
    ), flush=True)

certificate = {
    "schema_version": ZZ(1),
    "scope": "high-precision CM-jet bridge for the audited canonical three-torsion candidate",
    "precision_bits": PREC_BITS,
    "theta_max": THETA_MAX,
    "CM_point": "t=i",
    "L": L,
    "q0": q0,
    "certified_q_max": CERTIFIED_Q_MAX,
    "working_range": (working_min, working_max),
    "unsafe_buffer_ignored": (CERTIFIED_Q_MAX + 1, working_max),
    "characters": characters,
    "basis_labels6": basis_labels6,
    "basis_labels5": basis_labels5,
    "relation_vector": relation_vector,
    "relation_coefficient_checks": relation_checks,
    "relation_jet_residuals": relation_jet_residuals,
    "cutoffs": cutoffs,
    "convergence": convergence,
    "max_last_step": max_last_step,
    "jets": jets,
    "observable_pairs": pairs,
    "observable_labels15": labels15,
    "sector_pairs": sector_pairs,
    "sector_tensors": sector_tensors,
    "leading_channel_data": {
        "E2_i": E2_i,
        "E4_i": E4_i,
        "E6_i": E6_i,
        "B0": B0,
        "B_jets": B_jets,
        "Ycm": Ycm,
        "alpha0": alpha0,
        "beta0": beta0,
    },
    "next_step": "Cell 37B degree-three linear curvature response and bare-selectivity matrix",
}
save(certificate, OUT_SOBJ)

summary_lines = [
    "CELL 37A: DEGREE-THREE TORSION-SERIES CM-JET EXTRACTION",
    "precision bits: {}".format(PREC_BITS),
    "CM point: t=i, q0=exp(-2*pi)",
    "certified q cutoff used: {}".format(CERTIFIED_Q_MAX),
    "unsafe working buffer ignored: q^{} through q^{}".format(CERTIFIED_Q_MAX + 1, working_max),
    "theta jets: 0 through {}".format(THETA_MAX),
    "all three exact coefficientwise relation checks: {}".format(all(relation_checks.values())),
    "maximum final-cutoff jet change: {}".format(max_last_step),
    "reduced one-sided basis dimension: 5",
    "symmetric observable tensors per parked sector: 15",
    "leading B-channel jets prepared through theta^{}".format(THETA_MAX + 1),
    "next: Cell 37B degree-three linear curvature response",
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
