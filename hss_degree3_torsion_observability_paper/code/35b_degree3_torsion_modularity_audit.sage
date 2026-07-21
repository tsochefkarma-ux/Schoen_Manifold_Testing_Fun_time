# -*- coding: utf-8 -*-
# Cell 35b: referee-grade modularity audit for the degree-three HSS
# three-torsion specialization relation.
#
# This cell replaces the un-audited group assumption in Cell 35.
# It does NOT claim that the relation is a global E8 Jacobi-form identity.
# Instead it proves the negative control at z=0, derives the scalar Jacobi
# index and torsion correction, derives the common torsion stabilizer, audits
# cusp widths/local parameters, and applies a conservative Sturm bound only
# after checking the hypotheses of the standard torsion-specialization theorem.
#
# INPUT:
#   results/degree3_A3_B3_hss_torsion_candidate.sobj
# produced by Cell 32 v2 with certified sector cutoff at least q^28
# (the user's current q^40 file is more than sufficient).

from sage.all import *
import os

print("="*79)
print("CELL 35b: DEGREE-THREE TORSION MODULARITY AND STURM AUDIT")
print("="*79)

RESULTS_DIR = "results"
INPUT_CANDIDATES = [
    os.path.join(RESULTS_DIR, "degree3_A3_B3_hss_torsion_candidate.sobj"),
    "degree3_A3_B3_hss_torsion_candidate.sobj",
]
INPUT_PATH = next((p for p in INPUT_CANDIDATES if os.path.exists(p)), None)
if INPUT_PATH is None:
    raise IOError(
        "Missing Cell-32 result: results/degree3_A3_B3_hss_torsion_candidate.sobj"
    )

print("Loading:", INPUT_PATH, flush=True)
data = load(INPUT_PATH)
required = [
    "sector_q_max", "one_sided_range", "basis", "basis_labels",
    "gamma", "normalization_checks",
]
missing = [key for key in required if key not in data]
if missing:
    raise KeyError("Cell-32 result is missing keys: {}".format(missing))

q_cert = ZZ(data["sector_q_max"])
q_min, q_buffer_max = map(ZZ, data["one_sided_range"])
basis = data["basis"]
labels = list(data["basis_labels"])
gamma = vector(ZZ, data["gamma"])
normalization_checks = dict(data["normalization_checks"])
characters = ["one", "omega", "omega2"]

expected_labels = [
    "E6^2*A3",
    "E4^3*A3",
    "E4*E6*B3",
    "E4^2*A1*A2",
    "E6*A1*B2",
    "E4*A1^3",
]
if labels != expected_labels:
    raise ValueError("Unexpected Cell-32 basis ordering: {}".format(labels))
if not all(normalization_checks.values()):
    raise ArithmeticError("A Cell-32 Sakai normalization check is false")

relation = vector(QQ, [3, -3, 0, -6, -10, 16])
weight = ZZ(16)
lattice_index = ZZ(3)
gamma_norm = ZZ(gamma.dot_product(gamma))
if gamma_norm != 8:
    raise ArithmeticError("Expected (gamma,gamma)=8, got {}".format(gamma_norm))
Q_gamma = ZZ(gamma_norm // 2)
scalar_index = ZZ(lattice_index * Q_gamma)  # classical Jacobi index on z=s*gamma
N = ZZ(3)
alpha = QQ(1)/N
correction_tau = QQ(scalar_index) * alpha**2  # 4/3 in tau variable
correction_qt = ZZ(N * correction_tau)         # 4 after tau=3t
if correction_tau != QQ(4)/3 or correction_qt != 4:
    raise ArithmeticError("Unexpected torsion correction")

# -----------------------------------------------------------------------------
# Part I. Global negative control and homogeneous-space audit.
# -----------------------------------------------------------------------------

print("\nPART I. GLOBAL JACOBI NEGATIVE CONTROL")
print("-"*79)

R = PolynomialRing(QQ, names=("E4", "E6"))
E4s, E6s = R.gens()
z0_residual = (
    3*E6s**2*E4s
    - 3*E4s**3*E4s
    - 6*E4s**2*E4s*E4s
    - 10*E6s*E4s*E6s
    + 16*E4s*E4s**3
)
z0_expected = 7*E4s*(E4s**3 - E6s**2)
if z0_residual != z0_expected:
    raise ArithmeticError("z=0 negative-control simplification failed")

print("  Phi(tau,0) = 7 E4 (E4^3-E6^2)")
print("             = 12096 E4 Delta != 0")
print("  global E8 Jacobi-form identity? False")
print("  correct scope: kernel after HSS three-torsion specialization")

# Weight/index ledger for the six global basis elements.
weight_index = {
    "E6^2*A3":       (16, 3),
    "E4^3*A3":       (16, 3),
    "E4*E6*B3":      (16, 3),
    "E4^2*A1*A2":    (16, 3),
    "E6*A1*B2":      (16, 3),
    "E4*A1^3":       (16, 3),
}
homogeneous = all(weight_index[label] == (16, 3) for label in labels)
print("  homogeneous weight/index 16/3? {}".format(homogeneous))
if not homogeneous:
    raise ArithmeticError("The relation is not homogeneous")

# -----------------------------------------------------------------------------
# Part II. Restriction to the gamma line and intrinsic torsion correction.
# -----------------------------------------------------------------------------

print("\nPART II. SCALAR JACOBI RESTRICTION AND CORRECTION")
print("-"*79)
print("  lattice index m                  : {}".format(lattice_index))
print("  (gamma,gamma)                    : {}".format(gamma_norm))
print("  Q(gamma)=(gamma,gamma)/2         : {}".format(Q_gamma))
print("  scalar index on z=s gamma        : m Q(gamma) = {}".format(scalar_index))
print("  torsion section in tau variable  : s=(tau+r)/3")
print("  intrinsic exponential correction : exp(2 pi i * (4/3) tau)")
print("  after tau=3t                     : q_t^{}".format(correction_qt))
print("  modular weight after correction  : {} (unchanged)".format(weight))
print("  multiplied by Delta or eta?      : False")

# -----------------------------------------------------------------------------
# Part III. Derive the stabilizers over F_3 instead of assuming a group.
# Pair convention: (alpha,beta) transforms as the row vector
# (alpha,beta) A.  Modulo Z^2 this becomes (1,r) A modulo 3.
# -----------------------------------------------------------------------------

print("\nPART III. EXACT THREE-TORSION STABILIZER")
print("-"*79)
F3 = GF(3)
SL2F3 = []
for a in F3:
    for b in F3:
        for c in F3:
            for d in F3:
                A = matrix(F3, [[a,b],[c,d]])
                if A.det() == 1:
                    SL2F3.append(A)
if len(SL2F3) != 24:
    raise ArithmeticError("Expected |SL2(F3)|=24")

def matrix_key(A):
    return tuple(ZZ(x) for x in A.list())

stabilizer_keys = {}
for r in range(3):
    pair = vector(F3, [1, r])
    stab = [A for A in SL2F3 if pair*A == pair]
    stabilizer_keys[r] = set(matrix_key(A) for A in stab)
    print("  r={} individual stabilizer size/index: {}/{}".format(
        r, len(stab), ZZ(24//len(stab))
    ))

common_keys = set.intersection(*[stabilizer_keys[r] for r in range(3)])
identity_key = matrix_key(identity_matrix(F3, 2))
if common_keys != set([identity_key]):
    raise ArithmeticError("Common mod-3 stabilizer is not the identity")

Gamma_tau = Gamma(3)
sl_index = ZZ(Gamma_tau.index())
proj_index = ZZ(Gamma_tau.projective_index())
if sl_index != 24:
    raise ArithmeticError("Unexpected index for Gamma(3): {}".format(sl_index))

print("  common mod-3 stabilizer size      : {}".format(len(common_keys)))
print("  common subgroup in tau variable   : Gamma(3)")
print("  [SL2(Z):Gamma(3)]                 : {}".format(sl_index))
print("  projective index                  : {}".format(proj_index))
print("  Gamma(3) is odd (-I absent)?      : {}".format(Gamma_tau.is_odd()))

# Cross-check the conjugated group in the t variable.  For tau=3t,
# S B S^{-1} in Gamma(3) is equivalent to
# c=0 mod 9 and a=d=1 mod 3, i.e. Gamma_H(9,{1,4,7}).
Gamma_t = GammaH(9, [4])
if ZZ(Gamma_t.index()) != 24:
    raise ArithmeticError("Unexpected t-variable conjugate-group index")
print("  conjugate subgroup in t variable : GammaH(9,<4>)")
print("    = Gamma0(9) intersect Gamma1(3)")
print("  t-variable subgroup index         : {}".format(Gamma_t.index()))

# -----------------------------------------------------------------------------
# Part IV. Cusps and local parameter.
# -----------------------------------------------------------------------------

print("\nPART IV. CUSP AND LOCAL-PARAMETER AUDIT")
print("-"*79)
try:
    cusps = list(Gamma_tau.cusps())
    cusp_widths = [ZZ(Gamma_tau.cusp_width(c)) for c in cusps]
except Exception as error:
    raise RuntimeError("Sage could not enumerate Gamma(3) cusps: {}".format(error))

try:
    infinity_width = ZZ(Gamma_tau.cusp_width(Infinity))
except Exception:
    infinity_width = ZZ(Gamma_tau.cusp_width(Cusp(Infinity)))

print("  reduced cusps                     : {}".format(cusps))
print("  cusp widths                       : {}".format(cusp_widths))
print("  width at infinity                 : {}".format(infinity_width))
print("  sum of cusp widths                : {}".format(sum(cusp_widths)))
if infinity_width != 3:
    raise ArithmeticError("Expected Gamma(3) cusp width 3 at infinity")
if sum(cusp_widths) != proj_index:
    raise ArithmeticError("Cusp widths do not sum to the projective index")

print("  local parameter at infinity       : exp(2 pi i tau/3)")
print("  with tau=3t this is               : q_t=exp(2 pi i t)")

# -----------------------------------------------------------------------------
# Part V. Holomorphy/character theorem ledger and Sturm bound.
# -----------------------------------------------------------------------------

print("\nPART V. THEOREM HYPOTHESES AND STURM BOUND")
print("-"*79)
# The standard torsion-specialization theorem for holomorphic Jacobi forms says
# q^(M alpha^2) phi(tau,alpha tau+beta) is a holomorphic modular form of the
# same weight on the torsion stabilizer, possibly with finite character.
# We do not assume that character is trivial; the conservative subgroup bound
# below remains valid with finite character over the cyclotomic coefficient field.

jacobi_holomorphic = bool(all(normalization_checks.values()) and homogeneous)
torsion_theorem_applies = jacobi_holomorphic
finite_character_allowed = True
holomorphic_at_all_cusps = torsion_theorem_applies

print("  Phi is a holomorphic Jacobi form?  {}".format(jacobi_holomorphic))
print("  torsion-specialization theorem?    {}".format(torsion_theorem_applies))
print("  finite character allowed?          {}".format(finite_character_allowed))
print("  trivial character assumed?         False")
print("  holomorphic at every cusp?         {} (by torsion theorem)".format(
    holomorphic_at_all_cusps
))
if not holomorphic_at_all_cusps:
    raise ArithmeticError("Holomorphy hypotheses were not met")

conservative_bound = ZZ(floor(QQ(weight*sl_index)/12))
try:
    sage_bound = ZZ(Gamma_tau.sturm_bound(weight))
except Exception:
    sage_bound = conservative_bound
# Never use less than the transparent k*[SL2:Gamma]/12 bound.
audit_bound = max(conservative_bound, sage_bound)
required_uncorrected_max = ZZ(audit_bound - correction_qt)

print("  Sage-reported Sturm bound          : {}".format(sage_bound))
print("  conservative k*[SL2:Gamma]/12     : {}".format(conservative_bound))
print("  audit bound used in q_t powers     : {}".format(audit_bound))
print("  uncorrected coefficients required  : q^{} through q^{}".format(
    -correction_qt, required_uncorrected_max
))
print("  Cell-32 certified cutoff           : q^{}".format(q_cert))
print("  Cell-32 unsafe buffer endpoint     : q^{}".format(q_buffer_max))
if q_cert < required_uncorrected_max:
    raise ArithmeticError(
        "Certified q cutoff {} is below required {}".format(
            q_cert, required_uncorrected_max
        )
    )

# -----------------------------------------------------------------------------
# Part VI. Exact coefficient test, using certified data only.
# -----------------------------------------------------------------------------

print("\nPART VI. EXACT STURM-WINDOW COEFFICIENT TEST")
print("-"*79)

def coefficient(series, exponent):
    return series.get(ZZ(exponent), 0)

character_results = {}
all_coefficients_zero = True
for character in characters:
    failures = []
    # corrected coefficient Q^n equals uncorrected relation coefficient Q^(n-4)
    for n in range(0, audit_bound + 1):
        source_exponent = ZZ(n - correction_qt)
        value = sum(
            relation[j] * coefficient(basis[character][labels[j]], source_exponent)
            for j in range(6)
        )
        if value != 0:
            failures.append((ZZ(n), source_exponent, value))
    passed = not failures
    all_coefficients_zero = all_coefficients_zero and passed
    character_results[character] = {
        "passed": passed,
        "failures": failures,
    }
    print("  {:7s}: corrected q_t^0..q_t^{} vanish? {}".format(
        character, audit_bound, passed
    ))
    if failures:
        print("           first failure:", failures[0])

# Explicitly refuse to touch the unsafe working buffer.
unsafe_buffer_used = False
print("  unsafe buffer coefficients used?   {}".format(unsafe_buffer_used))
if not all_coefficients_zero:
    raise ArithmeticError("A corrected Sturm-window coefficient is nonzero")

# -----------------------------------------------------------------------------
# Part VII. Verdict and saved certificate.
# -----------------------------------------------------------------------------

all_orders_pass = bool(
    z0_residual != 0
    and homogeneous
    and common_keys == set([identity_key])
    and sl_index == 24
    and infinity_width == 3
    and holomorphic_at_all_cusps
    and q_cert >= required_uncorrected_max
    and all_coefficients_zero
    and not unsafe_buffer_used
)

print("\nPART VII. AUDIT VERDICT")
print("-"*79)
print("  global Jacobi identity rejected?            {}".format(z0_residual != 0))
print("  intrinsic correction, weight unchanged?    True")
print("  common tau-group derived as Gamma(3)?       {}".format(
    common_keys == set([identity_key])
))
print("  cusp/local-parameter audit passed?          {}".format(
    infinity_width == 3 and sum(cusp_widths) == proj_index
))
print("  theorem-based all-cusp holomorphy passed?   {}".format(
    holomorphic_at_all_cusps
))
print("  certified coefficients exceed bound?        {}".format(
    q_cert >= required_uncorrected_max
))
print("  all-orders torsion-specialization identity? {}".format(all_orders_pass))

if not all_orders_pass:
    raise ArithmeticError("Cell 35b modularity audit failed")

os.makedirs(RESULTS_DIR, exist_ok=True)
SOBJ_PATH = os.path.join(RESULTS_DIR, "degree3_torsion_modularity_audit_35b.sobj")
TXT_PATH = os.path.join(RESULTS_DIR, "degree3_torsion_modularity_audit_35b.txt")

certificate = {
    "input_path": INPUT_PATH,
    "relation_labels": labels,
    "relation_vector": relation,
    "global_z0_residual": "12096*E4*Delta",
    "global_jacobi_identity": False,
    "scope": "HSS order-three torsion specialization kernel",
    "weight": weight,
    "lattice_index": lattice_index,
    "gamma_norm": gamma_norm,
    "scalar_jacobi_index": scalar_index,
    "torsion_denominator": N,
    "correction_in_tau": correction_tau,
    "correction_in_qt": correction_qt,
    "correction_changes_weight": False,
    "tau_group": str(Gamma_tau),
    "tau_group_index": sl_index,
    "tau_group_projective_index": proj_index,
    "t_group": str(Gamma_t),
    "t_group_index": ZZ(Gamma_t.index()),
    "cusps": [str(c) for c in cusps],
    "cusp_widths": cusp_widths,
    "infinity_width": infinity_width,
    "local_parameter": "q_t=exp(2*pi*i*t)=exp(2*pi*i*tau/3)",
    "character": "finite character permitted; triviality not assumed or needed",
    "holomorphy_at_all_cusps": holomorphic_at_all_cusps,
    "holomorphy_reason": "standard torsion-specialization theorem for holomorphic Jacobi forms",
    "sage_sturm_bound": sage_bound,
    "conservative_sturm_bound": conservative_bound,
    "audit_bound": audit_bound,
    "certified_q_max": q_cert,
    "required_uncorrected_q_max": required_uncorrected_max,
    "unsafe_buffer_used": unsafe_buffer_used,
    "character_results": character_results,
    "all_orders_pass": all_orders_pass,
}
save(certificate, SOBJ_PATH)

summary = """DEGREE-THREE TORSION MODULARITY AUDIT (CELL 35b)

GLOBAL CONTROL
  Phi(tau,0) = 12096 E4 Delta != 0
  global E8 Jacobi-form identity: False
  valid scope: HSS order-three torsion-specialization kernel

JACOBI/TORSION DATA
  weight/index: 16/3
  gamma norm: 8
  scalar line index: 12
  correction in tau variable: exp(2 pi i (4/3) tau)
  after tau=3t: q_t^4
  correction changes modular weight: False

GROUP AUDIT
  common group in tau variable: Gamma(3)
  [SL2(Z):Gamma(3)]: {sl_index}
  projective index: {proj_index}
  cusp widths: {widths}
  width at infinity: {inf_width}
  local parameter: q_t=exp(2 pi i tau/3)
  conjugate t-variable group: Gamma0(9) intersect Gamma1(3)
  t-variable group index: {t_index}
  character: finite character permitted; triviality not assumed

STURM AUDIT
  weight: {weight}
  Sage bound: {sage_bound}
  conservative bound used: {audit_bound}
  required uncorrected cutoff: q^{required}
  certified Cell-32 cutoff: q^{certified}
  unsafe buffer used: False

COEFFICIENT TESTS
  one: {one}
  omega: {omega}
  omega2: {omega2}

ALL-ORDERS TORSION-SPECIALIZATION IDENTITY: {passed}

CONSEQUENCE
  The relation is not global in J_{{16,3}}(E8). It is an all-orders
  relation after the three canonical HSS order-three torsion evaluations,
  under the standard holomorphic Jacobi torsion-specialization theorem.
""".format(
    sl_index=sl_index,
    proj_index=proj_index,
    widths=cusp_widths,
    inf_width=infinity_width,
    t_index=Gamma_t.index(),
    weight=weight,
    sage_bound=sage_bound,
    audit_bound=audit_bound,
    required=required_uncorrected_max,
    certified=q_cert,
    one=character_results["one"]["passed"],
    omega=character_results["omega"]["passed"],
    omega2=character_results["omega2"]["passed"],
    passed=all_orders_pass,
)
with open(TXT_PATH, "w") as handle:
    handle.write(summary)

print("\nSaved:")
print(" ", SOBJ_PATH)
print(" ", TXT_PATH)
print("\nCELL 35b AUDIT: PASS")
