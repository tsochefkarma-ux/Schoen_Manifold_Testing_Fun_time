from sage.all import *
from itertools import combinations
import os

# ============================================================================
# CHANNEL-GAUGE REGULARIZATION AND FALSE-DIVISOR CONTROLS
# ============================================================================
#
# This cell implements the two controls selected after cell 25:
#
#   A. Determine the U=3X+1 valuations of the two degree-one response
#      channels and of the fixed-geometry degree-two defect.  Use those
#      valuations to pass to pole-regularized channels and determine how much
#      of the U^{-3},U^{-2},U^{-1} connection support is normalization/gauge.
#
#   B. Repeat the unique nine-coefficient Laurent closure test with matched
#      false divisor variables
#
#          W = 2X+1,   T = 3X+2,   Q = X+1,   R = X,
#
#      using exactly the successful U-support pattern and independent held-out
#      radial points.  U is the positive control.
#
# RESTART BEHAVIOUR
# -----------------
# The exact classical-frame section runs in a fresh kernel.
# The numerical controls require the live geometry namespace from cells 13--24.
# If that namespace is absent, the cell exits cleanly after the exact section
# and tells you what must be rerun.
#
# Recommended order in a live session:
#
#     cells 13--24, then this cell.
#
# Cell 25 is not required for the numerical controls, although its checkpoint
# is loaded when available and recorded in the output checkpoint.
# ============================================================================

print("="*78, flush=True)
print("CHANNEL-GAUGE REGULARIZATION / FALSE-DIVISOR CONTROL TEST", flush=True)
print("Separating intrinsic U-divisor structure from channel normalization", flush=True)
print("="*78, flush=True)
print(flush=True)

# ----------------------------------------------------------------------------
# PART I. Exact classical symmetric-frame certificate.
# ----------------------------------------------------------------------------

print("PART I. EXACT CLASSICAL SYMMETRIC-FRAME CERTIFICATE", flush=True)
print("-"*78, flush=True)

P26 = PolynomialRing(QQ, names=("x26", "y26", "z26"))
x26, y26, z26 = P26.gens()
K26 = FractionField(P26)

V26 = (
    9*x26*y26*z26
    + QQ(3)/2*y26*y26*z26
    + QQ(3)/2*y26*z26*z26
)
vars26 = [x26, y26, z26]

g26 = matrix(K26, 3, 3, [
    K26(V26.derivative(vars26[i])*V26.derivative(vars26[j])/V26**2
        - V26.derivative(vars26[i]).derivative(vars26[j])/V26)
    for i in range(3) for j in range(3)
])
M26 = matrix(K26, 3, 3, [
    K26(V26.derivative(vars26[i]).derivative(vars26[j]))
    for i in range(3) for j in range(3)
])

subs26 = {y26: QQ(1), z26: QQ(1)}
gline26 = matrix(K26, [[g26[i,j].subs(subs26) for j in range(3)] for i in range(3)])
Mline26 = matrix(K26, [[M26[i,j].subs(subs26) for j in range(3)] for i in range(3)])

U26_symbol = 3*x26 + 1
e26 = vector(K26, [-(2*x26+1), 1, 1])
o26 = vector(K26, [0, 1, -1])
r26 = vector(K26, [x26, 1, 1])
frame26 = [("e", e26), ("o", o26), ("r", r26)]


def quad26(A, left, right):
    left_matrix26 = matrix(K26, 1, 3, list(left))
    right_matrix26 = matrix(K26, 3, 1, list(right))
    return K26((left_matrix26*A*right_matrix26)[0,0]).factor()


def exact_equal26(left, right):
    return K26(left-right) == 0

expected_g26 = {
    ("e","e"): K26(6),
    ("o","o"): K26(2),
    ("r","r"): K26(3),
    ("e","o"): K26(0),
    ("o","e"): K26(0),
    ("e","r"): K26(0),
    ("r","e"): K26(0),
    ("o","r"): K26(0),
    ("r","o"): K26(0),
}

print("Classical metric frame products:", flush=True)
for name_left, left in frame26:
    for name_right, right in frame26:
        value = quad26(gline26, left, right)
        print("  g0({},{}) = {}".format(name_left, name_right, value), flush=True)
        if not exact_equal26(value, expected_g26[(name_left,name_right)]):
            raise ArithmeticError("Exact classical frame identity failed.")

det_g26 = K26(gline26.det()).factor()
print("  det(g0) = {}".format(det_g26), flush=True)
if not exact_equal26(det_g26, K26(9)/U26_symbol**2):
    raise ArithmeticError("Exact determinant identity failed.")

print("Classical Hessian M frame norms:", flush=True)
expected_M26 = {
    "e": -18*U26_symbol,
    "o": -6*U26_symbol,
    "r": 18*U26_symbol,
}
for name, direction in frame26:
    value = quad26(Mline26, direction, direction)
    print("  M({},{}) = {}".format(name, name, value), flush=True)
    if not exact_equal26(value, expected_M26[name]):
        raise ArithmeticError("Exact M-norm identity failed for {}.".format(name))

gradVline26 = vector(K26, [V26.derivative(v).subs(subs26) for v in vars26])
print("Primitive directional derivatives:", flush=True)
for name, direction in [("e",e26),("o",o26)]:
    value = K26(sum(gradVline26[i]*direction[i] for i in range(3))).factor()
    print("  dV({}) = {}".format(name, value), flush=True)
    if not exact_equal26(value, 0):
        raise ArithmeticError("Primitive directional derivative did not vanish.")

print("Exact frame certificate passed.", flush=True)
print("  intrinsic classical degeneration divisor: U=3X+1=0", flush=True)
print("  W=2X+1 occurs in the frame vector but not in det(g0) or M norms", flush=True)
print(flush=True)

# Optional restart-safe affine checkpoint from cell 25.
CHECKPOINT25_26 = None
if os.path.exists("connection_affine_checkpoint.sobj"):
    try:
        CHECKPOINT25_26 = load("connection_affine_checkpoint.sobj")
        print("Loaded optional cell-25 affine checkpoint.", flush=True)
    except Exception as error26:
        print("Warning: could not load cell-25 checkpoint: {}".format(error26), flush=True)

# ----------------------------------------------------------------------------
# Detect whether the live degree-two geometry pipeline is available.
# ----------------------------------------------------------------------------

required_live26 = [
    "RF", "alpha0", "beta0", "component_data_at_X",
    "BASE70_ROWS_24", "BASE70_RHS_24", "ALL_POINTS_24",
    "FIT_COUNT_24", "TOTAL_COUNT_24", "CONNECTION_NAMES_24",
    "solve_connection_24", "rank_pair_23", "canonical_solution_23",
    "residuals_23", "max_abs_23", "compact_rank_profile_23",
]
missing_live26 = [name for name in required_live26 if name not in globals()]
LIVE26 = not missing_live26

if not LIVE26:
    print("="*78, flush=True)
    print("EXACT SECTION COMPLETE; LIVE NUMERICAL CONTROLS NOT RUN", flush=True)
    print("Missing live objects: {}".format(missing_live26), flush=True)
    print("Rerun cells 13--24 in one Sage session, then rerun this cell.", flush=True)
    print("The cell-25 checkpoint alone does not contain the radial response rows", flush=True)
    print("needed for divisor valuations and false-divisor controls.", flush=True)
    print("="*78, flush=True)
else:
    # ========================================================================
    # LIVE NUMERICAL SECTION
    # ========================================================================

    print("PART II. LIVE CHANNEL VALUATIONS AND DIVISOR CONTROLS", flush=True)
    print("-"*78, flush=True)
    print("Working real precision: {} bits".format(RF.precision()), flush=True)

    TOL26 = RF(10)**(-70)
    DISPLAY_ZERO26 = RF(10)**(-60)

    GEOMETRIES26 = [
        ("primitive_product", "primitive product F_A", (QQ(1), QQ(0), QQ(0))),
        ("primitive_plus_double_cover", "F_A + F_B/8", (QQ(1), QQ(1)/8, QQ(0))),
        ("Z2_product", "Z_2(y)Z_2(z)", (QQ(1), QQ(1)/64, QQ(1)/8)),
    ]

    # New positive points, disjoint from the 70 points through cell 24.
    VALIDATION_X26 = [
        QQ(1)/3,
        QQ(7)/20,
        QQ(11)/25,
        QQ(14)/25,
        QQ(17)/25,
        QQ(19)/25,
        QQ(6)/7,
        QQ(24)/25,
        QQ(26)/25,
        QQ(23)/18,
        QQ(20)/13,
        QQ(27)/16,
        QQ(29)/12,
        QQ(23)/6,
        QQ(21)/4,
        QQ(10),
    ]

    TRAIN_POINTS26 = list(ALL_POINTS_24)
    ALL_POINTS26 = list(TRAIN_POINTS26)
    for value26 in VALIDATION_X26:
        if value26 not in ALL_POINTS26:
            ALL_POINTS26.append(value26)

    TRAIN_COUNT26 = len(TRAIN_POINTS26)
    TOTAL_COUNT26 = len(ALL_POINTS26)
    HELDOUT_COUNT26 = TOTAL_COUNT26-TRAIN_COUNT26

    print("Training radial equations inherited from cell 24: {}".format(TRAIN_COUNT26), flush=True)
    print("New held-out radial equations: {}".format(HELDOUT_COUNT26), flush=True)
    print("Total radial equations: {}".format(TOTAL_COUNT26), flush=True)
    print(flush=True)

    def max_abs26(values):
        vals26 = list(values)
        if not vals26:
            return RF(0)
        return max(abs(RF(v26)) for v26 in vals26)

    def base_row_at_X26(X_value26):
        base26 = component_data_at_X(X_value26, None)
        colA26 = component_data_at_X(X_value26, "A")
        colB26 = component_data_at_X(X_value26, "B")
        colC26 = component_data_at_X(X_value26, "C")

        def corrected26(total26):
            return RF(
                total26["Scurv2"] - base26["Scurv2"]
                - alpha0*total26["sigma2"] - beta0*total26["D2"]
            )

        return (
            [
                corrected26(colA26),
                corrected26(colB26),
                corrected26(colC26),
                RF(-base26["sigma1"]),
                RF(-base26["D1"]),
            ],
            RF(-base26["Scurv2"]),
            {
                "sigma1": RF(base26["sigma1"]),
                "D1": RF(base26["D1"]),
                "Rquad": RF(base26["Scurv2"]),
                "rA": corrected26(colA26),
                "rB": corrected26(colB26),
                "rC": corrected26(colC26),
            },
        )

    # Extend the response rows from 70 to the new held-out points.
    BASE_ALL_ROWS26 = [list(row26) for row26 in BASE70_ROWS_24]
    BASE_ALL_RHS26 = [RF(v26) for v26 in BASE70_RHS_24]
    # Recover the 70 inherited component dictionaries algebraically from the
    # cached cell-24 rows.  This avoids 280 expensive component_data_at_X calls.
    # Row convention is [rA,rB,rC,-sigma1,-D1], rhs=-Rquad.
    COMPONENTS_ALL26 = [
        {
            "sigma1": RF(-row26[3]),
            "D1": RF(-row26[4]),
            "Rquad": RF(-rhs26),
            "rA": RF(row26[0]),
            "rB": RF(row26[1]),
            "rC": RF(row26[2]),
        }
        for row26, rhs26 in zip(BASE70_ROWS_24, BASE70_RHS_24)
    ]
    print("Reused the 70 cached cell-24 response rows without recomputation.", flush=True)

    print("Building response rows at new held-out points...", flush=True)
    for index26, X_value26 in enumerate(VALIDATION_X26, start=1):
        row26, rhs26, comp26 = base_row_at_X26(X_value26)
        BASE_ALL_ROWS26.append(row26)
        BASE_ALL_RHS26.append(rhs26)
        COMPONENTS_ALL26.append(comp26)
        print("  completed {}/{} at X={}".format(index26, HELDOUT_COUNT26, X_value26), flush=True)

    if len(BASE_ALL_ROWS26) != TOTAL_COUNT26 or len(COMPONENTS_ALL26) != TOTAL_COUNT26:
        raise ArithmeticError("Live response cache has inconsistent length.")

    def fixed_rhs26(coefficients26, row_count26=TOTAL_COUNT26):
        a26, b26, c26 = [RF(v26) for v26 in coefficients26]
        return [
            RF(BASE_ALL_RHS26[i26])
            - a26*RF(BASE_ALL_ROWS26[i26][0])
            - b26*RF(BASE_ALL_ROWS26[i26][1])
            - c26*RF(BASE_ALL_ROWS26[i26][2])
            for i26 in range(row_count26)
        ]

    def defect_values26(coefficients26, row_count26=TOTAL_COUNT26):
        # fixed_rhs = -F, where F=Rquad+a*rA+b*rB+c*rC.
        return [-v26 for v26 in fixed_rhs26(coefficients26, row_count26)]

    # ------------------------------------------------------------------------
    # PART II-A. Near-divisor valuation estimates.
    # ------------------------------------------------------------------------

    print(flush=True)
    print("PART II-A. U=0 CHANNEL AND DEFECT VALUATIONS", flush=True)
    print("-"*78, flush=True)

    # Exact rational U samples approaching 0 from the V>0 side.
    VALUATION_POWERS26 = [6, 8, 10, 12, 14, 16, 18, 20]
    VALUATION_U26 = [QQ(1)/(QQ(2)**n26) for n26 in VALUATION_POWERS26]
    VALUATION_X26 = [(u26-1)/3 for u26 in VALUATION_U26]

    NEAR_DATA26 = []
    print("Evaluating the degree-two response near U=0...", flush=True)
    for index26, (u26, X_value26) in enumerate(zip(VALUATION_U26, VALUATION_X26), start=1):
        _, _, comp26 = base_row_at_X26(X_value26)
        NEAR_DATA26.append((RF(u26), comp26))
        print("  completed {}/{} at U=2^-{}".format(
            index26, len(VALUATION_U26), VALUATION_POWERS26[index26-1]
        ), flush=True)

    def median26(values26):
        ordered26 = sorted(RF(v26) for v26 in values26)
        n26 = len(ordered26)
        if n26 == 0:
            raise ValueError("median of empty list")
        if n26 % 2:
            return ordered26[n26//2]
        return (ordered26[n26//2-1]+ordered26[n26//2])/2

    def estimate_U_order26(U_values26, function_values26):
        slopes26 = []
        for i26 in range(len(U_values26)-1):
            u_left26 = RF(U_values26[i26])
            u_right26 = RF(U_values26[i26+1])
            f_left26 = RF(function_values26[i26])
            f_right26 = RF(function_values26[i26+1])
            if f_left26 == 0 or f_right26 == 0:
                continue
            slopes26.append(
                (abs(f_right26/f_left26).log())
                /(abs(u_right26/u_left26).log())
            )
        if len(slopes26) < 3:
            return None
        tail26 = slopes26[-min(4, len(slopes26)):]
        estimate26 = median26(tail26)
        rounded26 = ZZ(estimate26.round())
        error26 = max_abs26([s26-RF(rounded26) for s26 in tail26])
        pole_order26 = max(ZZ(0), -rounded26)
        scaled_tail26 = [
            RF(U_values26[i26])**(-rounded26)*RF(function_values26[i26])
            for i26 in range(max(0,len(U_values26)-3), len(U_values26))
        ]
        scale_spread26 = max_abs26([
            scaled_tail26[i26]/scaled_tail26[-1]-1
            for i26 in range(len(scaled_tail26)-1)
            if scaled_tail26[-1] != 0
        ])
        return {
            "order": rounded26,
            "pole_order": pole_order26,
            "estimate": estimate26,
            "tail_error": error26,
            "scaled_tail_spread": scale_spread26,
            "slopes": slopes26,
        }

    sigma_near26 = [comp26["sigma1"] for _, comp26 in NEAR_DATA26]
    D_near26 = [comp26["D1"] for _, comp26 in NEAR_DATA26]
    sigma_val26 = estimate_U_order26([u26 for u26,_ in NEAR_DATA26], sigma_near26)
    D_val26 = estimate_U_order26([u26 for u26,_ in NEAR_DATA26], D_near26)

    print("Degree-one channel valuations:", flush=True)
    for label26, result26 in [("sigma_1",sigma_val26),("D_1",D_val26)]:
        print("  {} ~ U^{} (pole order {}, tail slope error {})".format(
            label26, result26["order"], result26["pole_order"], result26["tail_error"]
        ), flush=True)

    DEFECT_VALUATIONS26 = {}
    for key26, description26, coefficients26 in GEOMETRIES26:
        a26,b26,c26 = [RF(v26) for v26 in coefficients26]
        values26 = []
        for _, comp26 in NEAR_DATA26:
            values26.append(
                comp26["Rquad"]
                + a26*comp26["rA"]
                + b26*comp26["rB"]
                + c26*comp26["rC"]
            )
        val26 = estimate_U_order26([u26 for u26,_ in NEAR_DATA26], values26)
        DEFECT_VALUATIONS26[key26] = val26
        print("  defect {} ~ U^{} (pole order {}, tail slope error {})".format(
            key26, val26["order"], val26["pole_order"], val26["tail_error"]
        ), flush=True)

    valuation_stable26 = all(
        result26 is not None and result26["tail_error"] < RF("0.05")
        for result26 in [sigma_val26, D_val26] + list(DEFECT_VALUATIONS26.values())
    )
    print("Integer valuation estimates stable to 0.05? {}".format(valuation_stable26), flush=True)

    # ------------------------------------------------------------------------
    # PART II-B. Transform and re-solve in the pole-regularized channel gauge.
    # ------------------------------------------------------------------------

    print(flush=True)
    print("PART II-B. POLE-REGULARIZED CONNECTION GAUGE", flush=True)
    print("-"*78, flush=True)

    RAW_ALPHA_EXP26 = [0, -3, -2, -1, 1, 2]
    RAW_BETA_EXP26 = [0, -3, -2]

    def raw_connection_dicts26(solution26):
        return (
            {
                0: RF(solution26[0]),
                -3: RF(solution26[2]),
                -2: RF(solution26[4]),
                -1: RF(solution26[6]),
                1: RF(solution26[7]),
                2: RF(solution26[8]),
            },
            {
                0: RF(solution26[1]),
                -3: RF(solution26[3]),
                -2: RF(solution26[5]),
            },
        )

    def shifted_dict26(data26, shift26):
        out26 = {}
        for exponent26, coefficient26 in data26.items():
            new_exp26 = ZZ(exponent26)+ZZ(shift26)
            out26[new_exp26] = out26.get(new_exp26, RF(0))+RF(coefficient26)
        return {e26:c26 for e26,c26 in out26.items() if abs(c26)>DISPLAY_ZERO26}

    def eval_laurent26(data26, U_value26):
        Ureal26 = RF(U_value26)
        return sum(RF(c26)*Ureal26**ZZ(e26) for e26,c26 in data26.items())

    def regularized_rows26(alpha_exponents26, beta_exponents26,
                           nu_sigma26, nu_D26, row_count26=TOTAL_COUNT26):
        rows26 = []
        for i26 in range(row_count26):
            Ureal26 = RF(3)*RF(ALL_POINTS26[i26])+RF(1)
            sigma_hat26 = Ureal26**nu_sigma26*COMPONENTS_ALL26[i26]["sigma1"]
            D_hat26 = Ureal26**nu_D26*COMPONENTS_ALL26[i26]["D1"]
            rows26.append(
                [Ureal26**e26*sigma_hat26 for e26 in alpha_exponents26]
                + [Ureal26**e26*D_hat26 for e26 in beta_exponents26]
            )
        return rows26

    GAUGE_RESULTS26 = {}

    for key26, description26, coefficients26 in GEOMETRIES26:
        print(flush=True)
        print("GEOMETRY: {}".format(description26), flush=True)

        raw_solution26 = solve_connection_24(coefficients26, TOTAL_COUNT_24)
        alpha_raw26, beta_raw26 = raw_connection_dicts26(raw_solution26)

        nu_sigma26 = sigma_val26["pole_order"]
        nu_D26 = D_val26["pole_order"]
        nu_F26 = DEFECT_VALUATIONS26[key26]["pole_order"]
        alpha_shift26 = nu_F26-nu_sigma26
        beta_shift26 = nu_F26-nu_D26

        alpha_hat26 = shifted_dict26(alpha_raw26, alpha_shift26)
        beta_hat26 = shifted_dict26(beta_raw26, beta_shift26)
        alpha_exponents26 = sorted(alpha_hat26)
        beta_exponents26 = sorted(beta_hat26)

        raw_poles26 = sum(1 for e26 in alpha_raw26 if e26 < 0) + sum(1 for e26 in beta_raw26 if e26 < 0)
        reg_poles26 = sum(1 for e26 in alpha_hat26 if e26 < 0) + sum(1 for e26 in beta_hat26 if e26 < 0)

        print("  regularization powers (nu_F,nu_sigma,nu_D) = ({},{},{})".format(
            nu_F26, nu_sigma26, nu_D26
        ), flush=True)
        print("  coefficient exponent shifts (alpha,beta) = ({},{})".format(
            alpha_shift26, beta_shift26
        ), flush=True)
        print("  raw alpha support: {}".format(sorted(alpha_raw26)), flush=True)
        print("  raw beta support:  {}".format(sorted(beta_raw26)), flush=True)
        print("  regularized alpha support: {}".format(alpha_exponents26), flush=True)
        print("  regularized beta support:  {}".format(beta_exponents26), flush=True)
        print("  negative-power terms before/after: {}/{}".format(raw_poles26, reg_poles26), flush=True)

        defect26 = defect_values26(coefficients26, TOTAL_COUNT26)
        scaled_rhs26 = []
        direct_residuals26 = []
        for i26 in range(TOTAL_COUNT26):
            Ureal26 = RF(3)*RF(ALL_POINTS26[i26])+RF(1)
            sigma_hat26 = Ureal26**nu_sigma26*COMPONENTS_ALL26[i26]["sigma1"]
            D_hat26 = Ureal26**nu_D26*COMPONENTS_ALL26[i26]["D1"]
            F_hat26 = Ureal26**nu_F26*defect26[i26]
            scaled_rhs26.append(F_hat26)
            direct_residuals26.append(
                F_hat26
                - eval_laurent26(alpha_hat26,Ureal26)*sigma_hat26
                - eval_laurent26(beta_hat26,Ureal26)*D_hat26
            )

        Areg26 = regularized_rows26(
            alpha_exponents26, beta_exponents26,
            nu_sigma26, nu_D26, TOTAL_COUNT26
        )
        pair_reg26 = rank_pair_23(Areg26, scaled_rhs26)
        solve_reg26 = canonical_solution_23(Areg26[:TRAIN_COUNT26], scaled_rhs26[:TRAIN_COUNT26])
        heldout_reg26 = None
        if solve_reg26 is not None:
            heldout_reg26 = max_abs26(residuals_23(
                Areg26[TRAIN_COUNT26:],
                scaled_rhs26[TRAIN_COUNT26:],
                solve_reg26["solution"],
            ))

        expected_solution26 = [alpha_hat26[e26] for e26 in alpha_exponents26] + [beta_hat26[e26] for e26 in beta_exponents26]
        expected_residual26 = max_abs26(direct_residuals26)
        expected_matrix_residual26 = max_abs26(residuals_23(Areg26, scaled_rhs26, expected_solution26))

        print("  regularized all-point rank(A)={} rank([A|b])={}".format(*pair_reg26), flush=True)
        print("  transformed-connection residual: {}".format(expected_residual26), flush=True)
        print("  transformed matrix residual: {}".format(expected_matrix_residual26), flush=True)
        print("  70-fit / held-out regularized residual: {}".format(heldout_reg26), flush=True)
        print("  all polar terms removed? {}".format(reg_poles26 == 0), flush=True)

        GAUGE_RESULTS26[key26] = {
            "nu_F":nu_F26, "nu_sigma":nu_sigma26, "nu_D":nu_D26,
            "alpha_shift":alpha_shift26, "beta_shift":beta_shift26,
            "alpha_support":alpha_exponents26, "beta_support":beta_exponents26,
            "raw_poles":raw_poles26, "regularized_poles":reg_poles26,
            "pair":pair_reg26, "heldout_residual":heldout_reg26,
            "direct_residual":expected_residual26,
        }

    # ------------------------------------------------------------------------
    # PART II-C. Matched false-divisor controls.
    # ------------------------------------------------------------------------

    print(flush=True)
    print("PART II-C. MATCHED FALSE-DIVISOR LAURENT CONTROLS", flush=True)
    print("-"*78, flush=True)

    FACTORS26 = [
        ("U", "3X+1 (classical-volume divisor; positive control)", lambda X26: RF(3)*RF(X26)+RF(1)),
        ("W", "2X+1 (frame-component control)", lambda X26: RF(2)*RF(X26)+RF(1)),
        ("T", "3X+2 (shifted-volume false divisor)", lambda X26: RF(3)*RF(X26)+RF(2)),
        ("Q", "X+1 (affine false divisor)", lambda X26: RF(X26)+RF(1)),
        ("R", "X (radial-coordinate false divisor)", lambda X26: RF(X26)),
    ]

    def factor_rows26(factor_function26, row_count26=TOTAL_COUNT26):
        rows26 = []
        for i26 in range(row_count26):
            t26 = RF(factor_function26(ALL_POINTS26[i26]))
            if t26 == 0:
                raise ZeroDivisionError("Control factor vanished at a sampled point.")
            sigma26 = COMPONENTS_ALL26[i26]["sigma1"]
            D26 = COMPONENTS_ALL26[i26]["D1"]
            rows26.append([
                sigma26,
                D26,
                t26**(-3)*sigma26,
                t26**(-3)*D26,
                t26**(-2)*sigma26,
                t26**(-2)*D26,
                t26**(-1)*sigma26,
                t26*sigma26,
                t26**2*sigma26,
            ])
        return rows26

    FACTOR_ROWS26 = {key26:factor_rows26(function26) for key26,_,function26 in FACTORS26}
    FALSE_RESULTS26 = {}

    for key_geom26, description_geom26, coefficients26 in GEOMETRIES26:
        print(flush=True)
        print("GEOMETRY: {}".format(description_geom26), flush=True)
        rhs_all26 = defect_values26(coefficients26, TOTAL_COUNT26)
        FALSE_RESULTS26[key_geom26] = {}

        for key_factor26, description_factor26, _ in FACTORS26:
            Aall26 = FACTOR_ROWS26[key_factor26]
            Atrain26 = Aall26[:TRAIN_COUNT26]
            rhs_train26 = rhs_all26[:TRAIN_COUNT26]
            pair_train26 = rank_pair_23(Atrain26, rhs_train26)
            pair_all26 = rank_pair_23(Aall26, rhs_all26)
            closes26 = pair_all26[0] == pair_all26[1]
            saturated26 = pair_all26[0] == TOTAL_COUNT26
            structural26 = closes26 and not saturated26

            train_solution26 = canonical_solution_23(Atrain26, rhs_train26)
            heldout26 = None
            all_residual26 = None
            if train_solution26 is not None:
                heldout26 = max_abs26(residuals_23(
                    Aall26[TRAIN_COUNT26:], rhs_all26[TRAIN_COUNT26:],
                    train_solution26["solution"]
                ))
                all_residual26 = max_abs26(residuals_23(
                    Aall26, rhs_all26, train_solution26["solution"]
                ))

            print("  {}: rank(A)={} rank([A|b])={} structural={} held-out={}".format(
                key_factor26, pair_all26[0], pair_all26[1], structural26, heldout26
            ), flush=True)

            FALSE_RESULTS26[key_geom26][key_factor26] = {
                "description":description_factor26,
                "train_pair":pair_train26,
                "all_pair":pair_all26,
                "closes":closes26,
                "saturated":saturated26,
                "structural":structural26,
                "heldout_residual":heldout26,
                "all_residual":all_residual26,
            }

    positive_control26 = all(FALSE_RESULTS26[g26]["U"]["structural"] for g26,_,_ in GEOMETRIES26)
    negative_controls26 = all(
        not FALSE_RESULTS26[g26][f26]["closes"]
        for g26,_,_ in GEOMETRIES26
        for f26 in ["W","T","Q","R"]
    )

    # ------------------------------------------------------------------------
    # Save compact checkpoint.
    # ------------------------------------------------------------------------

    checkpoint26 = {
        "exact_frame": {
            "g_norms": {"e":6,"o":2,"r":3},
            "det_g": "9/(3X+1)^2",
            "M_norms": {"e":"-18U","o":"-6U","r":"18U"},
        },
        "valuation_U_powers": VALUATION_POWERS26,
        "sigma_valuation": sigma_val26,
        "D_valuation": D_val26,
        "defect_valuations": DEFECT_VALUATIONS26,
        "gauge_results": GAUGE_RESULTS26,
        "false_divisor_results": FALSE_RESULTS26,
        "positive_control_passed": positive_control26,
        "all_negative_controls_failed": negative_controls26,
        "cell25_checkpoint_loaded": CHECKPOINT25_26 is not None,
    }
    save(checkpoint26, "channel_gauge_false_divisor_checkpoint.sobj")

    # ------------------------------------------------------------------------
    # Final summary.
    # ------------------------------------------------------------------------

    print(flush=True)
    print("="*78, flush=True)
    print("CHANNEL-GAUGE / FALSE-DIVISOR SUMMARY", flush=True)
    print("  exact classical frame certified? True", flush=True)
    print("  intrinsic classical divisor: U=3X+1", flush=True)
    print("  channel valuation sigma_1: order {} / pole {}".format(
        sigma_val26["order"], sigma_val26["pole_order"]
    ), flush=True)
    print("  channel valuation D_1:     order {} / pole {}".format(
        D_val26["order"], D_val26["pole_order"]
    ), flush=True)
    print("  integer valuation estimates stable? {}".format(valuation_stable26), flush=True)
    for key26,_,_ in GEOMETRIES26:
        gauge26 = GAUGE_RESULTS26[key26]
        print("  {}: defect pole={} ; polar terms raw->regularized {}->{}".format(
            key26, gauge26["nu_F"], gauge26["raw_poles"], gauge26["regularized_poles"]
        ), flush=True)
    print("  U positive control closes every fixed geometry? {}".format(positive_control26), flush=True)
    print("  W,T,Q,R controls all fail every fixed geometry? {}".format(negative_controls26), flush=True)
    print("  checkpoint written: channel_gauge_false_divisor_checkpoint.sobj", flush=True)
    print("="*78, flush=True)
    print("SUCCESS", flush=True)
    print("The pole-regularized gauge and matched false-divisor controls completed.", flush=True)
    print("Copy the valuation block, each gauge support, and the final summary back.", flush=True)
    print("="*78, flush=True)
