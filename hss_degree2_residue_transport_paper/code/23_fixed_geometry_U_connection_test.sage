from sage.all import *
from itertools import combinations

# ============================================================================
# FIXED DEGREE-TWO GEOMETRY / U-CONNECTION REDUCTION TEST
# ============================================================================
#
# Run in the SAME SageMath session after successful cells 13--22.
#
# PURPOSE
# -------
# Cell 22 found a non-row-saturated closure after promoting alpha_1,beta_1
# from constants to Laurent functions of
#
#     U = 3 X + 1 = V(X,1,1)/3.
#
# In that experiment the three degree-two factorized coefficients (a,b,c)
# were still free.  This cell removes that flexibility and tests the three
# concrete geometric prescriptions already studied:
#
#   primitive product:            (a,b,c) = (1, 0,    0)
#   primitive + double cover:     (a,b,c) = (1, 1/8,  0)
#   Z_2(y) Z_2(z):                (a,b,c) = (1, 1/64, 1/8)
#
# For each fixed geometry we test:
#
#   (i)  constant alpha_1,beta_1 only;
#   (ii) a COMMON U-connection
#
#            alpha_1(U) = A_0 + alpha_0 gamma(U),
#            beta_1(U)  = B_0 + beta_0  gamma(U),
#
#        gamma(U)=sum_{k=+-1,+-2,+-3} g_k U^k;
#
#   (iii) independent U-Laurent functions alpha_1(U), beta_1(U).
#
# The cell also performs the previously missing FREE-GEOMETRY / COMMON-U
# benchmark, because cell 22 tested a common U/W module but not common U alone.
#
# ROBUSTNESS
# ----------
# * Reuses the 33 points from cell 22 as a training set.
# * Adds 17 completely new positive rational X-values.
# * Tests all 50 points simultaneously.
# * Fits a training solution separately and reports its held-out residual.
# * Searches all common-U supports and a quotient basis of the independent-U
#   module for the smallest support closing all 50 equations.
# * No tested system has enough columns to row-saturate 50 equations.
# ============================================================================

print("="*78, flush=True)
print("FIXED DEGREE-TWO GEOMETRY / U-CONNECTION REDUCTION TEST", flush=True)
print("Locking (a,b,c) and reducing the evolving connection to U=3X+1", flush=True)
print("="*78, flush=True)
print(flush=True)

required_names = [
    "RF", "alpha0", "beta0", "component_data_at_X",
    "BASE_ROWS", "BASE_RHS", "ALL_CONNECTION_POINTS",
    "numerical_rref", "equilibrate_system",
]
missing = [name for name in required_names if name not in globals()]
if missing:
    raise RuntimeError(
        "Run cells 13 through 22 first in the same Sage session. Missing: {}"
        .format(missing)
    )

RANK_TOL_EXPONENT = 70
RANK_TOL = RF(10)**(-RANK_TOL_EXPONENT)
LAURENT_EXPONENTS_23 = [-3, -2, -1, 1, 2, 3]

# Disjoint from the 33 points used through cell 22.
FIXED_GEOMETRY_VALIDATION_X = [
    QQ(2)/5,
    QQ(4)/7,
    QQ(5)/7,
    QQ(7)/8,
    QQ(19)/20,
    QQ(21)/20,
    QQ(9)/8,
    QQ(11)/8,
    QQ(10)/7,
    QQ(12)/7,
    QQ(15)/8,
    QQ(19)/8,
    QQ(14)/5,
    QQ(16)/5,
    QQ(11)/3,
    QQ(13)/3,
    QQ(11)/2,
]

TRAIN_POINTS_23 = list(ALL_CONNECTION_POINTS)
ALL_POINTS_23 = list(TRAIN_POINTS_23)
for value in FIXED_GEOMETRY_VALIDATION_X:
    if value not in ALL_POINTS_23:
        ALL_POINTS_23.append(value)

TRAIN_COUNT_23 = len(TRAIN_POINTS_23)
TOTAL_COUNT_23 = len(ALL_POINTS_23)
VALIDATION_COUNT_23 = TOTAL_COUNT_23 - TRAIN_COUNT_23

print("Working real precision: {} bits".format(RF.precision()), flush=True)
print("Training equations inherited from cell 22: {}".format(TRAIN_COUNT_23), flush=True)
print("New held-out validation equations: {}".format(VALIDATION_COUNT_23), flush=True)
print("Total radial equations: {}".format(TOTAL_COUNT_23), flush=True)
print("U-Laurent powers: {}".format(LAURENT_EXPONENTS_23), flush=True)
print(flush=True)


def max_abs_23(values):
    values = list(values)
    if not values:
        return RF(0)
    return max(abs(RF(value)) for value in values)


def rank_pair_23(A_rows, rhs_values, tolerance=RANK_TOL):
    A_eq, b_eq, _ = equilibrate_system(A_rows, rhs_values)
    augmented = [row + [rhs] for row, rhs in zip(A_eq, b_eq)]
    return (
        len(numerical_rref(A_eq, tolerance)[1]),
        len(numerical_rref(augmented, tolerance)[1]),
    )


def compact_rank_profile_23(A_rows, rhs_values):
    decimal_digits = max(20, int(RF.precision()*0.3010299956639812))
    exponents = [30, 70, 120, 192, min(345, decimal_digits-10)]
    exponents = sorted(set(e for e in exponents if 10 <= e <= decimal_digits-10))
    return [
        (exponent,) + rank_pair_23(A_rows, rhs_values, RF(10)**(-exponent))
        for exponent in exponents
    ]


def canonical_solution_23(A_rows, rhs_values, tolerance=RANK_TOL):
    """Canonical RREF solution: every free equilibrated variable is set to zero."""
    A_eq, b_eq, column_scales = equilibrate_system(A_rows, rhs_values)
    augmented = [row + [rhs] for row, rhs in zip(A_eq, b_eq)]
    R_A, pivots_A = numerical_rref(A_eq, tolerance)
    R_aug, pivots_aug = numerical_rref(augmented, tolerance)

    if len(pivots_A) != len(pivots_aug):
        return None

    n_columns = len(A_rows[0])
    z = [RF(0)]*n_columns
    for row_index, pivot in enumerate(pivots_A):
        z[pivot] = R_aug[row_index][-1]

    solution = [z[index]/column_scales[index] for index in range(n_columns)]
    return {
        "solution": solution,
        "rank": len(pivots_A),
        "nullity": n_columns-len(pivots_A),
        "pivots": list(pivots_A),
    }


def residuals_23(A_rows, rhs_values, solution):
    return [
        sum(RF(a)*RF(x) for a, x in zip(row, solution)) - RF(rhs)
        for row, rhs in zip(A_rows, rhs_values)
    ]


def validation_residual_23(A_all, rhs_all, training_solution):
    if training_solution is None:
        return None
    solution = training_solution["solution"]
    return max_abs_23(residuals_23(
        A_all[TRAIN_COUNT_23:], rhs_all[TRAIN_COUNT_23:], solution
    ))


# ----------------------------------------------------------------------------
# Extend the five-column response system from 33 to 50 points.
# Column convention:
#
#     [rA, rB, rC, -sigma_1, -D_1],     rhs = -Rquad.
# ----------------------------------------------------------------------------

BASE50_ROWS = [list(row) for row in BASE_ROWS]
BASE50_RHS = [RF(value) for value in BASE_RHS]

if len(BASE50_ROWS) != TRAIN_COUNT_23:
    raise ArithmeticError(
        "Cell-22 BASE_ROWS length {} does not match its point count {}."
        .format(len(BASE50_ROWS), TRAIN_COUNT_23)
    )

print("Building the old response system at 17 new validation points...", flush=True)
for index, X_value in enumerate(FIXED_GEOMETRY_VALIDATION_X, start=1):
    base = component_data_at_X(X_value, None)
    colA = component_data_at_X(X_value, "A")
    colB = component_data_at_X(X_value, "B")
    colC = component_data_at_X(X_value, "C")

    def corrected(total):
        return RF(
            total["Scurv2"] - base["Scurv2"]
            - alpha0*total["sigma2"] - beta0*total["D2"]
        )

    BASE50_ROWS.append([
        corrected(colA),
        corrected(colB),
        corrected(colC),
        RF(-base["sigma1"]),
        RF(-base["D1"]),
    ])
    BASE50_RHS.append(RF(-base["Scurv2"]))

    print("  completed {}/{} at X={}".format(
        index, len(FIXED_GEOMETRY_VALIDATION_X), X_value
    ), flush=True)

if len(BASE50_ROWS) != TOTAL_COUNT_23:
    raise ArithmeticError("Expanded base-row count is inconsistent.")

base_train_pair_23 = rank_pair_23(
    BASE50_ROWS[:TRAIN_COUNT_23], BASE50_RHS[:TRAIN_COUNT_23]
)
base_all_pair_23 = rank_pair_23(BASE50_ROWS, BASE50_RHS)

print(flush=True)
print("Baseline free-geometry checks:", flush=True)
print("  33-point training rank(A)={} rank([A|b])={}".format(
    *base_train_pair_23
), flush=True)
print("  50-point total rank(A)={} rank([A|b])={}".format(
    *base_all_pair_23
), flush=True)
print(flush=True)

if base_train_pair_23 != (4, 5) or base_all_pair_23 != (4, 5):
    raise ArithmeticError("The old response system did not retain rank 4/5.")


# ----------------------------------------------------------------------------
# Construct U-dependent connection columns on all 50 points.
# ----------------------------------------------------------------------------

COMMON_U_NAMES_23 = ["gamma_U_{}".format(k) for k in LAURENT_EXPONENTS_23]
INDEPENDENT_U_NAMES_23 = []
for exponent in LAURENT_EXPONENTS_23:
    INDEPENDENT_U_NAMES_23.extend([
        "alpha_U_{}".format(exponent),
        "beta_U_{}".format(exponent),
    ])

COMMON_U_COLUMNS_23 = {}
INDEPENDENT_U_COLUMNS_23 = {}

for exponent in LAURENT_EXPONENTS_23:
    common_name = "gamma_U_{}".format(exponent)
    alpha_name = "alpha_U_{}".format(exponent)
    beta_name = "beta_U_{}".format(exponent)

    common_column = []
    alpha_column = []
    beta_column = []

    for X_value, row in zip(ALL_POINTS_23, BASE50_ROWS):
        U = RF(3)*RF(X_value) + RF(1)
        multiplier = U**exponent
        minus_sigma = RF(row[3])
        minus_D = RF(row[4])

        alpha_column.append(multiplier*minus_sigma)
        beta_column.append(multiplier*minus_D)
        common_column.append(
            multiplier*(alpha0*minus_sigma + beta0*minus_D)
        )

    COMMON_U_COLUMNS_23[common_name] = common_column
    INDEPENDENT_U_COLUMNS_23[alpha_name] = alpha_column
    INDEPENDENT_U_COLUMNS_23[beta_name] = beta_column


# ----------------------------------------------------------------------------
# Generic assemblers.
# ----------------------------------------------------------------------------

GEOMETRIES_23 = [
    (
        "primitive_product",
        "primitive product F_A",
        (QQ(1), QQ(0), QQ(0)),
    ),
    (
        "primitive_plus_double_cover",
        "F_A + F_B/8",
        (QQ(1), QQ(1)/8, QQ(0)),
    ),
    (
        "Z2_product",
        "Z_2(y) Z_2(z) = F_A + F_C/8 + F_B/64",
        (QQ(1), QQ(1)/64, QQ(1)/8),
    ),
]


def fixed_rhs_23(coefficients, row_count=TOTAL_COUNT_23):
    a, b, c = [RF(value) for value in coefficients]
    return [
        RF(BASE50_RHS[row])
        - a*RF(BASE50_ROWS[row][0])
        - b*RF(BASE50_ROWS[row][1])
        - c*RF(BASE50_ROWS[row][2])
        for row in range(row_count)
    ]


def fixed_constant_rows_23(row_count=TOTAL_COUNT_23):
    return [
        [RF(BASE50_ROWS[row][3]), RF(BASE50_ROWS[row][4])]
        for row in range(row_count)
    ]


def fixed_common_rows_23(names, row_count=TOTAL_COUNT_23):
    return [
        [RF(BASE50_ROWS[row][3]), RF(BASE50_ROWS[row][4])]
        + [COMMON_U_COLUMNS_23[name][row] for name in names]
        for row in range(row_count)
    ]


def fixed_independent_rows_23(names, row_count=TOTAL_COUNT_23):
    return [
        [RF(BASE50_ROWS[row][3]), RF(BASE50_ROWS[row][4])]
        + [INDEPENDENT_U_COLUMNS_23[name][row] for name in names]
        for row in range(row_count)
    ]


def free_common_rows_23(names, row_count=TOTAL_COUNT_23):
    return [
        list(BASE50_ROWS[row])
        + [COMMON_U_COLUMNS_23[name][row] for name in names]
        for row in range(row_count)
    ]


def free_independent_rows_23(names, row_count=TOTAL_COUNT_23):
    return [
        list(BASE50_ROWS[row])
        + [INDEPENDENT_U_COLUMNS_23[name][row] for name in names]
        for row in range(row_count)
    ]


def test_family_23(label, A_train, rhs_train, A_all, rhs_all, verbose=True):
    train_pair = rank_pair_23(A_train, rhs_train)
    all_pair = rank_pair_23(A_all, rhs_all)
    closes = bool(all_pair[0] == all_pair[1])
    saturated = bool(all_pair[0] == TOTAL_COUNT_23)
    structural = bool(closes and not saturated)

    training_solution = canonical_solution_23(A_train, rhs_train)
    heldout_residual = validation_residual_23(A_all, rhs_all, training_solution)

    all_solution = canonical_solution_23(A_all, rhs_all) if closes else None
    all_residual = None
    if all_solution is not None:
        all_residual = max_abs_23(residuals_23(
            A_all, rhs_all, all_solution["solution"]
        ))

    if verbose:
        print("-"*78, flush=True)
        print(label, flush=True)
        print("  columns: {}".format(len(A_all[0])), flush=True)
        print("  training rank(A)={} rank([A|b])={}".format(*train_pair), flush=True)
        print("  all-point rank(A)={} rank([A|b])={}".format(*all_pair), flush=True)
        print("  closes all 50 equations? {}".format(closes), flush=True)
        print("  row-saturated interpolation? {}".format(saturated), flush=True)
        print("  non-saturated structural closure? {}".format(structural), flush=True)
        if training_solution is not None:
            print("  training canonical nullity: {}".format(
                training_solution["nullity"]
            ), flush=True)
            print("  held-out residual of training canonical solution: {}".format(
                heldout_residual
            ), flush=True)
        if all_solution is not None:
            print("  all-point canonical nullity: {}".format(
                all_solution["nullity"]
            ), flush=True)
            print("  all-point canonical residual: {}".format(
                all_residual
            ), flush=True)
        print("  stable all-point rank profile:", flush=True)
        for exponent, rA, rAug in compact_rank_profile_23(A_all, rhs_all):
            print("    1e-{:<4d}: rank(A)={} rank([A|b])={}".format(
                exponent, rA, rAug
            ), flush=True)

    return {
        "label": label,
        "A_train": A_train,
        "rhs_train": rhs_train,
        "A_all": A_all,
        "rhs_all": rhs_all,
        "train_pair": train_pair,
        "all_pair": all_pair,
        "closes": closes,
        "saturated": saturated,
        "structural": structural,
        "training_solution": training_solution,
        "heldout_residual": heldout_residual,
        "all_solution": all_solution,
        "all_residual": all_residual,
    }


def quotient_pivot_names_23(A_rows, rhs_values, base_column_count, candidate_names):
    A_eq, _, _ = equilibrate_system(A_rows, rhs_values)
    _, pivots = numerical_rref(A_eq, RANK_TOL)
    indices = [pivot-base_column_count for pivot in pivots if pivot >= base_column_count]
    return [candidate_names[index] for index in indices]


def minimal_common_support_23(rhs_all, rhs_train):
    closers = []
    for size in range(1, len(COMMON_U_NAMES_23)+1):
        current = []
        for subset in combinations(COMMON_U_NAMES_23, size):
            A_all = fixed_common_rows_23(subset, TOTAL_COUNT_23)
            pair_all = rank_pair_23(A_all, rhs_all)
            if pair_all[0] != pair_all[1] or pair_all[0] == TOTAL_COUNT_23:
                continue
            A_train = fixed_common_rows_23(subset, TRAIN_COUNT_23)
            pair_train = rank_pair_23(A_train, rhs_train)
            if pair_train[0] != pair_train[1]:
                continue
            current.append((subset, pair_all))
        if current:
            closers = current
            break
    return closers


def minimal_independent_support_23(rhs_all, rhs_train):
    A_full = fixed_independent_rows_23(INDEPENDENT_U_NAMES_23, TOTAL_COUNT_23)
    full_pair = rank_pair_23(A_full, rhs_all)
    if full_pair[0] != full_pair[1]:
        return [], []

    quotient_names = quotient_pivot_names_23(
        A_full, rhs_all, 2, INDEPENDENT_U_NAMES_23
    )

    closers = []
    for size in range(1, len(quotient_names)+1):
        current = []
        for subset in combinations(quotient_names, size):
            A_all = fixed_independent_rows_23(subset, TOTAL_COUNT_23)
            pair_all = rank_pair_23(A_all, rhs_all)
            if pair_all[0] != pair_all[1] or pair_all[0] == TOTAL_COUNT_23:
                continue
            A_train = fixed_independent_rows_23(subset, TRAIN_COUNT_23)
            pair_train = rank_pair_23(A_train, rhs_train)
            if pair_train[0] != pair_train[1]:
                continue
            current.append((subset, pair_all))
        if current:
            closers = current
            break

    return quotient_names, closers


def print_common_formula_23(result, names):
    if not result["closes"] or result["all_solution"] is None:
        return
    values = result["all_solution"]["solution"]
    alpha_const = values[0]
    beta_const = values[1]
    gamma_coefficients = values[2:]

    print("  canonical common-U formula:", flush=True)
    print("    alpha_1(U) = {} + alpha_0 * gamma(U)".format(
        alpha_const
    ), flush=True)
    print("    beta_1(U)  = {} + beta_0  * gamma(U)".format(
        beta_const
    ), flush=True)
    print("    gamma(U) nonzero terms:", flush=True)
    nonzero = False
    for name, coefficient in zip(names, gamma_coefficients):
        if abs(RF(coefficient)) <= RF(10)**(-60):
            continue
        exponent = ZZ(name.split("_")[-1])
        print("      {} * U^{}".format(coefficient, exponent), flush=True)
        nonzero = True
    if not nonzero:
        print("      none", flush=True)


def print_independent_formula_23(result, names):
    if not result["closes"] or result["all_solution"] is None:
        return
    values = result["all_solution"]["solution"]
    alpha_const = values[0]
    beta_const = values[1]

    print("  canonical independent-U formula:", flush=True)
    print("    alpha_1(U) constant = {}".format(alpha_const), flush=True)
    print("    beta_1(U) constant  = {}".format(beta_const), flush=True)
    print("    nonzero Laurent terms:", flush=True)
    nonzero = False
    for name, coefficient in zip(names, values[2:]):
        if abs(RF(coefficient)) <= RF(10)**(-60):
            continue
        print("      {}: {}".format(name, coefficient), flush=True)
        nonzero = True
    if not nonzero:
        print("      none", flush=True)


# ----------------------------------------------------------------------------
# Benchmark: free geometry with U-only connection modules.
# ----------------------------------------------------------------------------

print("Free-geometry U-only benchmarks:", flush=True)

FREE_COMMON_RESULT_23 = test_family_23(
    "FREE (a,b,c) + COMMON U-CONNECTION",
    free_common_rows_23(COMMON_U_NAMES_23, TRAIN_COUNT_23),
    BASE50_RHS[:TRAIN_COUNT_23],
    free_common_rows_23(COMMON_U_NAMES_23, TOTAL_COUNT_23),
    BASE50_RHS,
)

FREE_INDEPENDENT_RESULT_23 = test_family_23(
    "FREE (a,b,c) + INDEPENDENT U-CONNECTION",
    free_independent_rows_23(INDEPENDENT_U_NAMES_23, TRAIN_COUNT_23),
    BASE50_RHS[:TRAIN_COUNT_23],
    free_independent_rows_23(INDEPENDENT_U_NAMES_23, TOTAL_COUNT_23),
    BASE50_RHS,
)

if FREE_COMMON_RESULT_23["closes"]:
    # Free system has [a,b,c,A0,B0,gamma...].
    values = FREE_COMMON_RESULT_23["all_solution"]["solution"]
    print("  free common-U factorized coefficients (a,b,c) = {}".format(
        tuple(values[:3])
    ), flush=True)
    print("  free common-U constants (A0,B0) = {}".format(
        tuple(values[3:5])
    ), flush=True)
    print("  free common-U gamma terms:", flush=True)
    for name, coefficient in zip(COMMON_U_NAMES_23, values[5:]):
        if abs(RF(coefficient)) > RF(10)**(-60):
            print("    {}: {}".format(name, coefficient), flush=True)


# ----------------------------------------------------------------------------
# Fixed-geometry tests.
# ----------------------------------------------------------------------------

FIXED_RESULTS_23 = {}

print(flush=True)
print("Fixed-geometry connection tests:", flush=True)

for key, description, coefficients in GEOMETRIES_23:
    print(flush=True)
    print("="*78, flush=True)
    print("GEOMETRY: {}".format(description), flush=True)
    print("  fixed (a,b,c) = {}".format(coefficients), flush=True)
    print("="*78, flush=True)

    rhs_all = fixed_rhs_23(coefficients, TOTAL_COUNT_23)
    rhs_train = rhs_all[:TRAIN_COUNT_23]

    constant_result = test_family_23(
        "CONSTANT alpha_1,beta_1 ONLY",
        fixed_constant_rows_23(TRAIN_COUNT_23),
        rhs_train,
        fixed_constant_rows_23(TOTAL_COUNT_23),
        rhs_all,
    )

    common_result = test_family_23(
        "COMMON U-CONNECTION aligned with (alpha_0,beta_0)",
        fixed_common_rows_23(COMMON_U_NAMES_23, TRAIN_COUNT_23),
        rhs_train,
        fixed_common_rows_23(COMMON_U_NAMES_23, TOTAL_COUNT_23),
        rhs_all,
    )
    print_common_formula_23(common_result, COMMON_U_NAMES_23)

    independent_result = test_family_23(
        "INDEPENDENT U-LAURENT alpha_1(U), beta_1(U)",
        fixed_independent_rows_23(INDEPENDENT_U_NAMES_23, TRAIN_COUNT_23),
        rhs_train,
        fixed_independent_rows_23(INDEPENDENT_U_NAMES_23, TOTAL_COUNT_23),
        rhs_all,
    )
    print_independent_formula_23(independent_result, INDEPENDENT_U_NAMES_23)

    common_minimal = minimal_common_support_23(rhs_all, rhs_train)
    independent_basis, independent_minimal = minimal_independent_support_23(
        rhs_all, rhs_train
    )

    print(flush=True)
    print("Minimal common-U supports closing all 50 points:", flush=True)
    if common_minimal:
        for subset, pair in common_minimal:
            print("  {} -> rank(A)={} rank([A|b])={}".format(
                subset, pair[0], pair[1]
            ), flush=True)
    else:
        print("  none", flush=True)

    print("Independent-U quotient basis:", independent_basis, flush=True)
    print("Minimal independent-U supports closing all 50 points:", flush=True)
    if independent_minimal:
        for subset, pair in independent_minimal:
            print("  {} -> rank(A)={} rank([A|b])={}".format(
                subset, pair[0], pair[1]
            ), flush=True)
    else:
        print("  none", flush=True)

    # Print a formula for the first minimal common support, when present.
    minimal_common_result = None
    if common_minimal:
        minimal_names = list(common_minimal[0][0])
        minimal_common_result = test_family_23(
            "MINIMAL COMMON-U REPRESENTATIVE",
            fixed_common_rows_23(minimal_names, TRAIN_COUNT_23),
            rhs_train,
            fixed_common_rows_23(minimal_names, TOTAL_COUNT_23),
            rhs_all,
            verbose=False,
        )
        print_common_formula_23(minimal_common_result, minimal_names)
        print("  minimal common-U held-out residual: {}".format(
            minimal_common_result["heldout_residual"]
        ), flush=True)

    minimal_independent_result = None
    if independent_minimal:
        minimal_names = list(independent_minimal[0][0])
        minimal_independent_result = test_family_23(
            "MINIMAL INDEPENDENT-U REPRESENTATIVE",
            fixed_independent_rows_23(minimal_names, TRAIN_COUNT_23),
            rhs_train,
            fixed_independent_rows_23(minimal_names, TOTAL_COUNT_23),
            rhs_all,
            verbose=False,
        )
        print_independent_formula_23(minimal_independent_result, minimal_names)
        print("  minimal independent-U held-out residual: {}".format(
            minimal_independent_result["heldout_residual"]
        ), flush=True)

    FIXED_RESULTS_23[key] = {
        "description": description,
        "coefficients": coefficients,
        "constant": constant_result,
        "common": common_result,
        "independent": independent_result,
        "common_minimal": common_minimal,
        "independent_basis": independent_basis,
        "independent_minimal": independent_minimal,
        "minimal_common_result": minimal_common_result,
        "minimal_independent_result": minimal_independent_result,
    }


# ----------------------------------------------------------------------------
# Compact summary.
# ----------------------------------------------------------------------------

print(flush=True)
print("="*78, flush=True)
print("FIXED-GEOMETRY U-CONNECTION SUMMARY", flush=True)
print("  geometry                          constants      common-U      independent-U", flush=True)

for key, description, coefficients in GEOMETRIES_23:
    result = FIXED_RESULTS_23[key]
    constant_pair = result["constant"]["all_pair"]
    common_pair = result["common"]["all_pair"]
    independent_pair = result["independent"]["all_pair"]

    print("  {:<32s} {:>2d}/{:<2d}         {:>2d}/{:<2d}         {:>2d}/{:<2d}".format(
        key,
        constant_pair[0], constant_pair[1],
        common_pair[0], common_pair[1],
        independent_pair[0], independent_pair[1],
    ), flush=True)

print(flush=True)
print("Free-geometry benchmarks:", flush=True)
print("  common U only:      rank(A)={} rank([A|b])={} structural={}".format(
    FREE_COMMON_RESULT_23["all_pair"][0],
    FREE_COMMON_RESULT_23["all_pair"][1],
    FREE_COMMON_RESULT_23["structural"],
), flush=True)
print("  independent U:      rank(A)={} rank([A|b])={} structural={}".format(
    FREE_INDEPENDENT_RESULT_23["all_pair"][0],
    FREE_INDEPENDENT_RESULT_23["all_pair"][1],
    FREE_INDEPENDENT_RESULT_23["structural"],
), flush=True)

print(flush=True)
print("Interpretive diagnostics:", flush=True)
print("  pure U common connection works before locking geometry? {}".format(
    FREE_COMMON_RESULT_23["closes"]
), flush=True)
for key, _, _ in GEOMETRIES_23:
    result = FIXED_RESULTS_23[key]
    print("  {}: common-U closes? {} ; independent-U closes? {}".format(
        key,
        result["common"]["closes"],
        result["independent"]["closes"],
    ), flush=True)
    print("    minimal common supports: {}".format(
        [subset for subset, _ in result["common_minimal"]]
    ), flush=True)
    print("    minimal independent supports: {}".format(
        [subset for subset, _ in result["independent_minimal"]]
    ), flush=True)

print(flush=True)
print("="*78, flush=True)
print("SUCCESS", flush=True)
print("The U=3X+1 spectral connection was tested with free and fixed", flush=True)
print("degree-two geometry on 33 training plus 17 held-out radial points.", flush=True)
print("Copy the summary, every closing fixed geometry, and its minimal formula back.", flush=True)
print("="*78, flush=True)
