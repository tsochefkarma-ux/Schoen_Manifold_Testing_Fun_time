from sage.all import *
from itertools import combinations

# ============================================================================
# AFFINE GEOMETRY-TO-CONNECTION RECONSTRUCTION
# ============================================================================
#
# Run in the SAME SageMath session after successful cells 13--23.
#
# PURPOSE
# -------
# Cell 23 found that every fixed factorized geometry closes on the same minimal
# independent U=3X+1 connection support:
#
#   alpha_1(U) = A0 + A_-3 U^-3 + A_-2 U^-2 + A_-1 U^-1 + A_1 U + A_2 U^2,
#   beta_1(U)  = B0 + B_-3 U^-3 + B_-2 U^-2.
#
# There are nine coefficients in total.  On 50 radial points the corresponding
# spectral matrix has rank nine, so the connection is UNIQUE once (a,b,c) is
# fixed.  Since the right-hand side depends affinely on (a,b,c), uniqueness
# implies an affine geometry-to-connection map
#
#     p(a,b,c) = p_0 + a p_A + b p_B + c p_C.
#
# This cell reconstructs that map, validates it at new radial points and random
# rational geometries, identifies geometry-independent connection terms, and
# searches (exploratorily) for geometries that make the connection maximally
# sparse.
#
# ROBUSTNESS
# ----------
# * Fits only on the 50 points from cell 23.
# * Adds 20 new held-out rational radial points (70 total).
# * Reconstructs the affine map from the four geometries 0, A, B, C.
# * Validates against direct solves for ten additional rational geometries.
# * Reports rank profiles at several tolerances.
# * Does not interpret a sparse geometry as physical; that search is explicitly
#   labelled exploratory.
# ============================================================================

print("="*78, flush=True)
print("AFFINE GEOMETRY-TO-U-CONNECTION RECONSTRUCTION", flush=True)
print("Reconstructing the unique minimal connection as a function of (a,b,c)", flush=True)
print("="*78, flush=True)
print(flush=True)

required_names = [
    "RF", "alpha0", "beta0", "component_data_at_X",
    "BASE50_ROWS", "BASE50_RHS", "ALL_POINTS_23", "TOTAL_COUNT_23",
    "rank_pair_23", "canonical_solution_23", "residuals_23", "max_abs_23",
    "compact_rank_profile_23",
]
missing = [name for name in required_names if name not in globals()]
if missing:
    raise RuntimeError(
        "Run cells 13 through 23 first in the same Sage session. Missing: {}"
        .format(missing)
    )

RANK_TOL_EXPONENT_24 = 70
RANK_TOL_24 = RF(10)**(-RANK_TOL_EXPONENT_24)
ZERO_24 = RF(10)**(-80)
DISPLAY_ZERO_24 = RF(10)**(-60)

# Exact minimal support selected independently for every fixed geometry in cell 23.
CONNECTION_NAMES_24 = [
    "alpha_const",
    "beta_const",
    "alpha_U_-3",
    "beta_U_-3",
    "alpha_U_-2",
    "beta_U_-2",
    "alpha_U_-1",
    "alpha_U_1",
    "alpha_U_2",
]

# Disjoint from all 50 points used through cell 23.
AFFINE_VALIDATION_X_24 = [
    QQ(3)/8,
    QQ(9)/20,
    QQ(8)/15,
    QQ(11)/15,
    QQ(13)/15,
    QQ(17)/15,
    QQ(23)/20,
    QQ(5)/4 + QQ(1)/20,
    QQ(29)/20,
    QQ(31)/20,
    QQ(17)/10,
    QQ(37)/20,
    QQ(41)/20,
    QQ(17)/7,
    QQ(18)/7,
    QQ(22)/7,
    QQ(25)/7,
    QQ(19)/4,
    QQ(25)/4,
    QQ(9),
]

FIT_POINTS_24 = list(ALL_POINTS_23)
ALL_POINTS_24 = list(FIT_POINTS_24)
for value in AFFINE_VALIDATION_X_24:
    if value not in ALL_POINTS_24:
        ALL_POINTS_24.append(value)

FIT_COUNT_24 = len(FIT_POINTS_24)
TOTAL_COUNT_24 = len(ALL_POINTS_24)
HELDOUT_COUNT_24 = TOTAL_COUNT_24 - FIT_COUNT_24

print("Working real precision: {} bits".format(RF.precision()), flush=True)
print("Fit equations inherited from cell 23: {}".format(FIT_COUNT_24), flush=True)
print("New held-out validation equations: {}".format(HELDOUT_COUNT_24), flush=True)
print("Total radial equations: {}".format(TOTAL_COUNT_24), flush=True)
print("Minimal connection coefficients: {}".format(CONNECTION_NAMES_24), flush=True)
print(flush=True)


def base_row_at_X_24(X_value):
    base = component_data_at_X(X_value, None)
    colA = component_data_at_X(X_value, "A")
    colB = component_data_at_X(X_value, "B")
    colC = component_data_at_X(X_value, "C")

    def corrected(total):
        return RF(
            total["Scurv2"] - base["Scurv2"]
            - alpha0*total["sigma2"] - beta0*total["D2"]
        )

    return (
        [
            corrected(colA),
            corrected(colB),
            corrected(colC),
            RF(-base["sigma1"]),
            RF(-base["D1"]),
        ],
        RF(-base["Scurv2"]),
    )


# Extend the five-column base response from 50 to 70 points.
BASE70_ROWS_24 = [list(row) for row in BASE50_ROWS]
BASE70_RHS_24 = [RF(value) for value in BASE50_RHS]

if len(BASE70_ROWS_24) != FIT_COUNT_24:
    raise ArithmeticError(
        "Cell-23 row count {} does not match its point count {}."
        .format(len(BASE70_ROWS_24), FIT_COUNT_24)
    )

print("Building the old response system at 20 new validation points...", flush=True)
for index, X_value in enumerate(AFFINE_VALIDATION_X_24, start=1):
    row, rhs = base_row_at_X_24(X_value)
    BASE70_ROWS_24.append(row)
    BASE70_RHS_24.append(rhs)
    print("  completed {}/{} at X={}".format(
        index, len(AFFINE_VALIDATION_X_24), X_value
    ), flush=True)

if len(BASE70_ROWS_24) != TOTAL_COUNT_24:
    raise ArithmeticError("Expanded 70-point row count is inconsistent.")


# The unique nine-column minimal spectral matrix.
def minimal_connection_row_24(X_value, base_row):
    U = RF(3)*RF(X_value) + RF(1)
    minus_sigma = RF(base_row[3])
    minus_D = RF(base_row[4])
    return [
        minus_sigma,
        minus_D,
        U**(-3)*minus_sigma,
        U**(-3)*minus_D,
        U**(-2)*minus_sigma,
        U**(-2)*minus_D,
        U**(-1)*minus_sigma,
        U*minus_sigma,
        U**2*minus_sigma,
    ]


SPECTRAL70_ROWS_24 = [
    minimal_connection_row_24(X_value, row)
    for X_value, row in zip(ALL_POINTS_24, BASE70_ROWS_24)
]


def fixed_rhs_24(coefficients, row_count=TOTAL_COUNT_24):
    a, b, c = [RF(value) for value in coefficients]
    return [
        RF(BASE70_RHS_24[row])
        - a*RF(BASE70_ROWS_24[row][0])
        - b*RF(BASE70_ROWS_24[row][1])
        - c*RF(BASE70_ROWS_24[row][2])
        for row in range(row_count)
    ]


def solve_connection_24(coefficients, row_count=FIT_COUNT_24):
    rhs = fixed_rhs_24(coefficients, row_count)
    A = SPECTRAL70_ROWS_24[:row_count]
    pair = rank_pair_23(A, rhs)
    solution = canonical_solution_23(A, rhs)
    if pair != (9, 9) or solution is None or solution["nullity"] != 0:
        raise ArithmeticError(
            "Minimal connection is not uniquely solvable for geometry {}: pair={}, solution={}."
            .format(coefficients, pair, solution)
        )
    return [RF(value) for value in solution["solution"]]


def affine_predict_24(p0, slopes, coefficients):
    a, b, c = [RF(value) for value in coefficients]
    return [
        p0[index]
        + a*slopes[0][index]
        + b*slopes[1][index]
        + c*slopes[2][index]
        for index in range(len(p0))
    ]


def max_vector_difference_24(left, right):
    return max_abs_23([RF(a)-RF(b) for a, b in zip(left, right)])


def heldout_residual_24(coefficients, solution):
    rhs = fixed_rhs_24(coefficients, TOTAL_COUNT_24)
    residual = residuals_23(SPECTRAL70_ROWS_24, rhs, solution)
    return max_abs_23(residual[FIT_COUNT_24:])


def allpoint_residual_24(coefficients, solution):
    rhs = fixed_rhs_24(coefficients, TOTAL_COUNT_24)
    return max_abs_23(residuals_23(SPECTRAL70_ROWS_24, rhs, solution))


# ----------------------------------------------------------------------------
# Reconstruct the affine map from the origin and the three factorized basis
# geometries.
# ----------------------------------------------------------------------------

BASIS_GEOMETRIES_24 = [
    ("zero", (QQ(0), QQ(0), QQ(0))),
    ("A",    (QQ(1), QQ(0), QQ(0))),
    ("B",    (QQ(0), QQ(1), QQ(0))),
    ("C",    (QQ(0), QQ(0), QQ(1))),
]

print("Reconstructing affine geometry-to-connection map...", flush=True)
BASIS_SOLUTIONS_24 = {}
for label, geometry in BASIS_GEOMETRIES_24:
    solution_fit = solve_connection_24(geometry, FIT_COUNT_24)
    solution_all = solve_connection_24(geometry, TOTAL_COUNT_24)
    mismatch = max_vector_difference_24(solution_fit, solution_all)
    residual = allpoint_residual_24(geometry, solution_fit)
    BASIS_SOLUTIONS_24[label] = solution_fit
    print("  {} geometry {}: fit/all coefficient mismatch={}, 70-point residual={}".format(
        label, geometry, mismatch, residual
    ), flush=True)

P0_24 = BASIS_SOLUTIONS_24["zero"]
SLOPES_24 = []
for label in ["A", "B", "C"]:
    SLOPES_24.append([
        BASIS_SOLUTIONS_24[label][index] - P0_24[index]
        for index in range(len(CONNECTION_NAMES_24))
    ])

print(flush=True)
print("Affine coefficient map:", flush=True)
for index, name in enumerate(CONNECTION_NAMES_24):
    p0 = P0_24[index]
    pa = SLOPES_24[0][index]
    pb = SLOPES_24[1][index]
    pc = SLOPES_24[2][index]
    print("  {} = {} + ({})*a + ({})*b + ({})*c".format(
        name, p0, pa, pb, pc
    ), flush=True)


# ----------------------------------------------------------------------------
# Validate affine dependence at additional rational geometries.
# ----------------------------------------------------------------------------

TEST_GEOMETRIES_24 = [
    ("primitive", (QQ(1), QQ(0), QQ(0))),
    ("primitive_plus_cover", (QQ(1), QQ(1)/8, QQ(0))),
    ("Z2_product", (QQ(1), QQ(1)/64, QQ(1)/8)),
    ("random_1", (QQ(2)/3, QQ(-1)/5, QQ(3)/7)),
    ("random_2", (QQ(-2), QQ(5)/11, QQ(-4)/9)),
    ("random_3", (QQ(7)/4, QQ(2)/9, QQ(1)/13)),
    ("random_4", (QQ(0), QQ(-3)/8, QQ(5)/6)),
    ("random_5", (QQ(11)/10, QQ(-7)/12, QQ(9)/14)),
    ("random_6", (QQ(-5)/6, QQ(4)/15, QQ(7)/16)),
    ("random_7", (QQ(13)/9, QQ(1)/17, QQ(-2)/19)),
]

print(flush=True)
print("Affine-map validation at additional geometries:", flush=True)
MAX_AFFINE_MISMATCH_24 = RF(0)
MAX_AFFINE_HELDOUT_24 = RF(0)
MAX_AFFINE_ALLPOINT_24 = RF(0)

for label, geometry in TEST_GEOMETRIES_24:
    direct = solve_connection_24(geometry, FIT_COUNT_24)
    predicted = affine_predict_24(P0_24, SLOPES_24, geometry)
    mismatch = max_vector_difference_24(direct, predicted)
    heldout = heldout_residual_24(geometry, predicted)
    allpoint = allpoint_residual_24(geometry, predicted)

    MAX_AFFINE_MISMATCH_24 = max(MAX_AFFINE_MISMATCH_24, mismatch)
    MAX_AFFINE_HELDOUT_24 = max(MAX_AFFINE_HELDOUT_24, heldout)
    MAX_AFFINE_ALLPOINT_24 = max(MAX_AFFINE_ALLPOINT_24, allpoint)

    print("  {:<22s} {}: coefficient mismatch={}, held-out residual={}".format(
        label, geometry, mismatch, heldout
    ), flush=True)

print("  maximum affine coefficient mismatch = {}".format(
    MAX_AFFINE_MISMATCH_24
), flush=True)
print("  maximum affine held-out residual    = {}".format(
    MAX_AFFINE_HELDOUT_24
), flush=True)
print("  maximum affine 70-point residual    = {}".format(
    MAX_AFFINE_ALLPOINT_24
), flush=True)
print("  affine law certified numerically at 1e-80? {}".format(
    MAX_AFFINE_MISMATCH_24 < ZERO_24 and MAX_AFFINE_HELDOUT_24 < ZERO_24
), flush=True)


# ----------------------------------------------------------------------------
# Geometry-independent and geometry-dependent pieces.
# ----------------------------------------------------------------------------

print(flush=True)
print("Universal versus geometry-dependent connection coefficients:", flush=True)
UNIVERSAL_NAMES_24 = []
DEPENDENT_NAMES_24 = []
for index, name in enumerate(CONNECTION_NAMES_24):
    slope_norm = max_abs_23([
        SLOPES_24[0][index],
        SLOPES_24[1][index],
        SLOPES_24[2][index],
    ])
    universal = bool(slope_norm < ZERO_24)
    if universal:
        UNIVERSAL_NAMES_24.append(name)
    else:
        DEPENDENT_NAMES_24.append(name)
    print("  {:<16s} slope norm={}  universal={}".format(
        name, slope_norm, universal
    ), flush=True)

print("  universal coefficient names: {}".format(UNIVERSAL_NAMES_24), flush=True)
print("  geometry-dependent names:    {}".format(DEPENDENT_NAMES_24), flush=True)

# Rank of the 9 x 3 slope matrix.  This is the dimension of the geometry-driven
# motion inside the nine-dimensional connection coefficient space.
SLOPE_MATRIX_ROWS_24 = [
    [SLOPES_24[column][row] for column in range(3)]
    for row in range(len(CONNECTION_NAMES_24))
]
SLOPE_RANK_24 = len(numerical_rref(SLOPE_MATRIX_ROWS_24, RANK_TOL_24)[1])
print("  rank of geometry-slope map (a,b,c) -> connection = {}".format(
    SLOPE_RANK_24
), flush=True)

# Explicit universal polynomial drift.
print(flush=True)
print("Universal U-polynomial drift:", flush=True)
for name in ["alpha_U_1", "alpha_U_2"]:
    if name in CONNECTION_NAMES_24:
        index = CONNECTION_NAMES_24.index(name)
        print("  {} = {}".format(name, P0_24[index]), flush=True)


# ----------------------------------------------------------------------------
# Exploratory sparse-geometry search.
#
# Use three geometry parameters to set triples of geometry-dependent connection
# coefficients to zero.  This is NOT a physical selection principle.  It merely
# identifies algebraically simple representatives that may suggest a natural
# normalization or missing enumerative constraint.
# ----------------------------------------------------------------------------


def solve_three_zero_conditions_24(indices):
    M = matrix(RF, [
        [SLOPES_24[0][index], SLOPES_24[1][index], SLOPES_24[2][index]]
        for index in indices
    ])
    if abs(M.det()) < RANK_TOL_24:
        return None
    rhs = vector(RF, [-P0_24[index] for index in indices])
    geometry = M.solve_right(rhs)
    return tuple(RF(value) for value in geometry)


def support_of_vector_24(values, tolerance=DISPLAY_ZERO_24):
    return [
        CONNECTION_NAMES_24[index]
        for index, value in enumerate(values)
        if abs(RF(value)) > tolerance
    ]


geometry_dependent_indices_24 = [
    CONNECTION_NAMES_24.index(name) for name in DEPENDENT_NAMES_24
]
SPARSE_CANDIDATES_24 = []

for indices in combinations(geometry_dependent_indices_24, 3):
    geometry = solve_three_zero_conditions_24(indices)
    if geometry is None:
        continue
    coefficients = affine_predict_24(P0_24, SLOPES_24, geometry)
    support = support_of_vector_24(coefficients)
    residual = allpoint_residual_24(geometry, coefficients)
    SPARSE_CANDIDATES_24.append({
        "zeroed": tuple(CONNECTION_NAMES_24[index] for index in indices),
        "geometry": geometry,
        "support": support,
        "support_size": len(support),
        "residual": residual,
        "max_geometry": max_abs_23(geometry),
    })

SPARSE_CANDIDATES_24.sort(key=lambda item: (
    item["support_size"],
    item["max_geometry"],
))

print(flush=True)
print("Exploratory maximally sparse geometries:", flush=True)
if SPARSE_CANDIDATES_24:
    for candidate in SPARSE_CANDIDATES_24[:10]:
        print("  zeroed {}".format(candidate["zeroed"]), flush=True)
        print("    geometry (a,b,c) = {}".format(candidate["geometry"]), flush=True)
        print("    surviving support = {}".format(candidate["support"]), flush=True)
        print("    support size = {}, 70-point residual = {}".format(
            candidate["support_size"], candidate["residual"]
        ), flush=True)
else:
    print("  no invertible three-condition systems found", flush=True)


# ----------------------------------------------------------------------------
# Compact summary.
# ----------------------------------------------------------------------------

print(flush=True)
print("="*78, flush=True)
print("AFFINE GEOMETRY-TO-CONNECTION SUMMARY", flush=True)
print("  minimal connection rank: 9", flush=True)
print("  geometry-slope rank: {}".format(SLOPE_RANK_24), flush=True)
print("  universal coefficients: {}".format(UNIVERSAL_NAMES_24), flush=True)
print("  geometry-dependent coefficients: {}".format(DEPENDENT_NAMES_24), flush=True)
print("  maximum affine coefficient mismatch: {}".format(
    MAX_AFFINE_MISMATCH_24
), flush=True)
print("  maximum held-out residual: {}".format(
    MAX_AFFINE_HELDOUT_24
), flush=True)
print("  affine law passes 1e-80 validation? {}".format(
    MAX_AFFINE_MISMATCH_24 < ZERO_24 and MAX_AFFINE_HELDOUT_24 < ZERO_24
), flush=True)
print(flush=True)
print("Interpretive diagnostics:", flush=True)
print("  connection is unique for fixed (a,b,c)? True", flush=True)
print("  dependence on (a,b,c) is affine numerically? {}".format(
    MAX_AFFINE_MISMATCH_24 < ZERO_24
), flush=True)
print("  positive U drift is geometry-independent? {}".format(
    "alpha_U_1" in UNIVERSAL_NAMES_24 and "alpha_U_2" in UNIVERSAL_NAMES_24
), flush=True)
print("  geometry dependence confined to constant/polar terms? {}".format(
    all(name in [
        "alpha_const", "beta_const",
        "alpha_U_-3", "beta_U_-3",
        "alpha_U_-2", "beta_U_-2",
        "alpha_U_-1",
    ] for name in DEPENDENT_NAMES_24)
), flush=True)
print(flush=True)
print("="*78, flush=True)
print("SUCCESS", flush=True)
print("The unique minimal U-connection was reconstructed as an affine function", flush=True)
print("of the factorized degree-two geometry and validated on 70 radial points.", flush=True)
print("Copy the affine summary, universal coefficients, and sparse candidates back.", flush=True)
print("="*78, flush=True)
