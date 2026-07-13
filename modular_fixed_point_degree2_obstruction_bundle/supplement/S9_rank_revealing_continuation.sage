# ============================================================================
# RANK-REVEALING CONTINUATION FOR THE DEGREE-TWO SCHOEN GLUING TEST
# ============================================================================
#
# Run this in the SAME Sage session after the previous cell stopped at the
# tiny-determinant exception.  The expensive DATA dictionary has already been
# constructed, so this continuation does not rebuild the geometry.
#
# It replaces the unreliable determinant test by:
#   1. column and row equilibration;
#   2. high-precision pivoted RREF;
#   3. rank(A) versus rank([A|b]) over all 13 sample points;
#   4. extraction of hidden column relations;
#   5. validation of any consistent solution on every sample point;
#   6. separate tests of the three named gluing prescriptions.
#
# A determinant near 10^(-precision_digits) is usually evidence of a rank
# defect, not a reason to increase SERIES_PREC further.
# ============================================================================

print("="*78, flush=True)
print("RANK-REVEALING DEGREE-TWO GLUING SOLVE", flush=True)
print("Using the already-computed full p^2 geometry", flush=True)
print("="*78, flush=True)
print(flush=True)

required = [
    "DATA", "TRAIN_X", "VALIDATE_X", "RF",
]
missing = [name for name in required if name not in globals()]
if missing:
    raise RuntimeError(
        "Run the degree-two geometry cell first. Missing: {}".format(missing)
    )

# Keep the point ordering deterministic.
POINTS = []
for value in list(TRAIN_X) + list(VALIDATE_X):
    if value not in POINTS:
        POINTS.append(value)

N_UNKNOWNS = 5
UNKNOWN_NAMES = ["a", "b", "c", "alpha_1", "beta_1"]

# A*u=b, with u=(a,b,c,alpha_1,beta_1).
A_rows = []
b_values = []
for Xv in POINTS:
    row = DATA[Xv]
    A_rows.append([
        RF(row["rA"]),
        RF(row["rB"]),
        RF(row["rC"]),
        RF(-row["sigma1"]),
        RF(-row["D1"]),
    ])
    b_values.append(RF(-row["Rquad"]))

print("Equations:", len(A_rows), flush=True)
print("Unknowns :", N_UNKNOWNS, flush=True)

# Decimal precision estimate for diagnostics and tolerance selection.
decimal_digits = max(20, int(RF.precision()*0.3010299956639812))
print("Working precision: about {} decimal digits".format(decimal_digits), flush=True)
print(flush=True)


def max_abs(values):
    values = list(values)
    if not values:
        return RF(0)
    return max(abs(RF(value)) for value in values)


# ---------------------------------------------------------------------------
# Equilibrate A and b without changing ranks or consistency.
# ---------------------------------------------------------------------------

b_global = max_abs(b_values)
if b_global == 0:
    b_global = RF(1)

# First divide the entire system by a common RHS scale.
A0 = [[value/b_global for value in row] for row in A_rows]
b0 = [value/b_global for value in b_values]

# Then scale each A column to maximum absolute entry 1.
column_scales = []
for column in range(N_UNKNOWNS):
    scale = max_abs(A0[row][column] for row in range(len(A0)))
    if scale == 0:
        scale = RF(1)
    column_scales.append(scale)

A_scaled = [
    [A0[row][column]/column_scales[column]
     for column in range(N_UNKNOWNS)]
    for row in range(len(A0))
]

# Row scaling is applied equally to A and b.
A_equilibrated = []
b_equilibrated = []
for row, rhs in zip(A_scaled, b0):
    scale = max(RF(1), max_abs(row), abs(rhs))
    A_equilibrated.append([value/scale for value in row])
    b_equilibrated.append(rhs/scale)

print("Column scales after common RHS normalization:", flush=True)
for name, scale in zip(UNKNOWN_NAMES, column_scales):
    print("  {:>7s}: {}".format(name, scale), flush=True)
print(flush=True)


# ---------------------------------------------------------------------------
# Small pivoted numerical RREF, designed for high-precision RealField input.
# ---------------------------------------------------------------------------


def numerical_rref(rows, tolerance):
    """Return (RREF rows, pivot columns) using maximum-magnitude pivots."""
    M = [[RF(value) for value in row] for row in rows]
    if not M:
        return M, []

    m = len(M)
    n = len(M[0])
    pivot_row = 0
    pivots = []

    for column in range(n):
        if pivot_row >= m:
            break

        best_row = max(
            range(pivot_row, m),
            key=lambda index: abs(M[index][column]),
        )
        best_value = abs(M[best_row][column])

        if best_value <= tolerance:
            continue

        if best_row != pivot_row:
            M[pivot_row], M[best_row] = M[best_row], M[pivot_row]

        pivot = M[pivot_row][column]
        M[pivot_row] = [value/pivot for value in M[pivot_row]]

        for row in range(m):
            if row == pivot_row:
                continue
            factor_value = M[row][column]
            if abs(factor_value) <= tolerance:
                M[row][column] = RF(0)
                continue
            M[row] = [
                M[row][col] - factor_value*M[pivot_row][col]
                for col in range(n)
            ]

        # Clean only entries safely below the requested tolerance.
        for row in range(m):
            for col in range(n):
                if abs(M[row][col]) <= tolerance:
                    M[row][col] = RF(0)

        pivots.append(column)
        pivot_row += 1

    return M, pivots


def rank_at_tolerance(rows, tolerance):
    return len(numerical_rref(rows, tolerance)[1])


augmented_equilibrated = [
    A_equilibrated[row] + [b_equilibrated[row]]
    for row in range(len(A_equilibrated))
]

# Show whether the rank conclusion is stable over a broad tolerance range.
trial_exponents = sorted(set(
    exponent for exponent in [
        20, 30, 40, 60, 80, 100,
        decimal_digits//3,
        decimal_digits//2,
        (2*decimal_digits)//3,
        decimal_digits-40,
    ]
    if 10 <= exponent <= decimal_digits-10
))

print("Rank profile after equilibration:", flush=True)
print("  tolerance          rank(A)   rank([A|b])", flush=True)
rank_profile = []
for exponent in trial_exponents:
    tolerance = RF(10)**(-exponent)
    rank_A = rank_at_tolerance(A_equilibrated, tolerance)
    rank_aug = rank_at_tolerance(augmented_equilibrated, tolerance)
    rank_profile.append((exponent, rank_A, rank_aug))
    print("  1e-{:<4d}            {:d}          {:d}".format(
        exponent, rank_A, rank_aug
    ), flush=True)
print(flush=True)

# Use a conservative tolerance well above the arithmetic noise floor.
chosen_exponent = min(80, max(30, decimal_digits//4))
TOL = RF(10)**(-chosen_exponent)

R_A, pivots_A = numerical_rref(A_equilibrated, TOL)
R_aug, pivots_aug = numerical_rref(augmented_equilibrated, TOL)
rank_A = len(pivots_A)
rank_aug = len(pivots_aug)

print("Chosen structural tolerance: 1e-{}".format(chosen_exponent), flush=True)
print("  rank(A)     =", rank_A, flush=True)
print("  rank([A|b]) =", rank_aug, flush=True)
print("  pivot columns of A =", pivots_A, flush=True)
print(flush=True)


# ---------------------------------------------------------------------------
# Nullspace of the equilibrated coefficient matrix.
# ---------------------------------------------------------------------------


def nullspace_from_rref(rref_rows, pivots, n_columns):
    free_columns = [column for column in range(n_columns) if column not in pivots]
    basis = []
    pivot_to_row = {column: row for row, column in enumerate(pivots)}

    for free in free_columns:
        vector_value = [RF(0)]*n_columns
        vector_value[free] = RF(1)
        for pivot_column in pivots:
            row = pivot_to_row[pivot_column]
            vector_value[pivot_column] = -rref_rows[row][free]
        basis.append(vector_value)

    return basis


scaled_null_basis = nullspace_from_rref(R_A, pivots_A, N_UNKNOWNS)

if scaled_null_basis:
    print("Hidden column relation(s):", flush=True)
    for index, scaled_vector in enumerate(scaled_null_basis, start=1):
        # Scaled unknown v_j = column_scale_j * original unknown u_j.
        original_vector = [
            scaled_vector[column]/column_scales[column]
            for column in range(N_UNKNOWNS)
        ]
        normalization = max_abs(original_vector)
        if normalization != 0:
            original_vector = [value/normalization for value in original_vector]
        print("  relation {} in (a,b,c,alpha_1,beta_1):".format(index), flush=True)
        print("   ", [value for value in original_vector], flush=True)
    print(flush=True)
else:
    print("No structural column relation detected at the chosen tolerance.", flush=True)
    print(flush=True)


# ---------------------------------------------------------------------------
# Solve a consistent system, or certify incompatibility.
# ---------------------------------------------------------------------------


def residual_for_values(values, Xv):
    aa, bb, cc, al1, be1 = values
    row = DATA[Xv]
    return (
        row["Rquad"]
        + aa*row["rA"]
        + bb*row["rB"]
        + cc*row["rC"]
        - al1*row["sigma1"]
        - be1*row["D1"]
    )


span_preserves = False
solution = None

if rank_aug > rank_A:
    print("RESULT: the overdetermined system is inconsistent.", flush=True)
    print("No combination in the natural factorized span restores", flush=True)
    print("two-channel closure at order p^2.", flush=True)
    print(flush=True)
else:
    # The augmented RREF is consistent.  Set all free variables to zero to
    # obtain one representative solution in scaled variables.
    pivot_variables = [column for column in pivots_aug if column < N_UNKNOWNS]
    free_variables = [
        column for column in range(N_UNKNOWNS)
        if column not in pivot_variables
    ]

    scaled_solution = [RF(0)]*N_UNKNOWNS
    pivot_to_row = {
        column: row for row, column in enumerate(pivots_aug)
        if column < N_UNKNOWNS
    }

    for pivot_column in pivot_variables:
        row = pivot_to_row[pivot_column]
        scaled_solution[pivot_column] = R_aug[row][-1]

    # v_j = column_scale_j*u_j.
    solution = vector(RF, [
        scaled_solution[column]/column_scales[column]
        for column in range(N_UNKNOWNS)
    ])

    residuals = [residual_for_values(solution, Xv) for Xv in POINTS]
    max_residual = max_abs(residuals)

    # A dimensionless residual relative to the equation scale.
    equation_scale = max(
        RF(1),
        max_abs(DATA[Xv]["Rquad"] for Xv in POINTS),
    )
    relative_residual = max_residual/equation_scale

    print("RESULT: the system is consistent.", flush=True)
    print("  solution-family dimension =", N_UNKNOWNS-rank_A, flush=True)
    print("  representative solution:", flush=True)
    for name, value in zip(UNKNOWN_NAMES, solution):
        print("    {:>7s} = {}".format(name, value), flush=True)
    print("  max residual over all points =", max_residual, flush=True)
    print("  relative residual            =", relative_residual, flush=True)

    span_preserves = bool(relative_residual < RF(10)**(-chosen_exponent//2))
    print("  natural factorized span restores closure?", span_preserves, flush=True)
    print(flush=True)


# ---------------------------------------------------------------------------
# Named candidates, using a two-point alpha_1,beta_1 fit and all-point check.
# ---------------------------------------------------------------------------


def solve_two_by_two(a11, a12, a21, a22, y1, y2):
    determinant = a11*a22 - a12*a21
    scale = max(RF(1), abs(a11*a22), abs(a12*a21))
    if abs(determinant) <= TOL*scale:
        raise ArithmeticError("The alpha/beta interpolation system is singular.")
    return (
        (y1*a22-a12*y2)/determinant,
        (a11*y2-y1*a21)/determinant,
    )


def test_named_candidate(label, abc):
    aa, bb, cc = [RF(value) for value in abc]
    XA = QQ(5)/4
    XB = QQ(3)/2

    def source(Xv):
        row = DATA[Xv]
        return (
            row["Rquad"]
            + aa*row["rA"]
            + bb*row["rB"]
            + cc*row["rC"]
        )

    alpha_shift, beta_shift = solve_two_by_two(
        DATA[XA]["sigma1"], DATA[XA]["D1"],
        DATA[XB]["sigma1"], DATA[XB]["D1"],
        source(XA), source(XB),
    )

    values = vector(RF, [aa, bb, cc, alpha_shift, beta_shift])
    residuals = [residual_for_values(values, Xv) for Xv in POINTS]
    max_residual = max_abs(residuals)
    source_scale = max(RF(1), max_abs(source(Xv) for Xv in POINTS))
    relative_residual = max_residual/source_scale
    preserves = bool(relative_residual < RF(10)**(-chosen_exponent//2))

    print("-"*78, flush=True)
    print(label, flush=True)
    print("  (a,b,c) =", (aa, bb, cc), flush=True)
    print("  fitted alpha_1 =", alpha_shift, flush=True)
    print("  fitted beta_1  =", beta_shift, flush=True)
    print("  relative all-point residual =", relative_residual, flush=True)
    print("  preserves closure?", preserves, flush=True)

    return preserves


candidate_primitive = test_named_candidate(
    "Candidate 1: primitive product F_A",
    (1, 0, 0),
)

candidate_gv = test_named_candidate(
    "Candidate 2: F_A + F_B/8",
    (1, QQ(1)/8, 0),
)

candidate_uninserted = test_named_candidate(
    "Candidate 3: Z2(y)Z2(z)",
    (1, QQ(1)/64, QQ(1)/8),
)

print(flush=True)
print("="*78, flush=True)
print("SUCCESS", flush=True)
print("The determinant test was replaced by a stable rank/consistency test.", flush=True)
print("rank(A) = {}, rank([A|b]) = {}".format(rank_A, rank_aug), flush=True)
print("natural factorized span restores closure?", span_preserves, flush=True)
print("="*78, flush=True)
