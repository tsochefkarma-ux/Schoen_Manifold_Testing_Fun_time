from sage.all import *
import os, gc, time

# ============================================================================
# FAST PART-V RESUME FOR THE EXACT DEGREE-TWO CERTIFICATE
# ============================================================================
#
# Use this in the SAME CoCalc Sage kernel after interrupting the slow Part V of
# 27_exact_degree2_defect_certificate.sage.  All objects through Part IV remain
# in memory after a normal Jupyter interrupt.
#
# Instead of asking Sage to normalize one enormous fraction-field expression,
# this script clears denominators structurally and proves the resulting
# polynomial is zero by exact evaluation at degree_bound+1 distinct rational
# X-values.  This is an exact polynomial identity certificate, not a numerical
# test.
# ============================================================================

print("="*78)
print("FAST EXACT PART-V RESUME")
print("Checkpointing Part IV and using a degree/evaluation identity certificate")
print("="*78)

required_names = [
    "K", "RX", "FX", "X", "U", "sigma1", "D1",
    "defect_functions", "defect_names", "solution_matrix",
    "universal_a1", "universal_a2", "rho3", "AFFINE_CONNECTION_MAP",
]
missing = [name for name in required_names if name not in globals()]
if missing:
    raise RuntimeError(
        "Run the exact certificate through Part IV in this same kernel first. "
        "Missing: {}".format(missing)
    )

PART4_CHECKPOINT = "exact_degree2_part4_checkpoint.sobj"
FINAL_CHECKPOINT = "exact_degree2_fast_certificate.sobj"
TEXT_SUMMARY = "exact_degree2_fast_certificate_summary.txt"


def is_zero_K_fast(value):
    value = K(value)
    return value.numerator() == 0


def num_den_RX(value):
    value = FX(value)
    return RX(value.numerator()), RX(value.denominator())


def degree_or_minus_one(poly):
    poly = RX(poly)
    return -1 if poly == 0 else int(poly.degree())


def horner_eval(poly, point):
    """Evaluate RX polynomial at a rational/K point with explicit Horner."""
    poly = RX(poly)
    point = K(point)
    result = K(0)
    for coefficient in reversed(poly.list()):
        result = result*point + K(coefficient)
    return result


# ---------------------------------------------------------------------------
# Save a minimal standard-Sage checkpoint BEFORE doing any further work.
# ---------------------------------------------------------------------------

sigma_num, sigma_den = num_den_RX(sigma1)
D_num, D_den = num_den_RX(D1)
defect_num_den = [num_den_RX(item) for item in defect_functions]
solution_entries = [
    [K(solution_matrix[row, col]) for col in range(solution_matrix.ncols())]
    for row in range(solution_matrix.nrows())
]

PART4_DATA = {
    "sigma_num": sigma_num,
    "sigma_den": sigma_den,
    "D_num": D_num,
    "D_den": D_den,
    "defect_num_den": defect_num_den,
    "defect_names": list(defect_names),
    "solution_entries": solution_entries,
    "universal_a1": K(universal_a1),
    "universal_a2": K(universal_a2),
    "rho3": K(rho3),
    "connection_map": AFFINE_CONNECTION_MAP,
}

save(PART4_DATA, PART4_CHECKPOINT)
print("  Part-IV checkpoint written: {}".format(PART4_CHECKPOINT), flush=True)

# Compact rho3 verification/printing.  Avoid printing the expanded coercion.
L_local, E4_local = K.gens()[0], K.gens()[1]
rho3_compact = K(2*(L_local^2*E4_local + 108)/(3*(L_local^2*E4_local + 36)))
print("  rho3 compact recognition exact? {}".format(
    is_zero_K_fast(K(rho3)-rho3_compact)
), flush=True)
print("  rho3 = 2*(L^2*E4+108)/(3*(L^2*E4+36))", flush=True)

# ---------------------------------------------------------------------------
# Structured cleared-numerator identity.
#
# Write
#   alpha = Atilde/U^3, beta = Btilde/U^3,
# so the rational identity is equivalent to the polynomial identity
#
#   N(X) = nF*dS*dD*U^3
#          - dF*(Atilde*nS*dD + Btilde*nD*dS) = 0.
#
# We never construct N(X) globally.  We compute an exact degree upper bound D
# and prove N=0 by checking D+1 distinct rational X-values in K.
# ---------------------------------------------------------------------------

U_poly = RX(3*X + 1)
U3 = U_poly^3


def connection_polynomials(column_index):
    am3 = K(solution_entries[0][column_index])
    bm3 = K(solution_entries[1][column_index])
    am2 = K(solution_entries[2][column_index])
    bm2 = K(solution_entries[3][column_index])
    am1 = K(solution_entries[4][column_index])
    aconst = K(solution_entries[5][column_index])
    bconst = K(solution_entries[6][column_index])

    # U^3 * alpha(U), U^3 * beta(U).
    Atilde = RX(am3 + am2*U_poly + am1*U_poly^2 + aconst*U_poly^3)
    Btilde = RX(bm3 + bm2*U_poly + bconst*U_poly^3)
    if column_index == 0:
        Atilde += RX(K(universal_a1)*U_poly^4 + K(universal_a2)*U_poly^5)
    return RX(Atilde), RX(Btilde)


def identity_degree_bound(index, Atilde, Btilde):
    nF, dF = defect_num_den[index]
    term_degrees = [
        degree_or_minus_one(nF) + degree_or_minus_one(sigma_den)
        + degree_or_minus_one(D_den) + 3,
        degree_or_minus_one(dF) + degree_or_minus_one(Atilde)
        + degree_or_minus_one(sigma_num) + degree_or_minus_one(D_den),
        degree_or_minus_one(dF) + degree_or_minus_one(Btilde)
        + degree_or_minus_one(D_num) + degree_or_minus_one(sigma_den),
    ]
    return max(term_degrees), term_degrees


def cleared_value(index, Atilde, Btilde, point):
    nF, dF = defect_num_den[index]
    x0 = K(point)
    u0 = K(3*x0 + 1)

    nF0 = horner_eval(nF, x0)
    dF0 = horner_eval(dF, x0)
    nS0 = horner_eval(sigma_num, x0)
    dS0 = horner_eval(sigma_den, x0)
    nD0 = horner_eval(D_num, x0)
    dD0 = horner_eval(D_den, x0)
    A0 = horner_eval(Atilde, x0)
    B0 = horner_eval(Btilde, x0)

    return K(
        nF0*dS0*dD0*u0^3
        - dF0*(A0*nS0*dD0 + B0*nD0*dS0)
    )


def rational_test_points(count):
    """Deterministic distinct QQ points; no pole avoidance is needed for N(X)."""
    points = []
    shell = 0
    while len(points) < count:
        candidates = []
        if shell == 0:
            candidates = [QQ(0)]
        else:
            candidates = [QQ(shell), QQ(-shell), (QQ(2*shell-1)/QQ(2)), (QQ(-(2*shell-1))/QQ(2))]
        for value in candidates:
            if value not in points:
                points.append(value)
                if len(points) >= count:
                    break
        shell += 1
    return points


results = []
witnesses = []

for index, name in enumerate(defect_names):
    print("-"*78, flush=True)
    print("  verifying {}".format(name), flush=True)
    Atilde, Btilde = connection_polynomials(index)
    bound, term_bounds = identity_degree_bound(index, Atilde, Btilde)
    point_count = bound + 1
    points = rational_test_points(point_count)
    print("    degree bounds by term: {}".format(term_bounds), flush=True)
    print("    cleared numerator degree <= {}".format(bound), flush=True)
    print("    exact evaluations required: {}".format(point_count), flush=True)

    passed = True
    witness = None
    started = time.time()
    for j, point in enumerate(points):
        value = cleared_value(index, Atilde, Btilde, point)
        if not is_zero_K_fast(value):
            passed = False
            witness = {"point": point, "value": value}
            print("    NONZERO witness at X={}".format(point), flush=True)
            break
        if (j + 1) % 5 == 0 or j + 1 == point_count:
            print("    passed {}/{} exact evaluations".format(j+1, point_count), flush=True)
            gc.collect()

    elapsed = time.time() - started
    print("    exact identity certified? {}".format(passed), flush=True)
    print("    elapsed seconds: {:.2f}".format(elapsed), flush=True)
    results.append(passed)
    witnesses.append(witness)

FINAL_DATA = dict(PART4_DATA)
FINAL_DATA.update({
    "verification_method": "exact degree bound plus D+1 rational evaluations",
    "identity_results": results,
    "witnesses": witnesses,
})
save(FINAL_DATA, FINAL_CHECKPOINT)

with open(TEXT_SUMMARY, "w") as handle:
    handle.write("EXACT DEGREE-TWO FAST CERTIFICATE\n")
    handle.write("Method: exact cleared-polynomial degree bound and D+1 evaluations\n")
    for name, result, witness in zip(defect_names, results, witnesses):
        handle.write("{}: {}\n".format(name, result))
        if witness is not None:
            handle.write("  witness X={}\n".format(witness["point"]))
    handle.write("rho3 = 2*(L^2*E4+108)/(3*(L^2*E4+36))\n")

print()
print("="*78)
print("FAST EXACT CERTIFICATE SUMMARY")
for name, result in zip(defect_names, results):
    print("  {:18s}: {}".format(name, result))
print("  all identities certified exactly? {}".format(all(results)))
print("  final checkpoint: {}".format(FINAL_CHECKPOINT))
print("  text summary: {}".format(TEXT_SUMMARY))
print("="*78)

if all(results):
    print("SUCCESS")
    print("All four rational identities are exact over the formal CM field.")
else:
    print("PARTIAL SUCCESS")
    print("A nonzero exact witness identifies the first failed affine column.")
