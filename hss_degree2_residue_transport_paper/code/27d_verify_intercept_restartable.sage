from sage.all import *
import gc, os, time

# ============================================================================
# RESTARTABLE SINGLE-COLUMN EXACT VERIFIER
# ============================================================================
# Edit TARGET_INDEX before running:
#   0 = intercept Rquad
#   1 = slope rA
#   2 = slope rB
#   3 = slope rC
#
# The verifier loads exact_degree2_part4_checkpoint.sobj, proves only the
# selected identity, and saves after EVERY exact evaluation.  Re-running the
# same file resumes at the first unfinished point.
# ============================================================================

TARGET_INDEX = 0
PART4_CHECKPOINT = "exact_degree2_part4_checkpoint.sobj"

print("="*78)
print("RESTARTABLE SINGLE-COLUMN EXACT DEGREE-TWO VERIFIER")
print("="*78)

if not os.path.exists(PART4_CHECKPOINT):
    raise RuntimeError("Missing {} in the current project directory".format(PART4_CHECKPOINT))

DATA = load(PART4_CHECKPOINT)

defect_names = list(DATA["defect_names"])
if TARGET_INDEX < 0 or TARGET_INDEX >= len(defect_names):
    raise ValueError("TARGET_INDEX must be 0,1,2,or 3")

TARGET_NAME = str(defect_names[TARGET_INDEX])
SAFE_NAME = ["intercept_Rquad", "slope_rA", "slope_rB", "slope_rC"][TARGET_INDEX]
PROGRESS_FILE = "exact_degree2_{}_progress.sobj".format(SAFE_NAME)
RESULT_FILE = "exact_degree2_{}_result.sobj".format(SAFE_NAME)
TEXT_FILE = "exact_degree2_{}_result.txt".format(SAFE_NAME)

sigma_num = DATA["sigma_num"]
sigma_den = DATA["sigma_den"]
D_num = DATA["D_num"]
D_den = DATA["D_den"]
defect_num_den = DATA["defect_num_den"]
solution_entries = DATA["solution_entries"]
universal_a1 = DATA["universal_a1"]
universal_a2 = DATA["universal_a2"]

RX = sigma_num.parent()
K = RX.base_ring()
X = RX.gen()
U_poly = RX(3*X + 1)


def is_zero_K(value):
    value = K(value)
    return value.numerator() == 0


def degree_or_minus_one(poly):
    poly = RX(poly)
    return -1 if poly == 0 else int(poly.degree())


def horner_eval(poly, point):
    poly = RX(poly)
    point = K(point)
    ans = K(0)
    # Avoid constructing a reversed copy of a potentially large list.
    coeffs = poly.list()
    for i in range(len(coeffs)-1, -1, -1):
        ans = ans*point + K(coeffs[i])
    return ans


def connection_polynomials(column_index):
    am3,bm3,am2,bm2,am1,aconst,bconst = [
        K(solution_entries[row][column_index]) for row in range(7)
    ]
    Atilde = RX(am3 + am2*U_poly + am1*U_poly**2 + aconst*U_poly**3)
    Btilde = RX(bm3 + bm2*U_poly + bconst*U_poly**3)
    if column_index == 0:
        Atilde += RX(K(universal_a1)*U_poly**4 + K(universal_a2)*U_poly**5)
    return Atilde, Btilde


def degree_bound(index, Atilde, Btilde):
    nF,dF = defect_num_den[index]
    terms = [
        degree_or_minus_one(nF) + degree_or_minus_one(sigma_den) + degree_or_minus_one(D_den) + 3,
        degree_or_minus_one(dF) + degree_or_minus_one(Atilde) + degree_or_minus_one(sigma_num) + degree_or_minus_one(D_den),
        degree_or_minus_one(dF) + degree_or_minus_one(Btilde) + degree_or_minus_one(D_num) + degree_or_minus_one(sigma_den),
    ]
    return max(terms), terms


def cleared_value(index, Atilde, Btilde, point):
    nF,dF = defect_num_den[index]
    x0 = K(point)
    u0 = K(3*x0 + 1)

    # Evaluate one object at a time to reduce peak memory.
    nF0 = horner_eval(nF, x0)
    dF0 = horner_eval(dF, x0)
    nS0 = horner_eval(sigma_num, x0)
    dS0 = horner_eval(sigma_den, x0)
    nD0 = horner_eval(D_num, x0)
    dD0 = horner_eval(D_den, x0)
    A0 = horner_eval(Atilde, x0)
    B0 = horner_eval(Btilde, x0)

    left = nF0*dS0*dD0*u0**3
    right = dF0*(A0*nS0*dD0 + B0*nD0*dS0)
    return K(left-right)


Atilde, Btilde = connection_polynomials(TARGET_INDEX)
bound, term_bounds = degree_bound(TARGET_INDEX, Atilde, Btilde)

# Integer points minimize coefficient growth and are sufficient because the
# cleared object is a polynomial of degree <= bound.
test_points = [QQ(i) for i in range(bound+1)]

print("  target: {}".format(TARGET_NAME), flush=True)
print("  term degree bounds: {}".format(term_bounds), flush=True)
print("  cleared numerator degree <= {}".format(bound), flush=True)
print("  exact integer points: X=0,...,{}".format(bound), flush=True)
print("  progress file: {}".format(PROGRESS_FILE), flush=True)

if os.path.exists(PROGRESS_FILE):
    progress = load(PROGRESS_FILE)
    if int(progress.get("target_index", -1)) != TARGET_INDEX or int(progress.get("degree_bound", -2)) != bound:
        raise RuntimeError("Existing progress file belongs to a different target or degree bound")
    completed = set([int(i) for i in progress.get("completed_indices", [])])
    witness = progress.get("witness", None)
    print("  resumed with {}/{} evaluations already saved".format(len(completed), len(test_points)), flush=True)
else:
    completed = set()
    witness = None
    progress = {
        "target_index": TARGET_INDEX,
        "target_name": TARGET_NAME,
        "degree_bound": bound,
        "term_degree_bounds": term_bounds,
        "points": test_points,
        "completed_indices": [],
        "witness": None,
        "certified": False,
    }
    save(progress, PROGRESS_FILE)

if witness is not None:
    print("  existing exact nonzero witness at X={}".format(witness["point"]), flush=True)
else:
    for j, point in enumerate(test_points):
        if j in completed:
            continue

        print("  evaluating {}/{} at X={} ...".format(j+1, len(test_points), point), flush=True)
        started = time.time()
        value = cleared_value(TARGET_INDEX, Atilde, Btilde, point)
        elapsed = time.time() - started

        if not is_zero_K(value):
            witness = {"index": j, "point": point, "value": value}
            progress["witness"] = witness
            progress["certified"] = False
            save(progress, PROGRESS_FILE)
            print("  NONZERO exact witness at X={}".format(point), flush=True)
            print("  elapsed seconds: {:.2f}".format(elapsed), flush=True)
            break

        completed.add(j)
        progress["completed_indices"] = sorted(completed)
        progress["last_elapsed_seconds"] = elapsed
        progress["certified"] = (len(completed) == len(test_points))
        save(progress, PROGRESS_FILE)
        print("  exact zero; saved ({}/{}) in {:.2f} s".format(len(completed), len(test_points), elapsed), flush=True)

        # Release large temporary fraction-field elements before the next point.
        del value
        gc.collect()

certified = (witness is None and len(completed) == len(test_points))
result = {
    "target_index": TARGET_INDEX,
    "target_name": TARGET_NAME,
    "degree_bound": bound,
    "term_degree_bounds": term_bounds,
    "verification_method": "exact cleared-polynomial degree bound plus D+1 distinct integer evaluations",
    "points": test_points,
    "completed_indices": sorted(completed),
    "witness": witness,
    "certified": certified,
}
save(result, RESULT_FILE)

with open(TEXT_FILE, "w") as handle:
    handle.write("{}\n".format(TARGET_NAME))
    handle.write("degree_bound={}\n".format(bound))
    handle.write("completed={}/{}\n".format(len(completed), len(test_points)))
    handle.write("certified={}\n".format(certified))
    if witness is not None:
        handle.write("witness_X={}\n".format(witness["point"]))

print()
print("="*78)
print("SINGLE-COLUMN VERIFICATION SUMMARY")
print("  target: {}".format(TARGET_NAME))
print("  completed: {}/{}".format(len(completed), len(test_points)))
print("  certified exactly? {}".format(certified))
print("  result file: {}".format(RESULT_FILE))
print("="*78)
