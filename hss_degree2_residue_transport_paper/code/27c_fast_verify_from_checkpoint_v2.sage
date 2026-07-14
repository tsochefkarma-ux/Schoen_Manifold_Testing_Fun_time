from sage.all import *
import gc, time

# ============================================================================
# EXACT DEGREE-TWO FAST VERIFIER FROM A SAVED PART-IV CHECKPOINT
# ============================================================================
# Run after 27_exact_degree2_cocalc_checkpointed.sage has written
# exact_degree2_part4_checkpoint.sobj.  This file does not rebuild the geometry.
# ============================================================================

PART4_CHECKPOINT = "exact_degree2_part4_checkpoint.sobj"
FINAL_CHECKPOINT = "exact_degree2_fast_certificate.sobj"
TEXT_SUMMARY = "exact_degree2_fast_certificate_summary.txt"

print("="*78)
print("EXACT DEGREE-TWO FAST VERIFIER FROM CHECKPOINT")
print("="*78)

DATA = load(PART4_CHECKPOINT)

sigma_num = DATA["sigma_num"]
sigma_den = DATA["sigma_den"]
D_num = DATA["D_num"]
D_den = DATA["D_den"]
defect_num_den = DATA["defect_num_den"]
defect_names = DATA["defect_names"]
solution_entries = DATA["solution_entries"]
universal_a1 = DATA["universal_a1"]
universal_a2 = DATA["universal_a2"]
rho3 = DATA["rho3"]

RX = sigma_num.parent()
K = RX.base_ring()
X = RX.gen()
U_poly = RX(3*X+1)


def is_zero_K(value):
    value = K(value)
    return value.numerator() == 0


def degree_or_minus_one(poly):
    poly = RX(poly)
    return -1 if poly == 0 else int(poly.degree())


def horner_eval(poly, point):
    poly = RX(poly)
    point = K(point)
    result = K(0)
    for coefficient in reversed(poly.list()):
        result = result*point + K(coefficient)
    return result


def connection_polynomials(column_index):
    am3,bm3,am2,bm2,am1,aconst,bconst = [
        K(solution_entries[row][column_index]) for row in range(7)
    ]
    Atilde = RX(am3 + am2*U_poly + am1*U_poly^2 + aconst*U_poly^3)
    Btilde = RX(bm3 + bm2*U_poly + bconst*U_poly^3)
    if column_index == 0:
        Atilde += RX(K(universal_a1)*U_poly^4 + K(universal_a2)*U_poly^5)
    return Atilde,Btilde


def degree_bound(index,Atilde,Btilde):
    nF,dF = defect_num_den[index]
    terms = [
        degree_or_minus_one(nF)+degree_or_minus_one(sigma_den)+degree_or_minus_one(D_den)+3,
        degree_or_minus_one(dF)+degree_or_minus_one(Atilde)+degree_or_minus_one(sigma_num)+degree_or_minus_one(D_den),
        degree_or_minus_one(dF)+degree_or_minus_one(Btilde)+degree_or_minus_one(D_num)+degree_or_minus_one(sigma_den),
    ]
    return max(terms),terms


def cleared_value(index,Atilde,Btilde,point):
    nF,dF = defect_num_den[index]
    x0=K(point)
    u0=K(3*x0+1)
    nF0,dF0=horner_eval(nF,x0),horner_eval(dF,x0)
    nS0,dS0=horner_eval(sigma_num,x0),horner_eval(sigma_den,x0)
    nD0,dD0=horner_eval(D_num,x0),horner_eval(D_den,x0)
    A0,B0=horner_eval(Atilde,x0),horner_eval(Btilde,x0)
    return K(nF0*dS0*dD0*u0^3-dF0*(A0*nS0*dD0+B0*nD0*dS0))


def points(count):
    out=[]
    shell=0
    while len(out)<count:
        candidates=[QQ(0)] if shell==0 else [QQ(shell),QQ(-shell),(QQ(2*shell-1)/QQ(2)),(QQ(-(2*shell-1))/QQ(2))]
        for value in candidates:
            if value not in out:
                out.append(value)
                if len(out)>=count:
                    break
        shell+=1
    return out

results=[]
witnesses=[]
bounds=[]
for index,name in enumerate(defect_names):
    print("-"*78,flush=True)
    print("  verifying {}".format(name),flush=True)
    Atilde,Btilde=connection_polynomials(index)
    bound,term_bounds=degree_bound(index,Atilde,Btilde)
    test_points=points(bound+1)
    print("    term degree bounds: {}".format(term_bounds),flush=True)
    print("    exact evaluations: {}".format(bound+1),flush=True)
    passed=True
    witness=None
    started=time.time()
    for j,point in enumerate(test_points):
        value=cleared_value(index,Atilde,Btilde,point)
        if not is_zero_K(value):
            passed=False
            witness={"point":point,"value":value}
            print("    NONZERO witness at X={}".format(point),flush=True)
            break
        if (j+1)%5==0 or j+1==len(test_points):
            print("    passed {}/{}".format(j+1,len(test_points)),flush=True)
            gc.collect()
    print("    certified? {}".format(passed),flush=True)
    print("    elapsed seconds: {:.2f}".format(time.time()-started),flush=True)
    results.append(passed)
    witnesses.append(witness)
    bounds.append(bound)

FINAL=dict(DATA)
FINAL.update({
    "verification_method":"exact degree bound plus D+1 rational evaluations",
    "full_identity_results":results,
    "identity_witnesses":witnesses,
    "identity_degree_bounds":bounds,
})
save(FINAL,FINAL_CHECKPOINT)

with open(TEXT_SUMMARY,"w") as handle:
    handle.write("EXACT DEGREE-TWO FAST CERTIFICATE\n")
    for name,result,bound,witness in zip(defect_names,results,bounds,witnesses):
        handle.write("{}: {} (degree <= {})\n".format(name,result,bound))
        if witness is not None:
            handle.write("  witness X={}\n".format(witness["point"]))

print()
print("="*78)
print("FAST CHECKPOINT VERIFICATION SUMMARY")
for name,result,bound in zip(defect_names,results,bounds):
    print("  {:18s}: {} (degree <= {})".format(name,result,bound))
print("  all exact? {}".format(all(results)))
print("  final checkpoint: {}".format(FINAL_CHECKPOINT))
print("="*78)
