from sage.all import *
import os

print("="*78)
print("AGGREGATE EXACT DEGREE-TWO COLUMN CERTIFICATES")
print("="*78)

names = ["intercept_Rquad", "slope_rA", "slope_rB", "slope_rC"]
labels = ["intercept Rquad", "slope rA", "slope rB", "slope rC"]
results = []
missing = []

for safe,label in zip(names,labels):
    path = "exact_degree2_{}_result.sobj".format(safe)
    if not os.path.exists(path):
        missing.append(path)
        results.append(None)
        print("  {:18s}: MISSING".format(label))
        continue
    data = load(path)
    answer = bool(data.get("certified", False))
    results.append(answer)
    print("  {:18s}: {} (degree <= {}, {}/{} points)".format(
        label, answer, data.get("degree_bound"),
        len(data.get("completed_indices", [])),
        len(data.get("points", [])),
    ))
    if data.get("witness") is not None:
        print("    witness X={}".format(data["witness"]["point"]))

print()
print("  missing result files: {}".format(missing))
print("  all four present and exact? {}".format(not missing and all(results)))

summary = {
    "labels": labels,
    "results": results,
    "missing": missing,
    "all_exact": (not missing and all(results)),
}
save(summary, "exact_degree2_all_columns_summary.sobj")

with open("exact_degree2_all_columns_summary.txt", "w") as handle:
    for label,result in zip(labels,results):
        handle.write("{}: {}\n".format(label,result))
    handle.write("missing={}\n".format(missing))
    handle.write("all_exact={}\n".format(summary["all_exact"]))

print("="*78)
