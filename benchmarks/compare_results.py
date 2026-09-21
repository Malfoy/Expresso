#!/usr/bin/env python3
"""Compare every CSV byte and every dataset metric between two EXPRESSO runs."""
import argparse
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("baseline", type=Path)
p.add_argument("candidate", type=Path)
p.add_argument("--report", type=Path, required=True)
a = p.parse_args()
def files(directory):
    return {path.relative_to(directory): path for path in directory.rglob("*")
            if path.is_file() and (".csv." in path.name or path.suffix == ".csv")}
left, right = files(a.baseline), files(a.candidate)
assert left.keys() == right.keys(), "Output CSV sets differ"
checks = []
for name in sorted(left):
    digests = []
    for path in [left[name], right[name]]:
        with path.open("rb") as f:
            digests.append(hashlib.file_digest(f, "sha256").hexdigest())
    assert digests[0] == digests[1], f"CSV differs: {name}"
    checks.append({"file":str(name), "sha256":digests[0]})
lm = json.loads((a.baseline / "manifest.json").read_text())
rm = json.loads((a.candidate / "manifest.json").read_text())
assert lm["datasets"] == rm["datasets"], "Dataset order, modes, or exact metrics differ"
for key in ("k", "unitig_k", "rounding", "statistics"):
    assert lm[key] == rm[key], f"Counting semantics differ: {key}"
report = {"baseline":str(a.baseline), "candidate":str(a.candidate),
          "all_csv_files_byte_identical":True, "dataset_metrics_identical":True,
          "datasets":len(lm["datasets"]), "csv_files":checks}
a.report.write_text(json.dumps(report,indent=2)+"\n")
print(f"{len(checks)} CSV files byte-identical; {len(lm['datasets'])} dataset metrics identical")
