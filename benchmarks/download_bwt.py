#!/usr/bin/env python3
"""Resume the recorded BWT corpus with disjoint, balanced rsync streams."""
import argparse
import json
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--directory", type=Path, required=True)
p.add_argument("--streams", type=int, default=4)
a = p.parse_args()
selection = json.loads((a.directory / "selection.json").read_text())
shards = [[] for _ in range(a.streams)]
sizes = [0] * a.streams
for item in sorted(selection["files"], key=lambda item: item["bytes"], reverse=True):
    i = min(range(a.streams), key=sizes.__getitem__)
    shards[i].append(item["name"])
    sizes[i] += item["bytes"]
processes = []
handles = []
for i, names in enumerate(shards):
    listing = a.directory / f"rsync-shard-{i}.txt"
    listing.write_text("".join(name + "\n" for name in names))
    log = (a.directory / "logs" / f"download-shard-{i}.log").open("w")
    handles.append(log)
    processes.append(subprocess.Popen([
        "rsync", "-a", "--partial", "--append-verify", "--info=progress2",
        f"--files-from={listing}", "-e", "ssh -F /dev/null -o BatchMode=yes -o ConnectTimeout=20",
        "alimasse@bwt.lifl.fr:/DATA/alimasse/000/", str(a.directory / "data") + "/",
    ], stdout=log, stderr=subprocess.STDOUT))
while any(p.poll() is None for p in processes):
    downloaded = sum(min(f.stat().st_size, item["bytes"]) for item in selection["files"]
                     if (f := a.directory / "data" / item["name"]).exists())
    complete = sum((f := a.directory / "data" / item["name"]).exists() and f.stat().st_size == item["bytes"]
                   for item in selection["files"])
    print(f"{downloaded:,}/{selection['total_bytes']:,} bytes; {complete}/{len(selection['files'])} complete files", flush=True)
    time.sleep(15)
for h in handles:
    h.close()
if any(p.returncode for p in processes):
    raise SystemExit("A transfer failed; see logs/download-shard-*.log and rerun to resume")
print("Transfer complete", flush=True)
