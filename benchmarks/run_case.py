#!/usr/bin/env python3
"""Time one local EXPRESSO invocation and save CPU/RSS telemetry and logs."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--directory", type=Path, required=True)
p.add_argument("--label", required=True)
p.add_argument("command", nargs=argparse.REMAINDER)
a = p.parse_args()
command = a.command[1:] if a.command[:1] == ["--"] else a.command
assert command and "/" not in a.label
(a.directory / "logs").mkdir(parents=True, exist_ok=True)
(a.directory / "reports").mkdir(exist_ok=True)
started = datetime.now(timezone.utc).isoformat()
start = time.perf_counter()
with (a.directory / "logs" / (a.label + ".log")).open("w") as log, (a.directory / "reports" / (a.label + ".telemetry.jsonl")).open("w") as telemetry:
    timing = a.directory / "reports" / (a.label + ".time.txt")
    process = subprocess.Popen(["/usr/bin/time", "-v", "-o", str(timing), *command], stdout=log, stderr=subprocess.STDOUT, env={**os.environ, "LC_ALL":"C"})
    print(f"{a.label}: PID {process.pid}", flush=True)
    while process.poll() is None:
        try:
            children = Path(f"/proc/{process.pid}/task/{process.pid}/children").read_text().split()
            pid = children[0] if children else process.pid
            status = dict(line.split(":", 1) for line in Path(f"/proc/{pid}/status").read_text().splitlines() if ":" in line)
            stat = Path(f"/proc/{pid}/stat").read_text().split(") ", 1)[1].split()
            snapshot = {"elapsed":time.perf_counter()-start, "cpu_seconds":(int(stat[11])+int(stat[12]))/os.sysconf("SC_CLK_TCK"), "rss_kib":int(status.get("VmRSS", "0 kB").split()[0]), "threads":int(status["Threads"])}
            telemetry.write(json.dumps(snapshot)+"\n"); telemetry.flush()
        except (FileNotFoundError, ProcessLookupError):
            pass
        try:
            process.wait(timeout=1)
        except subprocess.TimeoutExpired:
            pass
elapsed = time.perf_counter()-start
usage = dict(line.strip().split(": ",1) for line in timing.read_text().splitlines() if ": " in line)
user, system = float(usage["User time (seconds)"]), float(usage["System time (seconds)"])
report = {"label":a.label, "command":command, "started_utc":started, "elapsed_seconds":elapsed, "user_seconds":user, "system_seconds":system, "average_cpu_cores":(user+system)/elapsed, "max_rss_kib":int(usage["Maximum resident set size (kbytes)"]), "major_faults":int(usage["Major (requiring I/O) page faults"]), "voluntary_context_switches":int(usage["Voluntary context switches"]), "involuntary_context_switches":int(usage["Involuntary context switches"]), "exit_code":process.returncode}
(a.directory / "reports" / (a.label + ".json")).write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2),flush=True)
raise SystemExit(process.returncode)
