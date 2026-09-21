#!/usr/bin/env python3
"""Wait for a logged EXPRESSO run, verify it, then compare Zstandard levels."""
import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--previous-results", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=64)
    parser.add_argument("--compact", action="store_true", help="Verify EAB results against exact CSVs")
    args = parser.parse_args()
    run = args.run_dir.resolve(strict=True)
    start = time.monotonic()
    state = dict(state="waiting_for_quantification", pid=os.getpid(),
                 started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())

    def update(**fields):
        state.update(fields, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     elapsed_seconds=time.monotonic() - start)
        tmp = run / "benchmark_status.json.tmp"
        tmp.write_text(json.dumps(state, indent=2) + "\n")
        tmp.replace(run / "benchmark_status.json")
        print(json.dumps(state), flush=True)

    try:
        update()
        while not (run / "summary.json").exists():
            if not (run / "status.json").exists():
                time.sleep(1)
                continue
            current = json.loads((run / "status.json").read_text())
            if current["state"] == "failed":
                raise RuntimeError(current.get("error", "Quantification failed"))
            if current["state"] != "finished":
                os.kill(current["runner_pid"], 0)
            time.sleep(5)
        quant = json.loads((run / "summary.json").read_text())
        if quant["state"] != "finished" or not quant["results_published"]:
            raise RuntimeError("Quantification did not publish successful results")
        if not args.compact:
            pilot = json.loads((run / "pilot/summary.json").read_text())
            if not pilot["all_requested_levels_finished"]:
                raise RuntimeError("Pilot incomplete")
            level3 = next(x for x in pilot["levels"] if x["level"] == 3)
            if not level3["all_compressed_bytes_identical_to_original"]:
                raise RuntimeError("Pilot level 3 does not reproduce original encoder output")
        update(state="verifying_quantification")
        with (run / "comparison.log").open("x") as log:
            verify = ([sys.executable, str(run / "verify_compact.py"),
                       "--exact", str(args.previous_results), "--compact", str(run / "results"),
                       "--report", str(run / "compact-verification.json"), "--jobs", str(min(32, args.jobs))]
                      if args.compact else
                      [sys.executable, str(run / "compare_results.py"),
                       str(args.previous_results), str(run / "results"),
                       "--report", str(run / "quantification-verification.json")])
            subprocess.run(verify,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        update(state="comparing_compression_levels")
        command = [str(run / "bin/zstd_levels"), "--results", str(run / "results"),
                   "--output", str(run / "compression"), "--levels", "1,3,6,9,12,15,19",
                   "--jobs", str(args.jobs)]
        state["command"] = command
        with (run / "compression.log").open("x") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        summary = json.loads((run / "compression/summary.json").read_text())
        if not summary["all_requested_levels_finished"] or not all(
                x.get("all_content_bytes_verified", x.get("all_csv_bytes_verified")) for x in summary["levels"]):
            raise RuntimeError("Incomplete compression verification")
        level3 = next(x for x in summary["levels"] if x["level"] == 3)
        if not level3["all_compressed_bytes_identical_to_original"]:
            raise RuntimeError("Level 3 differs from original output")
        update(state="finished", finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        return 0
    except BaseException as error:
        update(state="failed", error=str(error))
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
