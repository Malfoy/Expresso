#!/usr/bin/env python3
"""Reproducible synthetic throughput check; requires built EXPRESSO and GGCAT."""
import argparse
import hashlib
import json
import os
import pathlib
import platform
import random
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="target/release/expresso")
    parser.add_argument("--reads", type=int, default=1_000_000)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--report", type=pathlib.Path)
    args = parser.parse_args()
    if args.reads < 1 or any(t < 1 for t in args.threads):
        parser.error("reads and threads must be positive")
    binary = str(pathlib.Path(args.binary).resolve())
    rng = random.Random(20260917)
    exons = ["".join(rng.choices("ACGT", k=300)) for _ in range(1000)]
    with tempfile.TemporaryDirectory(prefix="expresso-benchmark-") as scratch:
        work = pathlib.Path(scratch)
        with (work / "exons.fa").open("w") as out:
            for i, exon in enumerate(exons):
                out.write(f">exon_{i}\n{exon}\n")
        with (work / "reads.fa").open("w") as out:
            for i in range(args.reads):
                sequence = exons[rng.randrange(len(exons))]
                start = rng.randrange(151)
                out.write(f">read_{i}\n{sequence[start:start+150]}\n")
        (work / "datasets.txt").write_text("reads.fa\n")
        started = time.perf_counter()
        subprocess.run(
            [binary, "build", "-e", "exons.fa", "-i", "index", "-t",
             str(max(args.threads))], cwd=work, check=True, capture_output=True,
        )
        report = {
            "platform": platform.platform(), "cpu": platform.processor(),
            "available_cpus": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count(),
            "seed": 20260917, "exons": len(exons), "exon_length": 300,
            "reads": args.reads, "read_length": 150, "k": 31,
            "input_bytes": (work / "reads.fa").stat().st_size,
            "build_seconds": time.perf_counter() - started,
            "query": [],
        }
        digests = set()
        for iteration, threads in enumerate(args.threads):
            destination = f"output-{iteration}"
            started = time.perf_counter()
            subprocess.run(
                [binary, "quantify", "-i", "index", "-f", "datasets.txt",
                 "-o", destination, "-t", str(threads), "--mode", "reads",
                 "--compression", "none", "--format", "csv"],
                cwd=work, check=True, capture_output=True,
            )
            elapsed = time.perf_counter() - started
            digest = hashlib.sha256((work / destination / "global.csv").read_bytes()).hexdigest()
            digests.add(digest)
            report["query"].append({
                "threads": threads, "seconds": elapsed,
                "million_kmers_per_second": args.reads * 120 / elapsed / 1e6,
                "global_csv_sha256": digest,
            })
        if len(digests) != 1:
            raise RuntimeError("Different thread counts produced different abundances")
        report["identical_results"] = True
        text = json.dumps(report, indent=2) + "\n"
        print(text, end="")
        if args.report:
            args.report.write_text(text)


if __name__ == "__main__":
    main()
