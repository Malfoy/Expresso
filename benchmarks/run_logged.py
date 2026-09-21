#!/usr/bin/env python3
"""Freeze a corpus manifest and log EXPRESSO completions without polling it.

status.json is an atomic snapshot as of updated_utc; progress.jsonl has one
event per completed dataset plus lifecycle events. Bytes mean compressed input
bytes of fully completed files. Dataset outputs are published only on success.
"""

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import resource
import socket
import subprocess
import time
import traceback


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def prepare(args):
    root = args.input.resolve(strict=True)
    run = args.run_dir.resolve()
    run.mkdir(parents=True, exist_ok=False)
    datasets, skipped = [], []
    folders = ([p for p in sorted(root.iterdir())
                if p.is_dir() and re.fullmatch(r"[0-9]{3}", p.name)]
               if args.numbered_folders else [root])
    if args.numbered_folders and len(folders) != 1000:
        raise ValueError(f"Expected folders 000–999, found {len(folders)}")
    folder_totals = {p.name: dict(datasets=0, compressed_bytes=0) for p in folders}
    paths = ((folder, path) for folder in folders for path in sorted(folder.rglob("*")))
    for folder, path in paths:
        if not path.is_file():
            continue
        name = path.name
        for suffix in (".gz", ".xz", ".zstd", ".zst"):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        if Path(name).suffix.lower() not in (".fa", ".fasta", ".fna", ".fq", ".fastq"):
            skipped.append(str(path))
            continue
        if any(c in str(path) for c in "\t\r\n"):
            raise ValueError(f"Unsupported characters in path: {path!r}")
        stat = path.stat()
        name = re.sub(r"[^A-Za-z0-9._-]", "_", Path(name).stem)
        prefix = f"{folder.name}_" if args.numbered_folders else ""
        datasets.append(dict(name=f"{prefix}{len(datasets)+1:06}_{name}", path=str(path),
                             folder=folder.name, compressed_bytes=stat.st_size,
                             mtime_ns=stat.st_mtime_ns))
        folder_totals[folder.name]["datasets"] += 1
        folder_totals[folder.name]["compressed_bytes"] += stat.st_size
    if not datasets:
        raise ValueError("No sequence datasets found")
    manifest = dict(created_utc=utc(), input_root=str(root), mode=args.mode,
                    planned_datasets=len(datasets),
                    planned_compressed_bytes=sum(x["compressed_bytes"] for x in datasets),
                    folders=folder_totals, excluded_files=skipped, datasets=datasets)
    save(run / "inputs.json", manifest)
    (run / "datasets.tsv").write_text("".join(
        f'{x["name"]}\t{x["path"]}\t{args.mode}\n' for x in datasets))
    print(json.dumps({k: v for k, v in manifest.items()
                      if k not in ("datasets", "folders", "excluded_files")}, indent=2))


def execute(args):
    run = args.run_dir.resolve(strict=True)
    manifest = json.loads((run / "inputs.json").read_text())
    datasets = {x["name"]: x for x in manifest["datasets"]}
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise ValueError("Missing command after --")
    # An exclusive file prevents accidentally launching the same manifest twice.
    with (run / "runner.pid").open("x") as handle:
        handle.write(f"{os.getpid()}\n")
    started = time.monotonic()
    state = dict(state="starting", started_utc=utc(), host=socket.gethostname(),
                 runner_pid=os.getpid(), command=command,
                 planned_datasets=manifest["planned_datasets"],
                 planned_compressed_bytes=manifest["planned_compressed_bytes"],
                 completed_datasets=0, completed_compressed_bytes=0,
                 failed_datasets=0, failed_compressed_bytes=0,
                 results_published=False)
    pattern = re.compile(r"^(.+): (\d+) records, (\d+) / (\d+) valid k-mers assigned$")
    completed = set()
    failed = set()
    child = None
    with (run / "progress.jsonl").open("x", buffering=1) as progress:
        def event(kind, **details):
            state.update(updated_utc=utc(), elapsed_seconds=time.monotonic() - started)
            state["completed_GB"] = state["completed_compressed_bytes"] / 1e9
            state["completed_GiB"] = state["completed_compressed_bytes"] / (1 << 30)
            entry = dict(event=kind, utc=state["updated_utc"],
                         elapsed_seconds=state["elapsed_seconds"],
                         completed_datasets=state["completed_datasets"],
                         completed_compressed_bytes=state["completed_compressed_bytes"],
                         **details)
            progress.write(json.dumps(entry) + "\n")
            save(run / "status.json", state)

        event("starting")
        try:
            for item in datasets.values():
                stat = Path(item["path"]).stat()
                if (stat.st_size, stat.st_mtime_ns) != (item["compressed_bytes"], item["mtime_ns"]):
                    raise ValueError(f'Input changed since manifest: {item["path"]}')
            with (run / "expresso.log").open("x", buffering=1) as log:
                child = subprocess.Popen(command, cwd=run, stdin=subprocess.DEVNULL,
                                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                         text=True, errors="replace", bufsize=1)
                state.update(state="running", expresso_pid=child.pid)
                event("started", expresso_pid=child.pid)
                for line in child.stdout:
                    log.write(line)
                    if line.startswith("DATASET_FAILED "):
                        entry = json.loads(line[len("DATASET_FAILED "):])
                        name = entry["name"]
                        if name not in datasets or name in completed or name in failed:
                            raise ValueError(f"Unexpected or repeated failure: {name}")
                        failed.add(name)
                        state["failed_datasets"] += 1
                        state["failed_compressed_bytes"] += datasets[name]["compressed_bytes"]
                        event("dataset_failed", name=name, error=entry["error"],
                              compressed_bytes=datasets[name]["compressed_bytes"])
                        continue
                    match = pattern.fullmatch(line.rstrip("\n"))
                    if not match:
                        continue
                    name, records, assigned, valid = match.groups()
                    if name not in datasets or name in completed or name in failed:
                        raise ValueError(f"Unexpected or repeated completion: {name}")
                    item = datasets[name]
                    completed.add(name)
                    state["completed_datasets"] += 1
                    state["completed_compressed_bytes"] += item["compressed_bytes"]
                    state["last_completion_elapsed_seconds"] = time.monotonic() - started
                    event("dataset_completed", name=name,
                          compressed_bytes=item["compressed_bytes"], records=int(records),
                          assigned_kmers=int(assigned), valid_kmers=int(valid))
                code = child.wait()
                state["exit_code"] = code
                if code != 0:
                    raise RuntimeError(f"EXPRESSO exited with code {code}")
                if len(completed) + len(failed) != len(datasets):
                    raise RuntimeError("Missing dataset completion events")
            state.update(state="finished_with_errors" if failed else "finished",
                         results_published=True, complete_coverage=not failed, finished_utc=utc())
            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            state.update(user_seconds=usage.ru_utime, system_seconds=usage.ru_stime,
                         max_rss_kib=usage.ru_maxrss)
            event(state["state"], exit_code=0, failed_datasets=len(failed))
            save(run / "summary.json", state)
            return 0
        except BaseException as error:
            if child is not None and child.poll() is None:
                child.terminate()
                child.wait()
            state.update(state="failed", error=str(error), finished_utc=utc())
            if child is not None:
                state["exit_code"] = child.returncode
            event("failed", error=str(error))
            save(run / "summary.json", state)
            traceback.print_exc()
            return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    prep = commands.add_parser("prepare")
    prep.add_argument("--input", type=Path, required=True)
    prep.add_argument("--run-dir", type=Path, required=True)
    prep.add_argument("--mode", choices=("reads", "unitigs", "auto"), default="auto")
    prep.add_argument("--numbered-folders", action="store_true",
                      help="Use only folders 000–999 directly under --input; require all 1,000")
    run = commands.add_parser("run")
    run.add_argument("--run-dir", type=Path, required=True)
    run.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args)
        return 0
    return execute(args)


if __name__ == "__main__":
    raise SystemExit(main())
