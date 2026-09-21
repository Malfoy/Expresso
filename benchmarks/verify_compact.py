#!/usr/bin/env python3
"""Independent EAB decoder: compare all vectors against an exact CSV result set."""
import argparse
import bisect
import concurrent.futures
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import time
import zlib


def open_zstd(path):
    return subprocess.Popen(["zstd", "-qdc", str(path)], stdout=subprocess.PIPE, text=True)


def verify_one(task):
    vector_path, csv_path, expected_hash, target_count, recorded = task
    data = subprocess.check_output(["zstd", "-qdc", str(vector_path)])
    magic, bits, encoding, reserved, n, size, payload_bytes, maximum, reference = struct.unpack("<8sBBHQIQQ32s", data[:72])
    assert magic == b"EXPRAB01" and 2 <= bits <= 16 and encoding == reserved == 0
    assert n == target_count and reference.hex() == expected_hash
    assert 1 <= size <= 1 << bits and payload_bytes == (n * bits + 7) // 8
    assert len(data) == 72 + size * 8 + payload_bytes + 4
    assert zlib.crc32(data[:-4]) == struct.unpack("<I", data[-4:])[0]
    levels = struct.unpack(f"<{size}Q", data[72:72 + 8 * size])
    assert levels[0] == 0 and levels[-1] == maximum
    assert all(a < b for a, b in zip(levels, levels[1:]))
    payload = memoryview(data)[72 + size * 8:-4]
    if n * bits % 8:
        assert payload[-1] >> (n * bits % 8) == 0
    last = (1 << bits) - 1
    prefix = maximum if maximum <= last else min(15, last // 4)
    assert list(levels[:prefix + 1]) == list(range(prefix + 1))
    exact = nonzero = count = observed_max = maximum_error = 0
    sum_error = sum_relative = maximum_relative = 0.0
    child = open_zstd(csv_path)
    try:
        rows = csv.reader(child.stdout)
        assert next(rows) == ["exon_id", "exon_name", "abundance"]
        for i, row in enumerate(rows):
            assert i < n and int(row[0]) == i + 1
            value = int(row[2])
            if bits == 8:
                code = payload[i]
            else:
                offset = i * bits
                word = int.from_bytes(payload[offset // 8: (offset + bits + 7) // 8], "little")
                code = (word >> (offset % 8)) & last
            assert code < len(levels)
            if value <= prefix:
                expected = value
            else:
                hi = bisect.bisect_left(levels, value)
                assert hi < len(levels)
                expected = hi - 1 if value - levels[hi - 1] <= levels[hi] - value else hi
            assert code == expected, f"Wrong quantization in {vector_path}, target {i+1}"
            decoded = levels[code]
            assert (decoded == 0) == (value == 0)
            error = abs(decoded - value)
            exact += error == 0
            count += 1
            observed_max = max(observed_max, value)
            maximum_error = max(maximum_error, error)
            sum_error += error
            if value:
                nonzero += 1
                relative = error / value
                maximum_relative = max(maximum_relative, relative)
                sum_relative += relative
        assert child.wait() == 0
    finally:
        child.stdout.close()
        if child.poll() is None:
            child.terminate()
            child.wait()
    assert count == n and observed_max == maximum
    measured = dict(targets=count, maximum=maximum, exact_values=exact,
                    nonzero_values=nonzero, maximum_absolute_error=maximum_error,
                    mean_absolute_error=sum_error / max(1, n),
                    maximum_relative_error_nonzero=maximum_relative,
                    mean_relative_error_nonzero=sum_relative / max(1, nonzero))
    for key, value in measured.items():
        assert math.isclose(value, recorded[key], rel_tol=1e-12, abs_tol=1e-12), (key, value, recorded[key])
    return dict(file=str(vector_path), **measured)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact", type=Path, required=True)
    parser.add_argument("--compact", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=32)
    args = parser.parse_args()
    started = time.monotonic()
    old = json.loads((args.exact / "manifest.json").read_text())
    new = json.loads((args.compact / "manifest.json").read_text())
    assert new["format"] == "compact" and new["compact_format_version"] == 1
    assert not old["stats"] and not new["stats"], "Use no-statistics results for this benchmark"
    assert len(old["datasets"]) == len(new["datasets"])
    for a, b in zip(old["datasets"], new["datasets"]):
        assert all(a[k] == b[k] for k in ["name", "path", "mode", "metrics"])
    targets = new["reference"]["targets"]
    digest = hashlib.sha256(b"EXPRESSO-reference-v1\0" + struct.pack("<Q", targets))
    process = open_zstd(args.compact / new["reference"]["file"])
    with (args.exact / "exons.csv").open() as old_table:
        left, right = csv.reader(old_table), csv.reader(process.stdout)
        assert next(left) == next(right) == ["exon_id", "exon_name", "length", "unique_kmers"]
        for i in range(targets):
            a, b = next(left), next(right)
            assert a == b and int(a[0]) == i + 1
            name = a[1].encode()
            digest.update(struct.pack("<Q", len(name)) + name + struct.pack("<QQ", int(a[2]), int(a[3])))
        assert next(left, None) is None and next(right, None) is None
    assert process.wait() == 0
    assert digest.hexdigest() == new["reference"]["sha256"]
    tasks = [(args.compact / b["output"], args.exact / a["output"], digest.hexdigest(), targets, b["quantization"])
             for a, b in zip(old["datasets"], new["datasets"])]
    tasks.append((args.compact / new["global"]["output"], args.exact / "global.csv.zst",
                  digest.hexdigest(), targets, new["global"]["quantization"]))
    records = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(verify_one, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            records.append(future.result())
            if len(records) % 50 == 0 or len(records) == len(tasks):
                print(f"Verified {len(records)}/{len(tasks)} vectors ({time.monotonic()-started:.1f}s)", flush=True)
    records.sort(key=lambda x: x["file"])
    nonzero = sum(r["nonzero_values"] for r in records)
    total = sum(r["targets"] for r in records)
    report = dict(all_vectors_verified=True, reference_table_identical=True,
                  exact_dataset_metrics_identical=True, zero_nonzero_status_preserved=True,
                  vectors=len(records), datasets=len(new["datasets"]), targets=targets,
                  values=total, exact_values=sum(r["exact_values"] for r in records),
                  nonzero_values=nonzero,
                  maximum_absolute_error=max(r["maximum_absolute_error"] for r in records),
                  maximum_relative_error_nonzero=max(r["maximum_relative_error_nonzero"] for r in records),
                  mean_relative_error_nonzero=sum(r["mean_relative_error_nonzero"]*r["nonzero_values"] for r in records)/max(1,nonzero),
                  mean_absolute_error=sum(r["mean_absolute_error"]*r["targets"] for r in records)/max(1,total),
                  elapsed_seconds=time.monotonic()-started, files=records)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k:v for k,v in report.items() if k != "files"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
