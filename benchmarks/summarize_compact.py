#!/usr/bin/env python3
"""Summarize a verified EAB run and its complete Zstandard comparison."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report_dir", type=Path)
    parser.add_argument("--csv-benchmark", type=Path, required=True)
    args = parser.parse_args()
    p = args.report_dir
    load = lambda name: json.loads((p / name).read_text())
    run, launch, accuracy = load("summary.json"), load("launch.json"), load("compact-verification.json")
    comp = load("compression-summary.json" if (p / "compression-summary.json").exists() else "compression/summary.json")
    old = json.loads(args.csv_benchmark.read_text())
    old_sizes = {r["level"]: r["total_compressed_bytes"] for r in old["levels"]}
    assert run["state"] == "finished" and run["results_published"]
    assert accuracy["all_vectors_verified"] and accuracy["reference_table_identical"]
    assert comp["all_requested_levels_finished"] and comp["format"] == "compact"
    assert comp["pilot_limit"] is None
    rows = sorted(comp["levels"], key=lambda r: r["level"])
    assert [r["level"] for r in rows] == [1,3,6,9,12,15,19]
    assert all(r["all_content_bytes_verified"] and r["files"] == accuracy["vectors"] + 1 for r in rows)
    baseline = next(r for r in rows if r["level"] == 3)
    assert baseline["all_compressed_bytes_identical_to_original"]
    with (p / "compression-levels.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    lines = ["# EAB 8-bit: full BWT 000 benchmark", "",
             f"**{run['completed_datasets']:,} datasets**, **{run['completed_compressed_bytes']/1e9:.6f} GB** "
             f"compressed input, quantified in **{run['elapsed_seconds']:.3f} seconds** with 64 workers / 16 concurrent datasets.", "",
             f"Reference: GENCODE v49 whole transcripts, k=31, {accuracy['targets']:,} targets. "
             "Same complete inputs and reference as the exact-CSV benchmark; unitig header abundances.", "",
             "## Storage", "",
             f"The {accuracy['vectors']:,} vectors hold {accuracy['values']:,} abundance codes: "
             f"**{accuracy['values']/1e6:.6f} MB** at 8 bits before codebooks, headers and compression.", "",
             f"Uncompressed EAB files plus shared reference and manifest: "
             f"**{(baseline['uncompressed_bytes']+baseline['manifest_bytes'])/1e6:.6f} MB**. "
             "The old uncompressed abundance CSVs occupied 53,026.407168 MB.", "",
             f"Default level 3 complete bundle: **{baseline['total_with_manifest_bytes']:,} bytes "
             f"({baseline['total_with_manifest_bytes']/1e6:.6f} MB)**:", "",
             f"- Dataset vectors: {baseline['dataset_compressed_bytes']:,} bytes.",
             f"- Global vector: {baseline['global_compressed_bytes']:,} bytes.",
             f"- Shared reference table, stored once: {baseline['shared_reference_compressed_bytes']:,} bytes.",
             f"- Plain manifest: {baseline['manifest_bytes']:,} bytes.", "",
             "Each total below includes all vectors, the shared reference compressed at that level, "
             "and the unchanged plain manifest. Benchmark logs and measurement JSON files are excluded.", "",
             "| Zstd level | Vectors MB | Shared reference MB | Complete bundle MB | Smaller than compact level 3 | Recompression + verification seconds | Old CSV GB |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        vectors = row["dataset_compressed_bytes"] + row["global_compressed_bytes"] + row["statistics_compressed_bytes"]
        saving = 100*(1-row["total_with_manifest_bytes"]/baseline["total_with_manifest_bytes"])
        lines.append(f"| {row['level']}{' (default)' if row['level']==3 else ''} | {vectors/1e6:.6f} "
                     f"| {row['shared_reference_compressed_bytes']/1e6:.6f} | {row['total_with_manifest_bytes']/1e6:.6f} "
                     f"| {saving:.2f}% | {row['elapsed_seconds_including_decode_and_verification']:.2f} "
                     f"| {old_sizes[row['level']]/1e9:.6f} |")
    lines += ["", f"At default level 3, the complete compact bundle is "
              f"**{old_sizes[3]/baseline['total_with_manifest_bytes']:.1f}× smaller** than the old level-3 abundance CSVs. "
              "The CSV baseline excludes its separate reference metadata, while the compact total includes its shared table and manifest.", "",
              "## Accuracy", "",
              "Quantization is lossy; the later Zstandard passes are lossless. An independent Python decoder "
              "checked every stored code against the nearest codebook representative for each original exact CSV count, "
              "validated CRCs and the reference SHA-256, and recomputed the per-vector error statistics.", "",
              f"- Verified: **{accuracy['values']:,} counts across {accuracy['vectors']:,} vectors**.",
              "- Zero/nonzero status is preserved for every count. Counts 0–15 are exact at 8 bits.",
              f"- Exactly reconstructed values (including zeros): **{100*accuracy['exact_values']/accuracy['values']:.6f}%**.",
              f"- Nonzero values: {accuracy['nonzero_values']:,}.",
              f"- Mean relative error over nonzero values: **{100*accuracy['mean_relative_error_nonzero']:.6f}%**.",
              f"- Maximum relative error over nonzero values: **{100*accuracy['maximum_relative_error_nonzero']:.6f}%**.",
              f"- Mean absolute error over all values: {accuracy['mean_absolute_error']:.6f}.",
              f"- Maximum absolute error: {accuracy['maximum_absolute_error']:,}.", "",
              "These measured errors include all dataset vectors and the global vector. Large counts can have "
              "large absolute errors even with small relative errors. Global is quantized from the exact total, "
              "so decoded dataset sums can differ from the decoded global. Original exact counts remain in the previous CSV run.", "",
              "## Measurement and files", "",
              f"Zstandard {comp['zstd_version']}, one independent streaming frame per file, no dictionary, "
              "64 file workers. Each level ran once. Timings include baseline decompression, recompression, writing "
              "and readback verification, separate from quantification. The host is shared and caches were not dropped. "
              "Every compressed file round-tripped byte-for-byte, and level 3 reproduced the baseline compressed files exactly.", "",
              f"BWT directory: `{launch['run_dir']}`.", "",
              "- `results/`: default compact result bundle.",
              "- `compression/level-XX/`: each recompressed bundle (includes a copy of the manifest).",
              "- `compact-verification.json`: independent verification and detailed per-vector error measurements.",
              "- `compression/summary.json`, `compression/level-XX/files.json`: exact byte totals and file measurements.", "",
              "MB and GB are decimal (10^6 and 10^9 bytes).", ""]
    (p / "REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines[:36]))


if __name__ == "__main__":
    main()
