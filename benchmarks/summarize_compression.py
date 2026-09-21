#!/usr/bin/env python3
"""Render the verified full-corpus compression measurements as Markdown and CSV."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report_dir", type=Path)
    args = parser.parse_args()
    p = args.report_dir
    quant = json.loads((p / "summary.json").read_text())
    launch = json.loads((p / "launch.json").read_text())
    baseline = json.loads((p / "baseline-sizes.json").read_text())
    compression_report = p / "compression-summary.json"
    if not compression_report.exists():
        compression_report = p / "compression/summary.json"
    comparison = json.loads(compression_report.read_text())
    verification = json.loads((p / "quantification-verification.json").read_text())
    assert quant["state"] == "finished" and quant["results_published"]
    assert comparison["all_requested_levels_finished"]
    assert verification["all_csv_files_byte_identical"] and verification["dataset_metrics_identical"]
    levels = sorted(comparison["levels"], key=lambda row: row["level"])
    assert all(row["files"] == quant["completed_datasets"] + 1 for row in levels)
    assert all(row.get("all_content_bytes_verified", row.get("all_csv_bytes_verified")) for row in levels)
    assert all(row["original_compressed_bytes"] == baseline["total_compressed_bytes"] for row in levels)
    assert len({row["uncompressed_bytes"] for row in levels}) == 1
    default = next(row for row in levels if row["level"] == 3)
    assert default["all_compressed_bytes_identical_to_original"]
    raw = default["uncompressed_bytes"]
    with (p / "compression-levels.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(levels[0]))
        writer.writeheader()
        writer.writerows(levels)
    lines = [
        "# BWT: full 000 abundance CSV compression comparison", "",
        f"EXPRESSO processed **{quant['completed_datasets']:,} datasets**, "
        f"**{quant['completed_compressed_bytes']/1e9:.6f} GB compressed input**, "
        f"in **{quant['elapsed_seconds']:.3f} seconds** with 64 query workers and 16 concurrent datasets.", "",
        "Reference: the existing GENCODE v49 whole-transcript index, k=31. "
        "This run uses transcript targets, not the subsequently generated exon/genomic-gene FASTAs.", "",
        f"Default Zstandard level 3: **{baseline['total_compressed_bytes']:,} bytes "
        f"({baseline['total_compressed_bytes']/1e9:.6f} GB)** total abundance CSV output.", "",
        f"- Dataset CSVs: {baseline['dataset_files']:,} files, {baseline['dataset_compressed_bytes']:,} bytes.",
        f"- Global CSV: {baseline['global_compressed_bytes']:,} bytes.",
        f"- Uncompressed CSV content across these files: {raw:,} bytes ({raw/1e9:.6f} GB), measured by decompression.",
        "- The reference metadata table `exons.csv`, manifests and logs are excluded from these sizes.", "",
        "| Zstandard level | Total bytes | Total GB | Smaller than level 3 | Recompression + verification (seconds) |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in levels:
        lines.append(f"| {row['level']}{' (default)' if row['level']==3 else ''} "
                     f"| {row['total_compressed_bytes']:,} | {row['total_compressed_bytes']/1e9:.6f} "
                     f"| {row['percent_smaller_than_original']:.2f}% "
                     f"| {row['elapsed_seconds_including_decode_and_verification']:.2f} |")
    lines += ["", "## Verification and measurement", "",
              f"All {len(verification['csv_files']):,} CSV files (including the reference metadata table) "
              f"and all {verification['datasets']:,} dataset metrics match the previous successful run. "
              "Every recompressed abundance CSV was decompressed and compared byte-for-byte with its source. "
              "Level 3 also reproduces all original compressed bytes exactly.", "",
              f"Each CSV was compressed separately with the same Zstandard {comparison['zstd_version']} "
              f"library and streaming settings as EXPRESSO, with {comparison['workers']} independent file workers. "
              "No dictionary, no frame checksum, unknown source size. These are full-corpus measurements, not projections.", "",
              "Level timings include reading/decompressing the baseline, compression, file writes, and readback "
              "verification. They are separate from the EXPRESSO quantification time. Each level ran once on "
              "a shared host; filesystem caches were not dropped. GB is decimal (10^9 bytes).", "",
              "## Files on BWT", "", f"Run directory: `{launch['run_dir']}`.", "",
              "- Original default output: `results/`.",
              "- Recompressed files: `compression/level-XX/`.",
              "- Detailed per-file measurements: `compression/level-XX/files.json`.",
              "- Source snapshot, manifests, logs and method description are saved in the run directory.", ""]
    (p / "REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines[:20]))


if __name__ == "__main__":
    main()
