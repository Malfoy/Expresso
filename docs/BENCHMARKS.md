# Benchmarks and measurement scope

These results measure three different workloads. All use GENCODE v49 **whole
transcripts** as reference records (533,740 targets), k=31, and unitig header
abundances. They do not benchmark the extracted exon or genomic-gene FASTAs.
Input sizes describe compressed sequence files; output sizes describe abundance
results, not recoverable copies of the sequence data. MB, GB, and TB are decimal.

## Full corpus: 1,000 folders

EXPRESSO attempted every eligible dataset in folders `000`–`999` on BWT. The
run used 64 query workers, 16 concurrent datasets, 8-bit EAB, zstd level 3, and
`--keep-going`. Optional statistics were disabled; the global sum is always
computed. The source machine exposes 64 cores from two Intel Xeon Gold 6430
processors and approximately 503 GiB RAM.

| Measurement | Result |
| --- | ---: |
| Started, UTC | 2026-09-19 00:11:08 |
| Finished, UTC | 2026-09-21 09:41:50 |
| Wall time | 207,042.27 s (57 h 30 min 42 s) |
| Attempted datasets | 822,141 |
| Successful datasets | 822,115 |
| Failed datasets | 26 |
| Attempted compressed bytes | 75,518,483,329,066 |
| Successful compressed bytes | 75,517,783,087,166 |
| Failed inputs, compressed bytes | 700,241,900 |
| Mean successful-input throughput | approximately 1.31 TB/hour |
| Peak child RSS | 10,088,764 KiB (9.62 GiB) |

The 26 failed inputs contain malformed abundance headers. They are excluded
entirely, including any counts accumulated before the error. The run published
results with `finished_with_errors` and `coverage.complete: false`. All files
were attempted; complete coverage is not claimed.

| Published output | Bytes |
| --- | ---: |
| 822,115 dataset vectors | 81,398,687,000 |
| Global sum vector | 360,626 |
| Shared reference table | 12,246,313 |
| Plain JSON manifest | 734,754,031 |
| **Complete bundle** | **82,146,047,970 (82.15 GB)** |

These totals exclude indexes, logs, source snapshots, and the separate small
validation runs. Global sums use exact rounded counts before independent 8-bit
quantization. Presence/absence is preserved; the decoded sum can differ from
sums of decoded dataset vectors.

The wall time covers the detached runner's input checks, index loading, all
query processing, per-dataset encoding/writing, and final aggregation/publishing.
The host/storage were shared; this was one complete run, not a controlled series
of repeated throughput trials. Peak RSS is the child-process statistic collected
by the wrapper, not total machine memory use. No independent exact-CSV audit of
all 822,115 output vectors was performed.

[Machine-readable final measurements](benchmarks/full-corpus.json).

To run the same configuration on an available local corpus:

```bash
expresso quantify --index transcript-index --fof datasets.tsv --output results \
  --threads 64 --jobs 16 --mode unitigs \
  --format compact --bits 8 --compression zstd --keep-going
```

The large input corpus, indexes, and output vectors are not distributed in this
repository. To use the background runner, follow the main README's
[long-run example](../README.md#long-runs-and-progress-logs).

## Compact format validation: folder 000

On BWT, **873 datasets / 83,665,109,503 compressed bytes** were quantified in
**193.641 seconds**, using 64 query workers and 16 concurrent datasets. All
inputs completed successfully. The 873 dataset vectors plus the global vector
contain **466,488,760 values**.

The exact CSV baseline and the compact run use the same inputs, index, weighting,
and row order. An independent Python decoder checked every compact code against
the nearest stored codebook representative for its original exact CSV count,
validated CRC checksums and the reference SHA-256, and recomputed error metrics.
Dataset counting metrics and shared reference metadata also matched.

| Quantization measurement | Result |
| --- | ---: |
| Checked values | 466,488,760 |
| Nonzero values | 45,165,013 |
| Exactly reconstructed, including zeros | 92.833961% |
| Mean relative error, nonzero values | 1.218569% |
| Maximum relative error, nonzero values | 5.555556% |
| Mean absolute error, all values | 17.448160 |
| Maximum absolute error | 202,442,575 |

Large values can have large absolute errors despite modest relative errors.
These errors include the independently quantized global vector and are specific
to this validation corpus. They are not universal bounds for 8-bit output.
Zero/nonzero status was preserved for every value; counts 0–15 are exact.

[Machine-readable accuracy measurements](benchmarks/compact-validation.json).

The bare 8-bit payloads occupy **466.49 MB**. The complete uncompressed EAB bundle,
including codebooks, headers, reference metadata, and manifest, occupies
**532.60 MB**. The original uncompressed abundance CSVs occupy **53.03 GB**,
largely because every row repeats target identifiers and names.

## Zstandard level comparison

Every level was applied to the complete 873-dataset compact result. Each total
includes dataset and global vectors, the shared reference compressed at that
level, and the unchanged 774,682-byte plain manifest. No optional statistics
vectors were present.

| Level | Compact bundle MB | Recompression + verification s | Previous abundance CSV GB |
| ---: | ---: | ---: | ---: |
| 1 | 97.612349 | 0.460 | 8.632689 |
| 3 | 99.277013 | 0.575 | 9.270625 |
| 6 | 93.808303 | 1.199 | 8.358948 |
| 9 | 90.531167 | 1.838 | 7.647872 |
| 12 | 90.282347 | 4.729 | 7.628324 |
| 15 | 89.693461 | 11.935 | 7.595507 |
| 19 | 81.445088 | 35.711 | 6.111308 |

At the default level 3, the complete compact bundle is **93.4× smaller** than
the previous compressed abundance CSVs. The CSV baseline excludes its separate
reference metadata, while the compact total includes its reference and manifest.
Compression sizes need not decrease at every adjacent level: level 1 happened
to be smaller than level 3 on these files.

The benchmark used Zstandard 1.5.7, one independent streaming frame per file,
no dictionary, and 64 file workers. Each level ran once. Timing includes source
decoding, recompression, writes, and readback verification; it does not include
quantification. Caches were not dropped and the machine was shared. Every
outer-compression round trip preserved all uncompressed bytes. Level 3 reproduced
the source compressed files byte for byte. The underlying abundance quantization
remains lossy.

[Exact level measurements](benchmarks/compact-zstd-levels.csv). The CSV's
`total_with_manifest_bytes` is the complete bundle size; its compression-ratio
and percentage fields use compressed files alone, excluding the plain manifest.

```bash
cargo build --release --locked --example zstd_levels
target/release/examples/zstd_levels \
  --results results --output compression-comparison \
  --levels 1,3,6,9,12,15,19 --jobs 64
```

The helper retains every level and writes `summary.json` plus per-file
measurements under `level-XX/files.json`. Use a new destination. For compact
results, the retained bundles can be passed directly to `expresso export`.

## Local 10 GB optimization benchmark

The earlier local benchmark used 97 files, totaling **10,000,533,099 compressed
bytes**, downloaded from folder 000. Files were selected by sorting eligible
names, shuffling with Python seed `20260918`, and taking whole files until at
least 10 GB. One selected file was empty and correctly produced a zero vector.

Machine: Intel Core i9-13950HX, 24 physical cores / 32 logical CPUs, approximately
61 GiB RAM, local NVMe storage. All comparisons used 32 query workers, exact
zstd CSV output, and no optional statistics.

| Configuration | Wall seconds | Peak query RSS GiB |
| --- | ---: | ---: |
| Original implementation, 1 dataset at a time | 348.95 | 1.24 |
| Original implementation, 8 concurrent datasets | 273.20 | 1.54 |
| Optimized implementation, 8 concurrent datasets | 48.94 | 1.30 |
| Optimized implementation, 8 datasets, repeat | 49.09 | 1.29 |
| Optimized implementation, 16 concurrent datasets | 54.64 | 1.35 |

The optimized 8-dataset median was **49.02 s**, about **7.12× faster** than the
original default and **5.57× faster** at matched concurrency. Improvements
included reusable buffers, bounded worker queues, streaming SSHash queries,
size-based scheduling, and concurrent output. Individual effects were not
isolated in a factorial benchmark.

All **99 CSVs** (97 datasets, global, and reference metadata) were byte-identical
to the baseline, and all dataset metrics matched. Timings include index loading,
input I/O, decompression, parsing, counting, aggregation, and writing. They exclude
index construction, downloads, and verification. Caches were not dropped: these
are warm-cache measurements, not predictions for remote storage.

[Machine-readable local summary](benchmarks/local-10gb.json).

The reference index had 186,672,793 canonical k-mers: 98,813,138 shared and
excluded, leaving 87,859,655 uniquely assigned. Of 533,740 transcript targets,
201,622 had no unique k-mers. Its original build took 19.01 s at approximately
4.12 GiB peak RSS; index build speed is not part of the query speedup claim.

For an exact CSV comparison with the current default-compact version, explicitly
pass `--format csv`:

```bash
expresso quantify --index transcript-index --fof datasets.tsv --output exact-results \
  --threads 32 --jobs 8 --mode unitigs --format csv --compression zstd
```

`benchmarks/run_case.py` records process telemetry, and
`benchmarks/compare_results.py` compares exact CSV outputs. The synthetic
benchmark in the main README is separately reproducible without this corpus.
