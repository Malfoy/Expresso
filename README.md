# EXPRESSO

**EXon-level RNA EXPRESsion quantificatiOn**

A parallel Rust tool for counting exon-specific k-mer observations in reads or
abundance-annotated unitigs. Build an index once, quantify many datasets, and
store one compact abundance vector per dataset plus a global sum.

EXPRESSO combines [**GGCAT simplitigs**](https://github.com/algbio/ggcat),
the [**Rust SSHash dictionary**](https://github.com/COMBINE-lab/sshash-rs), and
[**Helicase FASTA/FASTQ parsing**](https://github.com/imartayan/helicase).
Reference and query files can be uncompressed,
gzip, xz, or Zstandard. K-mers shared by multiple reference records are excluded.

## How it works

```mermaid
flowchart LR
    A[Reference FASTA] --> B[GGCAT simplitigs]
    B --> C[Rust SSHash index]
    A --> D[K-mer ownership]
    C --> E[Parallel queries]
    D --> E
    F[Reads or weighted unitigs] --> G[Helicase parsing]
    G --> E
    E --> H[Exact integer counts]
    H --> I[Per-dataset vectors]
    H --> J[Global sum]
    I --> K[Compact EAB or exact CSV]
    J --> K
```

## Installation

Requirements: **Rust 1.88+**, a C/C++ toolchain, and `pkg-config`.
**GGCAT is built into EXPRESSO** as a Rust dependency and installed automatically
with it; no separate GGCAT installation is required.

```bash
git clone https://github.com/Malfoy/Expresso.git
cd Expresso
cargo build --release --locked
export PATH="$PWD/target/release:$PATH"

expresso --help
```

Alternatively, install the executable into Cargo's binary directory:

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --path . --locked
```

The repository's Cargo configuration enables `target-cpu=native` for SIMD.
Build on the machine that will run EXPRESSO, or one with the same CPU instruction
set. `--features no-pdep` is available for older AMD CPUs.

## Quick start

Run the included small example:

```bash
expresso run \
  --exons examples/exons.fa --index example-index --k 7 \
  --fof examples/datasets.tsv --output example-results \
  --threads 4 --jobs 2 --stats

expresso export \
  --input example-results --output example-csv \
  --compression none --with-names
```

For your own data, build once and reuse the index:

```bash
expresso build \
  --exons exons.fa.gz --index exon-index \
  --k 31 --threads 16 --memory-gb 8

expresso quantify \
  --index exon-index --fof datasets.tsv --output results \
  --threads 16 --jobs 4 --mode unitigs

expresso inspect --index exon-index
```

Use `--mode reads` for FASTA/FASTQ reads, `--mode unitigs` for abundance-annotated
unitigs, or `--mode auto` to infer weighting record by record. Explicit modes in
the dataset list override the command-line mode.

Index, result, and export destinations must be **new directories**. Work is
staged beside the destination and published when complete. `--temp-dir /scratch`
relocates index-construction scratch files.

## Inputs

### Reference FASTA

Each FASTA **record** is one target: a header followed by its sequence.
One sequence line per exon is supported; wrapped sequences and CRLF line endings
are also accepted.

```fasta
>exon_1
ACGTTGCAACGTTGCAACGTTGCAACGTTGCA
>exon_2
GGTACCATGGTACCATGGTACCATGGTACCAT
```

Target IDs are 1-based in reference order. The first header token supplies the
name. Duplicate names or identical sequences in separate records remain
separate targets. Short targets and targets without usable unique k-mers remain
in every vector as zeros. A reference without any valid k-mers cannot be indexed.

Transcripts or genes can also be reference records; they are then the units
being quantified. Output metadata retains the names `exon_id` and `exon_name`.

### Dataset list

Supply one absolute or relative path per line:

```text
/data/sample_A.fastq.gz
reads/sample_B.fasta.xz
unitigs/sample_C.fa.zstd
```

Relative paths resolve from the **list's directory**, not the working directory.
Blank lines and lines starting with `#` are ignored. Paths may contain spaces.
Inferred names include a numeric prefix, such as `000001_sample_A.eab.zst`.
Every line is a dataset: listing the same file twice counts it twice in the sum.

For explicit names and mixed data types, use **actual tabs** between columns:

```text
sample_A	/data/sample_A.fastq.gz	reads
sample_B	reads/sample_B.fa.xz	reads
sample_C	unitigs/sample_C.fa.zstd	unitigs
```

The mode column is optional. Names must be unique and contain only ASCII
letters, digits, `_`, `-`, or `.`; `.` and `..` alone are not valid names.

### Compression and record formats

FASTA and conventional four-line FASTQ are supported. Compression is detected
from file contents, including concatenated gzip/xz/zstd streams. Files do not
need a particular extension for `quantify`; `.zst` and `.zstd` both work.

## Counting semantics

| Setting or rule | Behavior |
| --- | --- |
| k-mer length | Odd k from 3 to 63; default 31 |
| Strand | A k-mer and its reverse complement are equivalent |
| Ambiguity | Non-ACGT bases interrupt k-mers; flanking DNA is never joined |
| Shared sequence | K-mers present in multiple target records contribute to none |
| Repeated sequence | Repetitions within one target retain that target as owner |
| `--mode reads` | Every matching occurrence contributes 1 |
| `--mode unitigs` | Every matching occurrence contributes the header abundance |
| `--mode auto` | Use abundance when present; otherwise count as reads |

Lowercase DNA is accepted. Sequence formatting whitespace is removed.
Occurrences are counted independently, without read deduplication or paired-end
fragment correction. Abundance is a **sum of uniquely assigned k-mer
observations**, without length or library-size normalization: it is not a
transcript count, TPM, or RPKM.

Supported mean-abundance tags are `ka:f:`, `km:f:`, `ka:i:`, and `km:i:`.
A tag may be the first token, including a header without an identifier. Mean
tags take precedence over `KC:i:`. With only `KC:i:`, the weight is
`KC / (sequence_length - assembly_k + 1)`; set `--unitig-k` when the assembly k
differs from the index k. Missing abundance fails in `unitigs` mode; malformed
abundance fails in both `unitigs` and `auto` modes.

Unitig weights apply uniformly along each sequence. Mean coverage cannot recover
the original per-k-mer read counts. Generate GGCAT unitigs with abundance support.

Weights are accumulated as fixed-point millionths; extra decimal places round
half up at ingestion. Each target/dataset is rounded to the nearest integer,
ties up, **after summing**. Reads use exact `u64` counts. Unitig/auto mode uses
`u64` millionths, allowing approximately 18.4 trillion weighted observations per
target/dataset. Overflow is reported explicitly. The shared reference table
includes the number of distinct usable k-mers for downstream normalization.

## Outputs and CSV export

### Compact EAB: the default

```text
results/
├── datasets/
│   ├── sample_A.eab.zst
│   └── sample_B.eab.zst
├── global.eab.zst
├── exons.csv.zst
└── manifest.json
```

EAB v1 stores target names **once**, in `exons.csv.zst`. Each vector contains
one code per target in the same order, an integer codebook, a reference SHA-256,
and a CRC checksum. The manifest describes coverage, output files, counting
metrics, and measured approximation errors.

`--bits` accepts **2–16**, default **8**. Eight bits means one payload byte per
target before compression, plus a small header/codebook. Larger values use a
logarithmically spaced codebook fitted to each vector's observed maximum.

- Zero/nonzero status is preserved exactly at every supported width.
- At 8 bits, counts **0–15 are exact**. If the vector maximum is at most 255,
  every count in that vector is exact.
- Larger counts are approximate. More bits generally improve precision; error
  depends on the vector's range and is recorded in the manifest.
- Counts are not clipped to a fixed abundance ceiling.
- Quantization is lossy; the subsequent zstd/gzip/xz compression is lossless.

The design draws on the idea of logarithmic abundance discretization used in
[REINDEER2](https://github.com/Yohan-HernandezCourbevoie/REINDEER2), with its
own encoding and per-vector codebooks. Codes from different vectors must be
decoded before comparing or adding them.

### Global sum and optional statistics

**The global sum is always produced**; no extra flag is required. It sums exact
rounded dataset counts, then independently quantizes that total for compact
output. Consequently, adding decoded dataset values may not reproduce the
decoded global vector exactly.

`--stats` additionally writes mean, median, minimum, maximum, and detected-dataset
count under `statistics/`. Statistics include zero values and are computed from
exact counts before quantization. The detected-dataset count can therefore be
approximate, even though individual vectors preserve presence/absence exactly.
Exact statistics require temporary disk space proportional to targets × datasets.

### Export only what you need

```bash
# Every dataset and the global sum; one abundance column per dataset.
expresso export --input results --output csv --compression zstd

# Selected dataset, plain CSV; repeat --dataset to select more.
expresso export --input results --output selected \
  --dataset sample_A --compression none

# Global sum only, with IDs and names.
expresso export --input results --output global-csv --global-only --with-names
```

Without `--with-names`, dataset CSVs have a single `abundance` column and one
integer per target. Global CSVs also include any requested statistics. With
`--with-names`, `exon_id` and `exon_name` precede abundance. Export copies the
shared reference table and records source coverage in its manifest.

CSV export returns decoded representatives: it cannot restore values discarded
by quantization. To retain exact original counts from the start, quantify with
`--format csv`. This writes `exon_id,exon_name,abundance` in every dataset CSV
and an exact global CSV; repeating names makes these files much larger.

| `--compression` | Compact vectors | CSVs |
| --- | --- | --- |
| `zstd` (default, level 3) | `.eab.zst` | `.csv.zst` |
| `gz` | `.eab.gz` | `.csv.gz` |
| `xz` | `.eab.xz` | `.csv.xz` |
| `none` | `.eab` | `.csv` |

This option is independent of input compression. The shared compact reference
is always `exons.csv.zst`; exact CSV output uses `exons.csv`.

### Failed datasets

Normally, a dataset error stops the run. `--keep-going` records input/counting
failures and continues with the rest. Failed datasets contribute **no partial
counts** to the global sum. The final manifest lists their errors and records
`coverage.complete: false`. CSV export preserves this information.

Output/aggregation errors remain fatal. An all-failed run is not published.
`--keep-going` cannot be combined with `--stats`. Invalid dataset lists and
missing paths are rejected before processing. There is no automatic checkpoint
resume; completed vectors are staged until the overall run is published.

## Validation

```bash
cargo fmt --check
cargo test --locked
cargo clippy --all-targets --locked -- -D warnings
```

Rust tests cover the EAB codecs and bit widths, integer boundaries, integrity
checks, shared-k-mer exclusion, weighted counting, global sums, and failures.
Integration tests use the bundled GGCAT library with an empty `PATH` and compare
counts with an independent canonical-k-mer oracle. No separate tools or Python
scripts are needed to run the tests.

## Upstream projects

- [GGCAT](https://github.com/algbio/ggcat) — simplitig construction.
- [Rust SSHash](https://github.com/COMBINE-lab/sshash-rs) — compressed k-mer dictionary.
- [Helicase](https://github.com/imartayan/helicase) — SIMD FASTA/FASTQ parsing.
- [Logan](https://github.com/IndexThePlanet/Logan) — unitig datasets and abundance headers.
- [REINDEER2](https://github.com/Yohan-HernandezCourbevoie/REINDEER2) — logarithmic abundance discretization inspiration.

`Cargo.lock` pins dependency versions. GGCAT 2.1.0 is pinned to an upstream Git
revision and compiled into EXPRESSO; SSHash 0.7.1 and Helicase 0.2.0 are also
pinned explicitly in `Cargo.toml`. EXPRESSO handles xz through its shared liblzma
dependency before passing the stream to Helicase; gzip/zstd use Helicase's input
layer. EAB is EXPRESSO's own format and is not REINDEER2-compatible.
