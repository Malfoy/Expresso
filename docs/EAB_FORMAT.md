# EXPRESSO Abundance Binary (EAB), version 1

EAB stores a dataset's abundance vector without repeating target IDs or names.
The reference order is shared across the result set. Counts are accumulated and
rounded as `u64` exactly as in CSV mode, then quantized only when written.

The design draws on [REINDEER2's logarithmic abundance discretization](https://github.com/Yohan-HernandezCourbevoie/REINDEER2/blob/5bcc20f04f95ba94c801c2cdec9e00810a6f4540/src/reindeer2/mod.rs#L1586).
EAB defines its own encoding: an exact low-count region followed by a geometric
codebook, nearest-integer-representative encoding, and a maximum selected from
each vector's actual `u64` counts. No REINDEER2 file-format compatibility is implied.

## Result directory

```text
manifest.json
exons.csv.zst
datasets/sample_A.eab.zst
datasets/sample_B.eab.zst
global.eab.zst
statistics/mean.eab.zst                 # only with --stats
statistics/median.eab.zst
statistics/min.eab.zst
statistics/max.eab.zst
statistics/datasets_detected.eab.zst
```

The manifest records `format: "compact"`, `compact_format_version: 1`, the
reference filename, SHA-256 and target count, dataset filenames and metrics,
and quantization/error measurements for every vector. Target row IDs are
implicit, 1-based positions in the reference table. The reference table has
`exon_id,exon_name,length,unique_kmers` and is stored **once**, always with zstd.

`--compression zstd` (default), `gz`, `xz`, or `none` selects the outer lossless
codec for vectors: `.eab.zst`, `.eab.gz`, `.eab.xz`, or `.eab`. Compression is
detected from magic bytes on read. The manifest itself is plain JSON. A full
storage total must include vectors, the shared reference, and the manifest.

## Quantization

`--bits B` accepts 2–16 bits, default 8. Let `C = 2^B - 1` and let `M` be the
largest exact count in the vector.

* When `M <= C`, code `x` represents exactly `x`; store the table `0..M`.
* Otherwise preserve `0..P` exactly, where `P = min(15, floor(C/4))`.
* Set `A = P+1` and `b = (M/A)^(1/(C-A))`.
* For codes `c=A..C`, form `round(A*b^(c-A))`, constrain it to be at least the
  preceding representative plus one and at most `M-(C-c)`, and store it.
* The final representative is exactly the original integer `M`, avoiding
  floating-point loss at the endpoint, including `u64::MAX`.
* Encode an exact count using the nearest stored representative in **absolute
  count distance**, resolving ties toward the lower representative.

The integer codebook is persisted; decoding never evaluates floating-point
logarithms or exponentials. Stored representatives are strictly increasing.
The default preserves counts 0–15 exactly. Other counts may also be exact.
All-zero vectors have a one-entry codebook containing zero. **Zero and nonzero
status are always preserved. No abundance is clipped to a fixed ceiling.**

Each vector has its own base, range, and codebook. Codes from different vectors
cannot be added or compared as abundances without decoding them first.

Bit width controls storage and precision, not a universal error percentage.
The manifest records the actual maximum/mean absolute error, maximum/mean
relative error over **nonzero counts**, exact-value count, nonzero count,
codebook size, base and exact prefix. At wide dynamic ranges, fewer bits incur
larger errors. CSV export returns the integer representatives; original counts
cannot generally be recovered from a lossy vector.

## Uncompressed EAB byte layout

All integers are unsigned little-endian. Outer compression encloses the entire
following byte sequence.

| Offset | Bytes | Field |
| ---: | ---: | --- |
| 0 | 8 | ASCII magic/version `EXPRAB01` |
| 8 | 1 | Bit width, 2–16 |
| 9 | 1 | Encoding, currently 0 (integer codebook) |
| 10 | 2 | Reserved, must be zero |
| 12 | 8 | Number of targets `N` |
| 20 | 4 | Number of codebook entries `K`, 1..2^B |
| 24 | 8 | Packed payload bytes, `ceil(N*B/8)` |
| 32 | 8 | Original maximum count `M` |
| 40 | 32 | Shared-reference SHA-256 |
| 72 | 8*K | Integer codebook, starting with zero and ending with M |
| 72+8*K | ceil(N*B/8) | Packed codes in reference order |
| end-4 | 4 | CRC-32/ISO-HDLC of every preceding uncompressed byte |

Codes are packed least-significant bit first. Code `i` starts at bit offset
`i*B`. Unused high bits in the last byte must be zero. For example, the codes
0–7 at 3 bits pack to hexadecimal `88 c6 fa`.

The reference digest covers these concatenated bytes:

1. ASCII `EXPRESSO-reference-v1` followed by one zero byte.
2. Target count as a little-endian u64.
3. For every target in order: UTF-8 name length as u64, the exact UTF-8 name,
   sequence length as u64, and unique k-mer count as u64.

The reader validates the version, reserved fields, target count/reference hash,
codebook ordering and bounds, payload size, code ranges, padding, checksum and
end of stream. Parsing is bounded by the reference count and bit width.

An 8-bit vector uses `N` payload bytes and at most **2,124 bytes** of header,
codebook and checksum overhead before compression. At 533,740 targets and 874
vectors this is at most 468,345,136 bytes, plus the shared reference and manifest.

## Global and statistical semantics

Global sums are computed from the exact rounded per-dataset counts, then
quantized separately. Optional mean, median, min, max and datasets-detected
vectors are also computed from exact counts (including zeros), then separately
quantized at the selected width. Therefore decoded global/statistical values
can differ from aggregates recomputed from decoded dataset vectors. The
datasets-detected statistic can also be approximate; per-target presence in
each individual dataset remains exact.

With `--keep-going`, only successfully completed datasets contribute to the
global sum. Failed datasets contribute no partial counts and have no vector in
the `datasets` array. The manifest's `failed_datasets` array records their names,
paths, and errors; `coverage` reports requested/completed/failed counts and a
`complete` flag. An all-failed run is not published. This mode cannot be combined
with `--stats`. CSV export retains the original coverage and failure information
as `source_coverage` and `source_failed_datasets` in its own manifest.

## CSV export

```sh
expresso quantify --index index --fof datasets.tsv --output compact --bits 8
expresso export --input compact --output csv --compression zstd
expresso export --input compact --output selected --dataset sample_A --compression none
expresso export --input compact --output named --with-names
expresso export --input compact --output global-only --global-only
```

Default dataset CSVs contain only `abundance`, one integer per target in shared
reference order. Global CSVs also contain requested statistical columns.
`--with-names` adds `exon_id,exon_name` explicitly. Export copies the shared
reference table once and writes a manifest marking the decoded values as
potentially approximate. Selection and writing are transactional: the output
directory must be new and is published only after successful decoding.

`--format csv` on `quantify` or `run` retains the legacy exact CSV output with
names. Choose it when exact original counts must be retained.
