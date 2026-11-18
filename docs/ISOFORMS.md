# Infer transcript expression from exon/intersection counts

[README](../README.md) · [Viewer guide](VIEWER.md) · [CLI recipes](CLI.md)

This mode fits **aggregate k-mer feature abundances**. It does not assign reads
or fragments to transcripts. Build a dedicated `--isoforms` index, obtain exact
feature counts, then run `infer-isoforms` independently of sequence processing.
An existing ordinary exon/gene result is insufficient: those modes discarded
shared k-mers, and that missing evidence cannot be reconstructed from their totals.

## Run the exon-skipping example

The small [example](../examples/isoforms/annotation.gtf) has one gene and two
transcripts: `T_WITH_B = A–B–C` and `T_SKIP_B = A–C`. Its weighted unitigs supply
coverage 6 and 4 respectively. Run these commands from the repository root:

```bash
expresso build --isoforms \
  --gtf examples/isoforms/annotation.gtf \
  --genome examples/isoforms/genome.fa \
  --index isoform-demo-index --k 7 --threads 4 --memory-gb 1

expresso quantify \
  --index isoform-demo-index --fof examples/isoforms/datasets.tsv \
  --output isoform-demo-features --mode unitigs --threads 4

expresso infer-isoforms \
  --index isoform-demo-index --input isoform-demo-features \
  --output isoform-demo-results --threads 4 --compression none

cat isoform-demo-results/datasets/sample.csv
```

The inferred `abundance` values should be approximately **6** for `T_WITH_B`
and **4** for `T_SKIP_B`, with `within_gene_fraction` approximately **0.6** and
**0.4**. Both have `status=identifiable`. Floating-point tolerances affect the
last digits. [expected.csv](../examples/isoforms/expected.csv) records these values.
The example is synthetic; no network download is needed.

To explore those transcript results in the viewer:

```bash
expresso pack --input isoform-demo-results --output isoform-demo-results/transcripts.eai
expresso view
```

Choose the generated file, or open the ready-to-use
[transcript toy index](../examples/isoforms/toy.eai). Open **Heatmap** to compare
both transcripts in `sample`. **Transcript metric** switches between assigned k-mer
observations (approximately 272 and 618) and EM abundance (approximately 4 and 6).
The additional observable position in `T_WITH_B` is a junction occurrence of an
already indexed exonic k-mer, so its observable length is 103 rather than 102.
Use the **Transcript list** tab to summarize results across datasets. See the
[viewer guide](VIEWER.md#summarize-a-gene-exon-or-transcript).

For real data replace the annotation, matching genome, and dataset list, and
normally use `--k 31`. Genome/GTF and dataset compression support is the same as
for other modes. `run --isoforms` can combine index construction and exact
feature counting; EM remains a separate `infer-isoforms` invocation.

## Features and the exposure model

Isoform mode preserves **original distinct GTF exon intervals and each
transcript's original exon path**. It does not use the gene mode's merged exon
regions. Duplicate GTF rows are deduplicated; alternative exon boundaries remain
distinct. Paths follow genomic coordinates in strand order, with reverse
complementation on the negative strand. Every exon row must have exactly one
`gene_id` and `transcript_id`. Overlap within one transcript is rejected.

All distinct canonical k-mers occurring inside annotated exons are retained,
including shared and repeated sequence. Each k-mer belongs to exactly **one
disjoint feature** identified by its exon membership/multiplicities and its
transcript occurrence/multiplicities. This distinguishes exclusive sequence,
pairwise intersections, higher-order intersections, and repeated sequences;
observations are never duplicated across overlapping feature totals. Sharing
between genes is retained too, so coupled genes are fitted together.

For feature `f` and transcript `t`, exposure `A[f,t]` is the number of positions
in the annotated transcript producing a k-mer in that feature. It equals the
feature's distinct k-mer count times the transcript multiplicity. Full annotated
transcript sequences are scanned when building this matrix. Novel junction
k-mers are not added, but a junction occurrence matching an already indexed
exonic k-mer is included in its exposure. `--junctions` expansion and
`--level gene` cannot be combined with `--isoforms`.

Reference membership collection uses parallel SSHash queries, bounded input
queues, and a disk-backed sort with bounded merge fan-in. `--temp-dir` places
scratch files on a chosen disk; allow space for exon and transcript occurrence
runs. `--memory-gb` remains a construction hint, not a hard process memory cap.
The dictionary, annotations, and final sparse exposure model also require RAM.

The saved index contains:

| File | Contents |
|---|---|
| `isoform-model.json` | Sparse memberships, multiplicities, transcript paths, observable lengths |
| `features.csv.zst` | Feature sizes and exon/transcript memberships with aligned multiplicities |
| `original-exons.csv.zst` | Original interval coordinates and annotation exon IDs |
| `transcripts.csv.zst` | Transcript/gene IDs, lengths, and ordered exon paths |
| `reference-exons.fa`, `reference-transcripts.fa` | Reconstructed reference sequences |
| `metadata.json`, SSHash files, `owners.bin` | Model checksum and reusable feature-counting index |

CSV feature, exon, transcript, group, and component IDs are **1-based**. Internal
IDs in `isoform-model.json` and the EFC binary format are **0-based**. Feature IDs
are local to an index: use its tables, rather than transferring IDs between builds.

## Exact counts and count-only input

Feature counting defaults to compressed **EFC v1**, a sparse exact fixed-point
format with a model SHA-256 binding and CRC32 integrity check. Reads contribute
integer observations; weighted unitigs retain millionths without target rounding
or logarithmic quantization. Individual datasets use `u64` accumulation and
the global feature sum uses `u128`. `--bits` has no effect in this mode.
`--compression zstd|gz|xz|none` selects compression; zstd is the default.

Use `quantify --format csv` to write exact sparse feature CSV instead. The
header is `feature_id,abundance`, with 1-based IDs and omitted zeros. Results
include a manifest, `datasets/<name>.efc.zst` (or CSV), and `global.efc.zst`.
`--keep-going` excludes failed datasets completely and records exclusions.
`--stats` and `--viewer-index` are unavailable for isoform feature results;
ordinary EAB/EAI export and the exon/gene viewer use a different result format.

To fit counts produced elsewhere, provide a list of disjoint-feature CSV or
EFC files, using the feature IDs from this index:

```text
sample_A<TAB>/data/sample_A.features.csv.gz
sample_B<TAB>/data/sample_B.efc.zst
```

Here `<TAB>` means an actual tab. Relative paths resolve from the list's directory.

```bash
expresso infer-isoforms --index isoform-index --fof feature-files.tsv \
  --output isoform-results --threads 16 --jobs 4
```

CSV rows may occur in any order; missing rows mean zero. Duplicate IDs, IDs
outside the model, negative/nonfinite values, and incompatible binary files fail
without publishing results. CSV carries no model checksum: the caller must use
the correct index's feature IDs. **Inclusive exon totals or overlapping
intersection totals are not independent disjoint features** and cannot be
substituted directly. Normalized or log-quantized counts are not fitting input.

`infer-isoforms` reads only index metadata, the saved model, and aggregate count
files. It does not open the SSHash dictionary, reference sequences, or original
datasets. The sparse model and up to `--jobs` feature-count vectors are held in
memory, rather than the entire dataset collection. Fitting runs in parallel over
datasets and independent connected components within the total `--threads` budget.

## EM and result columns

For observed disjoint feature counts `y[f]`, the working mean model is
`lambda[f] = sum_t A[f,t] × theta[t]`. EM apportions each feature count in
proportion to `A[f,t] × theta[t]`, then updates each rate to its assigned
observations divided by its observable k-mer length `L[t] = sum_f A[f,t]`.
Lengths include repeated positions. No per-dataset CPM or length normalization
is applied to input; the exposure matrix handles length during fitting.

`theta`, reported as **abundance**, is an inferred k-mer coverage rate. It is not
a molecule or fragment count. `within_gene_fraction` divides identifiable
individual rates by the gene's total resolved group rate.
`relative_abundance_ppm` divides rates by the total observable group rate and
multiplies by one million; it is not fragment TPM. These fractions are blank
when their denominator cannot be resolved.

Identical exposure columns are collapsed to one group whose rate is the sum of
its members' rates. A numerical rank check then detects broader dependencies,
including ambiguity between distant alternative exon choices. The default
`--rank-limit 512` checks components with at most 512 groups; larger components
are conservatively marked `rank_unchecked`. Raise this limit to check them,
at potentially greater CPU/memory cost. Numerical rank uses maximum-scaled
columns and tolerance `1e-10`.

| Status | Meaning |
|---|---|
| `identifiable` | One transcript in a full-rank, converged component |
| `group_only` | Only the combined abundance of identical profiles is resolved |
| `unresolved` | Exposure columns remain rank deficient after grouping |
| `rank_unchecked` | Component exceeds the chosen rank-check limit |
| `unobservable` | Transcript has no indexed k-mer exposure |
| `zero_evidence` | Observable component has no counts |
| `non_converged` | EM did not satisfy the fixed-point tolerance |

Individual `abundance` is blank for grouped, unresolved, unchecked, unobservable,
or non-converged transcripts. `model_estimate` is a diagnostic solver value,
not an identified abundance; unresolved allocations can depend on initialization.
Group tables expose resolved combined rates. The default EM limit is 10,000
iterations with relative tolerance `1e-8`; adjust `--max-iterations` and
`--tolerance` if needed.

Each result directory contains transcript CSVs under `datasets/`, group CSVs
under `groups/datasets/`, and convergence/rank/residual diagnostics under
`diagnostics/datasets/`, plus corresponding `global.csv[.zst]` files and a
manifest. Compression defaults to zstd. The global result **sums independently
fitted dataset rates and assigned observations**; its fractions/ppm are computed
from those sums. Global diagnostics summarize the separate fits, rather than a
refit of pooled feature counts. Source dataset exclusions remain in the manifest.

The Poisson objective is a working model: adjacent k-mer observations are
correlated and unitig mean coverage does not recover original per-k-mer counts.
Uniform coverage and a sufficiently complete annotation are assumptions.
Diagnostics report relative L1 model residuals; there are no confidence intervals,
and full rank alone does not establish strong biological evidence for an isoform.
