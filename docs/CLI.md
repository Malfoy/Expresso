# Command-line recipes

[README](../README.md) · [Viewer guide](VIEWER.md) · [Isoforms](ISOFORMS.md)

Run commands from the repository root after installing EXPRESSO.
Use new output directory names for each analysis.

## Choose an input workflow

| Reference | Build options |
|---|---|
| Exon FASTA | `--exons exons.fa.gz` |
| Exons from annotation | `--level exon --gtf annotation.gtf.gz --genome genome.fa.gz` |
| Genes from annotation | `--level gene --gtf annotation.gtf.gz --genome genome.fa.gz` |
| Gene FASTA | `--level gene --genes genes.fa.gz` |
| Isoform features | `--isoforms --gtf annotation.gtf.gz --genome genome.fa.gz` |

Choose the genome assembly matching the annotation and use the same contig names.
`--k 31` is the default; the included toy examples use `--k 7`.

## Build once, quantify several collections

```bash
expresso build --exons exons.fa.gz --index exon-index \
  --k 31 --threads 16 --memory-gb 8

expresso quantify --index exon-index --fof batch-1.tsv \
  --output batch-1-results --threads 16 --jobs 4 --viewer-index

expresso quantify --index exon-index --fof batch-2.tsv \
  --output batch-2-results --threads 16 --jobs 4 --viewer-index
```

Inspect the reference and k-mer counts:

```bash
expresso inspect --index exon-index
```

## Make a dataset list

### Paths only

```bash
cat > reads.txt <<'LIST'
/data/sample_A.fastq.gz
/data/sample_B.fastq.xz
/data/sample_C.fasta.zstd
LIST

expresso quantify --index exon-index --fof reads.txt \
  --output read-results --mode reads --threads 16 --jobs 4
```

Names are generated from each path, for example `000001_sample_A`.
Relative paths resolve from the list's directory. Blank lines and `#` comments
are accepted.

### Named samples and mixed input types

`printf` writes the actual tab separators:

```bash
printf '%s\t%s\t%s\n' \
  sample_A /data/sample_A.fastq.gz reads \
  sample_B /data/sample_B.fasta.xz reads \
  sample_C /data/sample_C.fa.zstd unitigs > mixed.tsv

expresso quantify --index exon-index --fof mixed.tsv \
  --output mixed-results --threads 16 --jobs 4 --viewer-index
```

Names can contain letters, digits, `_`, `-`, and `.`.
A mode in the list overrides `--mode` for that sample.

## Quantify Logan or GGCAT unitigs

```bash
cat > unitigs.txt <<'LIST'
/data/sample_A.u.fa.zstd
/data/sample_B.u.fa.zstd
LIST

expresso quantify --index exon-index --fof unitigs.txt \
  --output unitig-results --mode unitigs \
  --threads 16 --jobs 4 --viewer-index
```

| Mode | Contribution of each matching k-mer occurrence |
|---|---|
| `reads` | 1 |
| `unitigs` | Header abundance |
| `auto` | Header abundance when present, otherwise 1 |

Supported abundance tags include `ka:f:`, `km:f:`, `ka:i:`, `km:i:`, and `KC:i:`.
For `KC:i:` unitigs assembled at a different k, add `--unitig-k`, for example:

```bash
expresso quantify --index exon-index --fof unitigs.txt \
  --output unitig-k63-results --mode unitigs --unitig-k 63 --threads 16
```

[Weight calculation and counting rules](IMPLEMENTATION.md#counting-rules).

## Extract an exon FASTA

```bash
expresso extract-exons --gtf annotation.gtf.gz --genome genome.fa.gz \
  --output exons.fa
```

For multiple genome FASTAs, repeat `--genome`:

```bash
expresso extract-exons --gtf annotation.gtf.gz \
  --genome chr1.fa.gz --genome chr2.fa.gz \
  --output exons.fa.zst --compression zstd
```

The FASTA has one header and one sequence line per exon.
You can also pass the GTF/genome directly to `build` or `run`.

## Build a gene reference

```bash
expresso build --level gene \
  --gtf annotation.gtf.gz --genome genome.fa.gz \
  --index gene-index --k 31 --threads 16
```

Genes use annotation `gene_id` values. Overlapping exon intervals within each
gene are merged. The index retains `reference-genes.fa` and its gene table.

For supplied gene FASTA, use one record per contiguous sequence; records sharing
the same first header token are grouped into one gene:

```bash
expresso build --level gene --genes genes.fa.gz \
  --index gene-fasta-index --k 31 --threads 16
```

[Add junctions](../README.md#include-exon-junctions-in-gene-quantification).

## Choose output precision and compression

```bash
# More compact quantized counts.
expresso quantify --index exon-index --fof reads.txt \
  --output compact-results --bits 8 --compression zstd

# Higher precision.
expresso quantify --index exon-index --fof reads.txt \
  --output precise-results --bits 16 --compression zstd

# Exact counts written directly as CSV.
expresso quantify --index exon-index --fof reads.txt \
  --output exact-csv-results --format csv --compression gz
```

| Option | Choices | Default |
|---|---|---|
| `--bits` | 2–16 | 8 |
| `--compression` | `zstd`, `gz`, `xz`, `none` | `zstd` |
| `--format` | `compact`, `csv` | `compact` |

Exporting a compact result decodes the stored quantized counts.
Choose exact CSV output during quantification to retain exact integer counts.

## Add statistics or continue past failed datasets

Add mean, median, min, max, and detected-dataset counts:

```bash
expresso quantify --index exon-index --fof reads.txt \
  --output statistics-results --stats --threads 16 --jobs 4
```

Record input errors and continue with successful samples:

```bash
expresso quantify --index exon-index --fof reads.txt \
  --output tolerant-results --keep-going --threads 16 --jobs 4
```

Use either `--stats` or `--keep-going` in a run. The manifest records failures;
failed datasets contribute no partial counts to the global sum.

## Export target and dataset lists

The commands here use the demo built in the [README](../README.md#run-your-first-analysis).

```bash
printf '%s\n' GENE_A GENE_B > selected-genes.txt
printf '%s\n' sample_A1 sample_B1 > selected-datasets.txt

expresso export --input viewer-demo-results --output list-selection \
  --target-list selected-genes.txt --dataset-list selected-datasets.txt \
  --min-value 100 --max-value 10000 --compression zstd
```

| Selection | Options |
|---|---|
| Target names | Repeat `--target NAME`; aliases: `--exon`, `--gene` |
| Target name list | `--target-list FILE` |
| Original numeric target IDs | Repeat `--target-id ID` |
| Dataset names | Repeat `--dataset NAME` |
| Dataset name list | `--dataset-list FILE` |
| Inclusive count range | `--min-value N --max-value N` |
| Global vector only | `--global-only` |

Lists contain one exact name per line. Repeated flags and lists are combined.
Target/value filters automatically include IDs and names in the CSV.
Value filters apply separately to each dataset's raw decoded counts.

### Export directly from the viewer file

```bash
expresso export --input examples/viewer/toy.eai --output packed-selection \
  --dataset sample_A1 --gene GENE_A --compression none
```

### Normalize exported values

```bash
expresso export --input viewer-demo-results --output length-adjusted \
  --dataset sample_A1 --normalization length-cpm \
  --with-names --compression none
```

`cpm` adjusts for full dataset k-mer abundance; `coverage` adjusts for target
length in k-mers; `length-cpm` does both. [Formulas](IMPLEMENTATION.md#normalization).

## Command help

```bash
expresso extract-exons --help
expresso build --help
expresso quantify --help
expresso pack --help
expresso view --help
expresso export --help
expresso infer-isoforms --help
```
