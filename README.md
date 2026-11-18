# EXPRESSO

**EXon-level RNA EXPRESsion quantificatiOn**

Count exon or gene expression across reads and abundance-annotated unitigs.
Explore the results in a WebAssembly viewer, select samples and targets,
and export the tables you need.

[Try the viewer](#try-the-webassembly-viewer) · [Install](#install-the-command-line-tool) · [First analysis](#run-your-first-analysis) · [Your data](#analyze-your-own-data) · [Export CSV](#export-results) · [Advanced analyses](#advanced-analyses)

## Try the WebAssembly viewer

### 1. Start the viewer — no compilation needed

```bash
git clone https://github.com/Malfoy/Expresso.git
cd Expresso
python3 -m http.server 8765 --bind 127.0.0.1 --directory viewer
```

Open **[http://127.0.0.1:8765](http://127.0.0.1:8765)**. Keep the terminal running.

Already installed EXPRESSO? Start the same viewer with:

```bash
expresso view --port 8765
```

### 2. Open the included toy dataset

Click **Choose index** and select **`examples/viewer/toy.eai`** in your checkout.
It contains **6 genes × 4 samples**; no dataset download or quantification is needed.

The **Global overview** opens first. Check these values:

| Total expression | Nonzero cells | Median |
|---:|---:|---:|
| 19,278 | 16 of 24 | 85 |

### 3. Choose what you want to analyze

| Tab | Use it to… |
|---|---|
| **Global overview** | Summarize the whole matrix: totals, detection, mean, median, and percentiles |
| **Dataset list** | Compare samples; click a sample name to inspect it |
| **Gene list** | Compare genes across samples; click a gene name to inspect it |
| **Heatmap** | Compare selected genes and samples using colors, counts, or fold changes |
| **Dataset focus** | See one sample's most expressed genes and their percentages |
| **Gene focus** | See one gene's expression across samples, distribution, and fold changes |

For exon or transcript indexes, **Gene list/focus** becomes **Exon** or **Transcript list/focus**.

**Try this:** open **Dataset list**, then click **sample_A1**.
In **Dataset focus**, `GENE_A` contributes **4,080 of 4,862**, or **83.92%** of the sample's gene signal.

### 4. Explore the heatmap

Click **Heatmap**. All six genes and four samples are already selected:

| Gene | sample_A1 | sample_A2 | sample_B1 | sample_B2 |
|---|---:|---:|---:|---:|
| GENE_A | 4,080 | 3,060 | 68 | 102 |
| GENE_B | 68 | 102 | 3,740 | 4,420 |
| GENE_C | 0 | 0 | 340 | 510 |
| GENE_D | 34 | 0 | 0 | 34 |
| HOUSEKEEPING | 680 | 680 | 680 | 680 |
| GENE_EMPTY | 0 | 0 | 0 | 0 |

Try these controls:

| Goal | What to select |
|---|---|
| Remove absent genes | **Hide all-zero rows**, then **Apply view options**; `GENE_EMPTY` disappears |
| Find the strongest signals | **Sort targets → Total abundance · highest first**; `GENE_B` moves to the top |
| Compare relative expression | **Cell values → Fold change** or **Log₂ fold change**, then **Update heatmap** |
| Adjust for dataset size and target length | **Normalization → CPM / target k-mer length**, then **Update heatmap** |
| Shorten large numbers | **Number format → Compact** or **Scientific** |
| Change the palette | **Main color**; blue is the default |

Update buttons light up when parameters need applying and become dimmed when
results are current. Computations lasting over 0.5 seconds show a progress bar.

### 5. Select genes and samples, then export

Paste this into **Paste target names or #IDs**, then click **Use target list**:

```text
GENE_A
GENE_B
```

Paste this into **Paste dataset names**, then click **Use dataset list**:

```text
sample_A1
sample_B1
```

Click **Update heatmap**, then **Visible CSV** to download the 2 × 2 comparison.
Use **Full selection CSV** to export a selection spanning multiple pages.

The [viewer guide](docs/VIEWER.md) has recipes for detection thresholds,
percentiles, top-gene plots, fold changes, and larger indexes.

## Install the command-line tool

Requires **Rust 1.89+**, a C/C++ toolchain, and `pkg-config`.
**GGCAT is installed with EXPRESSO.**

From the cloned repository:

```bash
cargo build --release --locked
export PATH="$PWD/target/release:$PATH"
expresso --help
```

To install into your Cargo binary directory instead:

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --path . --locked
```

## Run your first analysis

Recreate the viewer demo from the included FASTA files:

```bash
expresso run --level gene --genes examples/viewer/genes.fa \
  --fof examples/viewer/datasets.tsv \
  --index viewer-demo-index --output viewer-demo-results \
  --k 7 --threads 4 --jobs 2 --bits 16 --stats --viewer-index
```

Open **`viewer-demo-results/abundance.eai`** in the viewer to explore the results.
This reproduces the counts shown above. Use new index/output directory names
when repeating an analysis.

## Analyze your own data

### 1. Make a dataset list

This command creates a tab-separated file with **sample name, path, and mode**:

```bash
printf '%s\t%s\t%s\n' \
  sample_A /data/sample_A.fastq.gz reads \
  sample_B /data/sample_B.fastq.gz reads > datasets.tsv
```

Use `reads` for reads, `unitigs` for abundance-annotated Logan/GGCAT unitigs,
or `auto` to select weighting from each header. Input files can be plain,
gzip, xz, or Zstandard. Relative paths start from the dataset list's directory.

### 2. Build a gene index

Use a GTF annotation and its matching genome FASTA:

```bash
expresso build --level gene \
  --gtf annotation.gtf.gz --genome genome.fa.gz \
  --index gene-index --k 31 --threads 16
```

For an **exon index**, use `--level exon` and a different index directory.
For an existing exon FASTA, replace `--gtf` and `--genome` with `--exons exons.fa`.

### 3. Quantify your samples

```bash
expresso quantify --index gene-index --fof datasets.tsv \
  --output gene-results --threads 16 --jobs 4 --viewer-index
```

Reuse `gene-index` for later sample collections. This example shares 16 workers
across four concurrent datasets.

### 4. Open your results

In the viewer, choose **`gene-results/abundance.eai`**.

For results created without `--viewer-index`, package them for the viewer:

```bash
expresso pack --input gene-results --output gene-results/abundance.eai
```

### What gets written?

```text
gene-results/
├── abundance.eai           # Open in the viewer (--viewer-index)
├── datasets/               # One compact abundance vector per sample
├── global.eab.zst          # Sum across samples
├── genes.csv.zst           # Shared gene names and reference information
└── manifest.json           # Sample names, outputs, and counting metrics
```

Default output uses **8 bits per abundance + Zstandard**.
Choose `--bits 16` for more precision, or `--format csv` during quantification
for exact CSV counts. [More input and output recipes →](docs/CLI.md)

## Export results

The commands below use the runnable `viewer-demo-results` example.
Replace it with your own result directory or `.eai` file for other analyses.

### Export all samples and the global sum

```bash
expresso export --input viewer-demo-results --output viewer-demo-csv \
  --with-names --compression none
```

### Export selected genes from one sample

```bash
expresso export --input viewer-demo-results --output selected-genes \
  --dataset sample_A1 --target GENE_A --target GENE_B \
  --with-names --compression none
```

### Export normalized expression

```bash
expresso export --input viewer-demo-results --output normalized-genes \
  --dataset sample_A1 --normalization length-cpm \
  --with-names --compression none
```

| Normalization | What it does |
|---|---|
| `none` | Raw k-mer abundance |
| `cpm` | Divide by all input k-mer observations, then multiply by one million |
| `coverage` | Divide by target length in k-mer positions |
| `length-cpm` | Apply both dataset-size and target-length corrections |

### Export the global sum or a count range

```bash
expresso export --input viewer-demo-results --output global-counts \
  --global-only --with-names --compression none

expresso export --input viewer-demo-results --output filtered-counts \
  --dataset sample_B1 --min-value 100 --max-value 10000 \
  --with-names --compression none
```

[Export lists of targets/samples, compression choices, and statistics →](docs/CLI.md#export-target-and-dataset-lists)

## Advanced analyses

### Include exon junctions in gene quantification

Add junctions documented by the annotation:

```bash
expresso build --level gene \
  --gtf annotation.gtf.gz --genome genome.fa.gz \
  --junctions annotated --index annotated-gene-index --threads 16
```

Include all ordered exon pairs, in both directions, for genes with up to
100 merged exons:

```bash
expresso build --level gene \
  --gtf annotation.gtf.gz --genome genome.fa.gz \
  --junctions all --junction-max-exons 100 \
  --index all-junction-gene-index --threads 16
```

[Exon merging, threshold behavior, and junction definitions →](docs/IMPLEMENTATION.md#junction-construction)

### Infer transcript expression with EM

Try the included exon-skipping example:

```bash
expresso build --isoforms \
  --gtf examples/isoforms/annotation.gtf --genome examples/isoforms/genome.fa \
  --index isoform-demo-index --k 7 --threads 4 --memory-gb 1

expresso quantify --index isoform-demo-index \
  --fof examples/isoforms/datasets.tsv --output isoform-demo-features \
  --mode unitigs --threads 4

expresso infer-isoforms --index isoform-demo-index \
  --input isoform-demo-features --output isoform-demo-results \
  --threads 4 --compression none

expresso pack --input isoform-demo-results \
  --output isoform-demo-results/transcripts.eai
```

Expected EM abundances: **6 for `T_WITH_B`**, **4 for `T_SKIP_B`**.
Open `transcripts.eai`, then choose **EM abundance** under **Transcript metric**.

[Count-only input, output columns, and the EM model →](docs/ISOFORMS.md)

## More help

| Guide | Contents |
|---|---|
| [Viewer](docs/VIEWER.md) | Detailed walkthroughs for every tab, selection, plotting, normalization, and export |
| [CLI recipes](docs/CLI.md) | Reads/unitigs, exon extraction, sample lists, statistics, CSV filters, and compression |
| [Isoforms](docs/ISOFORMS.md) | Exon/intersection features, count-only fitting, EM results, and diagnostics |
| [Implementation reference](docs/IMPLEMENTATION.md) | Counting rules, reference construction, normalization formulas, and file semantics |
| [EAB format](docs/EAB_FORMAT.md) | Compact abundance encoding and binary layout |

```bash
expresso build --help
expresso quantify --help
expresso export --help
expresso infer-isoforms --help
```

Built with [GGCAT simplitigs](https://github.com/algbio/ggcat),
the [Rust SSHash dictionary](https://github.com/COMBINE-lab/sshash-rs), and
[Helicase FASTA/FASTQ parsing](https://github.com/imartayan/helicase).
