# WebAssembly viewer guide

[README](../README.md) · [CLI recipes](CLI.md) · [Isoforms](ISOFORMS.md)

## Start and open a file

From the repository root, choose either launcher:

```bash
# Use the bundled static application; Python 3 is enough.
python3 -m http.server 8765 --bind 127.0.0.1 --directory viewer
```

```bash
# Or use the installed EXPRESSO executable.
expresso view --port 8765
```

1. Open **http://127.0.0.1:8765** in Firefox or Chromium.
2. Click **Choose index**.
3. Select `examples/viewer/toy.eai`, or your own `.eai` file.

Keep the server running; Ctrl+C stops it. To use another port, change `8765`
in both the command and the browser URL.

## Start with the whole matrix

**Global overview** covers every target and real dataset in the file.
The toy opens with these results:

| Targets | Datasets | Cells | Total | Nonzero cells | Mean | Median |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 4 | 24 | 19,278 | 16 | 803.25 | 85 |

To count cells above a threshold, expand **Global summary options**, set
**Cells with expression above X**, and click **Update summary**.
Change **Percentiles** to choose the quantiles you want.

Sidebar selections affect the other views; the global overview keeps the
whole-matrix scope. Stored global/statistic vectors are excluded from summaries.

## Compare datasets

Open **Dataset list**. Each row summarizes one sample across targets.

1. Keep **Dataset scope → Selected real datasets** and **Target scope → Selected targets**.
2. Compare **Total**, **Expressed targets**, **Mean**, **Median**, and percentiles.
3. Set **Sort by → Total** and **Descending** to put the highest-signal sample first.
4. Click a sample name to open **Dataset focus**.

For the toy, `sample_A1` has total **4,862**, **4** expressed genes, and median **51**.

For a large index, **Update page** calculates 25 datasets at a time.
**Calculate all datasets** fills the selected scope and enables statistical
sorting across it. **Calculated rows CSV** exports completed rows.

## Plot the most expressed genes in one dataset

In **Dataset focus**, choose a dataset by its position, arrows, or exact name.

1. Open **Dataset summary options**.
2. Set **Top N** to `2`.
3. Keep **Rank candidates → All targets** and **Percentage denominator → All target signal**.
4. Click **Update dataset summary**.

For `sample_A1`, the top bars are:

| Gene | Expression | Share of all gene signal |
|---|---:|---:|
| GENE_A | 4,080 | 83.92% |
| HOUSEKEEPING | 680 | 13.99% |

**Show other targets as a bar** adds the remainder. Restricting candidates to
selected genes keeps the denominator over all targets.

Use **Bar scale → Fixed 0–100%** to compare composition across datasets.
Download **Ranking CSV**, **Dataset statistics CSV**, or **Chart SVG**.

## Summarize a gene, exon, or transcript

Open **Gene list**; exon/transcript indexes use the corresponding tab name.

1. Set **Datasets to summarize** to all datasets or your selection.
2. Set **Count above X** to `100`.
3. Choose the statistics and enter percentiles such as `10, 50, 90, 99`.
4. Click **Update summaries** when the button is active.
5. Sort by a statistic or click a target name for its focus view.

For `GENE_A` in the toy:

| Total | Datasets > 100 | % > 100 | Mean | Median | P10 |
|---:|---:|---:|---:|---:|---:|
| 7,310 | 3 | 75% | 1,827.5 | 1,581 | 78.2 |

**Summary CSV** exports all summarized targets in the chosen order.
Use **Values → Apply the view normalization** for normalized summaries.

**Gene focus** shows one target's statistics, top-dataset bar chart, and a paged
profile with expression, fold change, and log₂ fold change.
Enter its exact name or original `#ID`, then click **Update profile**.

## Select targets and datasets

To compare genes A and B, paste this into **Paste target names or #IDs**:

```text
GENE_A
GENE_B
```

Click **Use target list**. Then paste sample names into **Paste dataset names**:

```text
sample_A1
sample_B1
```

Click **Use dataset list**, then open **Heatmap** and **Update heatmap**.

| Control | Effect |
|---|---|
| Search box | Filter the available choices |
| Select matching targets/datasets | Replace the selection with all search matches |
| Use target/dataset list | Replace the selection with the pasted entries |
| Checkbox | Add or remove one entry |
| `#1`, `#2`, … | Select targets by their original numeric IDs |

## Color and filter the heatmap

| Task | Control |
|---|---|
| Hide absent targets | **Hide all-zero rows**, then **Apply view options** |
| Find high expression | **Sort targets → Total abundance · highest first** |
| Reorder samples | **Sort datasets** |
| Change colors | **Main color** or the preset selector |
| Compress the color range | **Color scale → Logarithmic** |
| Show log values in cells | **Cell values → log₁₀(count + 1)** |
| Shorten numbers | **Number format → Compact** (`4.08K`) or **Scientific** (`4.08e+3`) |
| Filter a count range | **Min count** / **Max count**; excluded cells become gray |

All-zero filtering checks the entire selected dataset scope.
Color intensity uses the current page; hover a cell for its full value.

## Normalize or show fold changes

### Dataset size and target length

Choose **Normalization**, then update the active view:

| Choice | Calculation |
|---|---|
| Raw counts | Decoded k-mer observations |
| Per million dataset k-mers · CPM | Count ÷ full dataset k-mer total × 1,000,000 |
| Per target k-mer position · coverage | Count ÷ target k-mer length |
| CPM / target k-mer length | Apply both corrections |

The normalized toy value for `GENE_A / sample_A1` is approximately **24,681.2**
with **CPM / target k-mer length**.

### Fold change across samples

In **Heatmap**, select **Cell values → Fold change** or **Log₂ fold change**,
then **Update heatmap**. Each target is compared with its mean across selected
real datasets. For `GENE_A / sample_A1`, FC is approximately **2.23256**.

Colors center on the mean. Change the dataset selection to change the reference
population. [Calculation details](IMPLEMENTATION.md#fold-change-definitions).

## Explore transcripts

Open `examples/isoforms/toy.eai`, or package your EM results:

```bash
expresso pack --input isoform-results --output isoform-results/transcripts.eai
```

Use **Transcript metric** to switch between **Assigned k-mer observations** and
**EM abundance**. The toy has EM abundances **6** and **4**.
Missing/unresolved estimates appear as missing cells and have their own counts
in summaries. [Run the complete example](ISOFORMS.md#run-the-exon-skipping-example).

## Export a table or figure

| Button | Download |
|---|---|
| Visible CSV | Current heatmap page as a target × dataset matrix |
| Full selection CSV | All selected pages, one target/dataset/value per row |
| Overview CSV | Whole-matrix statistics |
| Calculated rows CSV | Completed dataset-list rows |
| Summary CSV | All computed target summaries |
| Statistics CSV | The focused target's summary |
| Visible profile CSV | The focused target's current dataset page |
| Ranking CSV / Chart SVG | Ranked values/percentages or a standalone bar chart |

Heatmap CSVs use the applied normalization and fold-change mode.
Log₁₀ display and compact/scientific number formatting leave CSV precision intact.
After a full-selection download, **Clear temporary CSV** removes its browser copy.

You can also export from the command line:

```bash
expresso export --input examples/viewer/toy.eai --output toy-selection \
  --dataset sample_A1 --target GENE_A --target GENE_B --compression none
```

## Work with large indexes

- Opening a large index reads its catalogue first; choose **Update summary** for a global scan.
- Heatmaps page through **100 targets × 25 datasets**; selections can span any number of pages.
- Dataset lists and target profiles page through **25 datasets**.
- Target lists page through **100 targets**.
- The viewer reads requested dataset blocks and keeps summaries and visible pages.
- For bounded quantile memory, choose **Bounded sample**; choose **Exact distribution** for exact quantiles.
- Use **Cancel operation** to stop a running scan or export.

Update buttons are dimmed when the displayed result is current.
Color, formatting, and cached sorting changes redraw immediately.
A progress bar appears after 0.5 seconds of computation.

## Open your own results

Generate the viewer file during quantification:

```bash
expresso quantify --index gene-index --fof datasets.tsv \
  --output gene-results --threads 16 --jobs 4 --viewer-index
```

Or package an existing compact result directory:

```bash
expresso pack --input gene-results --output gene-results/abundance.eai
```

Select the resulting `.eai` through **Choose index**.

## Quick fixes

| Problem | Action |
|---|---|
| Port already in use | Start on `9000` and open `http://127.0.0.1:9000` |
| Opened `index.html` directly | Use one of the localhost launch commands above |
| No heatmap visible | Open **Heatmap**; check selections and click its active update button |
| Normalization metadata missing | Rebuild/requantify with the current EXPRESSO version |
| Full-selection browser export unavailable | Use **Visible CSV** or command-line export |
| Viewer still shows an older interface | Refresh after restarting/rebuilding the server |

[Storage, quantile algorithms, and binary formats](IMPLEMENTATION.md#browser-storage-and-file-formats).
