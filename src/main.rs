mod abundance;
mod compact;
mod eab;
mod exons;
mod export;
mod index;
mod input;
mod output;
mod portable;
mod quantify;
mod viewer;

#[derive(
    Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum Level {
    #[default]
    Exon,
    Gene,
}
impl Level {
    pub fn columns(self) -> [&'static str; 2] {
        match self {
            Self::Exon => ["exon_id", "exon_name"],
            Self::Gene => ["gene_id", "gene_name"],
        }
    }
    pub fn table(self) -> &'static str {
        match self {
            Self::Exon => "exons.csv",
            Self::Gene => "genes.csv",
        }
    }
}

use anyhow::{Result, ensure};
use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "expresso",
    version,
    about = "EXPRESSO — EXon-level RNA EXPRESsion quantificatiOn"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Extract strand-oriented exon FASTA from a GTF and reference genome FASTAs.
    #[command(visible_alias = "exons")]
    ExtractExons(ExtractExonsArgs),
    /// Build reusable GGCAT simplitigs, Rust SSHash, and exon or gene ownership tables.
    Build(BuildArgs),
    /// Count a file-of-files against a previously built index.
    Quantify(QuantifyArgs),
    /// Decode compact abundance vectors into CSVs (one abundance column by default).
    Export(ExportArgs),
    /// Package compact results into a single indexed file for the WebAssembly viewer.
    Pack {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Serve the bundled WebAssembly abundance viewer on localhost.
    View {
        #[arg(long, default_value_t = 8765)]
        port: u16,
    },
    /// Build an index and quantify in one invocation.
    Run {
        #[command(flatten)]
        build: BuildArgs,
        #[command(flatten)]
        query: QueryOptions,
    },
    /// Print index metadata as JSON.
    Inspect {
        #[arg(short, long)]
        index: PathBuf,
    },
}

#[derive(Args)]
pub struct BuildArgs {
    /// Quantification target: exon or gene. Gene FASTA records with the same ID are grouped.
    #[arg(long, value_enum, default_value = "exon")]
    level: Level,
    /// Reference FASTA; gene mode groups records with the same first header token.
    #[arg(short, long, visible_alias = "genes", required_unless_present = "gtf", conflicts_with_all = ["gtf", "genome"])]
    exons: Option<PathBuf>,
    /// Generate exons from this GTF instead of supplying --exons.
    #[arg(long, requires = "genome")]
    gtf: Option<PathBuf>,
    /// Reference genome FASTA; repeat for multiple files. Used with --gtf.
    #[arg(long, visible_alias = "reference", requires = "gtf")]
    genome: Vec<PathBuf>,
    /// New index directory. Existing paths are never overwritten.
    #[arg(short, long)]
    index: PathBuf,
    /// Odd k-mer length in 3..=63.
    #[arg(short, long, default_value_t = 31)]
    k: usize,
    /// SSHash minimizer length; default min(19, k-2).
    #[arg(short, long)]
    minimizer: Option<usize>,
    #[arg(short = 't', long, default_value_t = default_threads())]
    threads: usize,
    /// GGCAT memory hint and SSHash external-sort budget, in GiB (not a hard RSS cap).
    #[arg(long, default_value_t = 4)]
    memory_gb: usize,
    /// Parent directory for construction scratch files.
    #[arg(long)]
    temp_dir: Option<PathBuf>,
}

#[derive(Args)]
pub struct QuantifyArgs {
    #[arg(short, long)]
    index: PathBuf,
    #[arg(short = 't', long, default_value_t = default_threads())]
    threads: usize,
    #[command(flatten)]
    query: QueryOptions,
}

#[derive(Args)]
pub struct QueryOptions {
    /// One dataset path per line, relative to this list's directory, or absolute.
    /// Optional TSV columns: dataset_name<TAB>path[<TAB>reads|unitigs|auto].
    #[arg(short = 'f', long)]
    fof: PathBuf,
    /// New output directory with per-dataset vectors and a global vector.
    #[arg(short, long)]
    output: PathBuf,
    /// Auto uses header abundance when available, otherwise read counting.
    #[arg(long, value_enum, default_value_t = abundance::Mode::Auto)]
    mode: abundance::Mode,
    /// Number of simultaneous datasets, sharing the total thread budget.
    #[arg(short = 'j', long, default_value_t = 1)]
    jobs: usize,
    /// Approximate DNA bytes per work batch; long records are split.
    #[arg(long, default_value_t = 1048576)]
    batch_bases: usize,
    /// Compact bit-packed vectors, or legacy exact CSV output.
    #[arg(long, value_enum, default_value_t = output::Format::Compact)]
    format: output::Format,
    /// Also write abundance.eai for direct loading in the WebAssembly viewer (compact output only).
    #[arg(long)]
    viewer_index: bool,
    /// Bits per abundance in compact output; 8 means one byte before compression.
    #[arg(long, default_value_t = 8, value_parser = clap::value_parser!(u8).range(2..=16))]
    bits: u8,
    /// Compression for abundance vectors or CSVs.
    #[arg(long, value_enum, default_value_t = output::Compression::Zstd)]
    compression: output::Compression,
    /// Add mean, median, min, max, and detected-dataset count to global output.
    #[arg(long)]
    stats: bool,
    /// Record failed input datasets and continue; global sums include successful datasets only.
    #[arg(long, conflicts_with = "stats")]
    keep_going: bool,
    /// Assembly k for interpreting KC:i total counts; defaults to index k.
    #[arg(long)]
    unitig_k: Option<usize>,
}

#[derive(Args)]
pub struct ExportArgs {
    /// Compact EXPRESSO result directory or packed .eai viewer index.
    #[arg(short, long)]
    input: PathBuf,
    /// New destination directory for exported CSVs.
    #[arg(short, long)]
    output: PathBuf,
    /// Export only these datasets (repeat the flag); defaults to all plus global.
    #[arg(long, conflicts_with = "global_only")]
    dataset: Vec<String>,
    /// Dataset names, one per line; combined with --dataset.
    #[arg(long, conflicts_with = "global_only")]
    dataset_list: Option<PathBuf>,
    /// Target names to export; repeat to select more (also --exon or --gene).
    #[arg(long, visible_aliases = ["exon", "gene"])]
    target: Vec<String>,
    /// Target names, one per line; combined with --target.
    #[arg(long, visible_alias = "exon-list")]
    target_list: Option<PathBuf>,
    /// Explicit 1-based target IDs; repeat to select more.
    #[arg(long)]
    target_id: Vec<usize>,
    /// Inclusive lower bound on decoded abundance. Nonmatching rows are omitted.
    #[arg(long)]
    min_value: Option<u64>,
    /// Inclusive upper bound on decoded abundance. Nonmatching rows are omitted.
    #[arg(long)]
    max_value: Option<u64>,
    /// Export only the global vector and optional statistics.
    #[arg(long)]
    global_only: bool,
    /// Include target ID and name columns; default CSVs contain abundance only.
    #[arg(long)]
    with_names: bool,
    #[arg(long, value_enum, default_value_t = output::Compression::Zstd)]
    compression: output::Compression,
    #[arg(short = 't', long, default_value_t = default_threads())]
    threads: usize,
}

#[derive(Args)]
pub struct ExtractExonsArgs {
    /// GTF annotation containing exon features (plain, gzip, xz, or zstd).
    #[arg(long)]
    gtf: PathBuf,
    /// Reference genome FASTA; repeat for multiple files.
    #[arg(long, visible_alias = "reference", required = true)]
    genome: Vec<PathBuf>,
    /// New output FASTA file, with one unwrapped sequence line per exon.
    #[arg(short, long)]
    output: PathBuf,
    /// Output compression; explicitly choose a matching filename suffix.
    #[arg(long, value_enum, default_value_t = output::Compression::None)]
    compression: output::Compression,
}

fn default_threads() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::ExtractExons(args) => {
            exons::extract(&args.gtf, &args.genome, &args.output, args.compression)
        }
        Command::Build(args) => index::build(&args),
        Command::Quantify(args) => quantify::run(&args.index, args.threads, &args.query),
        Command::Export(args) => export::run(&args),
        Command::Pack { input, output } => portable::pack(&input, &output),
        Command::View { port } => viewer::serve(port),
        Command::Run { build, query } => {
            ensure!(
                !query.viewer_index || query.format == output::Format::Compact,
                "--viewer-index requires --format compact"
            );
            ensure!(
                !query.output.exists(),
                "output already exists: {}",
                query.output.display()
            );
            input::datasets(&query.fof, query.mode)?;
            index::build(&build)?;
            quantify::run(&build.index, build.threads, &query)
        }
        Command::Inspect { index } => {
            let meta = index::Metadata::load(&index)?;
            println!("{}", serde_json::to_string_pretty(&meta)?);
            Ok(())
        }
    }
}
