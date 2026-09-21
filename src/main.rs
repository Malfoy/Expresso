mod abundance;
mod compact;
mod export;
mod index;
mod input;
mod output;
mod quantify;

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
    /// Build reusable GGCAT simplitigs, Rust SSHash, and exon ownership tables.
    Build(BuildArgs),
    /// Count a file-of-files against a previously built index.
    Quantify(QuantifyArgs),
    /// Decode compact abundance vectors into CSVs (one abundance column by default).
    Export(ExportArgs),
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
    /// FASTA: one exon per record (wrapped sequences are accepted).
    #[arg(short, long)]
    exons: PathBuf,
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
    /// GGCAT executable path/name.
    #[arg(long, default_value = "ggcat")]
    ggcat: PathBuf,
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
    /// Compact EXPRESSO result directory, including its shared reference table.
    #[arg(short, long)]
    input: PathBuf,
    /// New destination directory for exported CSVs.
    #[arg(short, long)]
    output: PathBuf,
    /// Export only these datasets (repeat the flag); defaults to all plus global.
    #[arg(long, conflicts_with = "global_only")]
    dataset: Vec<String>,
    /// Export only the global vector and optional statistics.
    #[arg(long)]
    global_only: bool,
    /// Include exon_id and exon_name columns; default CSVs contain abundance only.
    #[arg(long)]
    with_names: bool,
    #[arg(long, value_enum, default_value_t = output::Compression::Zstd)]
    compression: output::Compression,
    #[arg(short = 't', long, default_value_t = default_threads())]
    threads: usize,
}

fn default_threads() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Build(args) => index::build(&args),
        Command::Quantify(args) => quantify::run(&args.index, args.threads, &args.query),
        Command::Export(args) => export::run(&args),
        Command::Run { build, query } => {
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
