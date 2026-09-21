//! Recompress CSV or EAB results using EXPRESSO's Zstandard library.
//! Verify every resulting file against the original uncompressed contents.
use anyhow::{Context, Result, ensure};
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Write},
    path::{Component, Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
    time::Instant,
};

#[derive(Parser)]
struct Args {
    /// Completed EXPRESSO results with Zstandard-compressed CSVs or EAB vectors.
    #[arg(long)]
    results: PathBuf,
    /// New directory for recompressed files and reports.
    #[arg(long)]
    output: PathBuf,
    #[arg(long, value_delimiter = ',', default_value = "1,3,6,9,12,15,19")]
    levels: Vec<i32>,
    /// Concurrent files; each encoder uses one thread.
    #[arg(long, default_value_t = 64)]
    jobs: usize,
    /// Limit file count for a pilot; omit for the full corpus.
    #[arg(long)]
    limit: Option<usize>,
}

#[derive(Serialize)]
struct FileResult {
    file: PathBuf,
    original_compressed_bytes: u64,
    uncompressed_bytes: u64,
    compressed_bytes: u64,
    content_bytes_verified: bool,
    compressed_bytes_identical_to_original: bool,
    elapsed_seconds: f64,
}

#[derive(Serialize)]
struct LevelResult {
    level: i32,
    files: usize,
    dataset_compressed_bytes: u64,
    global_compressed_bytes: u64,
    shared_reference_compressed_bytes: u64,
    statistics_compressed_bytes: u64,
    manifest_bytes: u64,
    total_with_manifest_bytes: u64,
    total_compressed_bytes: u64,
    original_compressed_bytes: u64,
    uncompressed_bytes: u64,
    percent_smaller_than_original: f64,
    compression_ratio: f64,
    elapsed_seconds_including_decode_and_verification: f64,
    all_content_bytes_verified: bool,
    all_compressed_bytes_identical_to_original: bool,
}

fn json(path: &Path, value: &impl Serialize) -> Result<()> {
    let temp = path.with_extension("json.tmp");
    let mut out = BufWriter::new(File::create(&temp)?);
    serde_json::to_writer_pretty(&mut out, value)?;
    writeln!(out)?;
    out.flush()?;
    fs::rename(temp, path)?;
    Ok(())
}

fn recompress(
    source: &Path,
    destination: &Path,
    relative: &Path,
    level: i32,
) -> Result<FileResult> {
    let started = Instant::now();
    let input = source.join(relative);
    let output = destination.join(relative);
    let original = fs::read(&input)?;
    let raw = zstd::stream::decode_all(original.as_slice())?;
    let writer = BufWriter::with_capacity(1 << 20, File::create(&output)?);
    let mut encoder = zstd::stream::write::Encoder::new(writer, level)?;
    // Match output::csv: unknown source size, default window, no checksum,
    // streaming writes, explicit flush followed by finish and file flush.
    for chunk in raw.chunks(8192) {
        encoder.write_all(chunk)?;
    }
    encoder.flush()?;
    encoder.finish()?.flush()?;
    let mut decoder = zstd::stream::read::Decoder::new(File::open(&output)?)?;
    let mut offset = 0;
    let mut buffer = [0u8; 65536];
    loop {
        let count = decoder.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        ensure!(
            raw.get(offset..offset + count) == Some(&buffer[..count]),
            "Contents differ: {}",
            relative.display()
        );
        offset += count;
    }
    ensure!(
        offset == raw.len(),
        "Truncated contents: {}",
        relative.display()
    );
    let compressed_bytes = fs::metadata(&output)?.len();
    let identical = compressed_bytes == original.len() as u64 && fs::read(&output)? == original;
    Ok(FileResult {
        file: relative.to_owned(),
        original_compressed_bytes: original.len() as u64,
        uncompressed_bytes: raw.len() as u64,
        compressed_bytes,
        content_bytes_verified: true,
        compressed_bytes_identical_to_original: identical,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    })
}

fn main() -> Result<()> {
    let args = Args::parse();
    ensure!(
        args.jobs > 0 && args.limit != Some(0),
        "jobs and limit must be positive"
    );
    ensure!(
        !args.levels.is_empty() && args.levels.iter().all(|n| (1..=22).contains(n)),
        "levels must be in 1..=22"
    );
    let mut levels = args.levels.clone();
    levels.sort_unstable();
    levels.dedup();
    ensure!(
        levels.len() == args.levels.len(),
        "Repeated compression level"
    );
    ensure!(!args.output.exists(), "Output already exists");
    let manifest: serde_json::Value = serde_json::from_reader(BufReader::new(File::open(
        args.results.join("manifest.json"),
    )?))?;
    ensure!(
        manifest["compression"] == "zstd",
        "Expected Zstandard source CSVs"
    );
    let mut files = manifest["datasets"]
        .as_array()
        .context("Missing datasets")?
        .iter()
        .map(|d| {
            d["output"]
                .as_str()
                .map(PathBuf::from)
                .context("Missing dataset output")
        })
        .collect::<Result<Vec<_>>>()?;
    let compact = manifest["format"] == "compact";
    let global_path = PathBuf::from(
        manifest["global"]["output"]
            .as_str()
            .unwrap_or("global.csv.zst"),
    );
    files.push(global_path.clone());
    let reference_path = if compact {
        let path = PathBuf::from(
            manifest["reference"]["file"]
                .as_str()
                .context("Missing shared reference")?,
        );
        files.push(path.clone());
        if let Some(columns) = manifest["global"]["statistics"].as_array() {
            for column in columns {
                files.push(PathBuf::from(
                    column["output"]
                        .as_str()
                        .context("Missing statistic path")?,
                ));
            }
        }
        Some(path)
    } else {
        None
    };
    let manifest_bytes = if compact {
        fs::metadata(args.results.join("manifest.json"))?.len()
    } else {
        0
    };
    if let Some(limit) = args.limit {
        files.truncate(limit);
    }
    ensure!(!files.is_empty(), "No abundance files");
    for path in &files {
        ensure!(
            path.components().all(|c| matches!(c, Component::Normal(_)))
                && path.extension().is_some_and(|ext| ext == "zst"),
            "Invalid abundance path"
        );
        ensure!(
            args.results.join(path).is_file(),
            "Missing input: {}",
            path.display()
        );
    }
    fs::create_dir_all(&args.output)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.jobs)
        .build()?;
    let mut results = Vec::new();
    for level in args.levels {
        let started = Instant::now();
        let out = args.output.join(format!("level-{level:02}"));
        fs::create_dir(&out)?;
        for path in &files {
            if let Some(parent) = path.parent() {
                fs::create_dir_all(out.join(parent))?;
            }
        }
        if compact {
            fs::copy(
                args.results.join("manifest.json"),
                out.join("manifest.json"),
            )?;
        }
        let completed = AtomicUsize::new(0);
        eprintln!(
            "Level {level}: {} files, {} workers",
            files.len(),
            args.jobs
        );
        let per_file: Vec<FileResult> = pool.install(|| {
            files
                .par_iter()
                .map(|path| {
                    let result = recompress(&args.results, &out, path, level)
                        .with_context(|| format!("level {level}: {}", path.display()))?;
                    let n = completed.fetch_add(1, Ordering::Relaxed) + 1;
                    if n.is_multiple_of(50) || n == files.len() {
                        eprintln!(
                            "Level {level}: {n}/{} complete ({:.1}s)",
                            files.len(),
                            started.elapsed().as_secs_f64()
                        );
                    }
                    Ok(result)
                })
                .collect::<Result<Vec<_>>>()
        })?;
        let total: u64 = per_file.iter().map(|f| f.compressed_bytes).sum();
        let original: u64 = per_file.iter().map(|f| f.original_compressed_bytes).sum();
        let raw: u64 = per_file.iter().map(|f| f.uncompressed_bytes).sum();
        let global: u64 = per_file
            .iter()
            .filter(|f| f.file == global_path)
            .map(|f| f.compressed_bytes)
            .sum();
        let reference: u64 = per_file
            .iter()
            .filter(|f| Some(&f.file) == reference_path.as_ref())
            .map(|f| f.compressed_bytes)
            .sum();
        let datasets: u64 = per_file
            .iter()
            .filter(|f| f.file.starts_with("datasets"))
            .map(|f| f.compressed_bytes)
            .sum();
        results.push(LevelResult {
            level,
            files: per_file.len(),
            dataset_compressed_bytes: datasets,
            global_compressed_bytes: global,
            shared_reference_compressed_bytes: reference,
            statistics_compressed_bytes: total - datasets - global - reference,
            manifest_bytes,
            total_with_manifest_bytes: total + manifest_bytes,
            total_compressed_bytes: total,
            original_compressed_bytes: original,
            uncompressed_bytes: raw,
            percent_smaller_than_original: 100.0 * (1.0 - total as f64 / original as f64),
            compression_ratio: raw as f64 / total as f64,
            elapsed_seconds_including_decode_and_verification: started.elapsed().as_secs_f64(),
            all_content_bytes_verified: per_file.iter().all(|f| f.content_bytes_verified),
            all_compressed_bytes_identical_to_original: per_file
                .iter()
                .all(|f| f.compressed_bytes_identical_to_original),
        });
        json(&out.join("files.json"), &per_file)?;
        json(
            &args.output.join("summary.json"),
            &serde_json::json!({
                "source": args.results, "zstd_version": zstd::zstd_safe::version_string(),
                "workers": args.jobs, "pilot_limit": args.limit,
                "format":if compact {"compact"} else {"csv"},
                "scope": if compact {"Dataset and global EAB vectors, optional statistics, shared reference table; plain manifest bytes reported separately and in total_with_manifest_bytes"} else {"Per-dataset and global abundance CSVs; excludes exons.csv and manifest.json"},
                "method": "Independent streaming compression per file, no dictionary, unknown content size, no frame checksum; lossless roundtrip verified for every file",
                "timing": "One pass per level, including baseline decompression, compression, output writes and readback verification; not EXPRESSO quantification time",
                "all_requested_levels_finished": results.len() == levels.len(), "levels": results
            }),
        )?;
        eprintln!(
            "Level {level}: {total} bytes ({:.1}% smaller than original), {:.1}s",
            100.0 * (1.0 - total as f64 / original as f64),
            started.elapsed().as_secs_f64()
        );
    }
    Ok(())
}
