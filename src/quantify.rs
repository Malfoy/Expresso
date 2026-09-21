use crate::{
    QueryOptions,
    abundance::round_ratio,
    compact,
    index::{self, Metadata, Owners},
    input::{self, Dataset},
    output,
};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use sshash_lib::{Dictionary, Kmer, KmerBits};
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
};

#[derive(Default, Serialize)]
struct Metrics {
    records: u64,
    valid_kmers: u64,
    assigned_kmers: u64,
}
struct Counts {
    values: Vec<u64>,
    metrics: Metrics,
}

struct OutputContext<'a> {
    meta: &'a Metadata,
    options: &'a QueryOptions,
    reference: [u8; 32],
}

fn add(a: &mut u64, b: u64) -> Result<()> {
    *a = a
        .checked_add(b)
        .context("count overflow (u64); split this dataset into smaller inputs")?;
    Ok(())
}

#[inline]
fn scan<const K: usize>(
    seq: &[u8],
    value: u64,
    owners: &Owners,
    counts: &mut Counts,
    query: &mut sshash_lib::streaming_query::StreamingQueryEngine<'_, K>,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    for run in seq.split(|b| !input::is_dna(*b)).filter(|s| s.len() >= K) {
        query.reset();
        add(&mut counts.metrics.valid_kmers, (run.len() - K + 1) as u64)?;
        for window in run.windows(K) {
            let hit = query.lookup(window);
            if hit.is_found()
                && let Some(exon) = owners.exon(hit.kmer_id as usize)
            {
                add(&mut counts.values[exon], value)?;
                add(&mut counts.metrics.assigned_kmers, 1)?;
            }
        }
    }
    Ok(())
}

fn count<const K: usize>(
    dataset: &Dataset,
    dict: &Dictionary,
    owners: &Owners,
    n: usize,
    pool: &rayon::ThreadPool,
    options: &QueryOptions,
) -> Result<Counts>
where
    Kmer<K>: KmerBits,
{
    if pool.current_num_threads() == 1 {
        // With one worker per dataset, query borrowed Helicase records directly.
        // Avoid copying, queue handoffs, and a separate decoder thread.
        let mut counts = Counts {
            values: vec![0; n],
            metrics: Metrics::default(),
        };
        let mut query = dict.create_streaming_query::<K>();
        input::records(&dataset.path, |header, seq| {
            let value = crate::abundance::weight(
                header,
                seq.len(),
                options.unitig_k.unwrap_or(K),
                dataset.mode,
            )?;
            add(&mut counts.metrics.records, 1)?;
            scan::<K>(seq, value, owners, &mut counts, &mut query)
        })
        .with_context(|| format!("counting dataset {}", dataset.name))?;
        return Ok(counts);
    }
    input::batches(
        &dataset.path,
        dataset.mode,
        options.unitig_k.unwrap_or(K),
        K,
        options.batch_bases,
        2 * pool.current_num_threads(),
        |rx| {
            // Let each worker block directly on the MPMC queue. par_bridge
            // serializes next() behind a mutex and generates excessive wakeups
            // when the decoder cannot keep every query worker busy.
            let stopped = AtomicBool::new(false);
            let partials = pool.broadcast(|_| {
                let result = (|| -> Result<Counts> {
                    let mut counts = Counts {
                        values: vec![0; n],
                        metrics: Metrics::default(),
                    };
                    let mut query = dict.create_streaming_query::<K>();
                    while !stopped.load(Ordering::Relaxed) {
                        let Ok(batch) = rx.recv() else { break };
                        add(&mut counts.metrics.records, batch.num_records())?;
                        for (record, seq) in batch.records() {
                            scan::<K>(seq, record.value, owners, &mut counts, &mut query)?;
                        }
                    }
                    Ok(counts)
                })();
                if result.is_err() {
                    stopped.store(true, Ordering::Relaxed);
                }
                result
            });
            let mut partials = partials.into_iter();
            let mut counts = partials.next().context("empty worker pool")??;
            for other in partials {
                let other = other?;
                for (a, b) in counts.values.iter_mut().zip(other.values) {
                    add(a, b)?;
                }
                add(&mut counts.metrics.records, other.metrics.records)?;
                add(&mut counts.metrics.valid_kmers, other.metrics.valid_kmers)?;
                add(
                    &mut counts.metrics.assigned_kmers,
                    other.metrics.assigned_kmers,
                )?;
            }
            Ok(counts)
        },
    )
    .with_context(|| format!("counting dataset {}", dataset.name))
}

pub fn run(index_dir: &Path, threads: usize, options: &QueryOptions) -> Result<()> {
    ensure!(
        threads > 0 && options.jobs > 0 && options.batch_bases > 0,
        "threads, jobs and batch-bases must be positive"
    );
    ensure!(options.unitig_k != Some(0), "unitig-k must be positive");
    let datasets = input::datasets(&options.fof, options.mode)?;
    let stage = index::staging(&options.output)?;
    let meta = Metadata::load(index_dir)?;
    let dict = Dictionary::load(index_dir.join("dictionary"))
        .map_err(|e| anyhow::anyhow!("loading SSHash: {e}"))?;
    ensure!(
        dict.k() == meta.k
            && dict.m() == meta.minimizer
            && index::num_kmers(&dict)? == meta.indexed_kmers,
        "SSHash and metadata disagree"
    );
    let owners = Owners::load(index_dir, &meta)?;
    let jobs = options.jobs.min(threads).min(datasets.len());
    let pools: Vec<_> = (0..jobs)
        .map(|i| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads / jobs + usize::from(i < threads % jobs))
                .build()
        })
        .collect::<std::result::Result<_, _>>()?;
    let n = meta.exons.len();
    let output_context = OutputContext {
        meta: &meta,
        options,
        reference: compact::reference_id(&meta.exons),
    };
    let mut total = vec![0u64; n];
    let scratch = tempfile::Builder::new()
        .prefix("counts-")
        .tempdir_in(stage.path())?;
    let dataset_dir = stage.path().join("datasets");
    fs::create_dir(&dataset_dir)?;
    let mut manifest = vec![serde_json::Value::Null; datasets.len()];
    let mut failed_datasets = Vec::new();
    // Schedule large inputs first to reduce the last-file tail. Preserve input
    // order in metadata and statistics, independently of completion order.
    let mut order: Vec<_> = (0..datasets.len()).collect();
    let sizes: Vec<_> = datasets
        .iter()
        .map(|d| fs::metadata(&d.path).map(|m| m.len()))
        .collect::<std::io::Result<_>>()?;
    order.sort_by_key(|i| std::cmp::Reverse(sizes[*i]));
    eprintln!(
        "Quantifying {} datasets with {threads} query workers ({jobs} concurrent datasets)",
        datasets.len()
    );
    let next = AtomicUsize::new(0);
    let cancelled = AtomicBool::new(false);
    std::thread::scope(|scope| -> Result<()> {
        let (tx, rx) = crossbeam_channel::bounded(jobs);
        let mut handles = Vec::new();
        for pool in &pools {
            let output_context = &output_context;
            let tx = tx.clone();
            let (next, cancelled, order, datasets, dict, owners, meta, dataset_dir, scratch) = (
                &next,
                &cancelled,
                &order,
                &datasets,
                &dict,
                &owners,
                &meta,
                &dataset_dir,
                scratch.path(),
            );
            handles.push(scope.spawn(move || {
                while !cancelled.load(Ordering::Relaxed) {
                    let position = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&id) = order.get(position) else { break };
                    let dataset = &datasets[id];
                    let counts = sshash_lib::dispatch_on_k!(meta.k, K => count::<K>(dataset, dict, owners, n, pool, options));
                    let skip_failed_input = counts.is_err() && options.keep_going;
                    let result = counts.and_then(|counts| finish_dataset(id, dataset, counts, output_context, dataset_dir, scratch));
                    // Output and aggregation errors remain fatal, including disk-full errors.
                    if result.is_err() && !skip_failed_input { cancelled.store(true, Ordering::Relaxed); }
                    if tx.send((id, result, skip_failed_input)).is_err() { break; }
                }
            }));
        }
        drop(tx);
        let mut failure = None;
        for (id, result, skip_failed_input) in rx {
            match result {
                Ok((counts, entry)) => {
                    for (sum, value) in total.iter_mut().zip(counts.values) {
                        if let Err(error) = add(sum, value) {
                            failure.get_or_insert(error);
                            cancelled.store(true, Ordering::Relaxed);
                            break;
                        }
                    }
                    manifest[id] = entry;
                }
                Err(error) if skip_failed_input => {
                    let dataset = &datasets[id];
                    let entry = serde_json::json!({"name":dataset.name, "path":dataset.path,
                        "error":format!("{error:#}")});
                    eprintln!("DATASET_FAILED {}", serde_json::to_string(&entry)?);
                    failed_datasets.push(entry);
                }
                Err(error) => {
                    failure.get_or_insert(error);
                }
            }
        }
        for handle in handles {
            if handle.join().is_err() {
                failure.get_or_insert_with(|| anyhow::anyhow!("dataset worker panicked"));
            }
        }
        if let Some(error) = failure {
            return Err(error);
        }
        ensure!(
            manifest.iter().filter(|m| !m.is_null()).count() + failed_datasets.len()
                == datasets.len(),
            "not all datasets completed"
        );
        ensure!(
            failed_datasets.len() < datasets.len(),
            "all input datasets failed"
        );
        Ok(())
    })?;
    manifest.retain(|entry| !entry.is_null());
    let mut statistics = options.stats.then(|| vec![vec![0u64; n]; 5]);
    if let Some(statistics) = statistics.as_mut() {
        // Exact medians using disk-backed count vectors, read in exon blocks.
        // At most ~32 MiB of matrix data, and one open scratch file at a time.
        let block_size = ((32 << 20) / (8 * datasets.len())).clamp(1, 4096);
        for start in (0..n).step_by(block_size) {
            let len = block_size.min(n - start);
            let mut matrix = vec![0u64; len * datasets.len()];
            for d in 0..datasets.len() {
                let mut file = BufReader::new(File::open(scratch.path().join(format!("{d}.bin")))?);
                file.seek(SeekFrom::Start(start as u64 * 8))?;
                for e in 0..len {
                    let mut b = [0; 8];
                    file.read_exact(&mut b)?;
                    matrix[e * datasets.len() + d] = u64::from_le_bytes(b);
                }
            }
            for e in 0..len {
                let id = start + e;
                let values = &mut matrix[e * datasets.len()..(e + 1) * datasets.len()];
                let min = *values.iter().min().unwrap();
                let max = *values.iter().max().unwrap();
                let detected = values.iter().filter(|v| **v > 0).count();
                let median = median(values)?;
                let mean = round_ratio(total[id] as u128, datasets.len() as u128)?;
                for (column, value) in
                    statistics
                        .iter_mut()
                        .zip([mean, median, min, max, detected as u64])
                {
                    column[id] = value;
                }
            }
        }
    }
    let global_entry;
    let (reference_file, reference_compression);
    if options.format == output::Format::Compact {
        let filename = format!("global{}", options.compression.vector_suffix());
        let info = compact::write_vector(
            &stage.path().join(&filename),
            options.compression,
            &total,
            options.bits,
            &output_context.reference,
        )?;
        let mut columns = Vec::new();
        if let Some(statistics) = &statistics {
            let dir = stage.path().join("statistics");
            fs::create_dir(&dir)?;
            for (name, values) in ["mean", "median", "min", "max", "datasets_detected"]
                .iter()
                .zip(statistics)
            {
                let path = format!("statistics/{name}{}", options.compression.vector_suffix());
                let quantization = compact::write_vector(
                    &stage.path().join(&path),
                    options.compression,
                    values,
                    options.bits,
                    &output_context.reference,
                )?;
                columns.push(
                    serde_json::json!({"name":name,"output":path,"quantization":quantization}),
                );
            }
        }
        global_entry =
            serde_json::json!({"output":filename,"quantization":info,"statistics":columns});
        reference_file = "exons.csv.zst";
        reference_compression = output::Compression::Zstd;
    } else {
        let filename = format!("global{}", options.compression.suffix());
        output::csv(&stage.path().join(&filename), options.compression, |csv| {
            if let Some(statistics) = &statistics {
                csv.write_record([
                    "exon_id",
                    "exon_name",
                    "abundance",
                    "mean",
                    "median",
                    "min",
                    "max",
                    "datasets_detected",
                ])?;
                for (id, e) in meta.exons.iter().enumerate() {
                    csv.serialize((
                        id + 1,
                        &e.name,
                        total[id],
                        statistics[0][id],
                        statistics[1][id],
                        statistics[2][id],
                        statistics[3][id],
                        statistics[4][id],
                    ))?;
                }
            } else {
                csv.write_record(["exon_id", "exon_name", "abundance"])?;
                for (id, (exon, value)) in meta.exons.iter().zip(&total).enumerate() {
                    csv.serialize((id + 1, &exon.name, value))?;
                }
            }
            Ok(())
        })?;
        global_entry = serde_json::json!({"output":filename});
        reference_file = "exons.csv";
        reference_compression = output::Compression::None;
    }
    output::csv(
        &stage.path().join(reference_file),
        reference_compression,
        |csv| {
            csv.write_record(["exon_id", "exon_name", "length", "unique_kmers"])?;
            for (i, e) in meta.exons.iter().enumerate() {
                csv.serialize((i + 1, &e.name, e.length, e.unique_kmers))?;
            }
            Ok(())
        },
    )?;
    let mut out = BufWriter::new(File::create(stage.path().join("manifest.json"))?);
    serde_json::to_writer_pretty(
        &mut out,
        &serde_json::json!({
            "expresso_version":env!("CARGO_PKG_VERSION"), "index":index_dir.canonicalize()?, "k":meta.k,
        "unitig_k":options.unitig_k.unwrap_or(meta.k), "threads":threads, "jobs":jobs,
        "compression":options.compression, "batch_bases":options.batch_bases, "stats":options.stats,
            "format": options.format, "compact_format_version": if options.format == output::Format::Compact {Some(1)} else {None},
            "reference":{"file":reference_file,"sha256":compact::hex(&output_context.reference),"targets":n},
            "global":global_entry,
            "coverage":{"requested_datasets":datasets.len(), "completed_datasets":manifest.len(),
                "failed_datasets":failed_datasets.len(), "complete":failed_datasets.is_empty(),
                "global_includes":"successfully completed datasets only; no partial counts from failed datasets"},
            "failed_datasets":failed_datasets,
            "rounding":"nearest integer, ties up; unitig weights use six decimal places",
            "statistics":if options.format == output::Format::Compact {
                "Global and statistics are calculated from exact rounded dataset counts, then independently quantized; their decoded values need not equal aggregates of decoded datasets. Zeros and nonzero status are preserved."
            } else { "across datasets, including zeros; global is sum of reported dataset abundances" },
            "datasets":manifest
        }),
    )?;
    out.flush()?;
    drop(out);
    drop(scratch);
    index::publish(stage, &options.output)?;
    eprintln!("Results saved: {}", options.output.display());
    Ok(())
}

fn finish_dataset(
    id: usize,
    dataset: &Dataset,
    mut counts: Counts,
    context: &OutputContext<'_>,
    dataset_dir: &Path,
    scratch: &Path,
) -> Result<(Counts, serde_json::Value)> {
    let (meta, options) = (context.meta, context.options);
    for value in &mut counts.values {
        *value = round_ratio(*value as u128, dataset.mode.scale() as u128)?;
    }
    let (filename, quantization) = if options.format == output::Format::Compact {
        let filename = format!("{}{}", dataset.name, options.compression.vector_suffix());
        let info = compact::write_vector(
            &dataset_dir.join(&filename),
            options.compression,
            &counts.values,
            options.bits,
            &context.reference,
        )?;
        (filename, serde_json::to_value(info)?)
    } else {
        let filename = format!("{}{}", dataset.name, options.compression.suffix());
        output::csv(&dataset_dir.join(&filename), options.compression, |csv| {
            csv.write_record(["exon_id", "exon_name", "abundance"])?;
            for (id, (exon, value)) in meta.exons.iter().zip(&counts.values).enumerate() {
                csv.serialize((id + 1, &exon.name, value))?;
            }
            Ok(())
        })?;
        (filename, serde_json::Value::Null)
    };
    if options.stats {
        let mut file = BufWriter::new(File::create(scratch.join(format!("{id}.bin")))?);
        for value in &counts.values {
            file.write_all(&value.to_le_bytes())?;
        }
        file.flush()?;
    }
    eprintln!(
        "{}: {} records, {} / {} valid k-mers assigned",
        dataset.name,
        counts.metrics.records,
        counts.metrics.assigned_kmers,
        counts.metrics.valid_kmers
    );
    let entry = serde_json::json!({"name":dataset.name, "path":dataset.path, "mode":dataset.mode,
        "output":format!("datasets/{filename}"), "metrics":counts.metrics,"quantization":quantization});
    Ok((counts, entry))
}

fn median(values: &mut [u64]) -> Result<u64> {
    let n = values.len();
    ensure!(n > 0, "median of empty collection");
    let (lower, mid, _) = values.select_nth_unstable(n / 2);
    if n % 2 == 1 {
        Ok(*mid)
    } else {
        round_ratio(*lower.iter().max().unwrap() as u128 + *mid as u128, 2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_medians_include_zeros_and_avoid_overflow() {
        assert_eq!(median(&mut [10, 0, 0]).unwrap(), 0);
        assert_eq!(median(&mut [4, 0, 1, 2]).unwrap(), 2);
        assert_eq!(median(&mut [u64::MAX, u64::MAX]).unwrap(), u64::MAX);
        assert!(median(&mut []).is_err());
    }
}
