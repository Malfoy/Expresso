//! Lossless, sparse feature counts. EFC v1 binds exact fixed-point counts to a model.
use super::{load, validate_query};
use crate::{
    QueryOptions, abundance,
    index::{self, Metadata, Owners},
    input, output,
};
use anyhow::{Context, Result, ensure};
use serde_json::json;
use sshash_lib::Dictionary;
use std::{
    fs,
    io::{BufRead, BufReader, Read},
    path::Path,
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
};

const MAGIC: &[u8; 8] = b"EXPREF01";
pub(super) fn write_counts<T: Copy + Into<u128>>(
    path: &Path,
    compression: output::Compression,
    format: output::Format,
    values: &[T],
    scale: u64,
    model: &[u8; 32],
) -> Result<()> {
    if format == output::Format::Csv {
        return output::csv(path, compression, |csv| {
            csv.write_record(["feature_id", "abundance"])?;
            for (id, &value) in values.iter().enumerate() {
                let value: u128 = value.into();
                if value > 0 {
                    csv.serialize((id + 1, decimal(value, scale)))?;
                }
            }
            Ok(())
        });
    }
    output::encoded(path, compression, |out| {
        let mut crc = crc32fast::Hasher::new();
        let mut write = |bytes: &[u8]| -> Result<()> {
            crc.update(bytes);
            out.write_all(bytes)?;
            Ok(())
        };
        write(MAGIC)?;
        write(model)?;
        write(&(values.len() as u64).to_le_bytes())?;
        write(&scale.to_le_bytes())?;
        write(&(values.iter().filter(|&&v| v.into() > 0u128).count() as u64).to_le_bytes())?;
        for (id, &value) in values.iter().enumerate() {
            let value: u128 = value.into();
            if value > 0 {
                write(&(id as u64).to_le_bytes())?;
                write(&value.to_le_bytes())?;
            }
        }
        out.write_all(&crc.finalize().to_le_bytes())?;
        Ok(())
    })
}
fn decimal(value: u128, scale: u64) -> String {
    let whole = value / scale as u128;
    let fraction = value % scale as u128;
    if fraction == 0 {
        whole.to_string()
    } else {
        format!("{whole}.{fraction:06}")
            .trim_end_matches('0')
            .to_string()
    }
}
pub(super) fn read_counts(path: &Path, features: usize, model: &[u8; 32]) -> Result<Vec<f64>> {
    let mut input = BufReader::new(output::reader(path)?);
    let mut values = vec![0.0; features];
    if input.fill_buf()?.starts_with(MAGIC) {
        let mut header = [0; 64];
        input.read_exact(&mut header)?;
        ensure!(
            &header[8..40] == model,
            "feature counts belong to a different isoform model"
        );
        ensure!(
            u64::from_le_bytes(header[40..48].try_into()?) == features as u64,
            "feature count dimension mismatch"
        );
        let scale = u64::from_le_bytes(header[48..56].try_into()?);
        let n = u64::from_le_bytes(header[56..64].try_into()?);
        ensure!(
            [1, abundance::SCALE].contains(&scale) && n <= features as u64,
            "invalid feature count header"
        );
        let mut crc = crc32fast::Hasher::new();
        crc.update(&header);
        let mut previous = None;
        for _ in 0..n {
            let mut row = [0; 24];
            input.read_exact(&mut row)?;
            crc.update(&row);
            let id = usize::try_from(u64::from_le_bytes(row[..8].try_into()?))?;
            let value = u128::from_le_bytes(row[8..].try_into()?);
            ensure!(
                id < features && previous.is_none_or(|old| id > old) && value > 0,
                "invalid sparse feature count entry"
            );
            values[id] = value as f64 / scale as f64;
            previous = Some(id);
        }
        let mut checksum = [0; 4];
        input.read_exact(&mut checksum)?;
        ensure!(
            crc.finalize() == u32::from_le_bytes(checksum),
            "feature count checksum mismatch"
        );
        ensure!(
            input.read(&mut [0; 1])? == 0,
            "trailing bytes in feature counts"
        );
    } else {
        let mut csv = csv::Reader::from_reader(input);
        ensure!(
            csv.headers()?.iter().eq(["feature_id", "abundance"]),
            "expected disjoint feature CSV columns feature_id,abundance; ordinary exon totals or inclusive intersections cannot be fitted as independent features"
        );
        let mut seen = std::collections::HashSet::new();
        for row in csv.records() {
            let row = row?;
            let id: usize = row[0].parse()?;
            let value: f64 = row[1].parse()?;
            ensure!(
                id > 0 && id <= features && seen.insert(id) && value.is_finite() && value >= 0.0,
                "invalid, repeated or nonfinite feature abundance"
            );
            values[id - 1] = value;
        }
    }
    Ok(values)
}
fn suffix(format: output::Format, compression: output::Compression) -> String {
    if format == output::Format::Csv {
        compression.suffix().into()
    } else {
        compression.vector_suffix().replace(".eab", ".efc")
    }
}
pub fn quantify(root: &Path, threads: usize, options: &QueryOptions) -> Result<()> {
    validate_query(options)?;
    ensure!(
        threads > 0 && options.jobs > 0 && options.batch_bases > 0 && options.unitig_k != Some(0),
        "invalid query worker, batch or unitig-k setting"
    );
    let (model, hash) = load(root)?;
    let n = model.features.len();
    drop(model);
    let datasets = input::datasets(&options.fof, options.mode)?;
    let meta = Metadata::load(root)?;
    let dict = Dictionary::load(root.join("dictionary"))
        .map_err(|e| anyhow::anyhow!("loading SSHash: {e}"))?;
    ensure!(
        dict.k() == meta.k
            && dict.m() == meta.minimizer
            && index::num_kmers(&dict)? == meta.indexed_kmers,
        "dictionary/model disagreement"
    );
    let owners = Owners::load(root, &meta)?;
    let stage = index::staging(&options.output)?;
    fs::create_dir(stage.path().join("datasets"))?;
    let jobs = options.jobs.min(threads).min(datasets.len());
    let pools = (0..jobs)
        .map(|id| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads / jobs + usize::from(id < threads % jobs))
                .build()
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let mut order: Vec<_> = (0..datasets.len()).collect();
    let sizes = datasets
        .iter()
        .map(|d| fs::metadata(&d.path).map(|m| m.len()))
        .collect::<std::io::Result<Vec<_>>>()?;
    order.sort_by_key(|&id| std::cmp::Reverse(sizes[id]));
    let next = AtomicUsize::new(0);
    let canceled = AtomicBool::new(false);
    let mut entries = vec![serde_json::Value::Null; datasets.len()];
    let mut failures = Vec::new();
    let mut global = vec![0u128; n];
    let mut total_kmers = 0u128;
    let extension = suffix(options.format, options.compression);
    eprintln!(
        "Quantifying {} datasets into exact disjoint exon/intersection features with {threads} workers ({jobs} concurrent datasets); no bit quantization",
        datasets.len()
    );
    std::thread::scope(|scope| -> Result<()> {
        let (tx, rx) = crossbeam_channel::bounded(jobs);
        for pool in &pools {
            let tx = tx.clone();
            let (datasets, order, next, canceled, dict, owners, meta) =
                (&datasets, &order, &next, &canceled, &dict, &owners, &meta);
            scope.spawn(move || {
                while !canceled.load(Ordering::Relaxed) {
                    let position = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&id) = order.get(position) else { break };
                    let result = sshash_lib::dispatch_on_k!(meta.k, K => crate::quantify::count::<K>(&datasets[id], dict, owners, n, pool, options));
                    if result.is_err() && !options.keep_going { canceled.store(true, Ordering::Relaxed); }
                    if tx.send((id, result)).is_err() { break; }
                }
            });
        }
        drop(tx);
        let mut error = None;
        for (id, result) in rx {
            let dataset = &datasets[id];
            match result {
                Ok(counts) if error.is_none() => {
                    let filename = format!("datasets/{}{extension}", dataset.name);
                    let publish = (|| -> Result<()> {
                        write_counts(
                            &stage.path().join(&filename),
                            options.compression,
                            options.format,
                            &counts.values,
                            dataset.mode.scale(),
                            &hash,
                        )?;
                        let factor = (abundance::SCALE / dataset.mode.scale()) as u128;
                        for (total, &value) in global.iter_mut().zip(&counts.values) {
                            *total = total
                                .checked_add(value as u128 * factor)
                                .context("global feature count overflow")?;
                        }
                        total_kmers = total_kmers
                            .checked_add(
                                counts
                                    .metrics
                                    .weighted_kmers
                                    .checked_mul(factor)
                                    .context("global k-mer total overflow")?,
                            )
                            .context("global k-mer total overflow")?;
                        entries[id] = json!({"name":dataset.name,"path":dataset.path,"mode":dataset.mode,"output":filename,"metrics":counts.metrics,
                            "total_kmers":{"numerator":counts.metrics.weighted_kmers.to_string(),"scale":dataset.mode.scale()}});
                        eprintln!("{}: exact exon/intersection features saved", dataset.name);
                        Ok(())
                    })();
                    if let Err(e) = publish {
                        error = Some(e);
                        canceled.store(true, Ordering::Relaxed);
                    }
                }
                Err(e) if options.keep_going && error.is_none() => {
                    let failure =
                        json!({"name":dataset.name,"path":dataset.path,"error":format!("{e:#}")});
                    eprintln!("DATASET_FAILED {failure}");
                    failures.push(failure);
                }
                Err(e) => {
                    error.get_or_insert(e);
                    canceled.store(true, Ordering::Relaxed);
                }
                _ => {}
            }
        }
        if let Some(error) = error {
            return Err(error);
        }
        Ok(())
    })?;
    entries.retain(|entry| !entry.is_null());
    ensure!(
        !entries.is_empty() && entries.len() + failures.len() == datasets.len(),
        "no successful datasets or incomplete feature quantification"
    );
    let global_name = format!("global{extension}");
    write_counts(
        &stage.path().join(&global_name),
        options.compression,
        options.format,
        &global,
        abundance::SCALE,
        &hash,
    )?;
    let manifest = json!({"version":1,"format":"isoform-features","exact":true,"model_sha256":crate::compact::hex(&hash),"index":root.canonicalize()?,"k":meta.k,"features":n,
        "compression":options.compression,"encoding":if options.format==output::Format::Compact {"EFC v1, sparse u128 fixed-point"} else {"disjoint-feature CSV"},
        "feature_semantics":"disjoint exon/intersection classes; all canonical k-mer occurrences counted once; shared sequence retained; unitig weights preserved before rounding",
        "datasets":entries,"global":{"output":global_name,"total_kmers":{"numerator":total_kmers.to_string(),"scale":abundance::SCALE}},
        "coverage":{"requested_datasets":datasets.len(),"completed_datasets":entries.len(),"failed_datasets":failures.len(),"complete":failures.is_empty()},"failed_datasets":failures});
    fs::write(
        stage.path().join("manifest.json"),
        serde_json::to_vec_pretty(&manifest)?,
    )?;
    index::publish(stage, &options.output)?;
    eprintln!("Exact feature results saved: {}", options.output.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_sparse_counts_all_codecs_and_corruption() {
        let dir = tempfile::tempdir().unwrap();
        for compression in [
            output::Compression::None,
            output::Compression::Gz,
            output::Compression::Xz,
            output::Compression::Zstd,
        ] {
            let path = dir.path().join("counts");
            write_counts(
                &path,
                compression,
                output::Format::Compact,
                &[0u128, 1, u128::MAX],
                abundance::SCALE,
                &[7; 32],
            )
            .unwrap();
            assert_eq!(
                read_counts(&path, 3, &[7; 32]).unwrap(),
                vec![0.0, 1e-6, u128::MAX as f64 / 1e6]
            );
            assert!(read_counts(&path, 3, &[8; 32]).is_err());
        }
        let path = dir.path().join("bad");
        write_counts(
            &path,
            output::Compression::None,
            output::Format::Compact,
            &[0u64, 7, 9],
            1,
            &[7; 32],
        )
        .unwrap();
        let mut raw = fs::read(&path).unwrap();
        raw[72] ^= 1;
        fs::write(&path, raw).unwrap();
        assert!(read_counts(&path, 3, &[7; 32]).is_err());
        assert_eq!(decimal(1, abundance::SCALE), "0.000001");
        assert_eq!(decimal(1250000, abundance::SCALE), "1.25");
    }
}
