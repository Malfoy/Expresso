use crate::{
    ExportArgs, compact,
    index::{self, Exon},
    output,
};
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use serde_json::Value;
use std::{
    collections::HashSet,
    fs,
    path::{Component, Path, PathBuf},
};

fn source(root: &Path, name: &str) -> Result<PathBuf> {
    let path = Path::new(name);
    ensure!(
        !name.is_empty() && path.components().all(|c| matches!(c, Component::Normal(_))),
        "Invalid compact file path"
    );
    let path = root.join(path).canonicalize()?;
    ensure!(
        path.starts_with(root),
        "Compact input escapes result directory"
    );
    Ok(path)
}

fn text<'a>(value: &'a Value, key: &str) -> Result<&'a str> {
    value[key]
        .as_str()
        .with_context(|| format!("Missing manifest field: {key}"))
}

pub fn run(args: &ExportArgs) -> Result<()> {
    ensure!(args.threads > 0, "threads must be positive");
    let input = args.input.canonicalize()?;
    let manifest: Value = serde_json::from_reader(output::reader(&input.join("manifest.json"))?)?;
    ensure!(
        manifest["format"] == "compact" && manifest["compact_format_version"] == 1,
        "Expected compact EXPRESSO v1 results"
    );
    let reference_path = source(&input, text(&manifest["reference"], "file")?)?;
    let mut reader = csv::Reader::from_reader(output::reader(&reference_path)?);
    ensure!(
        reader
            .headers()?
            .iter()
            .eq(["exon_id", "exon_name", "length", "unique_kmers"]),
        "Invalid reference table columns"
    );
    let mut exons = Vec::new();
    for row in reader.deserialize::<(usize, String, usize, u64)>() {
        let (id, name, length, unique_kmers) = row?;
        ensure!(
            id == exons.len() + 1,
            "Reference IDs must be consecutive and 1-based"
        );
        exons.push(Exon {
            name,
            length,
            unique_kmers,
        });
    }
    let reference = compact::reference_id(&exons);
    ensure!(
        manifest["reference"]["targets"].as_u64() == Some(exons.len() as u64)
            && manifest["reference"]["sha256"].as_str() == Some(&compact::hex(&reference)),
        "Reference table hash/length mismatch"
    );
    let entries = manifest["datasets"]
        .as_array()
        .context("Missing datasets")?;
    let wanted: HashSet<_> = args.dataset.iter().map(String::as_str).collect();
    let mut available = HashSet::new();
    let mut tasks = Vec::new();
    for entry in entries {
        let name = text(entry, "name")?;
        ensure!(
            !name.is_empty()
                && name != "."
                && name != ".."
                && name
                    .bytes()
                    .all(|c| c.is_ascii_alphanumeric() || b"._-".contains(&c))
                && available.insert(name),
            "Invalid/repeated dataset name"
        );
        if !args.global_only && (wanted.is_empty() || wanted.contains(name)) {
            tasks.push((
                entry,
                format!("datasets/{name}{}", args.compression.suffix()),
                false,
            ));
        }
    }
    ensure!(
        wanted.is_subset(&available),
        "An requested dataset is not in this result set"
    );
    if wanted.is_empty() {
        tasks.push((
            &manifest["global"],
            format!("global{}", args.compression.suffix()),
            true,
        ));
    }
    let stage = index::staging(&args.output)?;
    fs::create_dir(stage.path().join("datasets"))?;
    fs::copy(&reference_path, stage.path().join("exons.csv.zst"))?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    pool.install(|| {
        tasks
            .par_iter()
            .try_for_each(|(entry, output_name, global)| -> Result<()> {
                let values = compact::read_vector(
                    &source(&input, text(entry, "output")?)?,
                    exons.len(),
                    &reference,
                )?;
                let mut statistics = Vec::new();
                if *global && let Some(columns) = entry["statistics"].as_array() {
                    let mut names = HashSet::new();
                    for column in columns {
                        let name = text(column, "name")?;
                        ensure!(
                            ["mean", "median", "min", "max", "datasets_detected"].contains(&name)
                                && names.insert(name),
                            "Invalid statistic name"
                        );
                        let values = compact::read_vector(
                            &source(&input, text(column, "output")?)?,
                            exons.len(),
                            &reference,
                        )?;
                        statistics.push((name, values));
                    }
                }
                output::csv(&stage.path().join(output_name), args.compression, |csv| {
                    let mut headers = Vec::new();
                    if args.with_names {
                        headers.extend(["exon_id", "exon_name"]);
                    }
                    headers.push("abundance");
                    headers.extend(statistics.iter().map(|(name, _)| *name));
                    csv.write_record(headers)?;
                    for (i, &value) in values.iter().enumerate() {
                        if statistics.is_empty() {
                            if args.with_names {
                                csv.serialize((i + 1, &exons[i].name, value))?;
                            } else {
                                csv.serialize((value,))?;
                            }
                        } else {
                            let mut row = Vec::new();
                            if args.with_names {
                                row.extend([(i + 1).to_string(), exons[i].name.clone()]);
                            }
                            row.push(value.to_string());
                            row.extend(statistics.iter().map(|(_, values)| values[i].to_string()));
                            csv.write_record(row)?;
                        }
                    }
                    Ok(())
                })
            })
    })?;
    let report = serde_json::json!({
        "source":input, "source_format":"EAB v1", "format":"csv", "approximate":true,
        "with_names":args.with_names, "reference":{"file":"exons.csv.zst","sha256":compact::hex(&reference),"targets":exons.len()},
        "compression":args.compression, "statistics":manifest["statistics"],
        "source_coverage":manifest["coverage"], "source_failed_datasets":manifest["failed_datasets"],
        "files":tasks.iter().map(|(entry, output, global)| serde_json::json!({"output":output,"global":global,
            "name":entry["name"],"quantization":entry["quantization"],"statistics":entry["statistics"]})).collect::<Vec<_>>()
    });
    fs::write(
        stage.path().join("manifest.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    index::publish(stage, &args.output)?;
    eprintln!(
        "Exported {} CSV files: {}",
        tasks.len(),
        args.output.display()
    );
    Ok(())
}
