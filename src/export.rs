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

pub(crate) fn source(root: &Path, name: &str) -> Result<PathBuf> {
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

pub(crate) fn text<'a>(value: &'a Value, key: &str) -> Result<&'a str> {
    value[key]
        .as_str()
        .with_context(|| format!("Missing manifest field: {key}"))
}

pub(crate) struct Results {
    pub input: PathBuf,
    pub manifest: Value,
    pub exons: Vec<Exon>,
    pub reference: [u8; 32],
    pub level: crate::Level,
}
pub(crate) fn load(path: &Path) -> Result<Results> {
    let input = path.canonicalize()?;
    if input.is_file() {
        return crate::portable::load(&input);
    }
    let manifest: Value = serde_json::from_reader(output::reader(&input.join("manifest.json"))?)?;
    ensure!(
        manifest["format"] == "compact" && manifest["compact_format_version"] == 1,
        "Expected compact EXPRESSO v1 results"
    );
    let level: crate::Level = serde_json::from_value(
        manifest
            .get("level")
            .cloned()
            .unwrap_or(serde_json::json!("exon")),
    )?;
    let labels = level.columns();
    let reference_path = source(&input, text(&manifest["reference"], "file")?)?;
    let mut reader = csv::Reader::from_reader(output::reader(&reference_path)?);
    ensure!(
        reader
            .headers()?
            .iter()
            .eq([labels[0], labels[1], "length", "unique_kmers"]),
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
    Ok(Results {
        input,
        manifest,
        exons,
        reference,
        level,
    })
}

fn list(values: &[String], file: &Option<PathBuf>) -> Result<HashSet<String>> {
    let mut selected: HashSet<String> = values.iter().cloned().collect();
    if let Some(file) = file {
        use std::io::{BufRead, BufReader};
        for line in BufReader::new(output::reader(file)?).lines() {
            let line = line?;
            let value = line.trim();
            if !value.is_empty() && !value.starts_with('#') {
                selected.insert(value.to_owned());
            }
        }
    }
    Ok(selected)
}

pub fn run(args: &ExportArgs) -> Result<()> {
    ensure!(args.threads > 0, "threads must be positive");
    ensure!(
        args.min_value.unwrap_or(0) <= args.max_value.unwrap_or(u64::MAX),
        "min-value exceeds max-value"
    );
    let Results {
        input,
        manifest,
        exons,
        reference,
        level,
    } = load(&args.input)?;
    let labels = level.columns();
    let target_names = list(&args.target, &args.target_list)?;
    let ids: HashSet<usize> = args.target_id.iter().copied().collect();
    ensure!(
        ids.iter().all(|&id| id > 0 && id <= exons.len()),
        "target ID out of range"
    );
    let available_names: HashSet<_> = exons.iter().map(|e| e.name.clone()).collect();
    ensure!(
        target_names.is_subset(&available_names),
        "requested target name not found"
    );
    let subset = !target_names.is_empty() || !ids.is_empty();
    ensure!(args.target_list.is_none() || subset, "target list is empty");
    let selected: Vec<usize> = exons
        .iter()
        .enumerate()
        .filter_map(|(i, e)| {
            (!subset || target_names.contains(&e.name) || ids.contains(&(i + 1))).then_some(i)
        })
        .collect();
    let with_names =
        args.with_names || subset || args.min_value.is_some() || args.max_value.is_some();
    let entries = manifest["datasets"]
        .as_array()
        .context("Missing datasets")?;
    let wanted = list(&args.dataset, &args.dataset_list)?;
    ensure!(
        args.dataset_list.is_none() || !wanted.is_empty(),
        "dataset list is empty"
    );
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
                && available.insert(name.to_owned()),
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
        "A requested dataset is not in this result set"
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
    let reference_file = format!("{}.zst", level.table());
    output::csv(
        &stage.path().join(&reference_file),
        output::Compression::Zstd,
        |csv| {
            csv.write_record([labels[0], labels[1], "length", "unique_kmers"])?;
            for &i in &selected {
                let e = &exons[i];
                csv.serialize((i + 1, &e.name, e.length, e.unique_kmers))?;
            }
            Ok(())
        },
    )?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    pool.install(|| {
        tasks
            .par_iter()
            .try_for_each(|(entry, output_name, global)| -> Result<()> {
                let values =
                    crate::eab::decode(vector_reader(&input, entry)?, exons.len(), &reference)?;
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
                        let values = crate::eab::decode(
                            vector_reader(&input, column)?,
                            exons.len(),
                            &reference,
                        )?;
                        statistics.push((name, values));
                    }
                }
                output::csv(&stage.path().join(output_name), args.compression, |csv| {
                    let mut headers = Vec::new();
                    if with_names {
                        headers.extend(labels);
                    }
                    headers.push("abundance");
                    headers.extend(statistics.iter().map(|(name, _)| *name));
                    csv.write_record(headers)?;
                    for &i in &selected {
                        let value = values[i];
                        if value < args.min_value.unwrap_or(0)
                            || value > args.max_value.unwrap_or(u64::MAX)
                        {
                            continue;
                        }
                        if statistics.is_empty() {
                            if with_names {
                                csv.serialize((i + 1, &exons[i].name, value))?;
                            } else {
                                csv.serialize((value,))?;
                            }
                        } else {
                            let mut row = Vec::new();
                            if with_names {
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
        "level":level, "with_names":with_names,
        "filters":{"targets_selected":selected.len(), "min_value":args.min_value,"max_value":args.max_value},
        "reference":{"file":reference_file,"source_sha256":compact::hex(&reference),"targets":selected.len(), "ids":"original 1-based IDs; may be nonconsecutive"},
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

pub(crate) fn vector_reader(input: &Path, entry: &Value) -> Result<Box<dyn std::io::Read>> {
    if input.is_file() {
        crate::portable::reader(input, entry)
    } else {
        output::reader(&source(input, text(entry, "output")?)?)
    }
}
