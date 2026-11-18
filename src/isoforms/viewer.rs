//! Package inferred transcript values, retaining missing values and both metrics.
use crate::{efv, export, output};
use anyhow::{Context, Result, ensure};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::{BufWriter, Seek, SeekFrom, Write},
    path::Path,
};

type Values = (Vec<Value>, Vec<Option<f64>>, Vec<Option<f64>>);
fn read(path: &Path) -> Result<Values> {
    let mut csv = csv::Reader::from_reader(output::reader(path)?);
    let header = csv.headers()?.clone();
    let column = |name| {
        header
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("missing transcript column {name}"))
    };
    let id = column("transcript_id")?;
    let name = column("transcript_name")?;
    let gene = column("gene_name")?;
    let length = column("observable_kmers")?;
    let rate = column("abundance")?;
    let observations = column("assigned_observations")?;
    let mut targets = Vec::new();
    let mut rates = Vec::new();
    let mut counts = Vec::new();
    for row in csv.records() {
        let row = row?;
        ensure!(
            row[id].parse::<usize>()? == targets.len() + 1,
            "transcript IDs must be consecutive"
        );
        ensure!(
            !row[name].is_empty() && !row[gene].is_empty(),
            "empty transcript/gene name"
        );
        let exposure = row[length].parse::<u64>()?;
        targets.push(json!({"name":row[name],"gene":row[gene],"observable_kmers":exposure}));
        let value = |text: &str| -> Result<Option<f64>> {
            if text.is_empty() {
                return Ok(None);
            }
            let value: f64 = text.parse()?;
            ensure!(
                value.is_finite() && value >= 0.0,
                "invalid transcript abundance"
            );
            Ok(Some(value))
        };
        rates.push(value(&row[rate])?);
        counts.push(value(&row[observations])?);
    }
    ensure!(!targets.is_empty(), "no transcript targets");
    Ok((targets, rates, counts))
}

pub fn pack(input: &Path, destination: &Path, manifest: &Value) -> Result<()> {
    ensure!(!destination.exists(), "output already exists");
    ensure!(
        manifest["version"] == 1 && manifest["format"] == "isoform-abundance",
        "unsupported isoform results"
    );
    let root = input.canonicalize()?;
    let entries = manifest["datasets"]
        .as_array()
        .context("missing transcript datasets")?;
    ensure!(!entries.is_empty(), "no transcript datasets");
    let first = export::source(&root, export::text(&entries[0], "output")?)?;
    let (targets, _, _) = read(&first)?;
    let target_bytes = serde_json::to_vec(&targets)?;
    let reference: [u8; 32] = Sha256::digest(&target_bytes).into();
    let parent = destination
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    let stage = tempfile::NamedTempFile::new_in(parent)?;
    let mut out = BufWriter::new(stage.as_file());
    out.write_all(b"EXPREAI2")?;
    out.write_all(&[0; 16])?;
    let mut catalogue = Vec::new();
    let mut names = std::collections::HashSet::new();
    let mut vectors = entries
        .iter()
        .map(|entry| {
            Ok((
                export::text(entry, "name")?.to_owned(),
                export::text(entry, "output")?.to_owned(),
                "dataset",
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    vectors.push((
        "global sum".into(),
        manifest["global"]
            .as_str()
            .context("missing transcript global CSV")?
            .into(),
        "global",
    ));
    for (name, file, kind) in vectors {
        ensure!(
            kind != "dataset" || names.insert(name.clone()),
            "duplicate dataset name"
        );
        let (table, rates, counts) = read(&export::source(&root, &file)?)?;
        ensure!(table == targets, "inconsistent transcript table");
        let mut metrics = json!({});
        for (metric, values) in [("abundance", rates), ("assigned_observations", counts)] {
            let offset = out.stream_position()?;
            let mut encoder = zstd::stream::write::Encoder::new(&mut out, 3)?;
            efv::write(&mut encoder, &values, &reference)?;
            encoder.finish()?;
            let length = out.stream_position()? - offset;
            metrics[metric] = json!({"offset":offset,"length":length});
        }
        catalogue.push(json!({"name":name,"kind":kind,"global":kind=="global","metrics":metrics}));
    }
    let offset = out.stream_position()?;
    let coverage = manifest["source_coverage"]
        .get("coverage")
        .unwrap_or(&manifest["source_coverage"]);
    let catalogue = serde_json::to_vec(
        &json!({"version":2,"level":"transcript","encoding":"efv","reference":crate::compact::hex(&reference),
        "targets":targets,"datasets":catalogue,"metrics":["abundance","assigned_observations"],"coverage":coverage,
        "approximate":false,"model_sha256":manifest["model_sha256"],"failed_datasets":manifest["source_coverage"]["failed_datasets"],"global_semantics":"sum of independently fitted resolved transcript rates/assigned observations; missing stays missing"}),
    )?;
    ensure!(
        catalogue.len() <= 268_435_456,
        "transcript catalogue exceeds 256 MiB"
    );
    out.write_all(&catalogue)?;
    out.write_all(&Sha256::digest(&catalogue))?;
    out.seek(SeekFrom::Start(8))?;
    out.write_all(&offset.to_le_bytes())?;
    out.write_all(&(catalogue.len() as u64).to_le_bytes())?;
    out.flush()?;
    drop(out);
    stage.persist_noclobber(destination)?;
    eprintln!("Transcript viewer index saved: {}", destination.display());
    Ok(())
}
