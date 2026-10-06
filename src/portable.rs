//! EAI v1: independently compressed EAB vectors and a seekable JSON catalogue.
use crate::{compact, export};
use anyhow::{Context, Result, ensure};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::{BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
};

pub fn pack(input: &Path, destination: &Path) -> Result<()> {
    ensure!(
        !destination.exists(),
        "output already exists: {}",
        destination.display()
    );
    let results = export::load(input)?;
    let parent = destination
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    let stage = tempfile::NamedTempFile::new_in(parent)?;
    let mut out = BufWriter::new(stage.as_file());
    out.write_all(b"EXPREAI1")?;
    out.write_all(&[0; 16])?;
    let entries = results.manifest["datasets"]
        .as_array()
        .context("Missing datasets")?;
    let mut catalogue = Vec::new();
    let mut names = std::collections::HashSet::new();
    let mut vectors: Vec<(String, &serde_json::Value, &str)> = entries
        .iter()
        .map(|e| Ok((export::text(e, "name")?.to_owned(), e, "dataset")))
        .collect::<Result<_>>()?;
    // Global is identified by kind, so a dataset named "global" remains unambiguous.
    vectors.push(("global sum".into(), &results.manifest["global"], "global"));
    if let Some(statistics) = results.manifest["global"]["statistics"].as_array() {
        for statistic in statistics {
            vectors.push((
                export::text(statistic, "name")?.to_owned(),
                statistic,
                "statistic",
            ));
        }
    }
    for (name, entry, kind) in vectors {
        ensure!(
            kind != "dataset" || names.insert(name.clone()),
            "Repeated dataset name"
        );
        // Verify integrity before publishing; recompress encoded bytes, never quantize again.
        let limit = (results.exons.len() as u64)
            .checked_mul(2)
            .and_then(|n| n.checked_add(524_364))
            .context("Vector too large")?;
        let mut raw = Vec::new();
        export::vector_reader(&results.input, entry)?
            .take(limit + 1)
            .read_to_end(&mut raw)?;
        ensure!(raw.len() as u64 <= limit, "Encoded vector too large");
        crate::eab::decode(raw.as_slice(), results.exons.len(), &results.reference)?;
        let offset = out.stream_position()?;
        let mut encoder = zstd::stream::write::Encoder::new(&mut out, 3)?;
        encoder.write_all(&raw)?;
        encoder.finish()?;
        let length = out.stream_position()? - offset;
        catalogue
            .push(serde_json::json!({"name":name,"kind":kind,"global":kind == "global","offset":offset,"length":length}));
    }
    let offset = out.stream_position()?;
    let catalogue = serde_json::to_vec(&serde_json::json!({
        "version":1,"level":results.level,"reference":compact::hex(&results.reference),
        "targets":results.exons,"datasets":catalogue,"coverage":results.manifest["coverage"],
        "failed_datasets":results.manifest["failed_datasets"],
        "approximate":true,"global_semantics":"sum of exact rounded counts before independent quantization"
    }))?;
    ensure!(
        catalogue.len() <= 268_435_456,
        "EAI v1 catalogue exceeds 256 MiB; use the result directory for CLI export"
    );
    out.write_all(&catalogue)?;
    out.write_all(&Sha256::digest(&catalogue))?;
    out.seek(SeekFrom::Start(8))?;
    out.write_all(&offset.to_le_bytes())?;
    out.write_all(&(catalogue.len() as u64).to_le_bytes())?;
    out.flush()?;
    drop(out);
    stage.persist_noclobber(destination)?;
    eprintln!("Viewer index saved: {}", destination.display());
    Ok(())
}

pub(crate) fn load(path: &Path) -> Result<export::Results> {
    let mut file = fs::File::open(path)?;
    let size = file.metadata()?.len();
    let mut header = [0; 24];
    file.read_exact(&mut header)?;
    ensure!(&header[..8] == b"EXPREAI1", "Invalid EAI magic/version");
    let offset = u64::from_le_bytes(header[8..16].try_into()?);
    let length = u64::from_le_bytes(header[16..24].try_into()?);
    ensure!(
        offset >= 24
            && length <= 268_435_456
            && offset.checked_add(length).and_then(|n| n.checked_add(32)) == Some(size),
        "Invalid EAI catalogue size"
    );
    file.seek(SeekFrom::Start(offset))?;
    let mut raw = vec![0; usize::try_from(length)?];
    file.read_exact(&mut raw)?;
    let mut checksum = [0; 32];
    file.read_exact(&mut checksum)?;
    ensure!(
        Sha256::digest(&raw).as_slice() == checksum,
        "EAI catalogue checksum mismatch"
    );
    let catalogue: serde_json::Value = serde_json::from_slice(&raw)?;
    ensure!(catalogue["version"] == 1, "Unsupported EAI version");
    let level = serde_json::from_value(catalogue["level"].clone())?;
    let exons: Vec<crate::index::Exon> = serde_json::from_value(catalogue["targets"].clone())?;
    ensure!(!exons.is_empty(), "EAI contains no targets");
    let reference = compact::reference_id(&exons);
    ensure!(
        catalogue["reference"] == compact::hex(&reference),
        "EAI reference hash mismatch"
    );
    let mut datasets = Vec::new();
    let mut statistics = Vec::new();
    let mut global = None;
    let mut end = 24;
    for entry in catalogue["datasets"]
        .as_array()
        .context("Missing EAI datasets")?
    {
        let start = entry["offset"].as_u64().context("Invalid EAI offset")?;
        let length = entry["length"].as_u64().context("Invalid EAI length")?;
        ensure!(
            start == end && length > 0 && start.checked_add(length).is_some_and(|n| n <= offset),
            "Invalid EAI dataset bounds"
        );
        end = start + length;
        match export::text(entry, "kind")? {
            "dataset" => datasets.push(entry.clone()),
            "global" => {
                ensure!(global.is_none(), "Repeated EAI global sum");
                global = Some(entry.clone());
            }
            "statistic" => statistics.push(entry.clone()),
            _ => anyhow::bail!("Invalid EAI entry kind"),
        }
    }
    ensure!(end == offset, "Inconsistent EAI data size");
    let mut global = global.context("EAI global sum missing")?;
    global["statistics"] = statistics.into();
    let manifest = serde_json::json!({"format":"compact", "compact_format_version":1,"level":level,
        "global":global,"datasets":datasets,"coverage":catalogue["coverage"],"failed_datasets":catalogue["failed_datasets"],"statistics":catalogue["global_semantics"]});
    Ok(export::Results {
        input: path.to_owned(),
        manifest,
        exons,
        reference,
        level,
    })
}

pub(crate) fn reader(path: &Path, entry: &serde_json::Value) -> Result<Box<dyn Read>> {
    let mut file = fs::File::open(path)?;
    file.seek(SeekFrom::Start(
        entry["offset"].as_u64().context("Invalid EAI offset")?,
    ))?;
    let slice = file.take(entry["length"].as_u64().context("Invalid EAI length")?);
    Ok(Box::new(zstd::stream::read::Decoder::new(slice)?))
}
