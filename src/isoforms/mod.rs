//! Isoform inference from aggregate, disjoint exon/intersection abundance features.
mod counts;
mod em;
mod reference;
mod sort;
mod viewer;

use crate::{BuildArgs, QueryOptions, index, input, output};
use anyhow::{Context, Result, ensure};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use sshash_lib::{Dictionary, Kmer, KmerBits};
use std::{
    collections::HashMap,
    fs::{self, File},
    hash::{Hash, Hasher},
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

pub use counts::quantify;
pub use em::infer;
pub use viewer::pack;

#[derive(clap::Args)]
pub struct InferArgs {
    /// Index built with --isoforms. Inference reads the model, not the k-mer dictionary.
    #[arg(short, long)]
    index: PathBuf,
    /// Exact feature-result directory produced by quantify.
    #[arg(long, required_unless_present = "fof", conflicts_with = "fof")]
    input: Option<PathBuf>,
    /// List of exact disjoint-feature CSV/EFC files: name<TAB>path, or one path per line.
    #[arg(short = 'f', long, conflicts_with = "input")]
    fof: Option<PathBuf>,
    #[arg(short, long)]
    output: PathBuf,
    #[arg(short = 't', long, default_value_t = crate::default_threads())]
    threads: usize,
    /// Concurrent datasets; each shares the total thread budget with component fitting.
    #[arg(short = 'j', long, default_value_t = 1)]
    jobs: usize,
    #[arg(long, default_value_t = 10000)]
    max_iterations: usize,
    /// Relative rate-change tolerance, combined with an EM fixed-point check.
    #[arg(long, default_value_t = 1e-8)]
    tolerance: f64,
    /// Largest component checked for numerical rank. Larger ones are explicitly rank-unchecked.
    #[arg(long, default_value_t = 512)]
    rank_limit: usize,
    #[arg(long, value_enum, default_value_t = output::Compression::Zstd)]
    compression: output::Compression,
}

#[derive(Serialize, Deserialize)]
pub struct ModelRef {
    pub file: String,
    pub sha256: String,
    pub features: usize,
}
#[derive(Serialize, Deserialize)]
pub(super) struct Exon {
    name: String,
    contig: String,
    start: usize,
    end: usize,
    strand: char,
    annotation_ids: Vec<String>,
}
#[derive(Serialize, Deserialize)]
pub(super) struct Transcript {
    name: String,
    gene: String,
    contig: String,
    strand: char,
    exons: Vec<usize>,
    length: usize,
    observable_kmers: u64,
}
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(super) struct Occurrence {
    member: usize,
    multiplicity: u64,
}
#[derive(Serialize, Deserialize)]
pub(super) struct Feature {
    distinct_kmers: u64,
    exons: Vec<Occurrence>,
    transcripts: Vec<Occurrence>,
}
#[derive(Serialize, Deserialize)]
pub(super) struct Model {
    version: u32,
    k: usize,
    indexed_kmers: usize,
    exons: Vec<Exon>,
    transcripts: Vec<Transcript>,
    features: Vec<Feature>,
}

fn load(root: &Path) -> Result<(Model, [u8; 32])> {
    let meta = index::Metadata::load(root)?;
    let info = meta.isoforms.context("Expected an index built with --isoforms; ordinary exon/gene counts omit intersection evidence")?;
    ensure!(
        info.file == "isoform-model.json",
        "invalid isoform model path"
    );
    let raw = fs::read(root.join(&info.file))?;
    let hash: [u8; 32] = Sha256::digest(&raw).into();
    ensure!(
        crate::compact::hex(&hash) == info.sha256,
        "isoform model checksum mismatch"
    );
    let model: Model = serde_json::from_slice(&raw)?;
    ensure!(
        model.version == 1
            && model.k == meta.k
            && model.indexed_kmers == meta.indexed_kmers
            && model.features.len() == info.features
            && info.features == meta.exons.len(),
        "inconsistent isoform model/index"
    );
    model.validate()?;
    Ok((model, hash))
}
impl Model {
    fn validate(&self) -> Result<()> {
        ensure!(
            !self.features.is_empty() && !self.transcripts.is_empty() && !self.exons.is_empty(),
            "empty isoform model"
        );
        let mut totals = vec![0u64; self.transcripts.len()];
        let mut kmers = 0u64;
        for feature in &self.features {
            ensure!(
                feature.distinct_kmers > 0
                    && !feature.exons.is_empty()
                    && !feature.transcripts.is_empty(),
                "empty feature"
            );
            for (members, limit) in [
                (&feature.exons, self.exons.len()),
                (&feature.transcripts, self.transcripts.len()),
            ] {
                ensure!(
                    members
                        .iter()
                        .all(|m| m.member < limit && m.multiplicity > 0)
                        && members.windows(2).all(|p| p[0].member < p[1].member),
                    "invalid feature membership"
                );
            }
            kmers = kmers
                .checked_add(feature.distinct_kmers)
                .context("feature size overflow")?;
            for occurrence in &feature.transcripts {
                let length = occurrence
                    .multiplicity
                    .checked_mul(feature.distinct_kmers)
                    .context("exposure overflow")?;
                totals[occurrence.member] = totals[occurrence.member]
                    .checked_add(length)
                    .context("transcript exposure overflow")?;
            }
        }
        ensure!(
            kmers == self.indexed_kmers as u64,
            "feature k-mer total mismatch"
        );
        let mut names = std::collections::HashSet::new();
        for (id, transcript) in self.transcripts.iter().enumerate() {
            ensure!(
                !transcript.name.is_empty()
                    && !transcript.gene.is_empty()
                    && names.insert(&transcript.name)
                    && !transcript.exons.is_empty()
                    && transcript.exons.iter().all(|&e| e < self.exons.len())
                    && transcript.observable_kmers == totals[id],
                "invalid transcript metadata or observable length"
            );
        }
        Ok(())
    }
}

pub fn build(args: &BuildArgs) -> Result<()> {
    ensure!(
        args.level == crate::Level::Exon && args.gtf.is_some() && !args.genome.is_empty(),
        "--isoforms requires --gtf/--genome and the default exon level"
    );
    ensure!(
        args.junctions.junctions == crate::junctions::Mode::None,
        "--isoforms fits exon/intersection features; --junctions expansion is a separate gene mode"
    );
    ensure!(
        (3..=63).contains(&args.k) && args.k % 2 == 1 && args.threads > 0 && args.memory_gb > 0,
        "invalid k, threads or memory budget"
    );
    let stage = index::staging(&args.index)?;
    let (mut model, exon_path, transcript_path) = reference::extract(args, stage.path())?;
    let mut base = args.clone();
    base.isoforms = false;
    base.gtf = None;
    base.genome.clear();
    base.exons = Some(exon_path.clone());
    base.index = stage.path().join("base-index");
    index::build(&base)?;
    let mut meta = index::Metadata::load(&base.index)?;
    let dict = Dictionary::load(base.index.join("dictionary"))
        .map_err(|e| anyhow::anyhow!("loading SSHash: {e}"))?;
    model.indexed_kmers = meta.indexed_kmers;
    let scratch = if let Some(root) = &args.temp_dir {
        tempfile::tempdir_in(root)?
    } else {
        tempfile::tempdir_in(stage.path())?
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    // Half the hint is reserved for events; the dictionary/model also occupy RAM.
    let capacity = args
        .memory_gb
        .checked_mul(1 << 29)
        .context("memory budget overflow")?
        / std::mem::size_of::<sort::Event>();
    let mut sorter = sort::Sorter::new(scratch.path(), capacity, &pool)?;
    eprintln!("Collecting original exon memberships and transcript occurrence multiplicities");
    sshash_lib::dispatch_on_k!(args.k, K => collect::<K>(&exon_path, false, &dict, &mut sorter, &pool))?;
    sshash_lib::dispatch_on_k!(args.k, K => collect::<K>(&transcript_path, true, &dict, &mut sorter, &pool))?;
    let map_path = scratch.path().join("class-owners.bin");
    let mut out = BufWriter::new(File::create(&map_path)?);
    out.write_all(index::OWNER_MAGIC)?;
    out.write_all(&[32])?;
    let mut interns: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut exons = Vec::new();
    let mut transcripts = Vec::new();
    let mut current = 0u64;
    let mut seen = false;
    let mut emit = |exons: &mut Vec<Occurrence>, transcripts: &mut Vec<Occurrence>| -> Result<()> {
        ensure!(
            !exons.is_empty() && !transcripts.is_empty(),
            "dictionary k-mer lacks exon or transcript evidence"
        );
        let mut hash = std::collections::hash_map::DefaultHasher::new();
        exons.hash(&mut hash);
        transcripts.hash(&mut hash);
        let candidates = interns.entry(hash.finish()).or_default();
        let id = candidates.iter().copied().find(|&id| {
            model.features[id].exons == *exons && model.features[id].transcripts == *transcripts
        });
        let id = if let Some(id) = id {
            model.features[id].distinct_kmers += 1;
            exons.clear();
            transcripts.clear();
            id
        } else {
            let id = model.features.len();
            ensure!(id < u32::MAX as usize - 1, "too many intersection features");
            model.features.push(Feature {
                distinct_kmers: 1,
                exons: std::mem::take(exons),
                transcripts: std::mem::take(transcripts),
            });
            candidates.push(id);
            id
        };
        out.write_all(&(id as u32 + 1).to_le_bytes())?;
        Ok(())
    };
    sorter.finish(|event| {
        if seen && event.kmer != current {
            emit(&mut exons, &mut transcripts)?;
            ensure!(
                event.kmer == current + 1,
                "nonconsecutive dictionary k-mer IDs"
            );
            current = event.kmer;
        } else if !seen {
            ensure!(event.kmer == 0, "missing first dictionary k-mer");
            seen = true;
        }
        let target = if event.member & 1 == 0 {
            &mut exons
        } else {
            &mut transcripts
        };
        target.push(Occurrence {
            member: usize::try_from(event.member / 2)?,
            multiplicity: event.count,
        });
        Ok(())
    })?;
    emit(&mut exons, &mut transcripts)?;
    ensure!(
        seen && current + 1 == model.indexed_kmers as u64,
        "missing dictionary feature assignments"
    );
    out.flush()?;
    drop(out);
    drop(interns);
    drop(dict);
    meta.owner_bits = if model.features.len() < u16::MAX as usize {
        16
    } else {
        32
    };
    if meta.owner_bits == 16 {
        let mut input = BufReader::new(File::open(&map_path)?);
        let mut header = [0; 9];
        input.read_exact(&mut header)?;
        let mut out = BufWriter::new(File::create(base.index.join("owners.bin"))?);
        out.write_all(index::OWNER_MAGIC)?;
        out.write_all(&[16])?;
        for _ in 0..model.indexed_kmers {
            let mut value = [0; 4];
            input.read_exact(&mut value)?;
            out.write_all(&u16::try_from(u32::from_le_bytes(value))?.to_le_bytes())?;
        }
        out.flush()?;
    } else {
        fs::copy(&map_path, base.index.join("owners.bin"))?;
    }
    for feature in &model.features {
        for m in &feature.transcripts {
            let length = m
                .multiplicity
                .checked_mul(feature.distinct_kmers)
                .context("exposure overflow")?;
            let total = &mut model.transcripts[m.member].observable_kmers;
            *total = total
                .checked_add(length)
                .context("transcript exposure overflow")?;
        }
    }
    model.validate()?;
    let raw = serde_json::to_vec(&model)?;
    let hash: [u8; 32] = Sha256::digest(&raw).into();
    fs::write(stage.path().join("isoform-model.json"), raw)?;
    write_tables(&model, stage.path())?;
    meta.exons = model
        .features
        .iter()
        .enumerate()
        .map(|(id, f)| index::Exon {
            name: format!("feature_{:06}", id + 1),
            length: 0,
            unique_kmers: f.distinct_kmers,
            kmer_length: Some(f.distinct_kmers),
        })
        .collect();
    meta.duplicate_kmers = 0;
    meta.isoforms = Some(ModelRef {
        file: "isoform-model.json".into(),
        sha256: crate::compact::hex(&hash),
        features: model.features.len(),
    });
    fs::write(
        base.index.join("metadata.json"),
        serde_json::to_vec_pretty(&meta)?,
    )?;
    for file in fs::read_dir(&base.index)? {
        let file = file?;
        fs::rename(file.path(), stage.path().join(file.file_name()))?;
    }
    fs::remove_dir(&base.index)?;
    drop(scratch);
    index::publish(stage, &args.index)?;
    eprintln!(
        "Isoform index saved: {} transcripts, {} disjoint exon/intersection features; all {} canonical k-mers retained",
        model.transcripts.len(),
        model.features.len(),
        model.indexed_kmers
    );
    Ok(())
}
fn collect<const K: usize>(
    path: &Path,
    transcript: bool,
    dict: &Dictionary,
    sorter: &mut sort::Sorter<'_>,
    pool: &rayon::ThreadPool,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let capacity = (sorter.capacity() / pool.current_num_threads()).max(1);
    let root = sorter.root().join(if transcript {
        "transcript-runs"
    } else {
        "exon-runs"
    });
    let runs = input::batches_with(
        path,
        1 << 20,
        2 * pool.current_num_threads(),
        K - 1,
        |header, _| {
            let id: u64 = std::str::from_utf8(header)?.parse()?;
            id.checked_mul(2)
                .and_then(|id| id.checked_add(u64::from(transcript)))
                .context("reference ID overflow")
        },
        |rx| {
            let stopped = std::sync::atomic::AtomicBool::new(false);
            pool.broadcast(|worker| {
                let result = (|| -> Result<Vec<PathBuf>> {
                    let mut local =
                        sort::Sorter::new(&root.join(worker.index().to_string()), capacity, pool)?;
                    let mut query = dict.create_streaming_query::<K>();
                    while !stopped.load(std::sync::atomic::Ordering::Relaxed) {
                        let Ok(batch) = rx.recv() else { break };
                        for (record, seq) in batch.records() {
                            for run in seq
                                .split(|b| !input::is_dna(*b))
                                .filter(|run| run.len() >= K)
                            {
                                query.reset();
                                for window in run.windows(K) {
                                    let hit = query.lookup(window);
                                    if hit.is_found() {
                                        local.push(hit.kmer_id, record.value)?;
                                    } else {
                                        ensure!(transcript, "exon k-mer absent from dictionary");
                                    }
                                }
                            }
                        }
                    }
                    local.into_runs()
                })();
                if result.is_err() {
                    stopped.store(true, std::sync::atomic::Ordering::Relaxed);
                }
                result
            })
            .into_iter()
            .collect::<Result<Vec<_>>>()
        },
    )?;
    for runs in runs {
        sorter.add_runs(runs);
    }
    Ok(())
}
fn write_tables(model: &Model, root: &Path) -> Result<()> {
    output::csv(
        &root.join("features.csv.zst"),
        output::Compression::Zstd,
        |csv| {
            csv.write_record([
                "feature_id",
                "distinct_kmers",
                "exon_ids",
                "exon_multiplicities",
                "transcript_ids",
                "transcript_multiplicities",
            ])?;
            for (id, f) in model.features.iter().enumerate() {
                let ids = |members: &[Occurrence]| {
                    members
                        .iter()
                        .map(|m| (m.member + 1).to_string())
                        .collect::<Vec<_>>()
                        .join(";")
                };
                let mult = |members: &[Occurrence]| {
                    members
                        .iter()
                        .map(|m| m.multiplicity.to_string())
                        .collect::<Vec<_>>()
                        .join(";")
                };
                csv.serialize((
                    id + 1,
                    f.distinct_kmers,
                    ids(&f.exons),
                    mult(&f.exons),
                    ids(&f.transcripts),
                    mult(&f.transcripts),
                ))?;
            }
            Ok(())
        },
    )?;
    output::csv(
        &root.join("transcripts.csv.zst"),
        output::Compression::Zstd,
        |csv| {
            csv.write_record([
                "transcript_id",
                "transcript_name",
                "gene_name",
                "length",
                "observable_kmers",
                "exon_ids",
            ])?;
            for (id, t) in model.transcripts.iter().enumerate() {
                csv.serialize((
                    id + 1,
                    &t.name,
                    &t.gene,
                    t.length,
                    t.observable_kmers,
                    t.exons
                        .iter()
                        .map(|id| (id + 1).to_string())
                        .collect::<Vec<_>>()
                        .join(";"),
                ))?;
            }
            Ok(())
        },
    )?;
    output::csv(
        &root.join("original-exons.csv.zst"),
        output::Compression::Zstd,
        |csv| {
            csv.write_record([
                "exon_id",
                "exon_name",
                "contig",
                "start",
                "end",
                "strand",
                "annotation_ids",
            ])?;
            for (id, e) in model.exons.iter().enumerate() {
                csv.serialize((
                    id + 1,
                    &e.name,
                    &e.contig,
                    e.start,
                    e.end,
                    e.strand,
                    e.annotation_ids.join(";"),
                ))?;
            }
            Ok(())
        },
    )
}
pub fn validate_query(options: &QueryOptions) -> Result<()> {
    ensure!(
        !options.viewer_index && !options.stats,
        "isoform feature counts are exact EFC/CSV evidence; --viewer-index and --stats apply to ordinary exon/gene vectors"
    );
    Ok(())
}
