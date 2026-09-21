use crate::{BuildArgs, input};
use anyhow::{Context, Result, ensure};
use ggcat_api::{ExtraElaboration, GGCATConfig, GGCATInstance, GeneralSequenceBlockData};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sshash_lib::{BuildConfiguration, Dictionary, DictionaryBuilder, Kmer, KmerBits};
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Write},
    path::Path,
    sync::atomic::{AtomicU16, AtomicU32, Ordering},
};

const FORMAT_VERSION: u32 = 1;
const OWNER_MAGIC: &[u8; 8] = b"EXPRown1";

#[derive(Serialize, Deserialize)]
pub struct Exon {
    pub name: String,
    pub length: usize,
    pub unique_kmers: u64,
}

#[derive(Serialize, Deserialize)]
pub struct Metadata {
    pub format_version: u32,
    pub expresso_version: String,
    pub sshash_version: String,
    pub k: usize,
    pub minimizer: usize,
    pub owner_bits: u8,
    pub indexed_kmers: usize,
    pub duplicate_kmers: usize,
    pub exons: Vec<Exon>,
}

impl Metadata {
    pub fn load(dir: &Path) -> Result<Self> {
        let meta: Self =
            serde_json::from_reader(BufReader::new(File::open(dir.join("metadata.json"))?))?;
        ensure!(
            meta.format_version == FORMAT_VERSION && meta.sshash_version == "0.7.1",
            "unsupported index version; rebuild the index"
        );
        ensure!(
            (3..=63).contains(&meta.k) && meta.k % 2 == 1,
            "invalid index k"
        );
        ensure!(!meta.exons.is_empty(), "index contains no exons");
        Ok(meta)
    }
}

// One reserved value for unassigned, one for shared. No auxiliary hash map.
pub enum Owners {
    U16(Vec<u16>),
    U32(Vec<u32>),
}

impl Owners {
    #[inline]
    pub fn exon(&self, kmer: usize) -> Option<usize> {
        match self {
            Self::U16(v) => {
                let x = v[kmer];
                (x != 0 && x != u16::MAX).then(|| x as usize - 1)
            }
            Self::U32(v) => {
                let x = v[kmer];
                (x != 0 && x != u32::MAX).then(|| x as usize - 1)
            }
        }
    }
    fn write(&self, path: &Path) -> Result<()> {
        let mut out = BufWriter::new(File::create(path)?);
        out.write_all(OWNER_MAGIC)?;
        match self {
            Self::U16(v) => {
                out.write_all(&[16])?;
                for x in v {
                    out.write_all(&x.to_le_bytes())?;
                }
            }
            Self::U32(v) => {
                out.write_all(&[32])?;
                for x in v {
                    out.write_all(&x.to_le_bytes())?;
                }
            }
        }
        out.flush()?;
        Ok(())
    }
    pub fn load(dir: &Path, meta: &Metadata) -> Result<Self> {
        let file = File::open(dir.join("owners.bin"))?;
        ensure!(matches!(meta.owner_bits, 16 | 32), "invalid owner width");
        let expected = (meta.indexed_kmers as u64)
            .checked_mul(meta.owner_bits as u64 / 8)
            .and_then(|n| n.checked_add(9))
            .context("invalid owner table size")?;
        ensure!(
            file.metadata()?.len() == expected,
            "truncated or inconsistent owner table"
        );
        let mut reader = BufReader::new(file);
        let mut header = [0; 9];
        reader.read_exact(&mut header)?;
        ensure!(
            &header[..8] == OWNER_MAGIC && header[8] == meta.owner_bits,
            "invalid owner table header"
        );
        let owners = if meta.owner_bits == 16 {
            let mut v = Vec::with_capacity(meta.indexed_kmers);
            for _ in 0..meta.indexed_kmers {
                let mut b = [0; 2];
                reader.read_exact(&mut b)?;
                v.push(u16::from_le_bytes(b));
            }
            Self::U16(v)
        } else {
            let mut v = Vec::with_capacity(meta.indexed_kmers);
            for _ in 0..meta.indexed_kmers {
                let mut b = [0; 4];
                reader.read_exact(&mut b)?;
                v.push(u32::from_le_bytes(b));
            }
            Self::U32(v)
        };
        for i in 0..meta.indexed_kmers {
            ensure!(
                owners.exon(i).is_none_or(|e| e < meta.exons.len()),
                "invalid exon ID in owner table"
            );
        }
        Ok(owners)
    }
}

enum AtomicOwners {
    U16(Vec<AtomicU16>),
    U32(Vec<AtomicU32>),
}
impl AtomicOwners {
    fn new(kmers: usize, exons: usize) -> Result<Self> {
        ensure!(exons < u32::MAX as usize, "too many exons for 32-bit IDs");
        Ok(if exons < u16::MAX as usize {
            Self::U16((0..kmers).map(|_| AtomicU16::new(0)).collect())
        } else {
            Self::U32((0..kmers).map(|_| AtomicU32::new(0)).collect())
        })
    }
    fn assign(&self, kmer: usize, exon: u64) {
        macro_rules! assign {
            ($v:expr, $ty:ty) => {{
                let slot = &$v[kmer];
                let id = exon as $ty;
                let mut old = slot.load(Ordering::Relaxed);
                loop {
                    if old == id || old == <$ty>::MAX {
                        break;
                    }
                    let new = if old == 0 { id } else { <$ty>::MAX };
                    match slot.compare_exchange_weak(old, new, Ordering::Relaxed, Ordering::Relaxed)
                    {
                        Ok(_) => break,
                        Err(current) => old = current,
                    }
                }
            }};
        }
        match self {
            Self::U16(v) => assign!(v, u16),
            Self::U32(v) => assign!(v, u32),
        }
    }
    fn finish(self) -> Owners {
        match self {
            Self::U16(v) => Owners::U16(v.into_iter().map(AtomicU16::into_inner).collect()),
            Self::U32(v) => Owners::U32(v.into_iter().map(AtomicU32::into_inner).collect()),
        }
    }
}

pub fn num_kmers(dict: &Dictionary) -> Result<usize> {
    usize::try_from(dict.spss().total_bases() - dict.num_strings() * (dict.k() as u64 - 1))
        .context("index too large for this platform")
}

/// Staging in the destination parent permits a final atomic directory rename.
pub fn staging(destination: &Path) -> Result<tempfile::TempDir> {
    ensure!(
        !destination.exists(),
        "destination already exists: {}",
        destination.display()
    );
    let parent = destination
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    Ok(tempfile::Builder::new()
        .prefix(".expresso-")
        .tempdir_in(parent)?)
}

pub fn publish(stage: tempfile::TempDir, destination: &Path) -> Result<()> {
    ensure!(
        !destination.exists(),
        "destination appeared while working: {}",
        destination.display()
    );
    fs::rename(stage.path(), destination)?;
    Ok(())
}

pub fn build(args: &BuildArgs) -> Result<()> {
    ensure!(
        args.threads > 0 && args.memory_gb > 0,
        "threads and memory-gb must be positive"
    );
    let m = args.minimizer.unwrap_or(19.min(args.k.saturating_sub(2)));
    let mut config = BuildConfiguration::new(args.k, m).map_err(anyhow::Error::msg)?;
    let stage = staging(&args.index)?;
    let scratch = if let Some(dir) = &args.temp_dir {
        fs::create_dir_all(dir)?;
        tempfile::Builder::new()
            .prefix("expresso-build-")
            .tempdir_in(dir)?
    } else {
        tempfile::Builder::new()
            .prefix("scratch-")
            .tempdir_in(stage.path())?
    };
    let normalized = scratch.path().join("exons.fa");
    let mut out = BufWriter::new(File::create(&normalized)?);
    let mut exons = Vec::new();
    let mut runs = 0usize;
    input::records(&args.exons, |header, seq| {
        let name = std::str::from_utf8(header)?
            .split_ascii_whitespace()
            .next()
            .context("empty exon header")?
            .to_string();
        let id = exons.len() + 1;
        exons.push(Exon {
            name,
            length: seq.len(),
            unique_kmers: 0,
        });
        // Explicit splitting ensures GGCAT cannot bridge an ambiguous base.
        for run in seq
            .split(|b| !input::is_dna(*b))
            .filter(|s| s.len() >= args.k)
        {
            writeln!(out, ">{id}")?;
            out.write_all(&run.to_ascii_uppercase())?;
            out.write_all(b"\n")?;
            runs += 1;
        }
        Ok(())
    })?;
    out.flush()?;
    drop(out);
    ensure!(!exons.is_empty(), "reference contains no exon records");
    ensure!(runs > 0, "reference contains no valid {}-mers", args.k);
    ensure!(
        exons.len() < u32::MAX as usize,
        "too many exons for 32-bit IDs"
    );
    eprintln!(
        "Building simplitigs for {} exons (k={})",
        exons.len(),
        args.k
    );
    let simplitigs = stage.path().join("simplitigs.fa");
    // Link GGCAT into EXPRESSO so building an index needs no external executable.
    // Input has already been decompressed, parsed, and normalized by Helicase.
    let ggcat = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(scratch.path().join("ggcat")),
        memory: args.memory_gb as f64,
        prefer_memory: false,
        total_threads_count: args.threads,
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: Some(|_, message| eprintln!("GGCAT: {message}")),
    })
    .context("initializing bundled GGCAT")?;
    ggcat
        .build_graph(
            vec![GeneralSequenceBlockData::FASTA((normalized.clone(), None))],
            simplitigs.clone(),
            None,
            args.k,
            args.threads,
            false,
            None,
            false,
            1,
            ExtraElaboration::FastSimplitigs,
            None,
        )
        .context("building GGCAT simplitigs")?;
    let mut sequences = Vec::new();
    input::records(&simplitigs, |_, seq| {
        ensure!(
            seq.len() >= args.k && seq.iter().all(|b| input::is_dna(*b)),
            "invalid GGCAT simplitig"
        );
        sequences.push(String::from_utf8(seq.to_ascii_uppercase())?);
        Ok(())
    })?;
    ensure!(!sequences.is_empty(), "GGCAT produced no simplitigs");
    config.num_threads = args.threads;
    config.ram_limit_gib = args.memory_gb;
    config.tmp_dirname = scratch.path().join("sshash");
    fs::create_dir_all(&config.tmp_dirname)?;
    eprintln!(
        "Constructing Rust SSHash from {} simplitigs",
        sequences.len()
    );
    let dict = DictionaryBuilder::new(config)
        .map_err(anyhow::Error::msg)?
        .build_from_sequences(sequences)
        .map_err(anyhow::Error::msg)?;
    let count = num_kmers(&dict)?;
    let owners = AtomicOwners::new(count, exons.len())?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    eprintln!("Assigning {count} k-mers to exons");
    sshash_lib::dispatch_on_k!(args.k, K => assign::<K>(&normalized, &dict, &owners, &pool))?;
    let owners = owners.finish();
    let (bits, unassigned, duplicate) = match &owners {
        Owners::U16(v) => (
            16,
            v.iter().filter(|x| **x == 0).count(),
            v.iter().filter(|x| **x == u16::MAX).count(),
        ),
        Owners::U32(v) => (
            32,
            v.iter().filter(|x| **x == 0).count(),
            v.iter().filter(|x| **x == u32::MAX).count(),
        ),
    };
    ensure!(
        unassigned == 0,
        "{unassigned} simplitig k-mers have no exon; inconsistent GGCAT/SSHash index"
    );
    for i in 0..count {
        if let Some(e) = owners.exon(i) {
            exons[e].unique_kmers += 1;
        }
    }
    dict.save(stage.path().join("dictionary"))
        .map_err(|e| anyhow::anyhow!("saving SSHash: {e}"))?;
    owners.write(&stage.path().join("owners.bin"))?;
    let meta = Metadata {
        format_version: FORMAT_VERSION,
        expresso_version: env!("CARGO_PKG_VERSION").to_string(),
        sshash_version: "0.7.1".into(),
        k: args.k,
        minimizer: m,
        owner_bits: bits,
        indexed_kmers: count,
        duplicate_kmers: duplicate,
        exons,
    };
    let mut out = BufWriter::new(File::create(stage.path().join("metadata.json"))?);
    serde_json::to_writer_pretty(&mut out, &meta)?;
    out.flush()?;
    drop(scratch);
    publish(stage, &args.index)?;
    eprintln!(
        "Index saved: {} ({bits}-bit exon IDs, {duplicate} shared k-mers excluded)",
        args.index.display()
    );
    Ok(())
}

fn assign<const K: usize>(
    path: &Path,
    dict: &Dictionary,
    owners: &AtomicOwners,
    pool: &rayon::ThreadPool,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    input::batches_with(
        path,
        1 << 20,
        2 * pool.current_num_threads(),
        K - 1,
        |header, _| Ok(std::str::from_utf8(header)?.parse()?),
        |rx| {
            pool.install(|| {
                rx.into_iter().par_bridge().try_for_each_init(
                    || dict.create_streaming_query::<K>(),
                    |query, batch| -> Result<()> {
                        for (record, seq) in batch.records() {
                            query.reset();
                            for window in seq.windows(K) {
                                let hit = query.lookup(window);
                                ensure!(
                                    hit.is_found(),
                                    "reference k-mer absent from GGCAT/SSHash index"
                                );
                                owners.assign(hit.kmer_id as usize, record.value);
                            }
                        }
                        Ok(())
                    },
                )
            })
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn compact_ids_and_shared_ownership() {
        for n in [65_534, 65_535] {
            let table = AtomicOwners::new(3, n).unwrap();
            table.assign(0, 1);
            table.assign(0, 1);
            table.assign(1, 1);
            table.assign(1, 2);
            table.assign(1, 1);
            table.assign(2, n as u64);
            let table = table.finish();
            assert_eq!(table.exon(0), Some(0));
            assert_eq!(table.exon(1), None);
            assert_eq!(table.exon(2), Some(n - 1));
            assert_eq!(matches!(table, Owners::U16(_)), n == 65_534);
        }
    }
}
