use crate::abundance::{self, Mode};
use anyhow::{Context, Result, bail, ensure};
use crossbeam_channel::{Receiver, Sender, bounded};
use helicase::{Config, FastxParser, HelicaseParser, ParserOptions, input::*};
use serde::Serialize;
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

const CONFIG: Config = ParserOptions::default().compute_quality().config();

pub fn records(path: &Path, mut visit: impl FnMut(&[u8], &[u8]) -> Result<()>) -> Result<()> {
    // Helicase currently panics on decoder errors. Convert them into a failed run
    // with the source path; never publish partial results as a successful dataset.
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| -> Result<()> {
        let mut file = BufReader::new(File::open(path)?);
        // SSHash's needletail and Helicase's optional xz feature require
        // incompatible liblzma versions. Decode xz with the shared 0.3 version,
        // then pass the stream to Helicase for all FASTA/FASTQ parsing.
        let reader: Box<dyn std::io::Read + Send> = if file.fill_buf()?.starts_with(b"\xfd7zXZ\0") {
            Box::new(xz2::read::XzDecoder::new_multi_decoder(file))
        } else {
            Box::new(file)
        };
        let mut input = ReaderInput::new(reader);
        if input.first_byte() == 0 && input.next().is_none() {
            return Ok(());
        }
        let mut parser = FastxParser::<CONFIG>::from_input(input)?;
        let mut cleaned = Vec::new();
        while parser.next().is_some() {
            let header = parser.get_header().trim_ascii();
            let seq = parser.get_dna_string();
            // Helicase's string path retains CR bytes (and adjoining line
            // breaks) on CRLF-wrapped FASTA. Remove formatting whitespace
            // before identifying ambiguous bases or computing KC denominators.
            if seq.iter().any(u8::is_ascii_whitespace) {
                cleaned.clear();
                cleaned.extend(seq.iter().copied().filter(|b| !b.is_ascii_whitespace()));
                validate_quality(parser.get_quality(), cleaned.len())?;
                visit(header, &cleaned)?;
            } else {
                validate_quality(parser.get_quality(), seq.len())?;
                visit(header, seq)?;
            }
        }
        Ok(())
    }))
    .map_err(|p| {
        let message = p
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| p.downcast_ref::<&str>().copied())
            .unwrap_or("parser panic");
        anyhow::anyhow!("{message}")
    })?
    .with_context(|| format!("reading {}", path.display()))
}

fn validate_quality(quality: Option<&[u8]>, length: usize) -> Result<()> {
    if let Some(quality) = quality {
        ensure!(
            quality.trim_ascii_end().len() == length,
            "FASTQ sequence/quality lengths differ"
        );
    }
    Ok(())
}

#[inline]
pub fn is_dna(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

pub struct Record {
    pub start: usize,
    pub end: usize,
    pub value: u64,
}

#[derive(Default)]
struct Buffers {
    bases: Vec<u8>,
    records: Vec<Record>,
    records_seen: u64,
}

/// A contiguous sequence arena avoids millions of tiny cross-thread allocations.
/// Completed arenas go back to their producer through a bounded recycling queue.
pub struct Batch {
    buffers: Buffers,
    recycle: Sender<Buffers>,
}

impl Batch {
    fn new(buffers: Buffers, recycle: Sender<Buffers>) -> Self {
        Self { buffers, recycle }
    }

    pub fn records(&self) -> impl Iterator<Item = (&Record, &[u8])> {
        self.buffers
            .records
            .iter()
            .map(|record| (record, &self.buffers.bases[record.start..record.end]))
    }

    pub fn num_records(&self) -> u64 {
        self.buffers.records_seen
    }
}

impl Drop for Batch {
    fn drop(&mut self) {
        let mut buffers = std::mem::take(&mut self.buffers);
        buffers.bases.clear();
        buffers.records.clear();
        buffers.records_seen = 0;
        // Never block a worker while the producer is waiting for query results.
        let _ = self.recycle.try_send(buffers);
    }
}

/// One parser/decoder feeds a bounded queue; the consumer runs query workers.
/// The producer retains at most one batch plus Helicase's largest-record buffer.
pub fn batches<T: Send>(
    path: &Path,
    mode: Mode,
    k: usize,
    query_k: usize,
    batch_bases: usize,
    queue_size: usize,
    consume: impl FnOnce(Receiver<Batch>) -> Result<T> + Send,
) -> Result<T> {
    batches_with(
        path,
        batch_bases,
        queue_size,
        query_k - 1,
        |header, seq| abundance::weight(header, seq.len(), k, mode),
        consume,
    )
}

pub fn batches_with<T: Send>(
    path: &Path,
    batch_bases: usize,
    queue_size: usize,
    overlap: usize,
    mut value: impl FnMut(&[u8], &[u8]) -> Result<u64>,
    consume: impl FnOnce(Receiver<Batch>) -> Result<T> + Send,
) -> Result<T> {
    std::thread::scope(|scope| {
        let (tx, rx) = bounded(queue_size);
        let (recycle_tx, recycle_rx) = bounded(queue_size);
        let consumer = scope.spawn(move || consume(rx));
        let producer = (|| -> Result<()> {
            let mut batch = Batch::new(Buffers::default(), recycle_tx.clone());
            let mut bases = 0;
            records(path, |header, seq| {
                let value = value(header, seq)?;
                let mut start: usize = 0;
                loop {
                    let end = start
                        .saturating_add(batch_bases)
                        .saturating_add(overlap)
                        .min(seq.len());
                    if start == 0 {
                        batch.buffers.records_seen += 1;
                    }
                    bases += (end - start).max(1);
                    // Still validate headers and count records shorter than k,
                    // but avoid copying DNA that cannot contribute any k-mer.
                    if end - start > overlap {
                        let offset = batch.buffers.bases.len();
                        batch.buffers.bases.extend_from_slice(&seq[start..end]);
                        batch.buffers.bases[offset..].make_ascii_uppercase();
                        batch.buffers.records.push(Record {
                            start: offset,
                            end: batch.buffers.bases.len(),
                            value,
                        });
                    }
                    if bases >= batch_bases {
                        let next = Batch::new(
                            recycle_rx.try_recv().unwrap_or_default(),
                            recycle_tx.clone(),
                        );
                        tx.send(std::mem::replace(&mut batch, next))
                            .context("query workers stopped")?;
                        bases = 0;
                    }
                    if end == seq.len() {
                        break;
                    }
                    // Exactly k-1 bases overlap, so each k-mer is queried once.
                    start = end - overlap;
                }
                Ok(())
            })?;
            if !batch.buffers.records.is_empty() || batch.num_records() != 0 {
                tx.send(batch).context("query workers stopped")?;
            }
            Ok(())
        })();
        drop(tx);
        let result = consumer
            .join()
            .map_err(|_| anyhow::anyhow!("query worker panicked"))?;
        // Prefer the actual query error (overflow etc.) to a broken queue error.
        let result = result?;
        producer?;
        Ok(result)
    })
}

#[derive(Clone, Serialize)]
pub struct Dataset {
    pub name: String,
    pub path: PathBuf,
    pub mode: Mode,
}

pub fn datasets(fof: &Path, default: Mode) -> Result<Vec<Dataset>> {
    let base = fof.parent().unwrap_or(Path::new("."));
    let mut result = Vec::new();
    let mut names = HashSet::new();
    for (line_no, line) in
        BufReader::new(File::open(fof).with_context(|| format!("opening {}", fof.display()))?)
            .lines()
            .enumerate()
    {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').map(str::trim).collect();
        let (name, path, mode) = match fields.as_slice() {
            [path] => (
                format!("{:06}_{}", result.len() + 1, basename(path)),
                *path,
                default,
            ),
            [name, path] => (name.to_string(), *path, default),
            [name, path, mode] => (
                name.to_string(),
                *path,
                match *mode {
                    "reads" => Mode::Reads,
                    "unitigs" => Mode::Unitigs,
                    "auto" => Mode::Auto,
                    _ => bail!("invalid mode on FOF line {}: {mode}", line_no + 1),
                },
            ),
            _ => bail!(
                "expected path or name<TAB>path[<TAB>mode] on FOF line {}",
                line_no + 1
            ),
        };
        ensure!(
            !name.is_empty()
                && name != "."
                && name != ".."
                && name
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b)),
            "dataset name must contain only ASCII letters, digits, dot, dash, underscore: {name}"
        );
        ensure!(names.insert(name.clone()), "duplicate dataset name: {name}");
        ensure!(
            !path.is_empty(),
            "empty dataset path on FOF line {}",
            line_no + 1
        );
        let path = base
            .join(path)
            .canonicalize()
            .with_context(|| format!("resolving FOF line {}: {path}", line_no + 1))?;
        ensure!(
            path.is_file(),
            "dataset is not a regular file: {}",
            path.display()
        );
        result.push(Dataset { name, path, mode });
    }
    ensure!(!result.is_empty(), "file-of-files contains no datasets");
    Ok(result)
}

fn basename(path: &str) -> String {
    let mut name = Path::new(path)
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    for suffix in [
        ".gz", ".xz", ".zstd", ".zst", ".fasta", ".fastq", ".fna", ".fa", ".fq",
    ] {
        if name.ends_with(suffix) {
            name.truncate(name.len() - suffix.len());
        }
    }
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || ".-_".contains(c) {
                c
            } else {
                '_'
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn wrapped_crlf_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("wrapped.fa");
        std::fs::write(&path, b">one\r\nACGT\r\nTGCA\r\n>two\r\nAAAA\r\n").unwrap();
        let mut found = Vec::new();
        records(&path, |h, s| {
            found.push((h.to_vec(), s.to_vec()));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            found,
            vec![
                (b"one".to_vec(), b"ACGTTGCA".to_vec()),
                (b"two".to_vec(), b"AAAA".to_vec())
            ]
        );
    }
}
