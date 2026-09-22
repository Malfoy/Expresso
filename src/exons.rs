//! Strand-oriented exon extraction from GTF coordinates and genomic FASTA.
use crate::{input, output};
use anyhow::{Context, Result, bail, ensure};
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Locus {
    start: usize,
    end: usize,
    strand: u8,
}

#[derive(Default)]
struct Names {
    exons: BTreeSet<String>,
    genes: BTreeSet<String>,
}

type Annotation = BTreeMap<String, BTreeMap<Locus, Names>>;

// Parse quoted GTF values without treating semicolons inside quotes as separators.
// Repeated optional attributes (e.g. tag) are permitted. IDs are optional.
fn attributes(mut text: &str) -> Result<Names> {
    let mut names = Names::default();
    if text.trim() == "." {
        return Ok(names);
    }
    while !text.trim().is_empty() {
        text = text.trim_start();
        let end = text
            .find(char::is_whitespace)
            .context("GTF attribute missing value")?;
        let key = &text[..end];
        text = text[end..].trim_start();
        let value;
        if let Some(rest) = text.strip_prefix('"') {
            let mut parsed = String::new();
            let mut escaped = false;
            let mut closed = None;
            for (i, c) in rest.char_indices() {
                if escaped {
                    parsed.push(c);
                    escaped = false;
                } else if c == '\\' {
                    escaped = true;
                } else if c == '"' {
                    closed = Some(i + 1);
                    break;
                } else {
                    parsed.push(c);
                }
            }
            text = &rest[closed.context("unterminated quoted GTF attribute")?..];
            value = parsed;
        } else {
            let end = text.find(';').unwrap_or(text.len());
            value = text[..end].trim().to_string();
            ensure!(
                !value.contains(char::is_whitespace),
                "invalid unquoted GTF attribute"
            );
            text = &text[end..];
        }
        ensure!(!value.is_empty(), "empty GTF attribute {key}");
        text = text.trim_start();
        if !text.is_empty() {
            text = text
                .strip_prefix(';')
                .context("expected semicolon after GTF attribute")?;
        }
        match key {
            "exon_id" => {
                names.exons.insert(value);
            }
            "gene_id" => {
                names.genes.insert(value);
            }
            _ => {}
        }
    }
    Ok(names)
}

fn annotation(path: &Path) -> Result<Annotation> {
    let mut annotation: Annotation = BTreeMap::new();
    let mut reader = BufReader::new(output::reader(path)?);
    let mut line = String::new();
    let mut number = 0usize;
    let mut occurrences = 0usize;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        number += 1;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        (|| -> Result<()> {
            let fields: Vec<_> = line.split('\t').collect();
            ensure!(fields.len() == 9, "expected nine tab-separated GTF columns");
            if fields[2] != "exon" {
                return Ok(());
            }
            let contig = fields[0];
            ensure!(
                !contig.is_empty() && !contig.contains(char::is_whitespace),
                "invalid contig name"
            );
            let start: usize = fields[3].parse().context("invalid exon start")?;
            let end: usize = fields[4].parse().context("invalid exon end")?;
            ensure!(
                start > 0 && end >= start,
                "invalid 1-based inclusive exon coordinates"
            );
            ensure!(matches!(fields[6], "+" | "-"), "exon strand must be + or -");
            let locus = Locus {
                start,
                end,
                strand: fields[6].as_bytes()[0],
            };
            let names = attributes(fields[8])?;
            let entry = annotation
                .entry(contig.to_string())
                .or_default()
                .entry(locus)
                .or_default();
            entry.exons.extend(names.exons);
            entry.genes.extend(names.genes);
            occurrences += 1;
            Ok(())
        })()
        .with_context(|| format!("{}: GTF line {number}", path.display()))?;
    }
    ensure!(
        occurrences > 0,
        "GTF contains no exon features: {}",
        path.display()
    );
    let unique: usize = annotation.values().map(BTreeMap::len).sum();
    eprintln!("Extracting {unique} unique exon intervals from {occurrences} annotations");
    Ok(annotation)
}

fn complement(base: u8) -> Result<u8> {
    Ok(match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        _ => bail!("invalid genomic DNA byte 0x{base:02x}"),
    })
}

// Escape header separators, whitespace, and non-ASCII bytes unambiguously.
fn header_value(value: &str) -> String {
    let mut out = String::new();
    for b in value.bytes() {
        if b.is_ascii_alphanumeric() || b"._-".contains(&b) {
            out.push(b as char);
        } else {
            use std::fmt::Write;
            write!(out, "%{b:02X}").unwrap();
        }
    }
    out
}

fn require_fasta(path: &Path) -> Result<()> {
    let mut reader = BufReader::new(output::reader(path)?);
    let mut line = String::new();
    loop {
        line.clear();
        ensure!(
            reader.read_line(&mut line)? != 0,
            "empty genome FASTA: {}",
            path.display()
        );
        if !line.trim().is_empty() {
            ensure!(
                line.starts_with('>'),
                "expected genome FASTA: {}",
                path.display()
            );
            return Ok(());
        }
    }
}

pub fn extract(
    gtf: &Path,
    genomes: &[PathBuf],
    destination: &Path,
    compression: output::Compression,
) -> Result<()> {
    ensure!(
        !genomes.is_empty(),
        "provide at least one reference genome FASTA"
    );
    ensure!(
        !destination.exists(),
        "output already exists: {}",
        destination.display()
    );
    let mut annotation = annotation(gtf)?;
    let parent = destination
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    let stage = tempfile::NamedTempFile::new_in(parent)?;
    let mut seen = HashSet::new();
    let mut count = 0usize;
    let mut sequence = Vec::new();
    output::encoded(stage.path(), compression, |out| {
        for genome in genomes {
            require_fasta(genome)?;
            input::records(genome, |header, bases| {
                let contig = std::str::from_utf8(header)?
                    .split_ascii_whitespace()
                    .next()
                    .context("empty genome header")?;
                ensure!(
                    seen.insert(contig.to_string()),
                    "duplicate genome contig: {contig}"
                );
                let Some(exons) = annotation.remove(contig) else {
                    return Ok(());
                };
                for (locus, names) in exons {
                    ensure!(
                        locus.end <= bases.len(),
                        "exon {contig}:{}-{} exceeds contig length {}",
                        locus.start,
                        locus.end,
                        bases.len()
                    );
                    let fragment = &bases[locus.start - 1..locus.end];
                    sequence.clear();
                    for &base in fragment {
                        let paired = complement(base).with_context(|| {
                            format!("exon {contig}:{}-{}", locus.start, locus.end)
                        })?;
                        sequence.push(if locus.strand == b'-' {
                            paired
                        } else {
                            base.to_ascii_uppercase()
                        });
                    }
                    if locus.strand == b'-' {
                        sequence.reverse();
                    }
                    count += 1;
                    write!(
                        out,
                        ">exon_{count:06}|{}:{}-{}:{}",
                        header_value(contig),
                        locus.start,
                        locus.end,
                        locus.strand as char
                    )?;
                    for (label, ids) in [("exon_ids", names.exons), ("gene_ids", names.genes)] {
                        if !ids.is_empty() {
                            write!(
                                out,
                                " {label}={}",
                                ids.iter()
                                    .map(|id| header_value(id))
                                    .collect::<Vec<_>>()
                                    .join(",")
                            )?;
                        }
                    }
                    writeln!(out)?;
                    out.write_all(&sequence)?;
                    writeln!(out)?;
                }
                Ok(())
            })?;
        }
        ensure!(
            annotation.is_empty(),
            "annotated contigs missing from genome (names must match exactly): {}",
            annotation
                .keys()
                .take(10)
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        );
        Ok(())
    })?;
    stage
        .persist_noclobber(destination)
        .with_context(|| format!("publishing {}", destination.display()))?;
    eprintln!("Extracted {count} exons: {}", destination.display());
    Ok(())
}
