//! Gene-owned splice-junction contexts; every emitted k-mer crosses the boundary.
use crate::exons::{Locus, complement, header_value, merge_overlapping};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    io::Write,
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    #[default]
    None,
    Annotated,
    All,
}

#[derive(clap::Args)]
pub struct Options {
    /// Gene/GTF mode only: no junctions, annotated transcript junctions, or all ordered exon pairs.
    #[arg(long, value_enum, default_value = "none")]
    pub junctions: Mode,
    /// Expand all ordered pairs only for genes with at most this many merged exons.
    /// Larger genes retain annotated junctions. Both exon orders are included.
    #[arg(long, default_value_t = 50, value_parser = clap::value_parser!(u32).range(1..))]
    pub junction_max_exons: u32,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            junctions: Mode::None,
            junction_max_exons: 50,
        }
    }
}

impl Options {
    pub fn validate(&self, level: crate::Level, has_gtf: bool) -> Result<()> {
        ensure!(
            self.junctions == Mode::None || (level == crate::Level::Gene && has_gtf),
            "--junctions annotated/all requires --level gene, --gtf and --genome"
        );
        Ok(())
    }
}

#[derive(Default, Serialize, Deserialize)]
pub struct Stats {
    pub annotated_junctions: usize,
    pub additional_pairs: usize,
    pub sequence_records: usize,
    pub expanded_genes: usize,
    pub genes_above_limit: usize,
}

#[derive(Serialize, Deserialize)]
pub struct Metadata {
    pub mode: Mode,
    pub max_exons: u32,
    #[serde(default)]
    pub exon_counting: ExonCounting,
    pub stats: Stats,
}

#[derive(Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExonCounting {
    // Indices built before overlap merging had no exon_counting field.
    #[default]
    OriginalIntervals,
    MergedOverlaps,
}

pub struct Gene {
    contig: String,
    strand: u8,
    exons: BTreeSet<Locus>,
    transcripts: BTreeMap<String, BTreeSet<Locus>>,
}

#[derive(Default)]
pub struct Annotation {
    genes: BTreeMap<String, Gene>,
    transcript_genes: HashMap<String, String>,
}

impl Annotation {
    pub fn add(&mut self, contig: &str, gene: &str, transcript: &str, locus: Locus) -> Result<()> {
        let owner = self
            .transcript_genes
            .entry(transcript.to_string())
            .or_insert_with(|| gene.to_string());
        ensure!(
            owner == gene,
            "transcript {transcript} belongs to multiple genes"
        );
        let entry = self.genes.entry(gene.to_string()).or_insert_with(|| Gene {
            contig: contig.to_string(),
            strand: locus.strand,
            exons: BTreeSet::new(),
            transcripts: BTreeMap::new(),
        });
        ensure!(
            entry.contig == contig && entry.strand == locus.strand,
            "gene {gene} has inconsistent contig or strand"
        );
        entry.exons.insert(locus);
        entry
            .transcripts
            .entry(transcript.to_string())
            .or_default()
            .insert(locus);
        Ok(())
    }

    pub fn into_contigs(self) -> BTreeMap<String, BTreeMap<String, Gene>> {
        let mut contigs: BTreeMap<String, BTreeMap<String, Gene>> = BTreeMap::new();
        for (name, mut gene) in self.genes {
            // Use the same union as the gene reference for exhaustive pairs and
            // their limit. Original paths remain intact for annotated junctions.
            gene.exons = merge_overlapping(gene.exons.into_iter().collect())
                .into_iter()
                .collect();
            contigs
                .entry(gene.contig.clone())
                .or_default()
                .insert(name, gene);
        }
        contigs
    }
}

fn sequence(bases: &[u8], exon: Locus) -> Result<Vec<u8>> {
    ensure!(
        exon.end <= bases.len(),
        "junction exon exceeds contig length"
    );
    let seq = &bases[exon.start - 1..exon.end];
    if exon.strand == b'-' {
        seq.iter().rev().map(|&b| complement(b)).collect()
    } else {
        Ok(seq.to_ascii_uppercase())
    }
}

fn boundary(left: Locus, right: Locus) -> (usize, usize) {
    if left.strand == b'+' {
        (left.end, right.start)
    } else {
        (left.start, right.end)
    }
}

type Context = (usize, usize, Vec<u8>);

pub fn write_contig(
    genes: &BTreeMap<String, Gene>,
    bases: &[u8],
    k: usize,
    options: &Options,
    out: &mut dyn Write,
    stats: &mut Stats,
) -> Result<()> {
    for (name, gene) in genes {
        let mut contexts: HashSet<Context> = HashSet::new();
        let mut boundaries = BTreeSet::new();
        // Emit in transcript order, independently of GTF input line/exon_number order.
        // Full transcript context also covers microexons and k-mers crossing several junctions.
        for (transcript, exons) in &gene.transcripts {
            let mut path: Vec<_> = exons.iter().copied().collect();
            ensure!(
                path.windows(2).all(|p| p[0].end < p[1].start),
                "overlapping exon intervals within transcript {transcript}"
            );
            if gene.strand == b'-' {
                path.reverse();
            }
            let mut seq = Vec::new();
            let mut offsets = Vec::new();
            for &exon in &path {
                seq.extend(sequence(bases, exon)?);
                offsets.push(seq.len());
            }
            for (i, pair) in path.windows(2).enumerate() {
                let edge = boundary(pair[0], pair[1]);
                boundaries.insert(edge);
                let cut = offsets[i];
                let context = &seq[cut.saturating_sub(k - 1)..seq.len().min(cut + k - 1)];
                if context.len() >= k && contexts.insert((edge.0, edge.1, context.to_vec())) {
                    writeln!(
                        out,
                        ">{} kind=junction donor={} acceptor={} strand={} source=annotated",
                        header_value(name),
                        edge.0,
                        edge.1,
                        gene.strand as char
                    )?;
                    out.write_all(context)?;
                    writeln!(out)?;
                    stats.sequence_records += 1;
                }
            }
        }
        stats.annotated_junctions += boundaries.len();
        if options.junctions != Mode::All {
            continue;
        }
        if gene.exons.len() > options.junction_max_exons as usize {
            stats.genes_above_limit += 1;
            continue;
        }
        stats.expanded_genes += 1;
        // Cache only the two k-1 flanks per exon, rather than the full exonic sequences.
        let flanks: Vec<_> = gene
            .exons
            .iter()
            .map(|&exon| {
                let seq = sequence(bases, exon)?;
                let n = (k - 1).min(seq.len());
                Ok((exon, seq[..n].to_vec(), seq[seq.len() - n..].to_vec()))
            })
            .collect::<Result<_>>()?;
        for (left, _, tail) in &flanks {
            for (right, head, _) in &flanks {
                if !(left.end < right.start || right.end < left.start) {
                    continue;
                }
                let edge = boundary(*left, *right);
                if boundaries.insert(edge) {
                    stats.additional_pairs += 1;
                }
                let mut context = Vec::with_capacity(tail.len() + head.len());
                context.extend(tail);
                context.extend(head);
                if context.len() >= k && contexts.insert((edge.0, edge.1, context.clone())) {
                    writeln!(
                        out,
                        ">{} kind=junction donor={} acceptor={} strand={} source=all",
                        header_value(name),
                        edge.0,
                        edge.1,
                        gene.strand as char
                    )?;
                    out.write_all(&context)?;
                    writeln!(out)?;
                    stats.sequence_records += 1;
                }
            }
        }
    }
    Ok(())
}
