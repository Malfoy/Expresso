use super::{Exon, Model, Transcript};
use crate::{
    BuildArgs,
    exons::{self, Locus},
    input, output,
};
use anyhow::{Context, Result, ensure};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

struct Draft {
    gene: String,
    contig: String,
    strand: u8,
    exons: BTreeSet<Locus>,
}

pub(super) fn extract(args: &BuildArgs, root: &Path) -> Result<(Model, PathBuf, PathBuf)> {
    let mut drafts: BTreeMap<String, Draft> = BTreeMap::new();
    let mut loci: BTreeMap<(String, Locus), BTreeSet<String>> = BTreeMap::new();
    let gtf = args.gtf.as_ref().context("--isoforms requires --gtf")?;
    for (line_no, line) in BufReader::new(output::reader(gtf)?).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        (|| -> Result<()> {
            let fields: Vec<_> = line.split('\t').collect();
            ensure!(fields.len() == 9, "expected nine GTF columns");
            if fields[2] != "exon" {
                return Ok(());
            }
            ensure!(
                !fields[0].is_empty() && !fields[0].contains(char::is_whitespace),
                "invalid contig name"
            );
            let start = fields[3].parse()?;
            let end = fields[4].parse()?;
            ensure!(
                start > 0 && end >= start && matches!(fields[6], "+" | "-"),
                "invalid exon coordinates or strand"
            );
            let locus = Locus {
                start,
                end,
                strand: fields[6].as_bytes()[0],
            };
            let names = exons::attributes(fields[8])?;
            ensure!(
                names.genes.len() == 1 && names.transcripts.len() == 1,
                "isoform indexing requires exactly one gene_id and transcript_id on every exon"
            );
            let gene = names.genes.first().unwrap();
            let name = names.transcripts.first().unwrap();
            let draft = drafts.entry(name.clone()).or_insert_with(|| Draft {
                gene: gene.clone(),
                contig: fields[0].into(),
                strand: locus.strand,
                exons: BTreeSet::new(),
            });
            ensure!(
                draft.gene == *gene && draft.contig == fields[0] && draft.strand == locus.strand,
                "transcript {name} has inconsistent gene, contig or strand"
            );
            draft.exons.insert(locus);
            loci.entry((fields[0].into(), locus))
                .or_default()
                .extend(names.exons);
            Ok(())
        })()
        .with_context(|| format!("{}: GTF line {}", gtf.display(), line_no + 1))?;
    }
    ensure!(
        !drafts.is_empty(),
        "GTF contains no annotated transcripts with exons"
    );
    let mut model = Model {
        version: 1,
        k: args.k,
        indexed_kmers: 0,
        exons: Vec::new(),
        transcripts: Vec::new(),
        features: Vec::new(),
    };
    let mut exon_ids = BTreeMap::new();
    let mut exon_contigs: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for ((contig, locus), names) in loci {
        let id = model.exons.len();
        exon_ids.insert((contig.clone(), locus), id);
        exon_contigs.entry(contig.clone()).or_default().push(id);
        model.exons.push(Exon {
            name: format!("exon_{:06}", id + 1),
            contig,
            start: locus.start,
            end: locus.end,
            strand: locus.strand as char,
            annotation_ids: names.into_iter().collect(),
        });
    }
    let mut transcript_contigs: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    let mut gene_loci = BTreeMap::new();
    for (name, draft) in drafts {
        if let Some((contig, strand)) =
            gene_loci.insert(draft.gene.clone(), (draft.contig.clone(), draft.strand))
        {
            ensure!(
                contig == draft.contig && strand == draft.strand,
                "gene {} has inconsistent contig or strand",
                draft.gene
            );
        }
        let mut path: Vec<_> = draft.exons.into_iter().collect();
        ensure!(
            path.windows(2).all(|p| p[0].end < p[1].start),
            "overlapping exon intervals within transcript {name}"
        );
        if draft.strand == b'-' {
            path.reverse();
        }
        let exons = path
            .iter()
            .map(|locus| exon_ids[&(draft.contig.clone(), *locus)])
            .collect();
        let id = model.transcripts.len();
        transcript_contigs
            .entry(draft.contig.clone())
            .or_default()
            .push(id);
        model.transcripts.push(Transcript {
            name,
            gene: draft.gene,
            contig: draft.contig,
            strand: draft.strand as char,
            exons,
            length: 0,
            observable_kmers: 0,
        });
    }
    let exon_path = root.join("reference-exons.fa");
    let transcript_path = root.join("reference-transcripts.fa");
    let mut exon_out = BufWriter::new(File::create(&exon_path)?);
    let mut transcript_out = BufWriter::new(File::create(&transcript_path)?);
    let mut seen = BTreeSet::new();
    let mut sequence = Vec::new();
    for genome in &args.genome {
        exons::require_fasta(genome)?;
        input::records(genome, |header, bases| {
            let contig = std::str::from_utf8(header)?
                .split_ascii_whitespace()
                .next()
                .context("empty genome header")?;
            ensure!(
                seen.insert(contig.to_string()),
                "duplicate genome contig {contig}"
            );
            if let Some(exons) = exon_contigs.remove(contig) {
                for id in exons {
                    sequence.clear();
                    append(bases, &model.exons[id], &mut sequence)?;
                    writeln!(exon_out, ">{id}")?;
                    exon_out.write_all(&sequence)?;
                    writeln!(exon_out)?;
                }
            }
            if let Some(transcripts) = transcript_contigs.remove(contig) {
                for id in transcripts {
                    sequence.clear();
                    for &exon in &model.transcripts[id].exons {
                        append(bases, &model.exons[exon], &mut sequence)?;
                    }
                    model.transcripts[id].length = sequence.len();
                    writeln!(transcript_out, ">{id}")?;
                    transcript_out.write_all(&sequence)?;
                    writeln!(transcript_out)?;
                }
            }
            Ok(())
        })?;
    }
    ensure!(
        exon_contigs.is_empty() && transcript_contigs.is_empty(),
        "annotated contigs missing from reference genome"
    );
    exon_out.flush()?;
    transcript_out.flush()?;
    Ok((model, exon_path, transcript_path))
}
fn append(bases: &[u8], exon: &Exon, sequence: &mut Vec<u8>) -> Result<()> {
    ensure!(
        exon.start > 0 && exon.end <= bases.len(),
        "exon {}:{}-{} outside genome",
        exon.contig,
        exon.start,
        exon.end
    );
    let fragment = &bases[exon.start - 1..exon.end];
    if exon.strand == '-' {
        for &base in fragment.iter().rev() {
            sequence.push(exons::complement(base)?);
        }
    } else {
        for &base in fragment {
            exons::complement(base)?;
            sequence.push(base.to_ascii_uppercase());
        }
    }
    Ok(())
}
