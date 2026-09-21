#!/usr/bin/env python3
"""Extract strand-oriented, unwrapped exon and genomic-gene FASTAs from GENCODE.

Exons are unique (contig, start, end, strand) intervals, not unique sequences.
All exon/gene identifiers remain in companion TSVs. Genes include introns.
Optionally verify every annotated exon occurrence against a transcript FASTA.
Uses only Python's standard library; the largest genome contig is held in RAM.
"""

import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from dataclasses import dataclass, field
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time


COMPLEMENT = bytes.maketrans(b"ACGTRYMKBDHVN", b"TGCAYRKMVHDBN")
ATTRIBUTES = re.compile(r'(?:^|;\s*)(gene_id|gene_name|gene_type|exon_id|transcript_id|exon_number)\s+(?:"([^"]*)"|([^;\s]+))')


def attributes(text):
    return {key: quoted or plain for key, quoted, plain in ATTRIBUTES.findall(text)}


def open_input(path, mode):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def fasta(path):
    with open_input(path, "rb") as handle:
        name, pieces = None, []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(pieces).upper()
                name, pieces = line[1:].split()[0].decode("ascii"), []
            else:
                if name is None:
                    raise ValueError("FASTA sequence before first header")
                pieces.append(line)
        if name is not None:
            yield name, b"".join(pieces).upper()


@dataclass(slots=True)
class Exon:
    locus: tuple
    ids: set = field(default_factory=set)
    genes: set = field(default_factory=set)
    occurrences: int = 0
    digest: bytes = b""


def annotation(path):
    genes, exon_lookup, exons, transcripts = {}, {}, [], defaultdict(list)
    transcript_ids = set()
    with open_input(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GTF line {number}")
            chrom, _, kind, start, end, _, strand, _, text = fields
            if kind not in ("gene", "transcript", "exon"):
                continue
            start, end = int(start), int(end)
            if start < 1 or end < start or strand not in ("+", "-"):
                raise ValueError(f"Invalid coordinates/strand at line {number}")
            locus = (sys.intern(chrom), start, end, sys.intern(strand))
            attrs = attributes(text)
            gene = attrs["gene_id"]
            if kind == "gene":
                if gene in genes:
                    raise ValueError(f"Repeated gene identifier: {gene}")
                genes[gene] = (locus, attrs["gene_name"], attrs["gene_type"])
            elif kind == "transcript":
                tx = attrs["transcript_id"]
                if tx in transcript_ids:
                    raise ValueError(f"Repeated transcript identifier: {tx}")
                transcript_ids.add(tx)
            else:
                index = exon_lookup.get(locus)
                if index is None:
                    index = len(exons)
                    exon_lookup[locus] = index
                    exons.append(Exon(locus))
                exon = exons[index]
                exon.ids.add(attrs["exon_id"])
                exon.genes.add(gene)
                exon.occurrences += 1
                transcripts[attrs["transcript_id"]].append((int(attrs["exon_number"]), index))
    if not genes or not exons or set(transcripts) != transcript_ids:
        raise ValueError("Incomplete gene/transcript/exon annotation")
    for exon in exons:
        chrom, start, end, strand = exon.locus
        for gene in exon.genes:
            gc, gs, ge, gd = genes[gene][0]
            if chrom != gc or strand != gd or not (gs <= start <= end <= ge):
                raise ValueError(f"Exon lies outside gene {gene}")
    for tx, parts in transcripts.items():
        parts.sort()
        if [n for n, _ in parts] != list(range(1, len(parts) + 1)):
            raise ValueError(f"Invalid exon numbering for {tx}")
    return genes, exons, transcripts


def sequence(genome, locus):
    _, start, end, strand = locus
    if end > len(genome):
        raise ValueError(f"Feature outside genome contig: {locus}")
    seq = genome[start - 1:end]
    if seq.translate(None, b"ACGTRYSWKMBDHVN"):
        raise ValueError(f"Invalid DNA alphabet at {locus}")
    return seq if strand == "+" else seq.translate(COMPLEMENT)[::-1]


def locus_text(locus):
    chrom, start, end, strand = locus
    return f"{chrom}:{start}-{end}:{strand}"


def extract(genome_path, stage, prefix, genes, exons):
    by_chrom = defaultdict(lambda: [[], []])
    for index, exon in enumerate(exons):
        by_chrom[exon.locus[0]][0].append(index)
    for gene, (locus, _, _) in genes.items():
        by_chrom[locus[0]][1].append(gene)
    totals = {kind: dict(records=0, bases=0, ambiguous_bases=0) for kind in ("exons", "genes")}
    seen = set()
    with ExitStack() as stack:
        files = {kind: stack.enter_context((stage / f"{prefix}.{kind}.fa").open("wb", buffering=4 << 20))
                 for kind in totals}
        tables = {kind: stack.enter_context(gzip.open(stage / f"{prefix}.{kind}.tsv.gz", "wt"))
                  for kind in totals}
        tables["exons"].write("fasta_id\texon_ids\tgene_ids\tcontig\tstart_1based\tend_1based\tstrand\tlength\ttranscript_occurrences\n")
        tables["genes"].write("fasta_id\tgene_id\tgene_name\tgene_type\tcontig\tstart_1based\tend_1based\tstrand\tlength\n")
        for chrom, genome in fasta(genome_path):
            if chrom in seen:
                raise ValueError(f"Duplicate genome contig: {chrom}")
            seen.add(chrom)
            if chrom not in by_chrom:
                continue
            exon_indexes, gene_ids = by_chrom.pop(chrom)
            for kind, entries in (("exons", exon_indexes), ("genes", gene_ids)):
                for entry in entries:
                    if kind == "exons":
                        exon = exons[entry]
                        locus = exon.locus
                        ids = sorted(exon.ids)
                        feature_id = ids[0] + "|" + locus_text(locus)
                        seq = sequence(genome, locus)
                        exon.digest = hashlib.sha256(seq).digest()
                        fields = [feature_id, ",".join(ids), ",".join(sorted(exon.genes)),
                                  *locus, len(seq), exon.occurrences]
                        description = f"gene_ids={fields[2]} length={len(seq)}"
                    else:
                        locus, name, biotype = genes[entry]
                        feature_id = entry + "|" + locus_text(locus)
                        seq = sequence(genome, locus)
                        fields = [feature_id, entry, name, biotype, *locus, len(seq)]
                        description = f"gene_name={name} gene_type={biotype} length={len(seq)} includes_introns=true"
                    files[kind].write(f">{feature_id} {description}\n".encode())
                    files[kind].write(seq)
                    files[kind].write(b"\n")
                    tables[kind].write("\t".join(map(str, fields)) + "\n")
                    totals[kind]["records"] += 1
                    totals[kind]["bases"] += len(seq)
                    totals[kind]["ambiguous_bases"] += len(seq.translate(None, b"ACGT"))
            if chrom in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}:
                print(f"Extracted {chrom}: {len(exon_indexes)} exons, {len(gene_ids)} genes", flush=True)
    if by_chrom:
        raise ValueError(f"Annotated contigs missing from genome: {sorted(by_chrom)}")
    if totals["exons"]["records"] != len(exons) or totals["genes"]["records"] != len(genes):
        raise ValueError("Feature count mismatch")
    return totals


def validate_transcripts(path, transcripts, exons):
    seen, occurrences, bases = set(), 0, 0
    for header, seq in fasta(path):
        tx = header.split("|")[0]
        if tx in seen or tx not in transcripts:
            raise ValueError(f"Unexpected or repeated transcript: {tx}")
        seen.add(tx)
        offset = 0
        for number, index in transcripts[tx]:
            exon = exons[index]
            length = exon.locus[2] - exon.locus[1] + 1
            fragment = seq[offset:offset + length]
            if len(fragment) != length or hashlib.sha256(fragment).digest() != exon.digest:
                raise ValueError(f"Transcript mismatch: {tx}, exon {number}, {exon.locus}")
            offset += length
            occurrences += 1
        if offset != len(seq):
            raise ValueError(f"Transcript length mismatch: {tx}")
        bases += len(seq)
    if seen != set(transcripts):
        raise ValueError(f"Missing {len(set(transcripts) - seen)} reference transcripts")
    return dict(transcripts=len(seen), exon_occurrences=occurrences, bases=bases,
                all_sequences_match=True)


def checksum(path):
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 << 20):
            sha.update(block)
    return dict(bytes=path.stat().st_size, sha256=sha.hexdigest())


def run(args):
    start = time.monotonic()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".features-", dir=output.parent))
    try:
        genes, exons, transcripts = annotation(args.gtf)
        print(f"Annotation: {len(genes)} genes, {len(exons)} unique exon intervals, {len(transcripts)} transcripts", flush=True)
        totals = extract(args.genome, stage, args.prefix, genes, exons)
        validation = None
        if args.transcripts:
            validation = validate_transcripts(args.transcripts, transcripts, exons)
            print("Transcript validation: " + json.dumps(validation), flush=True)
        report = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      gtf=str(args.gtf.resolve()), genome=str(args.genome.resolve()),
                      transcript_reference=str(args.transcripts.resolve()) if args.transcripts else None,
                      exon_definition="Unique annotated (contig, start, end, strand) interval; repeated transcript occurrences collapsed; identical sequences at different loci retained",
                      gene_definition="Full genomic gene interval including introns; strand-oriented 5-prime to 3-prime",
                      fasta_format="Two lines per record: header then one unwrapped sequence line",
                      coordinates="1-based inclusive", features=totals,
                      annotated_exon_occurrences=sum(x.occurrences for x in exons),
                      unique_exon_ids=len({id for x in exons for id in x.ids}),
                      gene_biotypes=dict(sorted(Counter(v[2] for v in genes.values()).items())),
                      transcript_validation=validation,
                      files={p.name: checksum(p) for p in sorted(stage.iterdir())})
        report["elapsed_seconds"] = time.monotonic() - start
        (stage / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        (stage / "SHA256SUMS").write_text("".join(f'{v["sha256"]}  {k}\n' for k, v in report["files"].items()))
        stage.rename(output)
        print(json.dumps(report, indent=2), flush=True)
    except BaseException:
        shutil.rmtree(stage)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--transcripts", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prefix", default="gencode.v49")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.prefix):
        parser.error("prefix must be a simple filename component")
    run(args)


if __name__ == "__main__":
    main()
