import argparse
import gzip
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "gencode_features", Path(__file__).resolve().parents[1] / "scripts/gencode_features.py")
features = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = features
spec.loader.exec_module(features)


class FeatureExtractionTest(unittest.TestCase):
    def test_strands_deduplication_and_transcript_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genome = root / "genome.fa"
            genome.write_text(">chrTest\nAACCGGTT\nACGATTCG\n")
            rows = []

            def row(kind, start, end, strand, attributes):
                rows.append(f"chrTest\tTEST\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")

            for gene, start, end, strand in [("G1", 1, 8, "+"), ("G2", 9, 16, "-"), ("G3", 1, 4, "+")]:
                row("gene", start, end, strand, f'gene_id "{gene}"; gene_name "{gene}"; gene_type "test";')
            for gene, tx, strand, parts in [
                ("G1", "T1", "+", [(1, 4, "E1"), (7, 8, "E2")]),
                ("G1", "T2", "+", [(1, 4, "E1")]),
                ("G2", "T3", "-", [(13, 16, "E3"), (9, 10, "E4")]),
                ("G3", "T4", "+", [(1, 4, "E1bis")]),
            ]:
                attrs = f'gene_id "{gene}"; transcript_id "{tx}";'
                row("transcript", min(x[0] for x in parts), max(x[1] for x in parts), strand, attrs)
                for number, (start, end, exon) in enumerate(parts, 1):
                    row("exon", start, end, strand, attrs + f' exon_number {number}; exon_id "{exon}";')
            gtf = root / "annotation.gtf.gz"
            with gzip.open(gtf, "wt") as handle:
                handle.writelines(rows)
            transcripts = root / "transcripts.fa"
            transcripts.write_text(">T1|G1\nAACCTT\n>T2|G1\nAACC\n>T3|G2\nCGAAGT\n>T4|G3\nAACC\n")
            output = root / "features"
            features.run(argparse.Namespace(gtf=gtf, genome=genome, transcripts=transcripts,
                                           output=output, prefix="test"))
            exons = dict(features.fasta(output / "test.exons.fa"))
            genes = dict(features.fasta(output / "test.genes.fa"))
            self.assertEqual(exons, {
                "E1|chrTest:1-4:+": b"AACC", "E2|chrTest:7-8:+": b"TT",
                "E3|chrTest:13-16:-": b"CGAA", "E4|chrTest:9-10:-": b"GT"})
            self.assertEqual(genes["G2|chrTest:9-16:-"], b"CGAATCGT")
            self.assertEqual(genes["G1|chrTest:1-8:+"], b"AACCGGTT")
            report = json.loads((output / "summary.json").read_text())
            self.assertEqual(report["annotated_exon_occurrences"], 6)
            self.assertEqual(report["unique_exon_ids"], 5)
            self.assertEqual(report["transcript_validation"]["transcripts"], 4)
            with gzip.open(output / "test.exons.tsv.gz", "rt") as handle:
                self.assertIn("E1,E1bis\tG1,G3", handle.read())
            lines = (output / "test.exons.fa").read_text().splitlines()
            self.assertEqual(len(lines), 8)
            self.assertTrue(all(lines[i].startswith(">") for i in range(0, 8, 2)))
            parsed_genes, parsed_exons, parsed_transcripts = features.annotation(gtf)
            features.extract(genome, root, "bad", parsed_genes, parsed_exons)
            transcripts.write_text(">T1|G1\nAACCAA\n")
            with self.assertRaisesRegex(ValueError, "Transcript mismatch"):
                features.validate_transcripts(transcripts, parsed_transcripts, parsed_exons)

    def test_iupac_reverse_complement_and_bounds(self):
        self.assertEqual(features.sequence(b"ACGTRYMKBDHVNSW", ("chr", 1, 15, "-")), b"WSNBDHVMKRYACGT")
        with self.assertRaisesRegex(ValueError, "outside genome"):
            features.sequence(b"ACG", ("chr", 1, 4, "+"))


if __name__ == "__main__":
    unittest.main()
