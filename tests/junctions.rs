use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
    process::Command,
};

fn invoke(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_expresso"))
        .current_dir(dir)
        .args(args)
        .env("PATH", "")
        .output()
        .unwrap()
}

fn success(dir: &Path, args: &[&str]) {
    let result = invoke(dir, args);
    assert!(
        result.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
}

fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            _ => panic!(),
        })
        .collect()
}

fn kmers(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    seq.windows(7)
        .filter(|w| w.iter().all(|b| b"ACGT".contains(b)))
        .map(|w| w.to_vec().min(rc(w)))
}

type Interval = (usize, usize);
struct Gene {
    name: &'static str,
    minus: bool,
    paths: Vec<Vec<Interval>>,
}

fn transcript(genome: &[u8], path: &[Interval], minus: bool) -> Vec<u8> {
    let seq: Vec<u8> = path
        .iter()
        .flat_map(|&(start, end)| genome[start..end].iter().copied())
        .collect();
    if minus { rc(&seq) } else { seq }
}

#[test]
fn gene_junction_indexes_match_independent_transcript_and_pair_oracles() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let mut seed = 193751u64;
    let mut genome: Vec<u8> = (0..180)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed & 3) as usize]
        })
        .collect();
    genome[97] = b'N';
    let genes = [
        Gene {
            name: "g1",
            minus: false,
            paths: vec![
                vec![(0, 12), (20, 22), (35, 51)],
                vec![(0, 14), (35, 51)],
                vec![(10, 18)], // Overlap chain extends the merged exon to 18.
            ],
        },
        Gene {
            name: "g2",
            minus: true,
            paths: vec![
                vec![(90, 106), (115, 118), (135, 153)],
                vec![(90, 106), (135, 153)],
                vec![(104, 110)], // Minus-strand merged exon extends to 110.
            ],
        },
        Gene {
            name: "g3",
            minus: false,
            paths: vec![vec![(0, 12), (20, 22), (35, 51)]],
        },
    ];
    let mut rows = Vec::new();
    let mut reads = Vec::new();
    for gene in &genes {
        for (t, path) in gene.paths.iter().enumerate() {
            for &(start, end) in path {
                rows.push(format!("chr1\ttest\texon\t{}\t{end}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}t{t}\";\n", start+1, if gene.minus { "-" } else { "+" }, gene.name, gene.name));
            }
            reads.push(transcript(&genome, path, gene.minus));
        }
        let exons: BTreeSet<_> = gene.paths.iter().flatten().copied().collect();
        // Query every forward/reverse pair as well as the annotated microexon transcripts.
        for &left in &exons {
            for &right in &exons {
                if left.1 <= right.0 || right.1 <= left.0 {
                    let mut pair = transcript(&genome, &[left], gene.minus);
                    pair.extend(transcript(&genome, &[right], gene.minus));
                    reads.push(pair);
                }
            }
        }
    }
    rows.push(rows[0].clone()); // Duplicate annotations must not change ownership.
    rows.reverse(); // Do not depend on GTF line order or exon_number.
    fs::write(dir.join("a.gtf"), rows.concat()).unwrap();
    fs::write(
        dir.join("genome.fa"),
        format!(">chr1\n{}\n", String::from_utf8_lossy(&genome)),
    )
    .unwrap();
    let mut fasta = String::new();
    for (i, read) in reads.iter().enumerate() {
        fasta.push_str(&format!(
            ">r{i} ka:f:3\n{}\n",
            String::from_utf8_lossy(read)
        ));
    }
    fs::write(dir.join("reads.fa"), fasta).unwrap();
    fs::write(dir.join("files.tsv"), "sample\treads.fa\tunitigs\n").unwrap();
    let mut total_counts = Vec::new();
    for (mode, limit) in [
        ("none", "50"),
        ("annotated", "50"),
        ("all", "50"),
        ("all", "3"), // Original interval counts exceed 3; merged counts do not.
        ("all", "1"),
    ] {
        let index = format!("idx-{mode}-{limit}");
        let results = format!("results-{mode}-{limit}");
        success(
            dir,
            &[
                "run",
                "--gtf",
                "a.gtf",
                "--genome",
                "genome.fa",
                "--level",
                "gene",
                "--junctions",
                mode,
                "--junction-max-exons",
                limit,
                "-i",
                &index,
                "-k",
                "7",
                "-t",
                "2",
                "-f",
                "files.tsv",
                "-o",
                &results,
                "--format",
                "csv",
                "--compression",
                "none",
            ],
        );
        let mut owners: BTreeMap<Vec<u8>, BTreeSet<&str>> = BTreeMap::new();
        let mut lengths = BTreeMap::new();
        for gene in &genes {
            let exons: BTreeSet<_> = gene.paths.iter().flatten().copied().collect();
            let mut merged: Vec<Interval> = Vec::new();
            for &(start, end) in &exons {
                if let Some(last) = merged.last_mut()
                    && start < last.1
                {
                    last.1 = last.1.max(end);
                } else {
                    merged.push((start, end));
                }
            }
            lengths.insert(
                gene.name,
                merged
                    .iter()
                    .map(|&(start, end)| end - start)
                    .sum::<usize>(),
            );
            let mut add = |seq: Vec<u8>| {
                for kmer in kmers(&seq) {
                    owners.entry(kmer).or_default().insert(gene.name);
                }
            };
            for &interval in &merged {
                add(transcript(&genome, &[interval], gene.minus));
            }
            if mode != "none" {
                for path in &gene.paths {
                    add(transcript(&genome, path, gene.minus));
                }
            }
            if mode == "all" && merged.len() <= limit.parse().unwrap() {
                for &left in &merged {
                    for &right in &merged {
                        if left.1 <= right.0 || right.1 <= left.0 {
                            let left = transcript(&genome, &[left], gene.minus);
                            let right = transcript(&genome, &[right], gene.minus);
                            let mut context = left[left.len().saturating_sub(6)..].to_vec();
                            context.extend(&right[..right.len().min(6)]);
                            add(context);
                        }
                    }
                }
            }
        }
        let meta: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join(&index).join("metadata.json")).unwrap())
                .unwrap();
        assert_eq!(
            meta["indexed_kmers"].as_u64().unwrap() as usize,
            owners.len()
        );
        assert_eq!(
            meta["duplicate_kmers"].as_u64().unwrap() as usize,
            owners.values().filter(|v| v.len() > 1).count()
        );
        for target in meta["exons"].as_array().unwrap() {
            let name = target["name"].as_str().unwrap();
            assert_eq!(
                target["length"].as_u64().unwrap() as usize,
                lengths[name],
                "junctions must not inflate exonic length"
            );
            assert_eq!(
                target["unique_kmers"].as_u64().unwrap() as usize,
                owners
                    .values()
                    .filter(|v| v.len() == 1 && v.contains(name))
                    .count()
            );
        }
        let mut expected = BTreeMap::from([
            ("g1".to_string(), 0u64),
            ("g2".to_string(), 0),
            ("g3".to_string(), 0),
        ]);
        for read in &reads {
            for kmer in kmers(read) {
                if let Some(genes) = owners.get(&kmer)
                    && genes.len() == 1
                {
                    *expected.get_mut(*genes.first().unwrap()).unwrap() += 3;
                }
            }
        }
        let actual: BTreeMap<String, u64> =
            csv::Reader::from_path(dir.join(&results).join("datasets/sample.csv"))
                .unwrap()
                .deserialize::<(usize, String, u64)>()
                .map(|row| {
                    let (_, name, n) = row.unwrap();
                    (name, n)
                })
                .collect();
        assert_eq!(actual, expected);
        total_counts.push(owners.len());
        if mode != "none" {
            assert_eq!(meta["junctions"]["exon_counting"], "merged_overlaps");
        }
        if mode == "all" && limit == "3" {
            assert_eq!(meta["junctions"]["stats"]["expanded_genes"], 3);
            assert_eq!(meta["junctions"]["stats"]["genes_above_limit"], 0);
        }
        if limit == "1" {
            assert_eq!(meta["junctions"]["stats"]["genes_above_limit"], 3);
            assert_eq!(meta["junctions"]["stats"]["additional_pairs"], 0);
        }
    }
    assert!(total_counts[1] > total_counts[0]);
    assert!(total_counts[2] > total_counts[1]);
    assert_eq!(total_counts[3], total_counts[2]);
    assert_eq!(total_counts[4], total_counts[1]);
}

#[test]
fn invalid_junction_requests_never_publish_an_index() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    fs::write(dir.join("g.fa"), ">chr1\nACGTTGCATATCCGAG\n").unwrap();
    for attributes in [
        "gene_id \"g\";",
        "transcript_id \"t\";",
        "gene_id \"g\"; transcript_id \"t\"; transcript_id \"u\";",
    ] {
        fs::write(
            dir.join("a.gtf"),
            format!("chr1\ttest\texon\t1\t10\t.\t+\t.\t{attributes}\n"),
        )
        .unwrap();
        let out = invoke(
            dir,
            &[
                "build",
                "--level",
                "gene",
                "--gtf",
                "a.gtf",
                "--genome",
                "g.fa",
                "--junctions",
                "annotated",
                "-i",
                "bad",
                "-k",
                "7",
            ],
        );
        assert!(!out.status.success());
        assert!(!dir.join("bad").exists());
    }
    fs::write(
        dir.join("a.gtf"),
        "chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n",
    )
    .unwrap();
    for args in [
        vec![
            "build",
            "--gtf",
            "a.gtf",
            "--genome",
            "g.fa",
            "--junctions",
            "all",
            "-i",
            "bad",
        ],
        vec![
            "build",
            "--level",
            "gene",
            "--genes",
            "g.fa",
            "--junctions",
            "annotated",
            "-i",
            "bad",
        ],
        vec![
            "build",
            "--level",
            "gene",
            "--gtf",
            "a.gtf",
            "--genome",
            "g.fa",
            "--junctions",
            "all",
            "--junction-max-exons",
            "0",
            "-i",
            "bad",
        ],
    ] {
        assert!(!invoke(dir, &args).status.success());
        assert!(!dir.join("bad").exists());
    }
    for rows in [
        "chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\nchr1\ttest\texon\t5\t15\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n",
        "chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\nchr1\ttest\texon\t11\t15\t.\t-\t.\tgene_id \"g\"; transcript_id \"u\";\n",
    ] {
        fs::write(dir.join("a.gtf"), rows).unwrap();
        assert!(
            !invoke(
                dir,
                &[
                    "build",
                    "--level",
                    "gene",
                    "--gtf",
                    "a.gtf",
                    "--genome",
                    "g.fa",
                    "--junctions",
                    "annotated",
                    "-i",
                    "bad",
                    "-k",
                    "7"
                ]
            )
            .status
            .success()
        );
        assert!(!dir.join("bad").exists());
    }
}
