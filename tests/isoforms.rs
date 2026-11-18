//! Independent canonical-sequence oracle for exon/intersection membership,
//! transcript multiplicity, exact counts, and aggregate-only EM inference.
use serde_json::Value;
use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::{Read, Write},
    path::Path,
    process::{Command, Output},
};

fn invoke(root: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_expresso"))
        .current_dir(root)
        .args(args)
        .env("PATH", "")
        .output()
        .unwrap()
}
fn success(root: &Path, args: &[&str]) {
    let out = invoke(root, args);
    assert!(
        out.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}
fn canonical(seq: &[u8]) -> Vec<u8> {
    let rc = seq
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect::<Vec<_>>();
    seq.to_vec().min(rc)
}
fn kmers(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    seq.windows(7)
        .filter(|w| w.iter().all(|b| b"ACGT".contains(b)))
        .map(canonical)
}
fn fasta(path: &Path) -> Vec<(usize, Vec<u8>)> {
    let text = fs::read_to_string(path).unwrap();
    let lines: Vec<_> = text.lines().collect();
    lines
        .chunks(2)
        .map(|row| (row[0][1..].parse().unwrap(), row[1].as_bytes().to_vec()))
        .collect()
}
fn csv(path: &Path) -> Vec<BTreeMap<String, String>> {
    let file = fs::File::open(path).unwrap();
    let reader: Box<dyn Read> = if path.extension().unwrap() == "zst" {
        Box::new(zstd::stream::read::Decoder::new(file).unwrap())
    } else {
        Box::new(file)
    };
    let mut csv = csv::Reader::from_reader(reader);
    let header = csv.headers().unwrap().clone();
    csv.records()
        .map(|row| {
            header
                .iter()
                .zip(row.unwrap().iter())
                .map(|(k, v)| (k.to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}
type Members = Vec<(usize, u64)>;
type Signature = (Members, Members);
type ReferenceMembers = (BTreeMap<usize, u64>, BTreeMap<usize, u64>);
fn occurrences(value: &Value) -> Members {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            (
                v["member"].as_u64().unwrap() as usize,
                v["multiplicity"].as_u64().unwrap(),
            )
        })
        .collect()
}
fn fixture(root: &Path) -> Vec<(Vec<u8>, u64)> {
    let seqs =
        fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/viewer/genes.fa"))
            .unwrap();
    let genes: Vec<_> = seqs.lines().filter(|line| !line.starts_with('>')).collect();
    let a = genes[1];
    let b = genes[2];
    let c = genes[3];
    let reverse = genes[5];
    let mut genome = "N".repeat(400).into_bytes();
    for (start, seq) in [
        (0, a),
        (50, b),
        (100, c),
        (180, reverse),
        (250, "AAAAAAA"),
        (270, "AAAAAAA"),
        (290, "AAAAAAAAAA"),
        (350, "ACGTAC"),
    ] {
        genome[start..start + seq.len()].copy_from_slice(seq.as_bytes());
    }
    fs::write(
        root.join("genome.fa"),
        format!(">chr1\n{}\n", String::from_utf8(genome).unwrap()),
    )
    .unwrap();
    let paths = [
        ("g_skip", "T1", '+', vec![(1, 40), (51, 90), (101, 140)]),
        ("g_skip", "T2", '+', vec![(1, 40), (101, 140)]),
        ("g_skip", "T3", '+', vec![(1, 40), (101, 140)]),
        ("g_reverse", "R1", '-', vec![(181, 220)]),
        ("g_reverse", "R2", '-', vec![(181, 219)]),
        ("g_repeat", "Q1", '+', vec![(251, 257)]),
        ("g_repeat", "Q2", '+', vec![(271, 277)]),
        ("g_repeat", "Q3", '+', vec![(291, 300)]),
        ("g_repeat", "Q4", '+', vec![(251, 257), (271, 277)]),
        ("g_short", "Z1", '+', vec![(351, 356)]),
    ];
    let mut gtf = String::new();
    // Deliberately reverse GTF order and repeat a line: paths must be deduplicated and strand sorted.
    for (gene, transcript, strand, path) in paths.iter().rev() {
        for &(start, end) in path.iter().rev() {
            gtf.push_str(&format!("chr1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{transcript}\"; exon_id \"e_{start}_{end}\";\n"));
        }
    }
    gtf.push_str(gtf.lines().next().unwrap().to_string().as_str());
    gtf.push('\n');
    fs::write(root.join("reference.gtf"), gtf).unwrap();
    let records = vec![
        (format!("{a}{b}{c}").into_bytes(), 6_000_000),
        (format!("{a}{c}").into_bytes(), 4_000_000),
        (b"AAAAAAA".to_vec(), 1_250_000),
        (b"AAAAAAAAAAAAAA".to_vec(), 750_000),
        (b"CCCCCCCCCCCCCCC".to_vec(), 2_500_000),
    ];
    let mut data = String::new();
    for (id, (seq, weight)) in records.iter().enumerate() {
        data.push_str(&format!(
            ">u{id} ka:f:{}\n{}\n",
            *weight as f64 / 1e6,
            std::str::from_utf8(seq).unwrap()
        ));
    }
    fs::write(root.join("unitigs.fa"), &data).unwrap();
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    gz.write_all(data.as_bytes()).unwrap();
    fs::write(root.join("unitigs.fa.gz"), gz.finish().unwrap()).unwrap();
    fs::write(
        root.join("datasets.tsv"),
        "A\tunitigs.fa\tunitigs\nglobal\tunitigs.fa.gz\tunitigs\n",
    )
    .unwrap();
    records
}

#[test]
fn intersection_index_and_exact_counts_feed_em_without_sequences_or_dictionary() {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path();
    let records = fixture(root);
    success(
        root,
        &[
            "build",
            "--isoforms",
            "--gtf",
            "reference.gtf",
            "--genome",
            "genome.fa",
            "--index",
            "index",
            "--k",
            "7",
            "--threads",
            "2",
        ],
    );
    let model: Value =
        serde_json::from_slice(&fs::read(root.join("index/isoform-model.json")).unwrap()).unwrap();
    let metadata: Value =
        serde_json::from_slice(&fs::read(root.join("index/metadata.json")).unwrap()).unwrap();
    assert_eq!(metadata["duplicate_kmers"], 0);
    assert_eq!(metadata["owner_bits"], 16);
    assert_eq!(model["transcripts"].as_array().unwrap().len(), 10);
    let mut reference: BTreeMap<Vec<u8>, ReferenceMembers> = BTreeMap::new();
    for (exon, seq) in fasta(&root.join("index/reference-exons.fa")) {
        for kmer in kmers(&seq) {
            *reference
                .entry(kmer)
                .or_default()
                .0
                .entry(exon)
                .or_default() += 1;
        }
    }
    for (transcript, seq) in fasta(&root.join("index/reference-transcripts.fa")) {
        for kmer in kmers(&seq) {
            if let Some(members) = reference.get_mut(&kmer) {
                *members.1.entry(transcript).or_default() += 1;
            }
        }
    }
    let mut expected: HashMap<Signature, u64> = HashMap::new();
    for (exons, transcripts) in reference.values() {
        *expected
            .entry((
                exons.iter().map(|(&id, &n)| (id, n)).collect(),
                transcripts.iter().map(|(&id, &n)| (id, n)).collect(),
            ))
            .or_default() += 1;
    }
    let features = model["features"].as_array().unwrap();
    assert_eq!(features.len(), expected.len());
    assert_eq!(
        metadata["indexed_kmers"].as_u64().unwrap() as usize,
        reference.len()
    );
    let mut ids = HashMap::new();
    for (id, feature) in features.iter().enumerate() {
        let signature = (
            occurrences(&feature["exons"]),
            occurrences(&feature["transcripts"]),
        );
        assert_eq!(
            feature["distinct_kmers"].as_u64().unwrap(),
            expected[&signature]
        );
        ids.insert(signature, id);
    }
    assert!(
        features
            .iter()
            .any(|f| f["exons"].as_array().unwrap().len() >= 3),
        "retain triple intersections"
    );
    assert!(
        features
            .iter()
            .any(|f| occurrences(&f["exons"]).iter().any(|(_, n)| *n > 1)),
        "retain repeated positions within exons"
    );
    let mut counts = vec![0u64; features.len()];
    for (seq, weight) in &records {
        for kmer in kmers(seq) {
            if let Some((exons, transcripts)) = reference.get(&kmer) {
                let signature = (
                    exons.iter().map(|(&id, &n)| (id, n)).collect(),
                    transcripts.iter().map(|(&id, &n)| (id, n)).collect(),
                );
                counts[ids[&signature]] += weight;
            }
        }
    }
    for (codec, extension) in [
        ("none", "efc"),
        ("gz", "efc.gz"),
        ("xz", "efc.xz"),
        ("zstd", "efc.zst"),
    ] {
        let directory = format!("features-{codec}");
        success(
            root,
            &[
                "quantify",
                "--index",
                "index",
                "--fof",
                "datasets.tsv",
                "--output",
                &directory,
                "--threads",
                "4",
                "--jobs",
                "2",
                "--batch-bases",
                "16",
                "--compression",
                codec,
            ],
        );
        // Verify the native sparse EFC independent of EXPRESSO's reader.
        let path = root
            .join(&directory)
            .join(format!("datasets/A.{extension}"));
        let file = fs::File::open(path).unwrap();
        let mut reader: Box<dyn Read> = match codec {
            "gz" => Box::new(flate2::read::MultiGzDecoder::new(file)),
            "xz" => Box::new(xz2::read::XzDecoder::new_multi_decoder(file)),
            "zstd" => Box::new(zstd::stream::read::Decoder::new(file).unwrap()),
            _ => Box::new(file),
        };
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).unwrap();
        assert_eq!(&bytes[..8], b"EXPREF01");
        let mut actual = vec![0u128; counts.len()];
        for row in bytes[64..bytes.len() - 4].chunks_exact(24) {
            actual[u64::from_le_bytes(row[..8].try_into().unwrap()) as usize] =
                u128::from_le_bytes(row[8..].try_into().unwrap());
        }
        assert_eq!(
            actual,
            counts.iter().map(|&n| n as u128).collect::<Vec<_>>()
        );
        let destination = format!("inferred-{codec}");
        success(
            root,
            &[
                "infer-isoforms",
                "--index",
                "index",
                "--input",
                &directory,
                "--output",
                &destination,
                "--threads",
                "4",
                "--jobs",
                "2",
                "--compression",
                "none",
                "--tolerance",
                "0.0000000001",
            ],
        );
        let transcripts = csv(&root.join(&destination).join("datasets/A.csv"));
        let t1 = transcripts
            .iter()
            .find(|t| t["transcript_name"] == "T1")
            .unwrap();
        assert_eq!(t1["status"], "identifiable");
        assert!((t1["abundance"].parse::<f64>().unwrap() - 6.0).abs() < 1e-7);
        assert!((t1["within_gene_fraction"].parse::<f64>().unwrap() - 0.6).abs() < 1e-8);
        for name in ["T2", "T3"] {
            let t = transcripts
                .iter()
                .find(|t| t["transcript_name"] == name)
                .unwrap();
            assert_eq!(t["status"], "group_only");
            assert!(t["abundance"].is_empty());
            assert!(t["model_estimate"].is_empty());
        }
        for name in ["Q1", "Q2", "Q3", "Q4"] {
            let t = transcripts
                .iter()
                .find(|t| t["transcript_name"] == name)
                .unwrap();
            assert_eq!(t["status"], "unresolved");
            assert!(t["abundance"].is_empty());
        }
        assert_eq!(
            transcripts
                .iter()
                .find(|t| t["transcript_name"] == "Z1")
                .unwrap()["status"],
            "unobservable"
        );
        let groups = csv(&root.join(&destination).join("groups/datasets/A.csv"));
        let skip = groups
            .iter()
            .find(|g| g["transcript_names"].contains("T2"))
            .unwrap();
        assert!((skip["abundance"].parse::<f64>().unwrap() - 4.0).abs() < 1e-7);
        let sum: f64 = groups
            .iter()
            .map(|g| g["assigned_observations"].parse::<f64>().unwrap())
            .sum();
        assert!((sum - counts.iter().map(|&n| n as f64 / 1e6).sum::<f64>()).abs() < 1e-7);
        let global = csv(&root.join(&destination).join("global.csv"));
        assert!(
            (global
                .iter()
                .find(|t| t["transcript_name"] == "T1")
                .unwrap()["abundance"]
                .parse::<f64>()
                .unwrap()
                - 12.0)
                .abs()
                < 1e-7
        );
        assert!(
            root.join(&destination)
                .join("groups/datasets/global.csv")
                .exists(),
            "a dataset named global must not be overwritten by the sum"
        );
        assert!(root.join(&destination).join("groups/global.csv").exists());
    }
    let mut external = "feature_id,abundance\n".to_string();
    for (id, &n) in counts.iter().enumerate().rev() {
        external.push_str(&format!("{},{}\n", id + 1, n as f64 / 1e6));
    }
    fs::write(root.join("external.csv"), external).unwrap();
    fs::write(root.join("counts.tsv"), "manual\texternal.csv\n").unwrap();
    // Prove aggregate-only inference doesn't open FASTA, original input, GTF or dictionary.
    fs::remove_file(root.join("unitigs.fa")).unwrap();
    fs::remove_file(root.join("reference.gtf")).unwrap();
    fs::remove_file(root.join("genome.fa")).unwrap();
    fs::remove_file(root.join("index/dictionary.ssi")).unwrap();
    fs::remove_file(root.join("index/dictionary.ssi.mphf")).unwrap();
    fs::remove_file(root.join("index/owners.bin")).unwrap();
    fs::remove_file(root.join("index/reference-exons.fa")).unwrap();
    fs::remove_file(root.join("index/reference-transcripts.fa")).unwrap();
    success(
        root,
        &[
            "infer-isoforms",
            "--index",
            "index",
            "--fof",
            "counts.tsv",
            "--output",
            "manual",
            "--threads",
            "2",
            "--compression",
            "none",
        ],
    );
    let manual = csv(&root.join("manual/datasets/manual.csv"));
    assert!(
        (manual
            .iter()
            .find(|t| t["transcript_name"] == "T1")
            .unwrap()["abundance"]
            .parse::<f64>()
            .unwrap()
            - 6.0)
            .abs()
            < 1e-5
    );
    success(
        root,
        &[
            "infer-isoforms",
            "--index",
            "index",
            "--fof",
            "counts.tsv",
            "--output",
            "capped",
            "--rank-limit",
            "1",
            "--compression",
            "none",
        ],
    );
    assert_eq!(
        csv(&root.join("capped/datasets/manual.csv"))
            .iter()
            .find(|t| t["transcript_name"] == "T1")
            .unwrap()["status"],
        "rank_unchecked"
    );
    success(
        root,
        &[
            "infer-isoforms",
            "--index",
            "index",
            "--fof",
            "counts.tsv",
            "--output",
            "unfinished",
            "--max-iterations",
            "1",
            "--compression",
            "none",
        ],
    );
    assert_eq!(
        csv(&root.join("unfinished/datasets/manual.csv"))
            .iter()
            .find(|t| t["transcript_name"] == "T1")
            .unwrap()["status"],
        "non_converged"
    );
    fs::write(root.join("external.csv"), "feature_id,abundance\n1,NaN\n").unwrap();
    assert!(
        !invoke(
            root,
            &[
                "infer-isoforms",
                "--index",
                "index",
                "--fof",
                "counts.tsv",
                "--output",
                "invalid"
            ]
        )
        .status
        .success()
    );
    assert!(!root.join("invalid").exists());
    fs::write(root.join("external.csv"), "feature_id,abundance\n").unwrap();
    success(
        root,
        &[
            "infer-isoforms",
            "--index",
            "index",
            "--fof",
            "counts.tsv",
            "--output",
            "zeros",
            "--compression",
            "none",
        ],
    );
    assert_eq!(
        csv(&root.join("zeros/datasets/manual.csv"))
            .iter()
            .find(|t| t["transcript_name"] == "T1")
            .unwrap()["abundance"],
        "0.0"
    );
    // The model hash must protect the exposure matrix independently of count integrity.
    let mut raw = fs::read(root.join("index/isoform-model.json")).unwrap();
    raw.push(b' ');
    fs::write(root.join("index/isoform-model.json"), raw).unwrap();
    let error = invoke(
        root,
        &[
            "infer-isoforms",
            "--index",
            "index",
            "--fof",
            "counts.tsv",
            "--output",
            "bad-model",
        ],
    );
    assert!(!error.status.success());
    assert!(String::from_utf8_lossy(&error.stderr).contains("checksum"));
    assert!(!root.join("bad-model").exists());
}

#[test]
fn long_reference_batches_preserve_occurrences_across_chunk_boundaries() {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path();
    let mut sequence = vec![b'A'; (1 << 20) + 123];
    // An ambiguous base at the parser chunk boundary must split the DNA run.
    sequence[(1 << 20) - 1] = b'N';
    fs::write(
        root.join("genome.fa"),
        format!(">chr1\n{}\n", std::str::from_utf8(&sequence).unwrap()),
    )
    .unwrap();
    fs::write(
        root.join("annotation.gtf"),
        format!(
            "chr1\ttoy\texon\t1\t{}\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n",
            sequence.len()
        ),
    )
    .unwrap();
    success(
        root,
        &[
            "build",
            "--isoforms",
            "--gtf",
            "annotation.gtf",
            "--genome",
            "genome.fa",
            "--index",
            "index",
            "--k",
            "7",
            "--threads",
            "4",
            "--memory-gb",
            "1",
        ],
    );
    let model: Value =
        serde_json::from_slice(&fs::read(root.join("index/isoform-model.json")).unwrap()).unwrap();
    let expected = sequence
        .split(|&b| b == b'N')
        .map(|run| run.len().saturating_sub(6) as u64)
        .sum::<u64>();
    assert_eq!(model["indexed_kmers"], 1);
    assert_eq!(model["features"].as_array().unwrap().len(), 1);
    assert_eq!(model["features"][0]["exons"][0]["multiplicity"], expected);
    assert_eq!(
        model["features"][0]["transcripts"][0]["multiplicity"],
        expected
    );
    assert_eq!(model["transcripts"][0]["observable_kmers"], expected);
}

#[test]
fn documented_toy_recovers_expected_isoform_rates() {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path();
    let example = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/isoforms");
    for file in ["genome.fa", "annotation.gtf", "unitigs.fa", "datasets.tsv"] {
        fs::copy(example.join(file), root.join(file)).unwrap();
    }
    success(
        root,
        &[
            "run",
            "--isoforms",
            "--gtf",
            "annotation.gtf",
            "--genome",
            "genome.fa",
            "--index",
            "index",
            "--k",
            "7",
            "--threads",
            "4",
            "--memory-gb",
            "1",
            "--fof",
            "datasets.tsv",
            "--output",
            "features",
            "--mode",
            "unitigs",
            "--format",
            "csv",
        ],
    );
    success(
        root,
        &[
            "infer-isoforms",
            "--index",
            "index",
            "--input",
            "features",
            "--output",
            "results",
            "--threads",
            "4",
            "--compression",
            "none",
        ],
    );
    let inferred = csv(&root.join("results/datasets/sample.csv"));
    for expected in csv(&example.join("expected.csv")) {
        let actual = inferred
            .iter()
            .find(|row| row["transcript_name"] == expected["transcript_name"])
            .unwrap();
        assert_eq!(actual["status"], "identifiable");
        for column in ["abundance", "within_gene_fraction"] {
            let value = actual[column].parse::<f64>().unwrap();
            let target = expected[column].parse::<f64>().unwrap();
            assert!(
                (value - target).abs() < 1e-6,
                "{column}: {value} != {target}"
            );
        }
    }
}

#[test]
fn invalid_isoform_requests_never_publish_an_index() {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path();
    fixture(root);
    for extra in [vec!["--level", "gene"], vec!["--junctions", "all"]] {
        let mut args = vec![
            "build",
            "--isoforms",
            "--gtf",
            "reference.gtf",
            "--genome",
            "genome.fa",
            "--index",
            "bad",
            "--k",
            "7",
        ];
        args.extend(extra);
        assert!(!invoke(root, &args).status.success());
        assert!(!root.join("bad").exists());
    }
    fs::write(
        root.join("reference.gtf"),
        "chr1\tx\texon\t1\t20\t.\t+\t.\tgene_id \"g\";\n",
    )
    .unwrap();
    assert!(
        !invoke(
            root,
            &[
                "build",
                "--isoforms",
                "--gtf",
                "reference.gtf",
                "--genome",
                "genome.fa",
                "--index",
                "bad",
                "--k",
                "7"
            ]
        )
        .status
        .success()
    );
    assert!(!root.join("bad").exists());
    fs::write(root.join("reference.gtf"),"chr1\tx\texon\t1\t20\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\nchr1\tx\texon\t15\t25\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n").unwrap();
    assert!(
        !invoke(
            root,
            &[
                "build",
                "--isoforms",
                "--gtf",
                "reference.gtf",
                "--genome",
                "genome.fa",
                "--index",
                "bad",
                "--k",
                "7"
            ]
        )
        .status
        .success()
    );
    assert!(!root.join("bad").exists());
}
