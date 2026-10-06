use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{Read, Write},
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
fn canonical(seq: &[u8]) -> Vec<u8> {
    let rc: Vec<u8> = seq
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => unreachable!(),
        })
        .collect();
    seq.to_vec().min(rc)
}
fn counts(path: &Path) -> BTreeMap<String, u64> {
    csv::Reader::from_path(path)
        .unwrap()
        .deserialize::<(usize, String, u64)>()
        .map(|row| {
            let (_, name, value) = row.unwrap();
            (name, value)
        })
        .collect()
}

#[test]
fn genes_use_union_ownership_and_never_bridge_separated_exons() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let genome = b"ACGTTGCATATCCGAGGGCCATGATCTAGCCTTAACCGTTCGAGTTCATCGAG";
    let loci = [
        ("g1", 0, 12),
        ("g1", 5, 17),
        ("g1", 25, 34),
        ("g2", 0, 12),
        ("g2", 38, 52),
    ];
    let mut gtf = String::new();
    for (gene, start, end) in loci {
        gtf.push_str(&format!(
            "chr1\ttest\texon\t{}\t{end}\t.\t+\t.\tgene_id \"{gene}\";\n",
            start + 1
        ));
    }
    fs::write(dir.join("annotation.gtf"), gtf).unwrap();
    fs::write(
        dir.join("genome.fa"),
        format!(">chr1\n{}\n", std::str::from_utf8(genome).unwrap()),
    )
    .unwrap();
    fs::write(
        dir.join("reads.fa"),
        format!(">r\n{}\n", std::str::from_utf8(genome).unwrap()),
    )
    .unwrap();
    fs::write(
        dir.join("datasets.tsv"),
        "A\treads.fa\treads\nB\treads.fa\treads\n",
    )
    .unwrap();
    success(
        dir,
        &[
            "run",
            "--gtf",
            "annotation.gtf",
            "--genome",
            "genome.fa",
            "--level",
            "gene",
            "--index",
            "idx",
            "--k",
            "7",
            "--threads",
            "2",
            "--fof",
            "datasets.tsv",
            "--output",
            "results",
            "--compression",
            "none",
            "--stats",
            "--viewer-index",
        ],
    );
    assert!(dir.join("results/abundance.eai").exists());
    let mut owners: BTreeMap<Vec<u8>, BTreeSet<&str>> = BTreeMap::new();
    // Independent oracle uses overlapping-exon unions, preserving disconnected boundaries.
    for (gene, start, end) in [("g1", 0, 17), ("g1", 25, 34), ("g2", 0, 12), ("g2", 38, 52)] {
        for kmer in genome[start..end].windows(7) {
            owners.entry(canonical(kmer)).or_default().insert(gene);
        }
    }
    let mut expected = BTreeMap::from([("g1".to_string(), 0), ("g2".to_string(), 0)]);
    for kmer in genome.windows(7) {
        if let Some(genes) = owners.get(&canonical(kmer))
            && genes.len() == 1
        {
            *expected.get_mut(*genes.first().unwrap()).unwrap() += 1;
        }
    }
    assert!(expected.values().all(|v| *v > 0));
    success(
        dir,
        &[
            "export",
            "--input",
            "results",
            "--output",
            "csv",
            "--with-names",
            "--compression",
            "none",
        ],
    );
    assert_eq!(counts(&dir.join("csv/datasets/A.csv")), expected);
    type GlobalRow = (usize, String, u64, u64, u64, u64, u64, u64);
    let global: Vec<GlobalRow> = csv::Reader::from_path(dir.join("csv/global.csv"))
        .unwrap()
        .deserialize()
        .map(Result::unwrap)
        .collect();
    for (_, gene, sum, mean, median, min, max, detected) in global {
        assert_eq!(
            (sum, mean, median, min, max, detected),
            (
                2 * expected[&gene],
                expected[&gene],
                expected[&gene],
                expected[&gene],
                expected[&gene],
                2
            )
        );
    }
    let header = fs::read_to_string(dir.join("csv/datasets/A.csv")).unwrap();
    assert!(header.starts_with("gene_id,gene_name,abundance"));
    assert!(dir.join("results/genes.csv.zst").exists());
    fs::write(dir.join("targets.txt"), "# comment\ng1\n").unwrap();
    fs::write(dir.join("datasets.txt"), "B\n").unwrap();
    let bound = expected["g1"].to_string();
    success(
        dir,
        &[
            "export",
            "--input",
            "results",
            "--output",
            "selected",
            "--target-list",
            "targets.txt",
            "--dataset-list",
            "datasets.txt",
            "--min-value",
            &bound,
            "--max-value",
            &bound,
            "--compression",
            "none",
        ],
    );
    assert_eq!(
        counts(&dir.join("selected/datasets/B.csv")),
        BTreeMap::from([("g1".to_string(), expected["g1"])])
    );
    assert!(!dir.join("selected/datasets/A.csv").exists());
    assert!(!dir.join("selected/global.csv").exists());
    let min = (expected["g1"] + 1).to_string();
    success(
        dir,
        &[
            "export",
            "--input",
            "results",
            "--output",
            "empty",
            "--dataset",
            "A",
            "--target",
            "g1",
            "--min-value",
            &min,
            "--compression",
            "none",
        ],
    );
    assert!(counts(&dir.join("empty/datasets/A.csv")).is_empty());
    for (flag, value) in [
        ("--target", "unknown"),
        ("--target-id", "0"),
        ("--dataset", "missing"),
    ] {
        assert!(
            !invoke(
                dir,
                &[
                    "export", "--input", "results", "--output", "bad", flag, value
                ]
            )
            .status
            .success()
        );
        assert!(!dir.join("bad").exists());
    }
    assert!(
        !invoke(
            dir,
            &[
                "export",
                "--input",
                "results",
                "--output",
                "bad",
                "--min-value",
                "2",
                "--max-value",
                "1"
            ]
        )
        .status
        .success()
    );
    success(
        dir,
        &["pack", "--input", "results", "--output", "genes.eai"],
    );
    assert!(
        !invoke(
            dir,
            &["pack", "--input", "results", "--output", "genes.eai"]
        )
        .status
        .success()
    );
    success(
        dir,
        &[
            "export",
            "--input",
            "genes.eai",
            "--output",
            "packed-csv",
            "--with-names",
            "--compression",
            "none",
        ],
    );
    for file in ["global.csv", "datasets/A.csv", "datasets/B.csv"] {
        assert_eq!(
            fs::read(dir.join("csv").join(file)).unwrap(),
            fs::read(dir.join("packed-csv").join(file)).unwrap()
        );
    }
    success(
        dir,
        &[
            "export",
            "--input",
            "genes.eai",
            "--output",
            "packed-subset",
            "--target-id",
            "2",
            "--dataset",
            "A",
            "--compression",
            "none",
        ],
    );
    let selected = counts(&dir.join("packed-subset/datasets/A.csv"));
    assert_eq!(selected.len(), 1);
    let metadata: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("idx/metadata.json")).unwrap()).unwrap();
    let second_name = metadata["exons"][1]["name"].as_str().unwrap();
    assert_eq!(selected[second_name], expected[second_name]);
    fs::write(dir.join("empty-list.txt"), "# no selections\n").unwrap();
    for flag in ["--target-list", "--dataset-list"] {
        assert!(
            !invoke(
                dir,
                &[
                    "export",
                    "--input",
                    "results",
                    "--output",
                    "empty-list-output",
                    flag,
                    "empty-list.txt"
                ]
            )
            .status
            .success()
        );
        assert!(!dir.join("empty-list-output").exists());
    }
    fs::write(
        dir.join("missing-gene.gtf"),
        "chr1\ttest\texon\t1\t12\t.\t+\t.\ttranscript_id \"t1\";\n",
    )
    .unwrap();
    assert!(
        !invoke(
            dir,
            &[
                "build",
                "--gtf",
                "missing-gene.gtf",
                "--genome",
                "genome.fa",
                "--level",
                "gene",
                "--index",
                "missing-gene-index",
                "--k",
                "7"
            ]
        )
        .status
        .success()
    );
    assert!(!dir.join("missing-gene-index").exists());
    let mut corrupted = fs::read(dir.join("genes.eai")).unwrap();
    *corrupted.last_mut().unwrap() ^= 1;
    fs::write(dir.join("corrupt.eai"), corrupted).unwrap();
    assert!(
        !invoke(
            dir,
            &[
                "export",
                "--input",
                "corrupt.eai",
                "--output",
                "bad-archive"
            ]
        )
        .status
        .success()
    );
    assert!(!dir.join("bad-archive").exists());
}

fn encode(raw: &[u8], codec: &str) -> Vec<u8> {
    match codec {
        "zstd" => zstd::stream::encode_all(raw, 3).unwrap(),
        "gz" => {
            let mut w = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            w.write_all(raw).unwrap();
            w.finish().unwrap()
        }
        "xz" => {
            let mut w = xz2::write::XzEncoder::new(Vec::new(), 1);
            w.write_all(raw).unwrap();
            w.finish().unwrap()
        }
        _ => raw.to_vec(),
    }
}
#[test]
fn portable_index_preserves_all_codecs_bits_and_full_u64_without_requantizing() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let mut hash = Sha256::new();
    hash.update(b"EXPRESSO-reference-v1\0");
    hash.update(2u64.to_le_bytes());
    for name in ["e1", "e2"] {
        hash.update(2u64.to_le_bytes());
        hash.update(name);
        hash.update(7u64.to_le_bytes());
        hash.update(1u64.to_le_bytes());
    }
    let reference: [u8; 32] = hash.finalize().into();
    let hex: String = reference.iter().map(|b| format!("{b:02x}")).collect();
    for bits in 2..=16u8 {
        let mut raw = b"EXPRAB01".to_vec();
        raw.extend([bits, 0, 0, 0]);
        raw.extend(2u64.to_le_bytes());
        raw.extend(2u32.to_le_bytes());
        let packed = (bits as usize * 2).div_ceil(8);
        raw.extend((packed as u64).to_le_bytes());
        raw.extend(u64::MAX.to_le_bytes());
        raw.extend(reference);
        raw.extend(0u64.to_le_bytes());
        raw.extend(u64::MAX.to_le_bytes());
        let payload = 1u64 << (bits as usize);
        raw.extend(&payload.to_le_bytes()[..packed]);
        let crc = crc32fast::hash(&raw);
        raw.extend(crc.to_le_bytes());
        for codec in ["none", "gz", "xz", "zstd"] {
            let source = format!("input-{bits}-{codec}");
            let source_dir = dir.join(&source);
            fs::create_dir(&source_dir).unwrap();
            fs::write(
                source_dir.join("ref"),
                b"exon_id,exon_name,length,unique_kmers\n1,e1,7,1\n2,e2,7,1\n",
            )
            .unwrap();
            fs::write(source_dir.join("vector"), encode(&raw, codec)).unwrap();
            let manifest = serde_json::json!({"format":"compact","compact_format_version":1,"reference":{"file":"ref","sha256":hex,"targets":2},"datasets":[{"name":"global","output":"vector"}],"global":{"output":"vector"}});
            fs::write(
                source_dir.join("manifest.json"),
                serde_json::to_vec(&manifest).unwrap(),
            )
            .unwrap();
            let destination = format!("{source}.eai");
            success(dir, &["pack", "--input", &source, "--output", &destination]);
            let bytes = fs::read(dir.join(destination)).unwrap();
            assert_eq!(&bytes[..8], b"EXPREAI1");
            let offset = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
            let length = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
            assert_eq!(offset + length + 32, bytes.len());
            assert_eq!(
                &bytes[offset + length..],
                Sha256::digest(&bytes[offset..offset + length]).as_slice()
            );
            let cat: serde_json::Value =
                serde_json::from_slice(&bytes[offset..offset + length]).unwrap();
            assert_eq!(cat["level"], "exon");
            assert_eq!(cat["datasets"][0]["global"], false);
            assert_eq!(cat["datasets"][1]["global"], true);
            for entry in cat["datasets"].as_array().unwrap() {
                let start = entry["offset"].as_u64().unwrap() as usize;
                let len = entry["length"].as_u64().unwrap() as usize;
                let mut decoded = Vec::new();
                zstd::stream::read::Decoder::new(&bytes[start..start + len])
                    .unwrap()
                    .read_to_end(&mut decoded)
                    .unwrap();
                assert_eq!(decoded, raw);
            }
            // Corruption must prevent publication of a partial archive.
            let mut bad = raw.clone();
            bad[72] ^= 1;
            fs::write(source_dir.join("vector"), encode(&bad, codec)).unwrap();
            assert!(
                !invoke(dir, &["pack", "--input", &source, "--output", "bad.eai"])
                    .status
                    .success()
            );
            assert!(!dir.join("bad.eai").exists());
        }
    }
}
