use std::{
    collections::HashMap,
    fs,
    io::{Read, Write},
    path::Path,
    process::{Command, Output},
};

// Serialize graph-building fixtures to bound their combined memory use.
static GGCAT_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn invoke(args: &[&str], cwd: &Path) -> Output {
    // These fixtures assert legacy exact CSV output; compact is tested below.
    let mut command = Command::new(env!("CARGO_BIN_EXE_expresso"));
    command.args(args);
    if matches!(args.first(), Some(&"quantify" | &"run")) && !args.contains(&"--format") {
        command.args(["--format", "csv"]);
    }
    // Index construction must work without a separately installed GGCAT.
    command.env("PATH", "").current_dir(cwd).output().unwrap()
}
fn success(args: &[&str], cwd: &Path) {
    let out = invoke(args, cwd);
    assert!(
        out.status.success(),
        "{args:?}\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
fn canonical(kmer: &[u8]) -> Vec<u8> {
    let fw: Vec<_> = kmer.iter().map(u8::to_ascii_uppercase).collect();
    let rc: Vec<_> = fw
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            _ => b'N',
        })
        .collect();
    fw.min(rc)
}
fn kmers(seq: &[u8], k: usize) -> impl Iterator<Item = Vec<u8>> + '_ {
    seq.windows(k)
        .filter(|s| s.iter().all(|b| b"ACGTacgt".contains(b)))
        .map(canonical)
}
fn oracle(exons: &[Vec<u8>], records: &[(Vec<u8>, u64)], k: usize) -> Vec<u64> {
    let mut owners: HashMap<Vec<u8>, Option<usize>> = HashMap::new();
    for (e, seq) in exons.iter().enumerate() {
        for key in kmers(seq, k) {
            owners
                .entry(key)
                .and_modify(|old| {
                    if *old != Some(e) {
                        *old = None
                    }
                })
                .or_insert(Some(e));
        }
    }
    let mut counts = vec![0; exons.len()];
    for (seq, weight) in records {
        for key in kmers(seq, k) {
            if let Some(Some(e)) = owners.get(&key) {
                counts[*e] += weight;
            }
        }
    }
    counts
        .into_iter()
        .map(|c| (c + 500_000) / 1_000_000)
        .collect()
}
fn compress(data: &[u8], kind: &str) -> Vec<u8> {
    match kind {
        "gz" => {
            let mut w = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            w.write_all(data).unwrap();
            w.finish().unwrap()
        }
        "xz" => {
            let mut w = xz2::write::XzEncoder::new(Vec::new(), 1);
            w.write_all(data).unwrap();
            w.finish().unwrap()
        }
        "zstd" => zstd::stream::encode_all(data, 1).unwrap(),
        _ => data.to_vec(),
    }
}
fn csv_rows(path: &Path, kind: &str) -> Vec<Vec<String>> {
    let file = fs::File::open(path).unwrap();
    let reader: Box<dyn Read> = match kind {
        "gz" => Box::new(flate2::read::MultiGzDecoder::new(file)),
        "xz" => Box::new(xz2::read::XzDecoder::new_multi_decoder(file)),
        "zstd" => Box::new(zstd::stream::read::Decoder::new(file).unwrap()),
        _ => Box::new(file),
    };
    csv::Reader::from_reader(reader)
        .records()
        .map(|r| r.unwrap().iter().map(str::to_string).collect())
        .collect()
}
fn values(path: &Path, kind: &str) -> Vec<u64> {
    csv_rows(path, kind)
        .iter()
        .map(|r| r[2].parse().unwrap())
        .collect()
}

#[test]
fn long_kmers_and_all_shared_exons() {
    let _guard = GGCAT_TEST_LOCK.lock().unwrap();
    let temp = tempfile::tempdir().unwrap();
    let dir = temp.path();
    let mut seed = 918_273u64;
    let mut random = || {
        (0..240)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                b"ACGT"[(seed & 3) as usize]
            })
            .collect::<Vec<_>>()
    };
    let mut exons = vec![random(), random(), random(), random()];
    let shared = exons[0][50..150].to_vec();
    exons[1][50..150].copy_from_slice(&shared);
    exons.push(exons[2].clone());
    let mut reference = Vec::new();
    for (i, seq) in exons.iter().enumerate() {
        writeln!(reference, ">e{i}").unwrap();
        reference.extend(seq);
        reference.push(b'\n');
    }
    fs::write(dir.join("exons.fa"), reference).unwrap();
    let reads: Vec<_> = exons
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut seq = if i % 2 == 0 { canonical(s) } else { s.clone() };
            seq[120] = b'N';
            (seq, 1_000_000)
        })
        .collect();
    let mut fasta = Vec::new();
    for (i, (seq, _)) in reads.iter().enumerate() {
        writeln!(fasta, ">r{i}").unwrap();
        fasta.extend(seq);
        fasta.push(b'\n');
    }
    fs::write(dir.join("reads.fa"), fasta).unwrap();
    // Exercise both absolute paths and inferred dataset naming.
    fs::write(
        dir.join("files.txt"),
        format!("{}\n", dir.join("reads.fa").display()),
    )
    .unwrap();
    for k in [31, 63] {
        let index = format!("index{k}");
        let out = format!("out{k}");
        success(
            &[
                "run",
                "-e",
                "exons.fa",
                "-i",
                &index,
                "-k",
                &k.to_string(),
                "-f",
                "files.txt",
                "-o",
                &out,
                "-t",
                "3",
                "--mode",
                "reads",
                "--batch-bases",
                "32",
            ],
            dir,
        );
        assert_eq!(
            values(&dir.join(format!("{out}/global.csv.zst")), "zstd"),
            oracle(&exons, &reads, k)
        );
    }
    fs::write(dir.join("identical.fa"), b">a\nAAAAAAA\n>b\nTTTTTTT\n").unwrap();
    success(
        &[
            "run",
            "-e",
            "identical.fa",
            "-i",
            "shared",
            "-k",
            "7",
            "-f",
            "files.txt",
            "-o",
            "shared-out",
            "-t",
            "2",
        ],
        dir,
    );
    assert_eq!(
        values(&dir.join("shared-out/global.csv.zst"), "zstd"),
        vec![0, 0]
    );
}

#[test]
fn bundled_ggcat_sshash_against_independent_oracle() {
    let _guard = GGCAT_TEST_LOCK.lock().unwrap();
    let temp = tempfile::tempdir().unwrap();
    let dir = temp.path();
    let exons = vec![
        b"AAACCGTAGCTTGACCGATCGATACGATCGTTAGCTACGATC".to_vec(),
        b"GGTACCTTTAGCGATCGATACGATCGTTAGCAAAACCGTAAA".to_vec(),
        b"AAAAAAAAAAAANacgtacgtacgtacgtNNNTTACGGCATGCA".to_vec(),
        b"NNNNNAC".to_vec(),
        b"CCCTTTGGGAAACCCGGGTTTCCCAAA".to_vec(),
    ];
    let mut reference = Vec::new();
    for (i, seq) in exons.iter().enumerate() {
        writeln!(reference, ">exon,{} description", i / 2).unwrap();
        for line in seq.chunks(13) {
            reference.extend_from_slice(line);
            reference.extend_from_slice(b"\r\n");
        }
    }
    fs::write(dir.join("exons.data"), compress(&reference, "xz")).unwrap();
    success(
        &[
            "build",
            "-e",
            "exons.data",
            "-i",
            "index",
            "-k",
            "7",
            "-t",
            "2",
            "--memory-gb",
            "1",
        ],
        dir,
    );
    let mut query_records: Vec<_> = exons.iter().map(|s| (s.clone(), 1_000_000)).collect();
    query_records.push((b"TGCAATCCGGTAGCTAGCTNNNAAAAAAA".to_vec(), 1_000_000));
    query_records.push((canonical(&exons[0]), 1_000_000));
    query_records.push((b"AAAAAAA".repeat(25_000), 1_000_000));
    let mut fasta = Vec::new();
    let mut fastq = Vec::new();
    for (i, (seq, _)) in query_records.iter().enumerate() {
        writeln!(fasta, ">r{i}").unwrap();
        fasta.extend(seq);
        fasta.push(b'\n');
        writeln!(fastq, "@r{i}").unwrap();
        fastq.extend(seq);
        fastq.extend(b"\n+\n");
        fastq.extend(vec![b'I'; seq.len()]);
        fastq.push(b'\n');
    }
    let expected = oracle(&exons, &query_records, 7);
    let mut fof = String::new();
    fs::create_dir(dir.join("inputs")).unwrap();
    for (i, kind) in ["none", "gz", "xz", "zstd"].iter().enumerate() {
        let data = if i % 2 == 0 { &fasta } else { &fastq };
        fs::write(
            dir.join(format!("inputs/{kind}.data")),
            compress(data, kind),
        )
        .unwrap();
        fof.push_str(&format!("{kind}\t{kind}.data\treads\n"));
    }
    for kind in ["gz", "xz", "zstd"] {
        let mut data = compress(&fasta, kind);
        data.extend(compress(&fasta, kind));
        fs::write(dir.join(format!("inputs/multi_{kind}.data")), data).unwrap();
        fof.push_str(&format!("multi_{kind}\tmulti_{kind}.data\treads\n"));
    }
    fs::write(dir.join("inputs/empty.data"), compress(b"", "zstd")).unwrap();
    fof.push_str("empty\tempty.data\treads\n");
    fs::write(dir.join("inputs/list.tsv"), fof).unwrap();
    success(
        &[
            "quantify",
            "-i",
            "index",
            "-f",
            "inputs/list.tsv",
            "-o",
            "out",
            "--stats",
            "-t",
            "4",
            "-j",
            "2",
            "--batch-bases",
            "512",
        ],
        dir,
    );
    for kind in ["none", "gz", "xz", "zstd"] {
        assert_eq!(
            values(&dir.join(format!("out/datasets/{kind}.csv.zst")), "zstd"),
            expected,
            "input {kind}"
        );
    }
    for kind in ["gz", "xz", "zstd"] {
        assert_eq!(
            values(
                &dir.join(format!("out/datasets/multi_{kind}.csv.zst")),
                "zstd"
            ),
            expected.iter().map(|v| v * 2).collect::<Vec<_>>(),
            "concatenated {kind}"
        );
    }
    assert_eq!(
        values(&dir.join("out/datasets/empty.csv.zst"), "zstd"),
        vec![0; exons.len()]
    );
    let global = csv_rows(&dir.join("out/global.csv.zst"), "zstd");
    for (e, row) in global.iter().enumerate() {
        assert_eq!(row[2].parse::<u64>().unwrap(), expected[e] * 10);
        assert_eq!(row[3].parse::<u64>().unwrap(), (expected[e] * 10 + 4) / 8);
        assert_eq!(row[4].parse::<u64>().unwrap(), expected[e]);
        assert_eq!(row[5], "0");
        assert_eq!(row[6].parse::<u64>().unwrap(), expected[e] * 2);
        assert_eq!(
            row[7].parse::<u64>().unwrap(),
            if expected[e] > 0 { 7 } else { 0 }
        );
    }
    let weights = [1_700_000, 2_500_000, 0, 3_000_000, 4_000_000];
    let weighted: Vec<_> = exons.iter().cloned().zip(weights).collect();
    let mut unitigs = Vec::new();
    for (i, (seq, weight)) in weighted.iter().enumerate() {
        writeln!(
            unitigs,
            ">u{i} {}:f:{}.{:06}",
            if i % 2 == 0 { "ka" } else { "km" },
            weight / 1_000_000,
            weight % 1_000_000
        )
        .unwrap();
        unitigs.extend(seq);
        unitigs.push(b'\n');
    }
    unitigs.pop();
    fs::write(dir.join("unitigs.fa.gz"), compress(&unitigs, "gz")).unwrap();
    fs::write(dir.join("unitigs.tsv"), "unitigs\tunitigs.fa.gz\tunitigs\n").unwrap();
    let weighted_expected = oracle(&exons, &weighted, 7);
    // The default format is compact; the small fixture fits in 8 bits exactly.
    let compact = Command::new(env!("CARGO_BIN_EXE_expresso"))
        .env("PATH", "")
        .args([
            "quantify",
            "-i",
            "index",
            "-f",
            "unitigs.tsv",
            "-o",
            "compact",
            "--stats",
            "-t",
            "3",
        ])
        .current_dir(dir)
        .output()
        .unwrap();
    assert!(
        compact.status.success(),
        "{}",
        String::from_utf8_lossy(&compact.stderr)
    );
    assert!(dir.join("compact/datasets/unitigs.eab.zst").exists());
    assert!(!dir.join("compact/datasets/unitigs.csv.zst").exists());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("compact/manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["format"], "compact");
    assert_eq!(manifest["datasets"][0]["quantization"]["bits"], 8);
    assert_eq!(
        manifest["datasets"][0]["quantization"]["exact_values"],
        exons.len()
    );
    success(
        &[
            "export",
            "-i",
            "compact",
            "-o",
            "decoded",
            "--compression",
            "none",
            "-t",
            "2",
        ],
        dir,
    );
    let decoded = csv_rows(&dir.join("decoded/datasets/unitigs.csv"), "none");
    assert!(decoded.iter().all(|row| row.len() == 1));
    assert_eq!(
        decoded
            .iter()
            .map(|row| row[0].parse::<u64>().unwrap())
            .collect::<Vec<_>>(),
        weighted_expected
    );
    let decoded_global = csv_rows(&dir.join("decoded/global.csv"), "none");
    for (i, row) in decoded_global.iter().enumerate() {
        let v = weighted_expected[i];
        assert_eq!(
            row.iter()
                .map(|s| s.parse::<u64>().unwrap())
                .collect::<Vec<_>>(),
            vec![v, v, v, v, v, u64::from(v > 0)]
        );
    }
    success(
        &[
            "export",
            "-i",
            "compact",
            "-o",
            "named",
            "--with-names",
            "--dataset",
            "unitigs",
            "--compression",
            "xz",
        ],
        dir,
    );
    assert_eq!(
        values(&dir.join("named/datasets/unitigs.csv.xz"), "xz"),
        weighted_expected
    );
    assert!(!dir.join("named/global.csv.xz").exists());
    assert!(
        !invoke(
            &[
                "export",
                "-i",
                "compact",
                "-o",
                "missing",
                "--dataset",
                "not_present"
            ],
            dir
        )
        .status
        .success()
    );
    assert!(!dir.join("missing").exists());
    success(
        &[
            "quantify",
            "-i",
            "index",
            "-f",
            "unitigs.tsv",
            "-o",
            "compact-5",
            "--format",
            "compact",
            "--bits",
            "5",
            "--compression",
            "none",
            "-t",
            "2",
        ],
        dir,
    );
    success(
        &[
            "export",
            "-i",
            "compact-5",
            "-o",
            "decoded-5",
            "--global-only",
            "--compression",
            "gz",
        ],
        dir,
    );
    assert!(dir.join("decoded-5/global.csv.gz").exists());
    let broken = dir.join("compact-5/global.eab");
    let mut bytes = fs::read(&broken).unwrap();
    let end = bytes.len() - 1;
    bytes[end] ^= 1;
    fs::write(&broken, bytes).unwrap();
    assert!(
        !invoke(
            &[
                "export",
                "-i",
                "compact-5",
                "-o",
                "bad-export",
                "--global-only"
            ],
            dir
        )
        .status
        .success()
    );
    assert!(!dir.join("bad-export").exists());
    for (i, (kind, ext)) in [
        ("none", "csv"),
        ("gz", "csv.gz"),
        ("xz", "csv.xz"),
        ("zstd", "csv.zst"),
    ]
    .iter()
    .enumerate()
    {
        let out = format!("weighted_{kind}");
        success(
            &[
                "quantify",
                "-i",
                "index",
                "-f",
                "unitigs.tsv",
                "-o",
                &out,
                "--compression",
                kind,
                "-t",
                if i == 0 { "1" } else { "3" },
                "--batch-bases",
                "8",
            ],
            dir,
        );
        assert_eq!(
            values(&dir.join(format!("{out}/datasets/unitigs.{ext}")), kind),
            weighted_expected
        );
    }
    fs::write(dir.join("bad.fa"), b">u km:f:NaN\nAAAAAAA\n").unwrap();
    fs::write(dir.join("bad.tsv"), "bad\tbad.fa\tunitigs\n").unwrap();
    assert!(
        !invoke(
            &["quantify", "-i", "index", "-f", "bad.tsv", "-o", "failed"],
            dir
        )
        .status
        .success()
    );
    assert!(!dir.join("failed").exists());
    fs::write(
        dir.join("mixed.tsv"),
        format!(
            "bad\tbad.fa\tunitigs\n{}",
            fs::read_to_string(dir.join("unitigs.tsv")).unwrap()
        ),
    )
    .unwrap();
    success(
        &[
            "quantify",
            "-i",
            "index",
            "-f",
            "mixed.tsv",
            "-o",
            "mixed",
            "--keep-going",
            "--format",
            "compact",
            "--jobs",
            "2",
            "-t",
            "4",
        ],
        dir,
    );
    let mixed: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("mixed/manifest.json")).unwrap()).unwrap();
    assert_eq!(mixed["coverage"]["complete"], false);
    assert_eq!(mixed["coverage"]["completed_datasets"], 1);
    assert_eq!(mixed["coverage"]["failed_datasets"], 1);
    assert_eq!(mixed["failed_datasets"][0]["name"], "bad");
    assert!(!dir.join("mixed/datasets/bad.eab.zst").exists());
    success(
        &[
            "export",
            "-i",
            "mixed",
            "-o",
            "mixed-csv",
            "--global-only",
            "--with-names",
            "--compression",
            "none",
        ],
        dir,
    );
    assert_eq!(
        values(&dir.join("mixed-csv/global.csv"), "none"),
        weighted_expected
    );
    assert!(
        !invoke(
            &[
                "quantify",
                "-i",
                "index",
                "-f",
                "bad.tsv",
                "-o",
                "all-failed",
                "--keep-going"
            ],
            dir
        )
        .status
        .success()
    );
    assert!(!dir.join("all-failed").exists());
    assert!(
        !invoke(
            &[
                "quantify",
                "-i",
                "index",
                "-f",
                "mixed.tsv",
                "-o",
                "invalid-options",
                "--keep-going",
                "--stats"
            ],
            dir
        )
        .status
        .success()
    );
    let mut corrupt = compress(&fasta, "gz");
    corrupt.truncate(corrupt.len() - 10);
    fs::write(dir.join("bad.fa"), corrupt).unwrap();
    fs::write(dir.join("bad.tsv"), "bad\tbad.fa\treads\n").unwrap();
    assert!(
        !invoke(
            &["quantify", "-i", "index", "-f", "bad.tsv", "-o", "corrupt"],
            dir
        )
        .status
        .success()
    );
    assert!(!dir.join("corrupt").exists());
    // A query-side error must wake/cancel all MPMC consumers and unblock the
    // bounded producer, rather than hanging or publishing partial counts.
    let mut overflowing = b">u km:f:18446744073709\n".to_vec();
    overflowing.extend(vec![b'A'; 100_000]);
    overflowing.push(b'\n');
    fs::write(dir.join("overflow.fa"), overflowing).unwrap();
    fs::write(dir.join("overflow.tsv"), "overflow\toverflow.fa\tunitigs\n").unwrap();
    let overflow = invoke(
        &[
            "quantify",
            "-i",
            "index",
            "-f",
            "overflow.tsv",
            "-o",
            "overflow-out",
            "-t",
            "4",
            "--batch-bases",
            "64",
        ],
        dir,
    );
    assert!(!overflow.status.success());
    assert!(String::from_utf8_lossy(&overflow.stderr).contains("count overflow"));
    assert!(!dir.join("overflow-out").exists());
    assert!(
        !invoke(
            &["build", "-e", "exons.data", "-i", "even", "-k", "30"],
            dir
        )
        .status
        .success()
    );
    assert!(!dir.join("even").exists());
    success(&["inspect", "-i", "index"], dir);
}
