use std::{
    fs,
    io::{Read, Write},
    path::Path,
    process::{Command, Output},
};

fn invoke(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_expresso"))
        .args(args)
        .current_dir(dir)
        .env("PATH", "")
        .output()
        .unwrap()
}
fn success(dir: &Path, args: &[&str]) {
    let out = invoke(dir, args);
    assert!(
        out.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}
fn compressed(data: &[u8], codec: &str) -> Vec<u8> {
    match codec {
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
fn decoded(path: &Path, codec: &str) -> String {
    let file = fs::File::open(path).unwrap();
    let mut reader: Box<dyn Read> = match codec {
        "gz" => Box::new(flate2::read::MultiGzDecoder::new(file)),
        "xz" => Box::new(xz2::read::XzDecoder::new_multi_decoder(file)),
        "zstd" => Box::new(zstd::stream::read::Decoder::new(file).unwrap()),
        _ => Box::new(file),
    };
    let mut s = String::new();
    reader.read_to_string(&mut s).unwrap();
    s
}
fn row(contig: &str, start: usize, end: usize, strand: &str, attributes: &str) -> String {
    format!("{contig}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
}

#[test]
fn strands_deduplication_attributes_and_all_codecs() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let rows = [
        row(
            "chr1",
            1,
            15,
            "+",
            "gene_id \"g2\"; exon_id \"e2\"; gene_name \"with; semicolon\"; tag \"a\"; tag \"b\";",
        ),
        row("chr1", 1, 15, "+", "exon_id \"e1\"; gene_id \"g1\";"),
        row("chr1", 1, 15, "-", "."),
        row("chr1", 3, 4, "+", "gene_id g3;"),
        row("chr1", 15, 15, "+", "."),
        row("chr2", 2, 6, "-", "exon_id \"e3\";"),
    ];
    let expected = concat!(
        ">exon_000001|chr2:2-6:- exon_ids=e3\nCGGTT\n",
        ">exon_000002|chr1:1-15:+ exon_ids=e1,e2 gene_ids=g1,g2\nACGTRYSWKMBDHVN\n",
        ">exon_000003|chr1:1-15:-\nNBDHVKMWSRYACGT\n",
        ">exon_000004|chr1:3-4:+ gene_ids=g3\nGT\n",
        ">exon_000005|chr1:15-15:+\nN\n"
    );
    for (i, codec) in ["none", "gz", "xz", "zstd"].iter().enumerate() {
        // Concatenated compression members must all be read; suffixes are irrelevant.
        let mut gtf = compressed(b"# annotation\r\n", codec);
        let text = if i % 2 == 0 {
            rows.concat()
        } else {
            rows.iter().rev().cloned().collect::<String>()
        };
        gtf.extend(compressed(text.as_bytes(), codec));
        fs::write(dir.join("annotation.data"), gtf).unwrap();
        fs::write(
            dir.join("one.data"),
            compressed(b">chr2 description\nAAACCGT\n", codec),
        )
        .unwrap();
        fs::write(
            dir.join("two.data"),
            compressed(
                b">unused\nACGT\n>chr1 wrapped\r\naCgTrYs\r\nWkMbDhVn\r\n",
                codec,
            ),
        )
        .unwrap();
        let output = format!("exons-{codec}");
        success(
            dir,
            &[
                "extract-exons",
                "--gtf",
                "annotation.data",
                "--genome",
                "one.data",
                "--reference",
                "two.data",
                "--output",
                &output,
                "--compression",
                codec,
            ],
        );
        assert_eq!(decoded(&dir.join(output), codec), expected);
    }
}

#[test]
fn extraction_errors_never_publish_partial_fasta() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let valid = row("chr1", 1, 4, "+", ".");
    let genome = ">chr1\nACGT\n";
    let cases = vec![
        (row("chr1", 0, 4, "+", "."), genome.to_string()),
        (row("chr1", 3, 2, "+", "."), genome.to_string()),
        (row("chr1", 1, 5, "+", "."), genome.to_string()),
        (row("missing", 1, 4, "+", "."), genome.to_string()),
        (row("chr1", 1, 4, ".", "."), genome.to_string()),
        (
            row("chr1", 1, 4, "+", "exon_id \"unterminated;"),
            genome.to_string(),
        ),
        ("bad\tgtf\n".to_string(), genome.to_string()),
        ("# no exons\n".to_string(), genome.to_string()),
        (valid.clone(), format!("{genome}{genome}")),
        (valid.clone(), ">chr1\nAC?T\n".to_string()),
        (valid.clone(), "@chr1\nACGT\n+\nIIII\n".to_string()),
        (valid.clone(), String::new()),
    ];
    for (gtf, fasta) in cases {
        fs::write(dir.join("a.gtf"), gtf).unwrap();
        fs::write(dir.join("g.fa"), fasta).unwrap();
        let out = invoke(
            dir,
            &[
                "exons", "--gtf", "a.gtf", "--genome", "g.fa", "-o", "out.fa",
            ],
        );
        assert!(!out.status.success());
        assert!(!dir.join("out.fa").exists());
        assert_eq!(
            fs::read_dir(dir).unwrap().count(),
            2,
            "temporary output was not cleaned up"
        );
    }
    fs::write(dir.join("a.gtf"), &valid).unwrap();
    fs::write(dir.join("g.fa"), genome).unwrap();
    let out = invoke(
        dir,
        &[
            "exons", "--gtf", "a.gtf", "--genome", "g.fa", "--genome", "g.fa", "-o", "out.fa",
        ],
    );
    assert!(!out.status.success());
    assert!(!dir.join("out.fa").exists());
    fs::write(dir.join("out.fa"), "keep me").unwrap();
    assert!(
        !invoke(
            dir,
            &[
                "exons", "--gtf", "a.gtf", "--genome", "g.fa", "-o", "out.fa"
            ]
        )
        .status
        .success()
    );
    assert_eq!(fs::read_to_string(dir.join("out.fa")).unwrap(), "keep me");
    let mut corrupt = compressed(valid.as_bytes(), "gz");
    corrupt.truncate(corrupt.len() - 5);
    fs::write(dir.join("a.gtf"), corrupt).unwrap();
    assert!(
        !invoke(
            dir,
            &[
                "exons",
                "--gtf",
                "a.gtf",
                "--genome",
                "g.fa",
                "-o",
                "corrupt.fa"
            ]
        )
        .status
        .success()
    );
    assert!(!dir.join("corrupt.fa").exists());
}

#[test]
fn build_and_run_accept_gtf_instead_of_exon_fasta() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let genome = ">chr1\nAAACCGTAGCTTGACCGATCGATACGATCGTTAGCTACGATC\n";
    fs::write(dir.join("g.fa"), genome).unwrap();
    fs::write(
        dir.join("a.gtf"),
        row("chr1", 3, 19, "+", ".") + &row("chr1", 22, 38, "-", "."),
    )
    .unwrap();
    success(
        dir,
        &[
            "extract-exons",
            "--gtf",
            "a.gtf",
            "--genome",
            "g.fa",
            "-o",
            "exons.fa",
        ],
    );
    success(
        dir,
        &[
            "build", "--gtf", "a.gtf", "--genome", "g.fa", "-i", "from-gtf", "-k", "7", "-t", "2",
        ],
    );
    assert_eq!(
        fs::read(dir.join("exons.fa")).unwrap(),
        fs::read(dir.join("from-gtf/reference-exons.fa")).unwrap()
    );
    success(
        dir,
        &[
            "build",
            "--exons",
            "exons.fa",
            "-i",
            "from-fasta",
            "-k",
            "7",
            "-t",
            "2",
        ],
    );
    assert_eq!(
        fs::read(dir.join("from-gtf/metadata.json")).unwrap(),
        fs::read(dir.join("from-fasta/metadata.json")).unwrap()
    );
    fs::write(dir.join("datasets.tsv"), "sample\tg.fa\treads\n").unwrap();
    success(
        dir,
        &[
            "run",
            "--gtf",
            "a.gtf",
            "--genome",
            "g.fa",
            "-i",
            "run-index",
            "-k",
            "7",
            "-t",
            "2",
            "-f",
            "datasets.tsv",
            "-o",
            "result",
        ],
    );
    assert!(dir.join("result/global.eab.zst").is_file());
    for args in [
        vec!["build", "--gtf", "a.gtf", "-i", "bad"],
        vec!["build", "--genome", "g.fa", "-i", "bad"],
        vec![
            "build", "--exons", "exons.fa", "--gtf", "a.gtf", "--genome", "g.fa", "-i", "bad",
        ],
        vec!["build", "-i", "bad"],
    ] {
        assert!(!invoke(dir, &args).status.success());
        assert!(!dir.join("bad").exists());
    }
}
