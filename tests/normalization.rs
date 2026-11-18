//! Verify normalizers against hand-counted weighted windows, before rounding.
use serde_json::{Value, json};
use std::{fs, path::Path, process::Command};

fn success(dir: &Path, args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_expresso"))
        .current_dir(dir)
        .args(args)
        .env("PATH", "")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn weighted_totals_survive_parallel_batches_rounding_and_packing() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let reference =
        fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/viewer/genes.fa"))
            .unwrap();
    let seq = reference.lines().nth(3).unwrap();
    assert!(!seq.contains("AAAAAAA") && !seq.contains("TTTTTTT"));
    fs::write(dir.join("reference.fa"), format!(">target\n{seq}\n")).unwrap();
    fs::write(
        dir.join("weighted.fa"),
        format!(
            ">w ka:f:1.25\n{seq}\n>unmatched ka:f:0.125\n{}\n>ambiguous ka:f:0.5\n{}N{}\n",
            "A".repeat(20),
            &seq[..20],
            &seq[20..]
        ),
    )
    .unwrap();
    fs::write(dir.join("reads.fa"), format!(">r\n{seq}\n")).unwrap();
    fs::write(dir.join("empty.fa"), ">empty\nNNNNNNNNNN\n").unwrap();
    fs::write(
        dir.join("datasets.tsv"),
        "weighted\tweighted.fa\tunitigs\nreads\treads.fa\treads\nempty\tempty.fa\treads\n",
    )
    .unwrap();
    success(
        dir,
        &[
            "build",
            "--exons",
            "reference.fa",
            "--index",
            "index",
            "--k",
            "7",
            "--threads",
            "2",
        ],
    );
    for (threads, output) in [("1", "single"), ("4", "parallel")] {
        success(
            dir,
            &[
                "quantify",
                "--index",
                "index",
                "--fof",
                "datasets.tsv",
                "--output",
                output,
                "--threads",
                threads,
                "--batch-bases",
                "16",
                "--bits",
                "16",
                "--viewer-index",
            ],
        );
        let manifest: Value =
            serde_json::from_slice(&fs::read(dir.join(output).join("manifest.json")).unwrap())
                .unwrap();
        assert_eq!(manifest["datasets"][0]["metrics"]["valid_kmers"], 76);
        assert_eq!(
            manifest["datasets"][0]["total_kmers"],
            json!({"numerator":"58250000","scale":1_000_000})
        );
        assert_eq!(
            manifest["datasets"][1]["total_kmers"],
            json!({"numerator":"34","scale":1})
        );
        assert_eq!(
            manifest["datasets"][2]["total_kmers"],
            json!({"numerator":"0","scale":1})
        );
        assert_eq!(
            manifest["global"]["total_kmers"],
            json!({"numerator":"92250000","scale":1_000_000})
        );
        let destination = format!("{output}-csv");
        success(
            dir,
            &[
                "export",
                "--input",
                &format!("{output}/abundance.eai"),
                "--output",
                &destination,
                "--normalization",
                "length-cpm",
                "--compression",
                "none",
            ],
        );
        for (dataset, expected) in [
            ("weighted", 57.0 / 34.0 / 58.25 * 1e6),
            ("reads", 1e6 / 34.0),
            ("empty", 0.0),
        ] {
            let mut reader = csv::Reader::from_path(
                dir.join(&destination)
                    .join(format!("datasets/{dataset}.csv")),
            )
            .unwrap();
            let row = reader.records().next().unwrap().unwrap();
            assert!((row[0].parse::<f64>().unwrap() - expected).abs() < 1e-8);
        }
    }
    // Old indexes still quantify. Missing lengths cannot silently normalize.
    let mut meta: Value =
        serde_json::from_slice(&fs::read(dir.join("index/metadata.json")).unwrap()).unwrap();
    meta["exons"][0]
        .as_object_mut()
        .unwrap()
        .remove("kmer_length");
    fs::write(
        dir.join("index/metadata.json"),
        serde_json::to_vec(&meta).unwrap(),
    )
    .unwrap();
    success(
        dir,
        &[
            "quantify",
            "--index",
            "index",
            "--fof",
            "datasets.tsv",
            "--output",
            "legacy",
            "--threads",
            "1",
        ],
    );
    let failed = Command::new(env!("CARGO_BIN_EXE_expresso"))
        .current_dir(dir)
        .args([
            "export",
            "--input",
            "legacy",
            "--output",
            "unavailable",
            "--normalization",
            "coverage",
        ])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("rebuild"));
    assert!(!dir.join("unavailable").exists());
    success(
        dir,
        &[
            "export",
            "--input",
            "legacy",
            "--output",
            "legacy-cpm",
            "--normalization",
            "cpm",
        ],
    );
}
