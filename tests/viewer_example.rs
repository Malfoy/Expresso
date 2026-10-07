//! Keep the ready-to-open viewer example consistent with its FASTA inputs and
//! the expected counts printed in the README.
use std::{collections::BTreeMap, fs, path::Path, process::Command};

fn success(args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_expresso"))
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

fn check_export(directory: &Path, expected: &Path) {
    let mut reader = csv::Reader::from_path(expected).unwrap();
    let headers = reader.headers().unwrap().clone();
    let rows: Vec<_> = reader.records().map(Result::unwrap).collect();
    for (column, name) in headers.iter().enumerate().skip(1) {
        let path = if name == "global" {
            directory.join("global.csv")
        } else {
            directory.join("datasets").join(format!("{name}.csv"))
        };
        let mut actual = csv::Reader::from_path(path).unwrap();
        let counts: BTreeMap<_, _> = actual
            .records()
            .map(|row| {
                let row = row.unwrap();
                (row[1].to_string(), row[2].parse::<u64>().unwrap())
            })
            .collect();
        assert_eq!(counts.len(), 6);
        for row in &rows {
            assert_eq!(counts[&row[0]], row[column].parse::<u64>().unwrap());
        }
    }
}

#[test]
fn bundled_and_regenerated_viewer_examples_match_documented_counts() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let example = root.join("examples/viewer");
    let tmp = tempfile::tempdir().unwrap();
    let index = tmp.path().join("index");
    let results = tmp.path().join("results");
    success(&[
        "run",
        "--level",
        "gene",
        "--genes",
        example.join("genes.fa").to_str().unwrap(),
        "--fof",
        example.join("datasets.tsv").to_str().unwrap(),
        "--index",
        index.to_str().unwrap(),
        "--output",
        results.to_str().unwrap(),
        "--k",
        "7",
        "--threads",
        "4",
        "--jobs",
        "2",
        "--bits",
        "16",
        "--stats",
        "--viewer-index",
    ]);
    let metadata: serde_json::Value =
        serde_json::from_slice(&fs::read(index.join("metadata.json")).unwrap()).unwrap();
    assert_eq!(metadata["indexed_kmers"], 204);
    assert_eq!(metadata["duplicate_kmers"], 0);
    for target in metadata["exons"].as_array().unwrap() {
        assert_eq!(target["unique_kmers"], 34);
    }
    for (label, input) in [
        ("bundled", example.join("toy.eai")),
        ("regenerated", results.join("abundance.eai")),
    ] {
        let destination = tmp.path().join(label);
        success(&[
            "export",
            "--input",
            input.to_str().unwrap(),
            "--output",
            destination.to_str().unwrap(),
            "--compression",
            "none",
            "--with-names",
        ]);
        check_export(&destination, &example.join("expected-counts.csv"));
    }
}
