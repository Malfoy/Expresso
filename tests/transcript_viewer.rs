use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{fs, io::Read, path::Path, process::Command};

fn pack(root: &Path, destination: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_expresso"))
        .current_dir(root)
        .args(["pack", "--input", "results", "--output", destination])
        .output()
        .unwrap()
}
#[test]
fn transcript_bundle_preserves_both_metrics_missing_values_and_source_coverage() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    fs::create_dir(root.join("results")).unwrap();
    let csv = "transcript_id,transcript_name,gene_name,observable_kmers,abundance,assigned_observations\n1,T1,G,10,1.25,12.5\n2,T2,G,20,,\n3,T3,H,5,0,0\n";
    fs::write(root.join("results/sample.csv"), csv).unwrap();
    fs::write(root.join("results/global.csv"), csv).unwrap();
    let manifest = json!({"version":1,"format":"isoform-abundance","model_sha256":"a".repeat(64),
        "datasets":[{"name":"global","output":"sample.csv"}],"global":"global.csv",
        "source_coverage":{"coverage":{"complete":false,"failed_datasets":1},"failed_datasets":[{"name":"bad"}]}});
    fs::write(
        root.join("results/manifest.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    let result = pack(root, "transcripts.eai");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = fs::read(root.join("transcripts.eai")).unwrap();
    assert_eq!(&bytes[..8], b"EXPREAI2");
    let offset = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let length = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
    let raw = &bytes[offset..offset + length];
    assert_eq!(Sha256::digest(raw).as_slice(), &bytes[offset + length..]);
    let catalogue: Value = serde_json::from_slice(raw).unwrap();
    assert_eq!(catalogue["level"], "transcript");
    assert_eq!(catalogue["coverage"]["complete"], false);
    assert_eq!(catalogue["datasets"][0]["name"], "global");
    assert_eq!(catalogue["datasets"][1]["kind"], "global");
    let reference = Sha256::digest(serde_json::to_vec(&catalogue["targets"]).unwrap());
    for dataset in catalogue["datasets"].as_array().unwrap() {
        for (metric, expected) in [
            ("abundance", [Some(1.25), None, Some(0.0)]),
            ("assigned_observations", [Some(12.5), None, Some(0.0)]),
        ] {
            let block = &dataset["metrics"][metric];
            let start = block["offset"].as_u64().unwrap() as usize;
            let n = block["length"].as_u64().unwrap() as usize;
            let mut raw = Vec::new();
            zstd::stream::read::Decoder::new(&bytes[start..start + n])
                .unwrap()
                .read_to_end(&mut raw)
                .unwrap();
            assert_eq!(&raw[..8], b"EXPREFV1");
            assert_eq!(&raw[8..40], reference.as_slice());
            assert_eq!(u64::from_le_bytes(raw[40..48].try_into().unwrap()), 3);
            let end = raw.len() - 4;
            assert_eq!(
                crc32fast::hash(&raw[..end]),
                u32::from_le_bytes(raw[end..].try_into().unwrap())
            );
            for (i, value) in expected.into_iter().enumerate() {
                let bits = u64::from_le_bytes(raw[48 + i * 8..56 + i * 8].try_into().unwrap());
                if let Some(value) = value {
                    assert_eq!(f64::from_bits(bits), value);
                } else {
                    assert_eq!(bits, 0x7ff8000000000000);
                }
            }
        }
    }
    // Packaging must reject mismatched tables and invalid values without publishing.
    fs::write(root.join("results/global.csv"), csv.replace("T1", "WRONG")).unwrap();
    assert!(!pack(root, "bad-table.eai").status.success());
    assert!(!root.join("bad-table.eai").exists());
    fs::write(root.join("results/global.csv"), csv.replace("1.25", "NaN")).unwrap();
    assert!(!pack(root, "bad-value.eai").status.success());
    assert!(!root.join("bad-value.eai").exists());
}
