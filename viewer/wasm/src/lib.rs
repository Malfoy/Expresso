//! Small raw WebAssembly ABI; no JavaScript build tool or bindings generator.
#[path = "../../../src/eab.rs"]
mod eab;
use anyhow::{Context, Result, ensure};
use serde::Deserialize;
use std::io::Read;

#[derive(Deserialize)]
struct Request {
    targets: usize,
    reference: String,
    selected: Vec<usize>,
}

fn decode(request: &[u8], vector: &[u8]) -> Result<Vec<String>> {
    let request: Request = serde_json::from_slice(request)?;
    ensure!(
        request.targets <= 100_000_000,
        "Reference exceeds browser decoder limit"
    );
    ensure!(
        request.reference.len() == 64 && request.reference.is_ascii(),
        "Invalid reference hash"
    );
    let mut reference = [0; 32];
    for (i, byte) in reference.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&request.reference[i * 2..i * 2 + 2], 16)?;
    }
    ensure!(
        request.selected.iter().all(|i| *i < request.targets),
        "Target ID out of range"
    );
    let limit = (request.targets as u64)
        .checked_mul(2)
        .and_then(|n| n.checked_add(524_364))
        .context("Vector too large")?;
    let decoder = ruzstd::decoding::StreamingDecoder::new(vector)?;
    let mut raw = Vec::new();
    decoder.take(limit + 1).read_to_end(&mut raw)?;
    ensure!(raw.len() as u64 <= limit, "Decompressed vector too large");
    let values = eab::decode(raw.as_slice(), request.targets, &reference)?;
    // Decimal strings preserve the entire u64 range in JavaScript and CSV.
    Ok(request
        .selected
        .iter()
        .map(|i| values[*i].to_string())
        .collect())
}

#[unsafe(no_mangle)]
pub extern "C" fn allocate(length: usize) -> *mut u8 {
    Box::into_raw(vec![0; length].into_boxed_slice()) as *mut u8
}

/// # Safety
/// The pointer/length must denote a live allocation returned by allocate or decode_vector.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn release(pointer: *mut u8, length: usize) {
    unsafe {
        drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
            pointer, length,
        )))
    };
}

/// # Safety
/// Both input ranges must point to live, initialized allocations in this WASM memory.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn decode_vector(
    request_ptr: *const u8,
    request_len: usize,
    vector_ptr: *const u8,
    vector_len: usize,
) -> u64 {
    let request = unsafe { std::slice::from_raw_parts(request_ptr, request_len) };
    let vector = unsafe { std::slice::from_raw_parts(vector_ptr, vector_len) };
    let result = match decode(request, vector) {
        Ok(values) => serde_json::json!({"values":values}),
        Err(error) => serde_json::json!({"error":format!("{error:#}")}),
    };
    let bytes = serde_json::to_vec(&result)
        .expect("JSON serialization")
        .into_boxed_slice();
    let length = bytes.len() as u64;
    let pointer = Box::into_raw(bytes) as *mut u8 as usize as u64;
    pointer | (length << 32)
}

#[cfg(test)]
mod tests {
    #[test]
    fn selected_counts_preserve_u64_and_non_byte_aligned_codes() {
        let request=br#"{"targets":5,"reference":"e2eddccdd8c0b2effc375a239e518c0717246d9200a7407a836b253718047501","selected":[3,0,1,3,4]}"#;
        let vector = include_bytes!("../tests/data/full-u64.eab.zst");
        assert_eq!(
            super::decode(request, vector).unwrap(),
            [
                "18446744073709551615",
                "0",
                "1",
                "18446744073709551615",
                "1"
            ]
        );
        let wrong=br#"{"targets":5,"reference":"0000000000000000000000000000000000000000000000000000000000000000","selected":[0]}"#;
        assert!(super::decode(wrong, vector).is_err());
        assert!(super::decode(request, &vector[..vector.len() - 1]).is_err());
    }
    #[test]
    fn invalid_requests_are_errors() {
        assert!(super::decode(b"{}", b"bad").is_err());
        assert!(
            super::decode(br#"{"targets":1,"reference":"bad","selected":[0]}"#, b"bad").is_err()
        );
    }
}
