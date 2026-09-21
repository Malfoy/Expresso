//! EAB v1: bit-packed log-quantized abundance vectors with a shared reference.
use crate::{index::Exon, output};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{io::Read, path::Path};

pub const MAGIC: &[u8; 8] = b"EXPRAB01";
const HEADER_SIZE: usize = 72;

pub fn reference_id(exons: &[Exon]) -> [u8; 32] {
    let mut hash = Sha256::new();
    hash.update(b"EXPRESSO-reference-v1\0");
    hash.update((exons.len() as u64).to_le_bytes());
    for exon in exons {
        hash.update((exon.name.len() as u64).to_le_bytes());
        hash.update(exon.name.as_bytes());
        hash.update((exon.length as u64).to_le_bytes());
        hash.update(exon.unique_kmers.to_le_bytes());
    }
    hash.finalize().into()
}

pub fn hex(id: &[u8; 32]) -> String {
    id.iter().map(|x| format!("{x:02x}")).collect()
}

pub struct Quantizer {
    pub levels: Vec<u64>,
    pub base: f64,
    pub exact_prefix: u64,
}

impl Quantizer {
    pub fn new(bits: u8, maximum: u64) -> Result<Self> {
        ensure!((2..=16).contains(&bits), "abundance bits must be in 2..=16");
        let last = (1u64 << bits) - 1;
        if maximum <= last {
            return Ok(Self {
                levels: (0..=maximum).collect(),
                base: 1.0,
                exact_prefix: maximum,
            });
        }
        // Preserve small counts before spending the remaining codes on a
        // geometric grid. Clamp rounded representatives to distinct integers;
        // leave enough room for every remaining code up to the exact maximum.
        let prefix = 15.min(last / 4);
        let first_log = prefix + 1;
        let steps = last - first_log;
        let log_base = ((maximum as f64).ln() - (first_log as f64).ln()) / steps as f64;
        let mut levels: Vec<u64> = (0..=prefix).collect();
        for code in first_log..=last {
            let value = if code == last {
                maximum
            } else {
                ((first_log as f64).ln() + (code - first_log) as f64 * log_base)
                    .exp()
                    .round() as u64
            };
            let lower = levels.last().unwrap() + 1;
            let upper = maximum - (last - code);
            levels.push(value.clamp(lower, upper));
        }
        Ok(Self {
            levels,
            base: log_base.exp(),
            exact_prefix: prefix,
        })
    }

    pub fn encode(&self, value: u64) -> usize {
        if value <= self.exact_prefix {
            return value as usize;
        }
        let hi = self.levels.partition_point(|&v| v < value);
        if hi == 0 {
            return 0;
        }
        if hi == self.levels.len() {
            return hi - 1;
        }
        if value - self.levels[hi - 1] <= self.levels[hi] - value {
            hi - 1
        } else {
            hi
        }
    }
}

#[derive(Debug, Serialize)]
pub struct VectorInfo {
    pub bits: u8,
    pub targets: usize,
    pub maximum: u64,
    pub codebook_entries: usize,
    pub logarithmic_base: f64,
    pub exact_prefix: u64,
    pub exact_values: usize,
    pub nonzero_values: usize,
    pub maximum_absolute_error: u64,
    pub mean_absolute_error: f64,
    pub maximum_relative_error_nonzero: f64,
    pub mean_relative_error_nonzero: f64,
    pub packed_bytes: usize,
    pub uncompressed_file_bytes: usize,
}

fn packed_len(n: usize, bits: u8) -> Result<usize> {
    n.checked_mul(bits as usize)
        .and_then(|n| n.checked_add(7))
        .map(|n| n / 8)
        .context("abundance vector size overflow")
}

pub fn write_vector(
    path: &Path,
    compression: output::Compression,
    values: &[u64],
    bits: u8,
    reference: &[u8; 32],
) -> Result<VectorInfo> {
    let maximum = values.iter().copied().max().unwrap_or(0);
    let quantizer = Quantizer::new(bits, maximum)?;
    let mut packed = vec![0u8; packed_len(values.len(), bits)?];
    let mut info = VectorInfo {
        bits,
        targets: values.len(),
        maximum,
        codebook_entries: quantizer.levels.len(),
        logarithmic_base: quantizer.base,
        exact_prefix: quantizer.exact_prefix,
        exact_values: 0,
        nonzero_values: 0,
        maximum_absolute_error: 0,
        mean_absolute_error: 0.0,
        maximum_relative_error_nonzero: 0.0,
        mean_relative_error_nonzero: 0.0,
        packed_bytes: packed.len(),
        uncompressed_file_bytes: HEADER_SIZE + quantizer.levels.len() * 8 + packed.len() + 4,
    };
    for (i, &value) in values.iter().enumerate() {
        let code = quantizer.encode(value);
        let offset = i * bits as usize;
        let word = (code as u32) << (offset % 8);
        let nbytes = (offset % 8 + bits as usize).div_ceil(8);
        for b in 0..nbytes {
            packed[offset / 8 + b] |= (word >> (8 * b)) as u8;
        }
        let decoded = quantizer.levels[code];
        let error = value.abs_diff(decoded);
        info.exact_values += usize::from(error == 0);
        info.maximum_absolute_error = info.maximum_absolute_error.max(error);
        info.mean_absolute_error += error as f64;
        if value > 0 {
            info.nonzero_values += 1;
            let relative = error as f64 / value as f64;
            info.mean_relative_error_nonzero += relative;
            info.maximum_relative_error_nonzero = info.maximum_relative_error_nonzero.max(relative);
        }
    }
    info.mean_absolute_error /= values.len().max(1) as f64;
    info.mean_relative_error_nonzero /= info.nonzero_values.max(1) as f64;
    let mut header = Vec::with_capacity(HEADER_SIZE);
    header.extend_from_slice(MAGIC);
    header.extend_from_slice(&[bits, 0, 0, 0]); // bit width, encoding=0, reserved=0
    header.extend_from_slice(&(values.len() as u64).to_le_bytes());
    header.extend_from_slice(&(quantizer.levels.len() as u32).to_le_bytes());
    header.extend_from_slice(&(packed.len() as u64).to_le_bytes());
    header.extend_from_slice(&maximum.to_le_bytes());
    header.extend_from_slice(reference);
    debug_assert_eq!(header.len(), HEADER_SIZE);
    output::encoded(path, compression, |out| {
        let mut crc = crc32fast::Hasher::new();
        out.write_all(&header)?;
        crc.update(&header);
        for level in &quantizer.levels {
            let bytes = level.to_le_bytes();
            out.write_all(&bytes)?;
            crc.update(&bytes);
        }
        out.write_all(&packed)?;
        crc.update(&packed);
        out.write_all(&crc.finalize().to_le_bytes())?;
        out.flush()?;
        Ok(())
    })?;
    Ok(info)
}

pub fn read_vector(path: &Path, targets: usize, reference: &[u8; 32]) -> Result<Vec<u64>> {
    let mut input = output::reader(path)?;
    let mut header = [0u8; HEADER_SIZE];
    input.read_exact(&mut header)?;
    ensure!(&header[..8] == MAGIC, "Invalid EAB magic/version");
    let bits = header[8];
    ensure!(
        (2..=16).contains(&bits) && header[9..12] == [0, 0, 0],
        "Unsupported EAB encoding"
    );
    let n = u64::from_le_bytes(header[12..20].try_into()?);
    let entries = u32::from_le_bytes(header[20..24].try_into()?) as usize;
    let bytes = u64::from_le_bytes(header[24..32].try_into()?);
    let maximum = u64::from_le_bytes(header[32..40].try_into()?);
    ensure!(
        n == targets as u64 && &header[40..72] == reference,
        "EAB reference mismatch"
    );
    ensure!(
        (1..=1usize << bits).contains(&entries),
        "Invalid EAB codebook size"
    );
    let packed_size = packed_len(targets, bits)?;
    ensure!(bytes == packed_size as u64, "Invalid EAB payload length");
    let mut table = vec![0u8; entries * 8];
    input.read_exact(&mut table)?;
    let levels: Vec<u64> = table
        .chunks_exact(8)
        .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
        .collect();
    ensure!(
        levels[0] == 0 && levels.last() == Some(&maximum) && levels.windows(2).all(|v| v[0] < v[1]),
        "Invalid EAB codebook"
    );
    let mut packed = vec![0u8; packed_size];
    input.read_exact(&mut packed)?;
    let mut checksum = [0u8; 4];
    input.read_exact(&mut checksum)?;
    let mut crc = crc32fast::Hasher::new();
    crc.update(&header);
    crc.update(&table);
    crc.update(&packed);
    ensure!(
        crc.finalize() == u32::from_le_bytes(checksum),
        "EAB checksum mismatch"
    );
    ensure!(
        input.read(&mut [0u8; 1])? == 0,
        "Trailing data after EAB vector"
    );
    let used_bits = targets * bits as usize % 8;
    if used_bits != 0 {
        ensure!(
            packed.last().unwrap() >> used_bits == 0,
            "Nonzero EAB padding"
        );
    }
    let mut values = Vec::with_capacity(targets);
    for i in 0..targets {
        let offset = i * bits as usize;
        let mut word = 0u32;
        for b in 0..(offset % 8 + bits as usize).div_ceil(8) {
            word |= u32::from(packed[offset / 8 + b]) << (8 * b);
        }
        let code = (word >> (offset % 8)) as usize & ((1usize << bits) - 1);
        values.push(*levels.get(code).context("EAB code outside codebook")?);
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roundtrip_all_widths_codecs_and_full_u64_range() {
        let dir = tempfile::tempdir().unwrap();
        let id = [17; 32];
        let values = [
            0,
            1,
            2,
            3,
            4,
            7,
            15,
            16,
            17,
            100,
            255,
            256,
            65535,
            1000000,
            u64::MAX,
        ];
        for bits in 2..=16 {
            for codec in [
                output::Compression::None,
                output::Compression::Gz,
                output::Compression::Xz,
                output::Compression::Zstd,
            ] {
                let path = dir.path().join("vector");
                let info = write_vector(&path, codec, &values, bits, &id).unwrap();
                let read = read_vector(&path, values.len(), &id).unwrap();
                assert_eq!(&read[..2], &[0, 1]);
                assert_eq!(read.last(), Some(&u64::MAX));
                assert!(read.windows(2).all(|w| w[0] <= w[1]));
                assert!(read.iter().skip(1).all(|v| *v > 0));
                let measured = values
                    .iter()
                    .zip(&read)
                    .map(|(x, y)| x.abs_diff(*y))
                    .max()
                    .unwrap();
                assert_eq!(info.maximum_absolute_error, measured);
                assert_eq!(
                    info.packed_bytes,
                    (values.len() * bits as usize).div_ceil(8)
                );
                assert!(read_vector(&path, values.len(), &[18; 32]).is_err());
            }
        }
    }
    #[test]
    fn exact_small_counts_and_corruption_detection() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("vector");
        let id = [0; 32];
        write_vector(
            &path,
            output::Compression::None,
            &(0..8).collect::<Vec<_>>(),
            3,
            &id,
        )
        .unwrap();
        let golden = std::fs::read(&path).unwrap();
        assert_eq!(&golden[72 + 8 * 8..72 + 8 * 8 + 3], &[0x88, 0xc6, 0xfa]);
        for values in [vec![], vec![0; 17], (0..=255).collect()] {
            write_vector(&path, output::Compression::None, &values, 8, &id).unwrap();
            assert_eq!(read_vector(&path, values.len(), &id).unwrap(), values);
        }
        let values: Vec<u64> = (0..10000).collect();
        let q = Quantizer::new(8, 9999).unwrap();
        for x in 0..=15 {
            assert_eq!(q.levels[q.encode(x)], x);
        }
        for x in values.iter().copied() {
            let decoded = q.levels[q.encode(x)];
            assert_eq!(
                decoded.abs_diff(x),
                q.levels.iter().map(|y| y.abs_diff(x)).min().unwrap()
            );
        }
        write_vector(&path, output::Compression::None, &values, 7, &id).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        for offset in [
            0,
            8,
            12,
            20,
            24,
            32,
            40,
            72,
            bytes.len() - 5,
            bytes.len() - 1,
        ] {
            let mut bad = bytes.clone();
            bad[offset] ^= 1;
            std::fs::write(&path, bad).unwrap();
            assert!(
                read_vector(&path, values.len(), &id).is_err(),
                "offset {offset}"
            );
        }
        std::fs::write(&path, &bytes[..bytes.len() - 1]).unwrap();
        assert!(read_vector(&path, values.len(), &id).is_err());
        let mut extra = bytes;
        extra.push(0);
        std::fs::write(&path, extra).unwrap();
        assert!(read_vector(&path, values.len(), &id).is_err());
    }
}
