//! Shared native/WebAssembly EAB validation and decoding.
use anyhow::{Context, Result, ensure};
use std::io::Read;
const MAGIC: &[u8; 8] = b"EXPRAB01";
const HEADER_SIZE: usize = 72;
fn packed_len(n: usize, bits: u8) -> Result<usize> {
    n.checked_mul(bits as usize)
        .and_then(|n| n.checked_add(7))
        .map(|n| n / 8)
        .context("abundance vector size overflow")
}

pub struct PackedVector {
    bits: u8,
    targets: usize,
    packed: Vec<u8>,
    levels: Vec<u64>,
}
impl PackedVector {
    pub fn read(mut input: impl Read, targets: usize, reference: &[u8; 32]) -> Result<Self> {
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
            levels[0] == 0
                && levels.last() == Some(&maximum)
                && levels.windows(2).all(|v| v[0] < v[1]),
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
        let vector = Self {
            bits,
            targets,
            packed,
            levels,
        };
        // Unselected codes must be valid as well. Full codebooks need no extra scan.
        if vector.levels.len() != 1usize << bits {
            for i in 0..targets {
                vector.get(i)?;
            }
        }
        Ok(vector)
    }

    pub fn get(&self, i: usize) -> Result<u64> {
        ensure!(i < self.targets, "Target ID out of range");
        let offset = i * self.bits as usize;
        let mut word = 0u32;
        for b in 0..(offset % 8 + self.bits as usize).div_ceil(8) {
            word |= u32::from(self.packed[offset / 8 + b]) << (8 * b);
        }
        let code = (word >> (offset % 8)) as usize & ((1usize << self.bits) - 1);
        self.levels
            .get(code)
            .copied()
            .context("EAB code outside codebook")
    }
}

pub fn decode(input: impl Read, targets: usize, reference: &[u8; 32]) -> Result<Vec<u64>> {
    let vector = PackedVector::read(input, targets, reference)?;
    (0..targets).map(|i| vector.get(i)).collect()
}
