//! EFV v1: lossless f64 observations/rates with an explicit missing-value sentinel.
//! Shared by native packaging and the WebAssembly decoder.
use anyhow::{Context, Result, ensure};
use std::io::{Read, Write};

const MAGIC: &[u8; 8] = b"EXPREFV1";
const MISSING: u64 = 0x7ff8_0000_0000_0000;

pub struct FloatVector {
    data: Vec<u8>,
}
impl FloatVector {
    pub fn read(mut input: impl Read, targets: usize, reference: &[u8; 32]) -> Result<Self> {
        let length = targets
            .checked_mul(8)
            .and_then(|n| n.checked_add(52))
            .context("vector size overflow")?;
        let mut data = Vec::new();
        input
            .by_ref()
            .take(length as u64 + 1)
            .read_to_end(&mut data)?;
        ensure!(data.len() == length, "invalid floating vector length");
        ensure!(
            &data[..8] == MAGIC && &data[8..40] == reference,
            "floating vector reference mismatch"
        );
        ensure!(
            u64::from_le_bytes(data[40..48].try_into()?) == targets as u64,
            "floating vector dimension mismatch"
        );
        let end = length - 4;
        ensure!(
            crc32fast::hash(&data[..end]) == u32::from_le_bytes(data[end..].try_into()?),
            "floating vector checksum mismatch"
        );
        let vector = Self { data };
        for i in 0..targets {
            vector.get(i)?;
        }
        Ok(vector)
    }
    pub fn get(&self, id: usize) -> Result<Option<f64>> {
        let targets = (self.data.len() - 52) / 8;
        ensure!(id < targets, "target ID out of range");
        let offset = 48 + id * 8;
        let bits = u64::from_le_bytes(self.data[offset..offset + 8].try_into()?);
        if bits == MISSING {
            return Ok(None);
        }
        let value = f64::from_bits(bits);
        ensure!(
            value.is_finite() && value >= 0.0,
            "invalid floating abundance"
        );
        Ok(Some(value))
    }
}

pub fn write(mut out: impl Write, values: &[Option<f64>], reference: &[u8; 32]) -> Result<()> {
    let mut crc = crc32fast::Hasher::new();
    let mut emit = |bytes: &[u8]| -> Result<()> {
        crc.update(bytes);
        out.write_all(bytes)?;
        Ok(())
    };
    emit(MAGIC)?;
    emit(reference)?;
    emit(&(values.len() as u64).to_le_bytes())?;
    for &value in values {
        let bits = if let Some(v) = value {
            ensure!(v.is_finite() && v >= 0.0, "invalid floating abundance");
            v.to_bits()
        } else {
            MISSING
        };
        emit(&bits.to_le_bytes())?;
    }
    out.write_all(&crc.finalize().to_le_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fractions_missing_and_integrity() {
        let values = [Some(0.0), Some(1.25), None, Some(f64::MAX)];
        let mut raw = Vec::new();
        write(&mut raw, &values, &[9; 32]).unwrap();
        let vector = FloatVector::read(raw.as_slice(), values.len(), &[9; 32]).unwrap();
        for (i, v) in values.into_iter().enumerate() {
            assert_eq!(vector.get(i).unwrap(), v);
        }
        assert!(vector.get(4).is_err());
        assert!(FloatVector::read(raw.as_slice(), 4, &[8; 32]).is_err());
        raw[60] ^= 1;
        assert!(FloatVector::read(raw.as_slice(), 4, &[9; 32]).is_err());
        assert!(write(Vec::new(), &[Some(f64::NAN)], &[9; 32]).is_err());
    }
}
