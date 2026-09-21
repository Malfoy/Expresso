use anyhow::{Context, Result, bail, ensure};
use clap::ValueEnum;
use serde::{Deserialize, Serialize};

// Fixed-point accumulation makes results independent of worker scheduling.
// Read mode uses scale=1, preserving the complete u64 integer count range.
pub const SCALE: u64 = 1_000_000;

#[derive(Clone, Copy, Debug, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    Reads,
    Unitigs,
    Auto,
}

impl Mode {
    pub fn scale(self) -> u64 {
        if matches!(self, Self::Reads) {
            1
        } else {
            SCALE
        }
    }
}

pub fn round_ratio(n: u128, d: u128) -> Result<u64> {
    ensure!(d != 0, "zero denominator");
    let rounded = n / d + u128::from(n % d >= d.div_ceil(2));
    u64::try_from(rounded).context("abundance exceeds u64")
}

/// Parse a nonnegative decimal into millionths without floating point.
/// More than six fractional digits are rounded half up once at ingestion.
pub fn decimal(s: &str) -> Result<u64> {
    let s = s.strip_prefix('+').unwrap_or(s);
    let (whole, fraction) = s.split_once('.').unwrap_or((s, ""));
    ensure!(!whole.is_empty() || !fraction.is_empty(), "empty abundance");
    ensure!(
        whole
            .bytes()
            .chain(fraction.bytes())
            .all(|b| b.is_ascii_digit()),
        "invalid nonnegative decimal abundance: {s}"
    );
    let whole: u64 = if whole.is_empty() { 0 } else { whole.parse()? };
    let mut frac = 0u64;
    for i in 0..6 {
        frac = frac * 10 + fraction.as_bytes().get(i).map_or(0, |b| (b - b'0') as u64);
    }
    frac += u64::from(fraction.as_bytes().get(6).is_some_and(|b| *b >= b'5'));
    whole
        .checked_mul(SCALE)
        .and_then(|n| n.checked_add(frac))
        .context("header abundance exceeds fixed-point range")
}

pub fn weight(header: &[u8], len: usize, assembly_k: usize, mode: Mode) -> Result<u64> {
    if matches!(mode, Mode::Reads) {
        return Ok(1);
    }
    let header = std::str::from_utf8(header).context("header is not UTF-8")?;
    let mut mean = None;
    let mut total = None;
    // Some Logan records contain only an abundance tag, with no identifier.
    for token in header.split_ascii_whitespace() {
        if let Some(v) = ["ka:f:", "km:f:", "ka:i:", "km:i:"]
            .iter()
            .find_map(|p| token.strip_prefix(p))
        {
            let v = decimal(v).with_context(|| format!("invalid abundance in {header}"))?;
            ensure!(
                mean.is_none_or(|old| old == v),
                "conflicting mean abundance tags in {header}"
            );
            mean = Some(v);
        } else if let Some(v) = token.strip_prefix("KC:i:") {
            let v: u64 = v
                .parse()
                .with_context(|| format!("invalid KC:i in {header}"))?;
            ensure!(total.replace(v).is_none(), "duplicate KC:i in {header}");
        }
    }
    if let Some(mean) = mean {
        return Ok(mean);
    }
    if let Some(total) = total {
        ensure!(
            assembly_k > 0 && len >= assembly_k,
            "KC:i requires sequence length >= assembly k ({assembly_k}): {header}"
        );
        return round_ratio(
            total as u128 * SCALE as u128,
            (len - assembly_k + 1) as u128,
        );
    }
    if matches!(mode, Mode::Unitigs) {
        bail!("unitig has no ka:f:, km:f:, or KC:i: abundance: {header}");
    }
    Ok(SCALE)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn headers_and_rounding() {
        for mode in [Mode::Unitigs, Mode::Auto] {
            assert_eq!(weight(b"ka:f:195.0", 34, 31, mode).unwrap(), 195_000_000);
            assert_eq!(weight(b"km:f:3.5", 34, 31, mode).unwrap(), 3_500_000);
            assert_eq!(weight(b"KC:i:14", 34, 31, mode).unwrap(), 3_500_000);
        }
        assert_eq!(
            weight(b"u ka:f:1.7", 34, 31, Mode::Unitigs).unwrap(),
            1_700_000
        );
        assert_eq!(
            weight(b"u KC:i:14", 34, 31, Mode::Unitigs).unwrap(),
            3_500_000
        );
        assert_eq!(
            weight(b"u KC:i:14 km:f:3.5", 34, 31, Mode::Auto).unwrap(),
            3_500_000
        );
        assert_eq!(weight(b"r", 34, 31, Mode::Reads).unwrap(), 1);
        assert!(weight(b"u", 34, 31, Mode::Unitigs).is_err());
        for bad in ["-1", "NaN", "inf", "1.2.3", ""] {
            assert!(decimal(bad).is_err());
        }
        assert_eq!(decimal("0.0000005").unwrap(), 1);
        assert_eq!(round_ratio(35, 10).unwrap(), 4);
        assert_eq!(round_ratio(u64::MAX as u128, 1).unwrap(), u64::MAX);
    }
}
