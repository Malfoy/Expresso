//! View/export transformations. Stored abundance vectors always retain raw counts.
use anyhow::{Context, Result, ensure};
use serde_json::Value;

#[derive(Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Mode {
    #[default]
    None,
    /// Observations per million weighted valid k-mer occurrences in the input.
    Cpm,
    /// Observations per valid exonic reference k-mer position.
    Coverage,
    /// CPM divided by the target's valid exonic k-mer length.
    LengthCpm,
}

impl Mode {
    pub fn dataset_total(self, entry: &Value) -> Result<Option<f64>> {
        if !matches!(self, Self::Cpm | Self::LengthCpm) {
            return Ok(None);
        }
        let total = &entry["total_kmers"];
        let numerator: u128 = total["numerator"]
            .as_str()
            .context("Dataset k-mer total unavailable; requantify with this EXPRESSO version")?
            .parse()?;
        let scale = total["scale"]
            .as_u64()
            .context("Invalid dataset k-mer total scale")?;
        ensure!(
            [1, crate::abundance::SCALE].contains(&scale),
            "Invalid dataset k-mer total scale"
        );
        Ok(Some(numerator as f64 / scale as f64))
    }

    pub fn apply(self, count: u64, length: Option<u64>, total: Option<f64>) -> Result<f64> {
        let mut value = count as f64;
        if matches!(self, Self::Coverage | Self::LengthCpm) {
            let length = length.context(
                "Target k-mer lengths unavailable; rebuild the reference index and requantify",
            )?;
            ensure!(
                length > 0 || count == 0,
                "Nonzero abundance has zero exonic k-mer length"
            );
            value = if length == 0 {
                0.0
            } else {
                value / length as f64
            };
        }
        if matches!(self, Self::Cpm | Self::LengthCpm) {
            let total = total.context("Dataset k-mer total unavailable")?;
            ensure!(
                total.is_finite() && total >= 0.0 && (total > 0.0 || count == 0),
                "Invalid or zero dataset k-mer total for nonzero abundance"
            );
            value = if total == 0.0 {
                0.0
            } else {
                value / total * 1_000_000.0
            };
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fractional_rates_and_missing_or_zero_denominators() {
        let entry = serde_json::json!({"total_kmers":{"numerator":"4000001","scale":1_000_000}});
        let total = Mode::LengthCpm.dataset_total(&entry).unwrap();
        assert_eq!(total, Some(4.000001));
        assert_eq!(Mode::Coverage.apply(1, Some(10), None).unwrap(), 0.1);
        assert!(
            (Mode::LengthCpm.apply(1, Some(10), total).unwrap() - 24999.993750001562).abs() < 1e-9
        );
        assert_eq!(Mode::LengthCpm.apply(0, Some(0), Some(0.0)).unwrap(), 0.0);
        assert!(Mode::Coverage.apply(1, Some(0), None).is_err());
        assert!(Mode::Coverage.apply(0, None, None).is_err());
        assert!(Mode::Cpm.dataset_total(&serde_json::json!({})).is_err());
        assert!(Mode::Cpm.apply(1, None, Some(0.0)).is_err());
        let large =
            serde_json::json!({"total_kmers":{"numerator":u128::MAX.to_string(),"scale":1}});
        assert!(
            Mode::Cpm
                .dataset_total(&large)
                .unwrap()
                .unwrap()
                .is_finite()
        );
    }
}
