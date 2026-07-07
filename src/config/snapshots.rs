use serde::{Deserialize, Serialize};

use std::fmt;

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct Snapshots {
    #[serde(default)]
    pub equilibration_interval: usize,
    #[serde(default)]
    pub measurement_interval: usize,
    #[serde(default)]
    pub compression_level: usize,
    #[serde(default = "default_save_dir")]
    pub save_directory: String,
}

fn default_save_dir() -> String {
    "snapshots".to_string()
}

impl Snapshots {
    pub fn validate(&self) -> anyhow::Result<()> {
        Ok(())
    }
}

impl fmt::Display for Snapshots {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\nSnapshots: Enable")?;
        writeln!(
            f,
            "  Equilibration Interval: {} steps",
            self.equilibration_interval
        )?;
        writeln!(
            f,
            "  Sampling Interval: {} steps",
            self.measurement_interval
        )?;
        writeln!(f, "  Compression Level: {}", self.compression_level)?;
        writeln!(f, "  Snapshot Directory: {}", self.save_directory)?;

        Ok(())
    }
}

#[cfg(test)]
#[path = "snapshots_tests.rs"]
mod tests;
