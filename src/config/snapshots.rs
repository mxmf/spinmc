use serde::{Deserialize, Serialize};

use std::fmt;
use std::fs::File;

use ndarray::Axis;
use ndarray_npy::NpzWriter;
use zip::CompressionMethod;
use zip::write::SimpleFileOptions;

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
        if self.compression_level > 9 {
            anyhow::bail!(
                "Snapshot compression_level must be between 0 and 9, got {}",
                self.compression_level
            );
        }
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

pub fn save_snapshots_to_npz(
    filename: &str,
    equil_data: &[ndarray::Array5<f64>],
    steps_data: &[ndarray::Array5<f64>],
    compression_level: usize,
) -> anyhow::Result<()> {
    fn stack_snapshots(
        arrays: &[ndarray::Array5<f64>],
        fallback_shape: Option<[usize; 5]>,
    ) -> anyhow::Result<ndarray::ArrayD<f64>> {
        if arrays.is_empty() {
            let [sub, x, y, z, components] = fallback_shape.unwrap_or([0, 0, 0, 0, 3]);
            return Ok(ndarray::ArrayD::zeros(vec![0, sub, x, y, z, components]));
        }
        let views: Vec<_> = arrays.iter().map(|a| a.view()).collect();
        Ok(ndarray::stack(Axis(0), &views)?.into_dyn())
    }

    let snapshot_shape = equil_data
        .first()
        .or_else(|| steps_data.first())
        .map(|array| {
            let shape = array.shape();
            [shape[0], shape[1], shape[2], shape[3], shape[4]]
        });

    let equil_stacked = stack_snapshots(equil_data, snapshot_shape)?;
    let steps_stacked = stack_snapshots(steps_data, snapshot_shape)?;

    let options = if compression_level == 0 {
        SimpleFileOptions::default().compression_method(CompressionMethod::Stored)
    } else {
        SimpleFileOptions::default()
            .compression_method(CompressionMethod::Deflated)
            .compression_level(Some(compression_level as i64))
    };
    let mut npz = NpzWriter::new_with_options(File::create(filename)?, options);
    npz.add_array("equil", &equil_stacked)?;
    npz.add_array("steps", &steps_stacked)?;
    npz.finish()?;

    Ok(())
}

#[cfg(test)]
#[path = "snapshots_tests.rs"]
mod tests;
