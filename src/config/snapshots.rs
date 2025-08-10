use crate::spin::SpinState;
use ndarray::Array4;
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

pub fn save_snapshots_to_hdf5<S: SpinState + hdf5_metno::H5Type>(
    filename: &str,
    equil_data: &[Array4<S>],
    steps_data: &[Array4<S>],
) -> hdf5_metno::Result<()> {
    let file = hdf5_metno::File::create(filename)?;
    let equil_views: Vec<_> = equil_data.iter().map(|a| a.view()).collect();
    let equil_stacked = ndarray::stack(ndarray::Axis(0), &equil_views)?;
    let steps_views: Vec<_> = steps_data.iter().map(|a| a.view()).collect();
    let steps_stacked = ndarray::stack(ndarray::Axis(0), &steps_views)?;

    let _equil_ds = file
        .new_dataset_builder()
        .with_data(&equil_stacked)
        .deflate(9)
        .create("snapshots/equil")?;
    let _steps_ds = file
        .new_dataset_builder()
        .with_data(&steps_stacked)
        .deflate(9)
        .create("snapshots/steps")?;
    Ok(())
}
