use std::fmt;

use crate::spin::SpinState;
use ndarray::Array4;
use serde::{Deserialize, Serialize};

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

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct Snapshots {
    pub snapshot_equil_interval: Option<usize>,
    pub snapshot_steps_interval: Option<usize>,
    pub snapshot_compression: Option<usize>,
    pub snapshot_dir: Option<String>,
}

impl Snapshots {
    pub fn parse(&self, n_equil: usize, n_steps: usize) -> ParsedSnapshots {
        ParsedSnapshots {
            snapshot_equil_interval: self.snapshot_equil_interval.unwrap_or(n_equil),
            snapshot_steps_interval: self.snapshot_steps_interval.unwrap_or(n_steps),
            snapshot_compression: self.snapshot_compression.unwrap_or_default(),
            snapshot_dir: self.snapshot_dir.clone().unwrap_or("snapshots".to_string()),
        }
    }
}

#[derive(Debug, Deserialize, Serialize, Clone, Default)]
pub struct ParsedSnapshots {
    pub snapshot_equil_interval: usize,
    pub snapshot_steps_interval: usize,
    pub snapshot_compression: usize,
    pub snapshot_dir: String,
}

impl fmt::Display for ParsedSnapshots {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\nSnapshots: Enable")?;
        writeln!(
            f,
            "  Equilibration Interval: {} steps",
            self.snapshot_equil_interval
        )?;
        writeln!(
            f,
            "  Sampling Interval: {} steps",
            self.snapshot_steps_interval
        )?;
        writeln!(f, "  Compression Level: {}", self.snapshot_compression)?;
        writeln!(f, "  Snapshot Directory: {}", self.snapshot_dir)?;

        Ok(())
    }
}
