use serde::{Deserialize, Serialize};
use std::fmt;

use crate::lattice::{self, Structure};

#[derive(Debug, Deserialize, Serialize)]
pub struct StructureConf {
    pub file: Option<String>,
    pub format: Option<String>,
    pub cell: Option<[[f64; 3]; 3]>,
    pub positions: Option<Vec<[f64; 3]>>,
    pub tolerance: Option<f64>,
    pub magnetic_indices: Option<Vec<u64>>,
}

impl StructureConf {
    pub fn parse(&self) -> anyhow::Result<Structure> {
        match (&self.file, &self.cell, &self.positions) {
            (Some(file), None, None) => lattice::load_from_file(file, self.format.clone()),
            (None, Some(cell), Some(positions)) => {
                if self.magnetic_indices.is_some() {
                    anyhow::bail!("`magnetic_indices` is only valid when loading from a `file`");
                }
                if self.format.is_some() {
                    anyhow::bail!("`format` is only valid when loading from a `file`");
                }
                Ok(Structure {
                    cell: *cell,
                    positions: positions.to_vec(),
                    tolerance: self.tolerance,
                    magnetic_indices: None,
                    frame: None,
                })
            }
            (None, None, None) => anyhow::bail!(
                "structure must be provided: specify either `file` or `cell` and `positions`"
            ),
            (Some(_), Some(_), _) => anyhow::bail!("cannot specify both `file` and `cell`"),
            (Some(_), _, Some(_)) => anyhow::bail!("cannot specify both `file` and `positions`"),
            (None, Some(_), None) => {
                anyhow::bail!("missing `positions`: required when `cell` is set")
            }
            (None, None, Some(_)) => {
                anyhow::bail!("missing `cell`: required when `positions` is set")
            }
        }
    }

    pub fn validate(&self, sublattices: usize) -> anyhow::Result<()> {
        if let (None, Some(_), Some(positions)) = (&self.file, &self.cell, &self.positions) {
            if positions.len() != sublattices {
                anyhow::bail!(
                    "number of provided positions ({}) does not match the number of sublattices ({sublattices})",
                    positions.len()
                );
            }
        }

        if let Some(indices) = &self.magnetic_indices {
            for &idx in indices {
                if idx as usize >= sublattices {
                    anyhow::bail!(
                        "magnetic_indices contains index {idx} which is out of range, must be less than sublattices count ({sublattices})"
                    );
                }
            }
            if indices.len() != sublattices {
                anyhow::bail!(
                    "number of magnetic_indices ({}) does not match the number of sublattices ({sublattices})",
                    indices.len()
                );
            }
        }

        Ok(())
    }
}

impl fmt::Display for StructureConf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\n Structure:")?;
        if let Some(file) = &self.file {
            writeln!(f, "  Mode: from file")?;
            writeln!(f, "    File: {file}")?;
            if let Some(format) = &self.format {
                writeln!(f, "    Format: {format}")?;
            }
        } else {
            writeln!(f, "  Mode: direct input")?;
        }

        match self.parse() {
            Ok(structure) => {
                writeln!(f, "  Cell:")?;
                let cell = structure.cell;
                writeln!(f, "    {}  {}  {}", cell[0][0], cell[0][1], cell[0][2])?;
                writeln!(f, "    {}  {}  {}", cell[1][0], cell[1][1], cell[1][2])?;
                writeln!(f, "    {}  {}  {}", cell[2][0], cell[2][1], cell[2][2])?;
                writeln!(f, "  Positions:")?;
                for position in &structure.positions {
                    writeln!(f, "    {} {} {}", position[0], position[1], position[2])?;
                }
            }
            Err(e) => {
                writeln!(f, "  [Failed to parse: {e}]")?;
            }
        }
        Ok(())
    }
}
