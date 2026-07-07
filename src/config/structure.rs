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
            (Some(file), None, None) => {
                let mut structure = lattice::load_from_file(file, self.format.clone())?;
                structure.tolerance = self.tolerance;
                if let Some(indices) = &self.magnetic_indices {
                    let indices: Vec<usize> = indices.iter().map(|&i| i as usize).collect();

                    let mut seen = std::collections::HashSet::new();
                    for &i in &indices {
                        if i >= structure.positions.len() {
                            anyhow::bail!(
                                "magnetic_indices contains index {i} which is out of range, file has only {} atoms",
                                structure.positions.len()
                            );
                        }
                        if !seen.insert(i) {
                            anyhow::bail!("magnetic_indices contains duplicate index {i}");
                        }
                    }

                    structure.positions = indices.iter().map(|&i| structure.positions[i]).collect();
                    structure.magnetic_indices = None;

                }
                Ok(structure)
            }
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
        if let (None, Some(_), Some(positions)) = (&self.file, &self.cell, &self.positions)
            && positions.len() != sublattices
        {
            anyhow::bail!(
                "number of provided positions ({}) does not match the number of sublattices ({sublattices})",
                positions.len()
            );
        }

        if let Some(indices) = &self.magnetic_indices
            && indices.len() != sublattices
        {
            anyhow::bail!(
                "number of magnetic_indices ({}) does not match the number of sublattices ({sublattices})",
                indices.len()
            );
        }

        // File mode without magnetic_indices: all atoms are magnetic,
        // so positions count must match sublattices.
        if let (Some(_), None) = (&self.file, &self.magnetic_indices) {
            let structure = self.parse()?;
            if structure.positions.len() != sublattices {
                anyhow::bail!(
                    "number of atoms in file ({}) does not match sublattices ({sublattices}); \
                     use `magnetic_indices` to select a subset if the file \
                     contains non-magnetic atoms",
                    structure.positions.len()
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
                writeln!(
                    f,
                    "    {:.8}  {:.8}  {:.8}",
                    cell[0][0], cell[0][1], cell[0][2]
                )?;
                writeln!(
                    f,
                    "    {:.8}  {:.8}  {:.8}",
                    cell[1][0], cell[1][1], cell[1][2]
                )?;
                writeln!(
                    f,
                    "    {:.8}  {:.8}  {:.8}",
                    cell[2][0], cell[2][1], cell[2][2]
                )?;
                writeln!(f, "  Positions:")?;
                for position in &structure.positions {
                    writeln!(
                        f,
                        "    {:.8} {:.8} {:.8}",
                        position[0], position[1], position[2]
                    )?;
                }
            }
            Err(e) => {
                writeln!(f, "  [Failed to parse: {e}]")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "structure_tests.rs"]
mod tests;
