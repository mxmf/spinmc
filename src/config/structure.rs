use serde::{Deserialize, Serialize};
use std::fmt;

use crate::lattice::{self, FullStructure, Structure, StructureAtom};

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
    fn validate_cell(cell: [[f64; 3]; 3]) -> anyhow::Result<()> {
        for (row_index, row) in cell.iter().enumerate() {
            for (component_index, component) in row.iter().enumerate() {
                if !component.is_finite() {
                    anyhow::bail!(
                        "cell[{row_index}][{component_index}] ({component}) must be finite"
                    );
                }
            }
        }

        let volume = cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
            - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
            + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0]);
        if volume.abs() <= f64::EPSILON {
            anyhow::bail!("cell vectors must define a non-zero volume");
        }

        Ok(())
    }

    fn validate_positions(positions: &[[f64; 3]]) -> anyhow::Result<()> {
        if positions.is_empty() {
            anyhow::bail!("positions must contain at least one position");
        }

        for (position_index, position) in positions.iter().enumerate() {
            for (component_index, component) in position.iter().enumerate() {
                if !component.is_finite() {
                    anyhow::bail!(
                        "positions[{position_index}][{component_index}] ({component}) must be finite"
                    );
                }
            }
        }

        Ok(())
    }

    fn validate_magnetic_indices(
        indices: &[u64],
        positions_len: Option<usize>,
    ) -> anyhow::Result<()> {
        let mut seen = std::collections::HashSet::new();
        for &index in indices {
            if !seen.insert(index) {
                anyhow::bail!("magnetic_indices contains duplicate index {index}");
            }
            if let Some(positions_len) = positions_len
                && index as usize >= positions_len
            {
                anyhow::bail!(
                    "magnetic_indices contains index {index} which is out of range, only {positions_len} atoms provided"
                );
            }
        }

        Ok(())
    }

    fn validate_format(format: &str) -> anyhow::Result<()> {
        if format.trim().is_empty() {
            anyhow::bail!("format must not be empty");
        }
        if !matches!(
            format.to_ascii_lowercase().as_str(),
            "poscar" | "vasp" | "contcar"
        ) {
            anyhow::bail!(
                "unsupported structure format `{format}`; only POSCAR/CONTCAR/VASP are supported"
            );
        }
        Ok(())
    }

    pub fn parse(&self) -> anyhow::Result<Structure> {
        match (&self.file, &self.cell, &self.positions) {
            (Some(file), None, None) => {
                let mut structure = lattice::load_from_file(file, self.format.clone())?;
                structure.tolerance = self.tolerance;
                if let Some(indices) = &self.magnetic_indices {
                    structure.magnetic_indices = indices.clone();
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
                }
                Ok(structure)
            }
            (None, Some(cell), Some(positions)) => {
                if self.format.is_some() {
                    anyhow::bail!("`format` is only valid when loading from a `file`");
                }
                let mut structure = Structure {
                    cell: *cell,
                    positions: positions.to_vec(),
                    tolerance: self.tolerance,
                    magnetic_indices: (0..positions.len() as u64).collect(),
                    full_structure: Some(FullStructure {
                        cell: *cell,
                        atoms: positions
                            .iter()
                            .copied()
                            .map(|position| StructureAtom {
                                element: "X".to_string(),
                                position,
                            })
                            .collect(),
                    }),
                };
                if let Some(indices) = &self.magnetic_indices {
                    structure.magnetic_indices = indices.clone();
                    let indices: Vec<usize> = indices.iter().map(|&i| i as usize).collect();

                    let mut seen = std::collections::HashSet::new();
                    for &i in &indices {
                        if i >= structure.positions.len() {
                            anyhow::bail!(
                                "magnetic_indices contains index {i} which is out of range, only {} atoms provided",
                                structure.positions.len()
                            );
                        }
                        if !seen.insert(i) {
                            anyhow::bail!("magnetic_indices contains duplicate index {i}");
                        }
                    }

                    structure.positions = indices.iter().map(|&i| structure.positions[i]).collect();
                }
                Ok(structure)
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
        if let Some(file) = &self.file
            && file.trim().is_empty()
        {
            anyhow::bail!("structure file path must not be empty");
        }

        if let Some(format) = &self.format {
            Self::validate_format(format)?;
            if self.file.is_none() {
                anyhow::bail!("`format` is only valid when loading from a `file`");
            }
        }

        if let Some(tolerance) = self.tolerance
            && (!tolerance.is_finite() || tolerance < 0.0)
        {
            anyhow::bail!("tolerance ({tolerance}) must be finite and non-negative");
        }

        if let Some(cell) = self.cell {
            Self::validate_cell(cell)?;
        }

        if let Some(positions) = &self.positions {
            Self::validate_positions(positions)?;
        }

        if let (Some(indices), Some(positions)) = (&self.magnetic_indices, &self.positions) {
            Self::validate_magnetic_indices(indices, Some(positions.len()))?;
        } else if let Some(indices) = &self.magnetic_indices {
            Self::validate_magnetic_indices(indices, None)?;
        }

        if let (None, Some(_), Some(positions)) = (&self.file, &self.cell, &self.positions)
            && self.magnetic_indices.is_none()
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

        // Without magnetic_indices: all atoms are magnetic,
        // so positions count must match sublattices.
        if self.magnetic_indices.is_none() {
            let structure = self.parse()?;
            if structure.positions.len() != sublattices {
                anyhow::bail!(
                    "number of atoms ({}) does not match sublattices ({sublattices}); \
                     use `magnetic_indices` to select a subset if there \
                     are non-magnetic atoms",
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
