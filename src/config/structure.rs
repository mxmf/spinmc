use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Deserialize, Serialize)]
pub struct Structure {
    pub cell: [[f64; 3]; 3],
    pub positions: Vec<[f64; 3]>,
    pub tolerance: Option<f64>,
}

impl Structure {
    pub fn validate(&self, sublattices: usize) -> anyhow::Result<()> {
        if self.positions.len() != sublattices {
            anyhow::bail!(
                "Number of provided positions ({}) does not match the number of sublattices ({sublattices})",
                self.positions.len()
            );
        }
        Ok(())
    }
}
impl fmt::Display for Structure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\n Structure:")?;
        writeln!(f, "  Cell:")?;
        writeln!(
            f,
            "    {}  {}  {}",
            self.cell[0][0], self.cell[0][1], self.cell[0][2]
        )?;
        writeln!(
            f,
            "    {}  {}  {}",
            self.cell[1][0], self.cell[1][1], self.cell[1][2]
        )?;
        writeln!(
            f,
            "    {}  {}  {}",
            self.cell[2][0], self.cell[2][1], self.cell[2][2]
        )?;

        writeln!(f, "  Positions:")?;
        for position in &self.positions {
            writeln!(f, "    {} {} {}", position[0], position[1], position[2])?;
        }
        Ok(())
    }
}
