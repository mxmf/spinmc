use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Deserialize, Serialize)]
pub struct Grid {
    pub dimensions: [usize; 3],
    pub sublattices: usize,
    pub spin_magnitudes: Vec<f64>,
    pub periodic_boundary: [bool; 3],
}

impl Grid {
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.sublattices != self.spin_magnitudes.len() {
            anyhow::bail!(
                "spin_magnitude length ({}) does not match sublattices ({})",
                self.spin_magnitudes.len(),
                self.sublattices
            );
        }
        Ok(())
    }
}
impl fmt::Display for Grid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\nGrid:")?;
        writeln!(
            f,
            "  Dimensions: {:?}\n  Periodic Boundary: {:?}",
            self.dimensions, self.periodic_boundary
        )?;
        writeln!(f, "  Sublattices: {}", self.sublattices)?;
        writeln!(f, "  Spin Magnitudes: {:?}", self.spin_magnitudes)?;
        Ok(())
    }
}
