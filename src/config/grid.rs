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
        for (axis, dimension) in self.dimensions.iter().enumerate() {
            if *dimension == 0 {
                anyhow::bail!("dimensions[{axis}] must be greater than zero");
            }
        }

        if self.sublattices == 0 {
            anyhow::bail!("sublattices must be greater than zero");
        }

        if self.sublattices != self.spin_magnitudes.len() {
            anyhow::bail!(
                "spin_magnitude length ({}) does not match sublattices ({})",
                self.spin_magnitudes.len(),
                self.sublattices
            );
        }

        for (index, magnitude) in self.spin_magnitudes.iter().enumerate() {
            if !magnitude.is_finite() || *magnitude <= 0.0 {
                anyhow::bail!(
                    "spin_magnitudes[{index}] ({magnitude}) must be finite and greater than zero"
                );
            }
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

#[cfg(test)]
#[path = "grid_tests.rs"]
mod tests;
