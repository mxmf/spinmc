use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Anisotropy {
    pub axis: Vec<[f64; 3]>,
    pub strength: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct ParsedAnisotropy {
    pub axis: [f64; 3],
    pub strength: f64,
}

impl Anisotropy {
    pub fn validate(&self, sublattices: usize) -> anyhow::Result<()> {
        if self.axis.len() != self.strength.len() {
            anyhow::bail!("anisotropy axis and strength arrays must have the same length");
        }
        if self.strength.len() != sublattices {
            anyhow::bail!(
                "anisotropy strength arrays must have the same length with {sublattices}"
            );
        }
        Ok(())
    }

    pub fn parse(&self) -> anyhow::Result<Vec<ParsedAnisotropy>> {
        let mut result = vec![];

        for (saxis, strength) in self.axis.iter().zip(self.strength.iter()) {
            let saxis_norm =
                (saxis[0] * saxis[0] + saxis[1] * saxis[1] + saxis[2] * saxis[2]).sqrt();
            if saxis_norm == 0.0 {
                anyhow::bail!("Anisotropy direction vector {:?} has zero length", saxis);
            }

            let ani = ParsedAnisotropy {
                axis: [
                    saxis[0] / saxis_norm,
                    saxis[1] / saxis_norm,
                    saxis[2] / saxis_norm,
                ],
                strength: *strength,
            };
            result.push(ani);
        }

        Ok(result)
    }
}

impl fmt::Display for ParsedAnisotropy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [x, y, z] = self.axis;
        writeln!(f, "{x:>4} {y:>4} {z:>4} | {:>8.12}", self.strength)?;

        Ok(())
    }
}
