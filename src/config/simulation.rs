use super::{Algorithm, InitialState, Model};
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Deserialize, Serialize)]
pub struct Simulation {
    pub initial_state: InitialState,
    pub model: Model,
    pub equilibration_steps: usize,
    pub measurement_steps: usize,

    #[serde(default)]
    pub temperatures: Vec<f64>,

    #[serde(default)]
    pub temperature_range: Vec<TemperatureRange>,

    pub num_threads: usize,
    pub algorithm: Algorithm,
    #[serde(default = "default_boltzmann_constant")]
    pub boltzmann_constant: f64,
}
fn default_boltzmann_constant() -> f64 {
    8.617333262145e-5 // eV/K
}

#[derive(Debug, Deserialize, Serialize)]
pub struct TemperatureRange {
    pub start: f64,
    pub end: f64,
    pub step: f64,
}

impl Simulation {
    pub fn validate(&mut self) -> anyhow::Result<()> {
        match (
            self.temperatures.is_empty(),
            self.temperature_range.is_empty(),
        ) {
            (true, true) => {
                anyhow::bail!("Either 'temperatures' or 'temperature_range' must be specified");
            }
            (false, false) => {
                anyhow::bail!(
                    "Only one of 'temperatures' or 'temperature_range' can be specified, not both"
                )
            }
            (true, false) => {
                for tem_range in &self.temperature_range {
                    let (start, end, step) = (tem_range.start, tem_range.end, tem_range.step);
                    let mut t = start;
                    while t <= end + 1e-8 {
                        self.temperatures.push(t);
                        t += step;
                    }
                }
                Ok(())
            }
            (false, true) => Ok(()),
        }
    }
}

impl fmt::Display for Simulation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\nSimulation Parameters:")?;
        writeln!(f, "  Initial State: {:?}", self.initial_state)?;
        writeln!(f, "  Model: {:?}", self.model)?;
        writeln!(f, "  Equilibration Steps: {}", self.equilibration_steps)?;
        writeln!(f, "  Simulation Steps: {}", self.measurement_steps)?;
        writeln!(
            f,
            "  Boltzmann Constant (kB): {} (eV/K)",
            self.boltzmann_constant
        )?;
        writeln!(f, "  Algorithm: {:?}", self.algorithm)?;
        writeln!(f, "  Threads: {}", self.num_threads)?;
        write!(f, "  Temperatures (K):\n  ")?;
        for t in &self.temperatures {
            write!(f, "{t:.4}   ")?;
        }
        writeln!(f)?;
        Ok(())
    }
}
