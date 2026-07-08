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

    #[serde(default = "default_pt_interval")]
    pub pt_interval: usize,

    pub algorithm: Algorithm,
    #[serde(default = "default_boltzmann_constant")]
    pub boltzmann_constant: f64,
}
fn default_pt_interval() -> usize {
    0
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
        if self
            .equilibration_steps
            .checked_add(self.measurement_steps)
            .is_none()
        {
            anyhow::bail!("equilibration_steps plus measurement_steps is too large");
        }
        if self.measurement_steps == 0 {
            anyhow::bail!("measurement_steps must be greater than zero");
        }
        if self.num_threads == 0 {
            anyhow::bail!("num_threads must be greater than zero");
        }
        if !self.boltzmann_constant.is_finite() || self.boltzmann_constant <= 0.0 {
            anyhow::bail!(
                "boltzmann_constant ({}) must be finite and greater than zero",
                self.boltzmann_constant
            );
        }

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
                    if !start.is_finite() {
                        anyhow::bail!("temperature_range start ({start}) must be finite");
                    }
                    if !end.is_finite() {
                        anyhow::bail!("temperature_range end ({end}) must be finite");
                    }
                    if !step.is_finite() {
                        anyhow::bail!("temperature_range step ({step}) must be finite");
                    }
                    if start <= 0.0 {
                        anyhow::bail!(
                            "temperature_range start ({start}) must be greater than zero"
                        );
                    }
                    if step <= 0.0 {
                        anyhow::bail!("temperature_range step ({step}) must be positive");
                    }
                    if start > end {
                        anyhow::bail!(
                            "temperature_range start ({start}) must be less than or equal to end ({end})"
                        );
                    }
                    let mut t = start;
                    while t <= end + 1e-8 {
                        self.temperatures.push(t);
                        t += step;
                    }
                }
                if self.pt_interval > 0 && self.temperatures.len() < 2 {
                    anyhow::bail!(
                        "parallel tempering requires at least two temperatures when pt_interval is greater than zero"
                    );
                }
                Ok(())
            }
            (false, true) => {
                for (index, temperature) in self.temperatures.iter().enumerate() {
                    if !temperature.is_finite() || *temperature <= 0.0 {
                        anyhow::bail!(
                            "temperatures[{index}] ({temperature}) must be finite and greater than zero"
                        );
                    }
                }
                if self.pt_interval > 0 && self.temperatures.len() < 2 {
                    anyhow::bail!(
                        "parallel tempering requires at least two temperatures when pt_interval is greater than zero"
                    );
                }
                Ok(())
            }
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
        if self.pt_interval > 0 {
            writeln!(
                f,
                "  PT (Parallel Tempering): enabled, swap every {} sweeps",
                self.pt_interval
            )?;
        } else {
            writeln!(f, "  PT (Parallel Tempering): disabled")?;
        }
        write!(f, "  Temperatures (K):\n  ")?;
        for t in &self.temperatures {
            write!(f, "{t:.4}   ")?;
        }
        writeln!(f)?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "simulation_tests.rs"]
mod tests;
